use anyhow::{bail, Context, Result};
use clap::{Parser, Subcommand};
use gmg_core::{
    downsample, guess, FracSketch, GenomeMeta, IndexBuilder, RefIndex, SketchFile, Thresholds,
    DEFAULT_SEED, RANKS,
};
use needletail::parse_fastx_file;
mod run;

use rayon::prelude::*;
use std::fs::File;
use std::io::{BufReader, BufWriter, Write};
use std::path::{Path, PathBuf};

#[derive(Parser)]
#[command(name = "gmg", version = concat!(env!("CARGO_PKG_VERSION"), " (", env!("GMG_GIT"), ")"), about = "Guess My Genome")]
struct Cli {
    /// worker threads (default: all cores)
    #[arg(long, global = true)]
    threads: Option<usize>,
    #[command(subcommand)]
    cmd: Cmd,
}

#[derive(Subcommand)]
enum Cmd {
    /// Sketch FASTA/FASTQ files (optionally gzipped) into .gsk files
    Sketch {
        #[arg(short, long, default_value_t = 21)]
        ksize: u32,
        #[arg(short, long, default_value_t = 1000)]
        scaled: u32,
        /// output directory; one <id>.gsk per input, existing ones are kept unless --force
        #[arg(short, long)]
        out_dir: PathBuf,
        /// track abundances (reads)
        #[arg(long)]
        abund: bool,
        /// stop after this many reads/records per input (subsample the start of a file)
        #[arg(long)]
        max_reads: Option<u64>,
        #[arg(long)]
        force: bool,
        /// file with one input path per line
        #[arg(long)]
        list: Option<PathBuf>,
        inputs: Vec<PathBuf>,
    },
    /// Build a reference index from a metadata table and a directory of .gsk files
    Build {
        /// TSV with header: id taxid phylum class order family genus species [cluster]
        #[arg(long)]
        meta: PathBuf,
        #[arg(long)]
        sketches: PathBuf,
        #[arg(short, long)]
        out: PathBuf,
        /// build at a coarser scaled than the sketches (sketches are downsampled)
        #[arg(long)]
        scaled: Option<u32>,
        #[arg(long, default_value_t = 64)]
        buckets: usize,
    },
    /// Guess the taxonomy of one or more queries (.gsk, FASTA or FASTQ)
    Guess {
        #[arg(long)]
        index: PathBuf,
        /// drop query hashes seen fewer than this many times (reads)
        #[arg(long, default_value_t = 1)]
        min_abund: u32,
        #[arg(long)]
        max_reads: Option<u64>,
        #[arg(long, default_value_t = 5)]
        top: usize,
        /// `known` (default), `novel`, or a JSON file (min_shared[6], min_ani[6], min_consistency[6], margin, top_n)
        #[arg(long)]
        thresholds: Option<String>,
        queries: Vec<PathBuf>,
    },
    /// Classify every genome / read set in a folder: writes results.tsv, hits.tsv and run.json
    Run {
        /// directories (searched recursively) and/or files: FASTA (.fa .fna .fasta .fas .ffn) or FASTQ (.fq .fastq), optionally .gz
        #[arg(required = true)]
        inputs: Vec<PathBuf>,
        #[arg(long)]
        index: PathBuf,
        /// output directory
        #[arg(short, long)]
        out: PathBuf,
        /// `known` (default), `novel` (queries that may belong to a genus absent from the reference), or a JSON file
        #[arg(long)]
        thresholds: Option<String>,
        /// auto (by extension), assembly or reads
        #[arg(long, default_value = "auto", value_parser = ["auto", "assembly", "reads"])]
        kind: String,
        /// drop hashes seen fewer than this many times (default 1 for assemblies, 2 for reads)
        #[arg(long)]
        min_abund: Option<u32>,
        /// reads: stop after this many reads per file
        #[arg(long)]
        max_reads: Option<u64>,
        /// reference hits per input in hits.tsv
        #[arg(long, default_value_t = 5)]
        top: usize,
        /// reuse finished results from a previous run in the same output directory
        #[arg(long)]
        resume: bool,
    },
    /// Summarise an index
    Info { index: PathBuf },
    /// Greedy ANI clustering of the genomes in an index (centroids picked in priority order)
    Cluster {
        #[arg(long)]
        index: PathBuf,
        #[arg(long)]
        sketches: PathBuf,
        /// a genome joins a centroid's cluster if its containment-ANI in the centroid is at least this
        #[arg(long, default_value_t = 0.95)]
        ani: f64,
        /// file with genome ids, highest priority first (e.g. BFD rule: BUSCO, N50); others follow in index order
        #[arg(long)]
        order: Option<PathBuf>,
        /// a genome joins a centroid only if the centroid has at most this many times its hashes
        /// (a much larger centroid is probably a mixed or contaminated assembly that "contains" its neighbours)
        #[arg(long, default_value_t = 1.6)]
        max_size_ratio: f64,
        #[arg(short, long)]
        out: PathBuf,
    },
    /// Rewrite an index in the compact v1 (delta + varint) format
    Convert {
        index: PathBuf,
        #[arg(short, long)]
        out: PathBuf,
    },
    /// Leave-out benchmark of genomes that are in the index
    Bench {
        #[arg(long)]
        index: PathBuf,
        #[arg(long)]
        sketches: PathBuf,
        /// loo = leave the genome out; los = leave its species out; log = leave its genus out
        #[arg(long, default_value = "loo")]
        mode: String,
        /// restrict queries to this phylum
        #[arg(long)]
        phylum: Option<String>,
        /// only genomes whose id is listed (one per line) may be references (simulates a smaller, representatives-only index)
        #[arg(long)]
        ref_list: Option<PathBuf>,
        /// with --ref-list: query only genomes that are NOT in the reference list
        #[arg(long)]
        queries_outside: bool,
        /// random sample of this many queries
        #[arg(long)]
        sample: Option<usize>,
        #[arg(long, default_value_t = 1)]
        seed: u64,
        /// split queries into this many shards and run shard --shard-index (for SLURM arrays)
        #[arg(long, default_value_t = 1)]
        shards: usize,
        #[arg(long, default_value_t = 0)]
        shard_index: usize,
        #[arg(short, long)]
        out: PathBuf,
    },
}

pub(crate) fn strip_ext(name: &str) -> String {
    for e in [
        ".fa.gz",
        ".fasta.gz",
        ".fna.gz",
        ".fastq.gz",
        ".fq.gz",
        ".fa",
        ".fasta",
        ".fna",
        ".fastq",
        ".fq",
        ".gsk",
    ] {
        if let Some(s) = name.strip_suffix(e) {
            return s.to_string();
        }
    }
    name.to_string()
}

fn id_of(path: &Path) -> String {
    strip_ext(
        &path
            .file_name()
            .map(|s| s.to_string_lossy().into_owned())
            .unwrap_or_default(),
    )
}

/// Per-input counts gathered while sketching (used for QC flags).
pub(crate) struct Stats {
    pub n_records: u64,
    pub total_len: u64,
    /// bases other than A/C/G/T (case-insensitive)
    pub non_acgt: u64,
    pub n50: u64,
}

/// Sketch one file and collect basic stats. Returns None if it has no sequence records.
pub(crate) fn sketch_with_stats(
    path: &Path,
    k: u32,
    scaled: u32,
    abund: bool,
    max_reads: Option<u64>,
) -> Result<Option<(SketchFile, Stats)>> {
    let mut r = parse_fastx_file(path).with_context(|| format!("opening {}", path.display()))?;
    let mut sk = FracSketch::new(k as usize, scaled, DEFAULT_SEED);
    let (mut n, mut total, mut non_acgt) = (0u64, 0u64, 0u64);
    let mut lens: Vec<u64> = Vec::new();
    while let Some(rec) = r.next() {
        let rec = rec.with_context(|| format!("parsing {}", path.display()))?;
        let seq = rec.seq();
        sk.add_sequence(&seq);
        n += 1;
        total += seq.len() as u64;
        non_acgt += seq
            .iter()
            .filter(|b| !matches!(b.to_ascii_uppercase(), b'A' | b'C' | b'G' | b'T'))
            .count() as u64;
        if !abund {
            lens.push(seq.len() as u64);
        }
        if max_reads.is_some_and(|m| n >= m) {
            break;
        }
    }
    if n == 0 {
        return Ok(None);
    }
    lens.sort_unstable_by(|a, b| b.cmp(a));
    let (mut acc, mut n50) = (0u64, 0u64);
    for l in &lens {
        acc += l;
        if acc * 2 >= total {
            n50 = *l;
            break;
        }
    }
    let (hashes, abunds) = if abund {
        let (h, a) = sk.finish_with_abundance();
        (h, Some(a))
    } else {
        (sk.finish(), None)
    };
    let sf = SketchFile {
        ksize: k,
        scaled,
        seed: DEFAULT_SEED,
        name: id_of(path),
        hashes,
        abunds,
    };
    Ok(Some((
        sf,
        Stats {
            n_records: n,
            total_len: total,
            non_acgt,
            n50,
        },
    )))
}

/// Sketch one file. Returns None if it has no sequence records.
fn sketch_file(
    path: &Path,
    k: u32,
    scaled: u32,
    abund: bool,
    max_reads: Option<u64>,
) -> Result<Option<SketchFile>> {
    Ok(sketch_with_stats(path, k, scaled, abund, max_reads)?.map(|(s, _)| s))
}

fn write_sketch(path: &Path, s: &SketchFile) -> Result<()> {
    let tmp = path.with_extension("gsk.tmp");
    let mut w = BufWriter::new(File::create(&tmp)?);
    s.write(&mut w)?;
    w.flush()?;
    std::fs::rename(&tmp, path)?;
    Ok(())
}

fn read_sketch(path: &Path) -> Result<SketchFile> {
    let f = File::open(path).with_context(|| format!("opening {}", path.display()))?;
    Ok(SketchFile::read(&mut BufReader::new(f))?)
}

pub(crate) fn load_index(path: &Path) -> Result<RefIndex> {
    let f = File::open(path).with_context(|| format!("opening {}", path.display()))?;
    Ok(RefIndex::read(&mut BufReader::with_capacity(1 << 20, f))?)
}

fn read_meta(path: &Path) -> Result<Vec<GenomeMeta>> {
    let text = std::fs::read_to_string(path)?;
    let mut out = Vec::new();
    for (i, line) in text.lines().enumerate() {
        if i == 0 || line.trim().is_empty() {
            continue;
        }
        let f: Vec<&str> = line.split('\t').collect();
        if f.len() != 8 && f.len() != 9 {
            bail!(
                "{}:{}: expected 8 tab-separated columns (+ optional cluster), got {}",
                path.display(),
                i + 1,
                f.len()
            );
        }
        let mut lineage: [String; 6] = Default::default();
        for r in 0..6 {
            lineage[r] = f[2 + r].to_string();
        }
        out.push(GenomeMeta {
            id: f[0].to_string(),
            taxid: f[1].parse().unwrap_or(0),
            lineage,
            n_hashes: 0,
            cluster: f.get(8).map(|c| c.to_string()).unwrap_or_default(),
        });
    }
    Ok(out)
}

/// `known` (default), `novel` (stricter: queries that may belong to a genus absent from the reference) or a JSON file.
pub(crate) fn load_thresholds(spec: Option<&str>) -> Result<(Thresholds, String)> {
    match spec {
        None | Some("known") => Ok((Thresholds::default(), "known".to_string())),
        Some("novel") => Ok((
            Thresholds::from_json_str(include_str!(
                "../../../data/thresholds/bfd_s10000.novel.json"
            ))
            .map_err(|e| anyhow::anyhow!(e))?,
            "novel".to_string(),
        )),
        Some(p) => Ok((read_thresholds(Path::new(p))?, p.to_string())),
    }
}

fn read_thresholds(path: &Path) -> Result<Thresholds> {
    Thresholds::from_json_str(&std::fs::read_to_string(path)?)
        .map_err(|e| anyhow::anyhow!("{}: {e}", path.display()))
}

fn main() -> Result<()> {
    let cli = Cli::parse();
    if let Some(t) = cli.threads {
        rayon::ThreadPoolBuilder::new()
            .num_threads(t)
            .build_global()?;
    }
    match cli.cmd {
        Cmd::Sketch {
            ksize,
            scaled,
            out_dir,
            abund,
            max_reads,
            force,
            list,
            mut inputs,
        } => {
            if let Some(l) = list {
                for line in std::fs::read_to_string(l)?
                    .lines()
                    .filter(|l| !l.trim().is_empty())
                {
                    inputs.push(PathBuf::from(line.trim()));
                }
            }
            std::fs::create_dir_all(&out_dir)?;
            let results: Vec<(String, &'static str, String)> = inputs
                .par_iter()
                .map(|p| {
                    let out = out_dir.join(format!("{}.gsk", id_of(p)));
                    if out.exists() && !force {
                        return (p.display().to_string(), "EXISTS", String::new());
                    }
                    match sketch_file(p, ksize, scaled, abund, max_reads) {
                        Ok(Some(s)) => match write_sketch(&out, &s) {
                            Ok(()) => (p.display().to_string(), "OK", s.hashes.len().to_string()),
                            Err(e) => (p.display().to_string(), "FAIL", format!("{e:#}")),
                        },
                        Ok(None) => (
                            p.display().to_string(),
                            "EMPTY",
                            "no sequence records".into(),
                        ),
                        Err(e) => (p.display().to_string(), "FAIL", format!("{e:#}")),
                    }
                })
                .collect();
            let mut counts = std::collections::BTreeMap::new();
            for (p, st, msg) in &results {
                *counts.entry(*st).or_insert(0usize) += 1;
                if *st == "EMPTY" || *st == "FAIL" {
                    eprintln!("{st}\t{p}\t{msg}");
                }
            }
            eprintln!("sketch summary: {counts:?}");
        }
        Cmd::Build {
            meta,
            sketches,
            out,
            scaled,
            buckets,
        } => {
            let metas = read_meta(&meta)?;
            let loaded: Vec<(GenomeMeta, Option<SketchFile>)> = metas
                .into_par_iter()
                .map(|m| {
                    let p = sketches.join(format!("{}.gsk", m.id));
                    let s = if p.exists() {
                        read_sketch(&p).ok()
                    } else {
                        None
                    };
                    (m, s)
                })
                .collect();
            let first = loaded
                .iter()
                .find_map(|(_, s)| s.as_ref())
                .context("no sketches found for any meta row")?;
            let (k, base) = (first.ksize, first.scaled);
            let target = scaled.unwrap_or(base);
            if target < base {
                bail!("--scaled {target} is finer than the sketches ({base}); re-sketch instead");
            }
            let mut b = IndexBuilder::new(k, target, DEFAULT_SEED);
            let (mut used, mut missing) = (0usize, 0usize);
            for (m, s) in loaded {
                match s {
                    Some(s) if s.ksize == k && s.scaled == base && s.seed == DEFAULT_SEED => {
                        let h = if target > base {
                            downsample(&s.hashes, target)
                        } else {
                            s.hashes
                        };
                        b.add(m, h);
                        used += 1;
                    }
                    Some(_) => {
                        eprintln!(
                            "MISMATCH\t{}\tk/scaled/seed differ from the first sketch",
                            m.id
                        );
                        missing += 1;
                    }
                    None => {
                        eprintln!("MISSING\t{}", m.id);
                        missing += 1;
                    }
                }
            }
            let ix = b.build(buckets);
            let mut w = BufWriter::new(File::create(&out)?);
            ix.write(&mut w)?;
            w.flush()?;
            eprintln!(
                "built {}: {used} genomes ({missing} skipped), k={k}, scaled={target}, {} unique hashes, {} postings",
                out.display(),
                ix.n_unique_hashes(),
                ix.n_postings()
            );
        }
        Cmd::Guess {
            index,
            min_abund,
            max_reads,
            top,
            thresholds,
            queries,
        } => {
            let ix = load_index(&index)?;
            let (th, _) = load_thresholds(thresholds.as_deref())?;
            for q in &queries {
                let hashes = query_hashes(&ix, q, min_abund, max_reads)?;
                let g = guess(&ix, &hashes, &|_| false, &th);
                println!(
                    "{}",
                    gmg_core::guess_json(&ix, &g, &q.display().to_string(), top)
                );
            }
        }
        Cmd::Cluster {
            index,
            sketches,
            ani,
            order,
            max_size_ratio,
            out,
        } => {
            let ix = load_index(&index)?;
            let n = ix.genomes.len();
            let by_id: std::collections::HashMap<&str, u32> = ix
                .genomes
                .iter()
                .enumerate()
                .map(|(i, g)| (g.id.as_str(), i as u32))
                .collect();
            let mut prio: Vec<u32> = Vec::with_capacity(n);
            let mut seen = vec![false; n];
            if let Some(p) = order {
                for l in std::fs::read_to_string(p)?.lines() {
                    if let Some(&g) = by_id.get(l.trim()) {
                        if !seen[g as usize] {
                            seen[g as usize] = true;
                            prio.push(g);
                        }
                    }
                }
            }
            for g in 0..n as u32 {
                if !seen[g as usize] {
                    prio.push(g);
                }
            }
            let k = ix.ksize as f64;
            let mut centroid: Vec<Option<(u32, f64)>> = vec![None; n]; // (centroid genome, member ANI)
            let mut n_clusters = 0usize;
            let mut n_too_big = 0usize;
            for &g in &prio {
                if centroid[g as usize].is_some() {
                    continue;
                }
                let m = &ix.genomes[g as usize];
                let sk = match read_sketch(&sketches.join(format!("{}.gsk", m.id))) {
                    Ok(s) => s,
                    Err(_) => {
                        centroid[g as usize] = Some((g, 1.0));
                        n_clusters += 1;
                        continue;
                    }
                };
                let q = if sk.scaled < ix.scaled {
                    downsample(&sk.hashes, ix.scaled)
                } else {
                    sk.hashes
                };
                centroid[g as usize] = Some((g, 1.0));
                n_clusters += 1;
                for h in ix.search(&q, &|x| centroid[x as usize].is_some()) {
                    let mem = &ix.genomes[h.genome as usize];
                    let a = (h.shared as f64 / mem.n_hashes.max(1) as f64).powf(1.0 / k);
                    if a >= ani {
                        if m.n_hashes as f64 > max_size_ratio * mem.n_hashes.max(1) as f64 {
                            n_too_big += 1;
                            continue;
                        }
                        centroid[h.genome as usize] = Some((g, a));
                    }
                }
            }
            let mut w = BufWriter::new(File::create(&out)?);
            writeln!(
                w,
                "id\tcluster\tis_centroid\tani_to_centroid\tspecies\tgenus\tn_hashes"
            )?;
            for (i, g) in ix.genomes.iter().enumerate() {
                let (c, a) = centroid[i].unwrap();
                writeln!(
                    w,
                    "{}\t{}\t{}\t{:.4}\t{}\t{}\t{}",
                    g.id,
                    ix.genomes[c as usize].id,
                    c as usize == i,
                    a,
                    g.lineage[5],
                    g.lineage[4],
                    g.n_hashes
                )?;
            }
            w.flush()?;
            eprintln!("{n} genomes -> {n_clusters} clusters at ANI >= {ani} ({}); {n_too_big} (centroid, genome) pairs above the ANI cutoff were refused by the size ratio", out.display());
        }
        Cmd::Convert { index, out } => {
            let ix = load_index(&index)?;
            let mut w = BufWriter::new(File::create(&out)?);
            ix.write_v1(&mut w)?;
            w.flush()?;
            let (a, b) = (
                std::fs::metadata(&index)?.len(),
                std::fs::metadata(&out)?.len(),
            );
            eprintln!(
                "{} -> {}: {} -> {} bytes ({:.2}x smaller)",
                index.display(),
                out.display(),
                a,
                b,
                a as f64 / b as f64
            );
        }
        Cmd::Run {
            inputs,
            index,
            out,
            thresholds,
            kind,
            min_abund,
            max_reads,
            top,
            resume,
        } => {
            run::run(run::RunArgs {
                inputs,
                index,
                out,
                thresholds,
                kind,
                min_abund,
                max_reads,
                top,
                resume,
            })?;
        }
        Cmd::Info { index } => {
            let ix = load_index(&index)?;
            let size = std::fs::metadata(&index)?.len();
            println!("file\t{}\t{} bytes", index.display(), size);
            println!("k\t{}\nscaled\t{}\nseed\t{}", ix.ksize, ix.scaled, ix.seed);
            println!("genomes\t{}", ix.genomes.len());
            println!(
                "unique_hashes\t{}\npostings\t{}",
                ix.n_unique_hashes(),
                ix.n_postings()
            );
            for (r, t) in RANKS.iter().zip(&ix.taxa) {
                println!("taxa_{r}\t{}", t.len() - 1);
            }
        }
        Cmd::Bench {
            index,
            sketches,
            mode,
            phylum,
            ref_list,
            queries_outside,
            sample,
            seed,
            shards,
            shard_index,
            out,
        } => {
            let ix = load_index(&index)?;
            if !["loo", "los", "log"].contains(&mode.as_str()) {
                bail!("--mode must be loo, los or log");
            }
            let allowed: Option<Vec<bool>> = match &ref_list {
                Some(p) => {
                    let ids: std::collections::HashSet<String> = std::fs::read_to_string(p)?
                        .lines()
                        .map(|l| l.trim().to_string())
                        .collect();
                    Some(ix.genomes.iter().map(|g| ids.contains(&g.id)).collect())
                }
                None => None,
            };
            let mut qs: Vec<u32> = (0..ix.genomes.len() as u32)
                .filter(|&g| {
                    phylum
                        .as_ref()
                        .is_none_or(|p| &ix.genomes[g as usize].lineage[0] == p)
                })
                .filter(|&g| !queries_outside || allowed.as_ref().is_none_or(|a| !a[g as usize]))
                .collect();
            // deterministic shuffle (xorshift) so samples and shards are reproducible
            let mut s = seed.wrapping_mul(0x9E3779B97F4A7C15) | 1;
            for i in (1..qs.len()).rev() {
                s ^= s << 13;
                s ^= s >> 7;
                s ^= s << 17;
                qs.swap(i, (s % (i as u64 + 1)) as usize);
            }
            if let Some(n) = sample {
                qs.truncate(n);
            }
            let qs: Vec<u32> = qs
                .into_iter()
                .enumerate()
                .filter(|(i, _)| i % shards == shard_index)
                .map(|(_, g)| g)
                .collect();
            let th = Thresholds::default();
            let rows: Vec<String> = qs
                .par_iter()
                .filter_map(|&gi| {
                    let m = &ix.genomes[gi as usize];
                    let sk = read_sketch(&sketches.join(format!("{}.gsk", m.id))).ok()?;
                    let q = if sk.scaled < ix.scaled {
                        downsample(&sk.hashes, ix.scaled)
                    } else {
                        sk.hashes
                    };
                    let mine = ix.lineage_ids[gi as usize];
                    let rk = match mode.as_str() {
                        "los" => 5,
                        "log" => 4,
                        _ => usize::MAX,
                    };
                    let skip = |g: u32| {
                        g == gi
                            || allowed.as_ref().is_some_and(|a| !a[g as usize])
                            || (rk != usize::MAX
                                && mine[rk] != 0
                                && ix.lineage_ids[g as usize][rk] == mine[rk])
                    };
                    let g = guess(&ix, &q, &skip, &th);
                    let (tid, tshared, tani, tcq, tcr) = match &g.top {
                        Some(t) => (
                            ix.genomes[t.genome as usize].id.clone(),
                            t.shared,
                            t.ani,
                            t.containment_query,
                            t.containment_ref,
                        ),
                        None => (String::new(), 0, 0.0, 0.0, 0.0),
                    };
                    let mut row = format!(
                        "{}\t{}\t{}\t{}\t{}\t{:.5}\t{:.5}\t{:.5}",
                        m.id,
                        g.query_hashes,
                        tid,
                        tshared,
                        tcq,
                        tcr,
                        tani,
                        ix.genomes.len()
                    );
                    let tl = g
                        .top
                        .as_ref()
                        .map(|t| &ix.genomes[t.genome as usize].lineage);
                    for r in 0..6 {
                        let truth = &m.lineage[r];
                        let top_taxon = tl.map(|l| l[r].as_str()).unwrap_or("");
                        let v = &g.votes[r];
                        row.push_str(&format!(
                            "\t{}\t{}\t{}\t{:.4}",
                            truth,
                            top_taxon,
                            v.taxon.clone().unwrap_or_default(),
                            v.consistency
                        ));
                    }
                    Some(row)
                })
                .collect();
            let mut w = BufWriter::new(File::create(&out)?);
            let mut header = String::from("id\tquery_hashes\ttop_id\ttop_shared\tcontainment_query\tcontainment_ref\tani\tn_ref");
            for r in RANKS {
                header.push_str(&format!("\t{r}_true\t{r}_top\t{r}_vote\t{r}_cons"));
            }
            writeln!(w, "{header}")?;
            for r in &rows {
                writeln!(w, "{r}")?;
            }
            w.flush()?;
            eprintln!(
                "bench {mode}: {} of {} queries scored -> {}",
                rows.len(),
                qs.len(),
                out.display()
            );
        }
    }
    Ok(())
}

/// Hashes of a query at the index's k/scaled, filtered by `min_abund`.
fn query_hashes(
    ix: &RefIndex,
    path: &Path,
    min_abund: u32,
    max_reads: Option<u64>,
) -> Result<Vec<u64>> {
    let sk = if path.extension().is_some_and(|e| e == "gsk") {
        read_sketch(path)?
    } else {
        sketch_file(path, ix.ksize, ix.scaled, min_abund > 1, max_reads)?
            .with_context(|| format!("{}: no sequence records", path.display()))?
    };
    query_hashes_from(ix, sk, min_abund)
}

/// Check a sketch against the index (k, seed, scaled), drop hashes seen fewer than `min_abund` times, downsample to the index.
pub(crate) fn query_hashes_from(ix: &RefIndex, sk: SketchFile, min_abund: u32) -> Result<Vec<u64>> {
    if sk.ksize != ix.ksize || sk.seed != ix.seed {
        bail!(
            "{}: k/seed ({}/{}) do not match the index ({}/{})",
            sk.name,
            sk.ksize,
            sk.seed,
            ix.ksize,
            ix.seed
        );
    }
    if sk.scaled > ix.scaled {
        bail!(
            "{}: query scaled {} is coarser than the index scaled {}",
            sk.name,
            sk.scaled,
            ix.scaled
        );
    }
    let mut pairs: Vec<(u64, u32)> = match &sk.abunds {
        Some(a) => sk.hashes.iter().copied().zip(a.iter().copied()).collect(),
        None => sk.hashes.iter().map(|&h| (h, 1)).collect(),
    };
    pairs.retain(|&(_, a)| a >= min_abund);
    let h: Vec<u64> = pairs.into_iter().map(|(h, _)| h).collect();
    Ok(if sk.scaled < ix.scaled {
        downsample(&h, ix.scaled)
    } else {
        h
    })
}
