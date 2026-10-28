use anyhow::{bail, Context, Result};
use clap::{Parser, Subcommand};
use gmg_core::{
    downsample, guess, FracSketch, GenomeMeta, IndexBuilder, RefIndex, SketchFile, Thresholds, DEFAULT_SEED, RANKS,
};
use needletail::parse_fastx_file;
use rayon::prelude::*;
use serde_json::json;
use std::fs::File;
use std::io::{BufReader, BufWriter, Write};
use std::path::{Path, PathBuf};

#[derive(Parser)]
#[command(name = "gmg", version, about = "Guess My Genome")]
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
        /// TSV with header: id taxid phylum class order family genus species
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
        /// JSON with min_shared[6], min_ani[6], min_consistency[6], margin, top_n (from scripts/calibrate.py)
        #[arg(long)]
        thresholds: Option<PathBuf>,
        queries: Vec<PathBuf>,
    },
    /// Summarise an index
    Info { index: PathBuf },
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

fn strip_ext(name: &str) -> String {
    for e in [".fa.gz", ".fasta.gz", ".fna.gz", ".fastq.gz", ".fq.gz", ".fa", ".fasta", ".fna", ".fastq", ".fq", ".gsk"] {
        if let Some(s) = name.strip_suffix(e) {
            return s.to_string();
        }
    }
    name.to_string()
}

fn id_of(path: &Path) -> String {
    strip_ext(&path.file_name().map(|s| s.to_string_lossy().into_owned()).unwrap_or_default())
}

/// Sketch one file. Returns None if it has no sequence records.
fn sketch_file(path: &Path, k: u32, scaled: u32, abund: bool, max_reads: Option<u64>) -> Result<Option<SketchFile>> {
    let mut r = parse_fastx_file(path).with_context(|| format!("opening {}", path.display()))?;
    let mut sk = FracSketch::new(k as usize, scaled, DEFAULT_SEED);
    let mut n = 0u64;
    while let Some(rec) = r.next() {
        let rec = rec.with_context(|| format!("parsing {}", path.display()))?;
        sk.add_sequence(&rec.seq());
        n += 1;
        if max_reads.is_some_and(|m| n >= m) {
            break;
        }
    }
    if n == 0 {
        return Ok(None);
    }
    let (hashes, abunds) = if abund {
        let (h, a) = sk.finish_with_abundance();
        (h, Some(a))
    } else {
        (sk.finish(), None)
    };
    Ok(Some(SketchFile { ksize: k, scaled, seed: DEFAULT_SEED, name: id_of(path), hashes, abunds }))
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

fn load_index(path: &Path) -> Result<RefIndex> {
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
        if f.len() != 8 {
            bail!("{}:{}: expected 8 tab-separated columns, got {}", path.display(), i + 1, f.len());
        }
        let mut lineage: [String; 6] = Default::default();
        for r in 0..6 {
            lineage[r] = f[2 + r].to_string();
        }
        out.push(GenomeMeta { id: f[0].to_string(), taxid: f[1].parse().unwrap_or(0), lineage, n_hashes: 0 });
    }
    Ok(out)
}

fn read_thresholds(path: &Path) -> Result<Thresholds> {
    let v: serde_json::Value = serde_json::from_str(&std::fs::read_to_string(path)?)?;
    let mut th = Thresholds::default();
    let arr = |k: &str| v.get(k).and_then(|a| a.as_array()).filter(|a| a.len() == 6);
    if let Some(a) = arr("min_shared") {
        for (i, x) in a.iter().enumerate() {
            th.min_shared[i] = x.as_u64().context("min_shared")? as u32;
        }
    }
    if let Some(a) = arr("min_ani") {
        for (i, x) in a.iter().enumerate() {
            th.min_ani[i] = x.as_f64().context("min_ani")?;
        }
    }
    if let Some(a) = arr("min_consistency") {
        for (i, x) in a.iter().enumerate() {
            th.min_consistency[i] = x.as_f64().context("min_consistency")?;
        }
    }
    if let Some(m) = v.get("margin").and_then(|m| m.as_f64()) {
        th.margin = m;
    }
    if let Some(n) = v.get("top_n").and_then(|n| n.as_u64()) {
        th.top_n = n as usize;
    }
    Ok(th)
}

fn main() -> Result<()> {
    let cli = Cli::parse();
    if let Some(t) = cli.threads {
        rayon::ThreadPoolBuilder::new().num_threads(t).build_global()?;
    }
    match cli.cmd {
        Cmd::Sketch { ksize, scaled, out_dir, abund, max_reads, force, list, mut inputs } => {
            if let Some(l) = list {
                for line in std::fs::read_to_string(l)?.lines().filter(|l| !l.trim().is_empty()) {
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
                        Ok(None) => (p.display().to_string(), "EMPTY", "no sequence records".into()),
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
        Cmd::Build { meta, sketches, out, scaled, buckets } => {
            let metas = read_meta(&meta)?;
            let loaded: Vec<(GenomeMeta, Option<SketchFile>)> = metas
                .into_par_iter()
                .map(|m| {
                    let p = sketches.join(format!("{}.gsk", m.id));
                    let s = if p.exists() { read_sketch(&p).ok() } else { None };
                    (m, s)
                })
                .collect();
            let first = loaded.iter().find_map(|(_, s)| s.as_ref()).context("no sketches found for any meta row")?;
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
                        let h = if target > base { downsample(&s.hashes, target) } else { s.hashes };
                        b.add(m, h);
                        used += 1;
                    }
                    Some(_) => {
                        eprintln!("MISMATCH\t{}\tk/scaled/seed differ from the first sketch", m.id);
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
        Cmd::Guess { index, min_abund, max_reads, top, thresholds, queries } => {
            let ix = load_index(&index)?;
            let th = match thresholds {
                Some(p) => read_thresholds(&p)?,
                None => Thresholds::default(),
            };
            for q in &queries {
                let hashes = query_hashes(&ix, q, min_abund, max_reads)?;
                let g = guess(&ix, &hashes, &|_| false, &th);
                let hits: Vec<_> = g
                    .hits
                    .iter()
                    .take(top)
                    .map(|&(gi, sh)| {
                        let m = &ix.genomes[gi as usize];
                        json!({"id": m.id, "taxid": m.taxid, "shared": sh, "lineage": m.lineage})
                    })
                    .collect();
                let calls: serde_json::Map<String, serde_json::Value> = RANKS
                    .iter()
                    .zip(&g.calls)
                    .zip(&g.votes)
                    .map(|((r, c), v)| (r.to_string(), json!({"taxon": c, "vote": v.taxon, "consistency": v.consistency})))
                    .collect();
                let top_json = g.top.as_ref().map(|t| {
                    json!({"id": ix.genomes[t.genome as usize].id, "shared": t.shared, "ani": t.ani,
                           "containment_query": t.containment_query, "containment_ref": t.containment_ref})
                });
                println!(
                    "{}",
                    json!({"query": q.display().to_string(), "query_hashes": g.query_hashes, "top": top_json,
                           "calls": calls, "hits": hits})
                );
            }
        }
        Cmd::Convert { index, out } => {
            let ix = load_index(&index)?;
            let mut w = BufWriter::new(File::create(&out)?);
            ix.write_v1(&mut w)?;
            w.flush()?;
            let (a, b) = (std::fs::metadata(&index)?.len(), std::fs::metadata(&out)?.len());
            eprintln!("{} -> {}: {} -> {} bytes ({:.2}x smaller)", index.display(), out.display(), a, b, a as f64 / b as f64);
        }
        Cmd::Info { index } => {
            let ix = load_index(&index)?;
            let size = std::fs::metadata(&index)?.len();
            println!("file\t{}\t{} bytes", index.display(), size);
            println!("k\t{}\nscaled\t{}\nseed\t{}", ix.ksize, ix.scaled, ix.seed);
            println!("genomes\t{}", ix.genomes.len());
            println!("unique_hashes\t{}\npostings\t{}", ix.n_unique_hashes(), ix.n_postings());
            for (r, t) in RANKS.iter().zip(&ix.taxa) {
                println!("taxa_{r}\t{}", t.len() - 1);
            }
        }
        Cmd::Bench { index, sketches, mode, phylum, sample, seed, shards, shard_index, out } => {
            let ix = load_index(&index)?;
            if !["loo", "los", "log"].contains(&mode.as_str()) {
                bail!("--mode must be loo, los or log");
            }
            let mut qs: Vec<u32> = (0..ix.genomes.len() as u32)
                .filter(|&g| phylum.as_ref().is_none_or(|p| &ix.genomes[g as usize].lineage[0] == p))
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
            let qs: Vec<u32> = qs.into_iter().enumerate().filter(|(i, _)| i % shards == shard_index).map(|(_, g)| g).collect();
            let th = Thresholds::default();
            let rows: Vec<String> = qs
                .par_iter()
                .filter_map(|&gi| {
                    let m = &ix.genomes[gi as usize];
                    let sk = read_sketch(&sketches.join(format!("{}.gsk", m.id))).ok()?;
                    let q = if sk.scaled < ix.scaled { downsample(&sk.hashes, ix.scaled) } else { sk.hashes };
                    let mine = ix.lineage_ids[gi as usize];
                    let rk = match mode.as_str() {
                        "los" => 5,
                        "log" => 4,
                        _ => usize::MAX,
                    };
                    let skip = |g: u32| g == gi || (rk != usize::MAX && mine[rk] != 0 && ix.lineage_ids[g as usize][rk] == mine[rk]);
                    let g = guess(&ix, &q, &skip, &th);
                    let (tid, tshared, tani, tcq, tcr) = match &g.top {
                        Some(t) => (ix.genomes[t.genome as usize].id.clone(), t.shared, t.ani, t.containment_query, t.containment_ref),
                        None => (String::new(), 0, 0.0, 0.0, 0.0),
                    };
                    let mut row = format!("{}\t{}\t{}\t{}\t{}\t{:.5}\t{:.5}\t{:.5}", m.id, g.query_hashes, tid, tshared, tcq, tcr, tani, ix.genomes.len());
                    let tl = g.top.as_ref().map(|t| &ix.genomes[t.genome as usize].lineage);
                    for r in 0..6 {
                        let truth = &m.lineage[r];
                        let top_taxon = tl.map(|l| l[r].as_str()).unwrap_or("");
                        let v = &g.votes[r];
                        row.push_str(&format!("\t{}\t{}\t{}\t{:.4}", truth, top_taxon, v.taxon.clone().unwrap_or_default(), v.consistency));
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
            eprintln!("bench {mode}: {} of {} queries scored -> {}", rows.len(), qs.len(), out.display());
        }
    }
    Ok(())
}

/// Hashes of a query at the index's k/scaled, filtered by `min_abund`.
fn query_hashes(ix: &RefIndex, path: &Path, min_abund: u32, max_reads: Option<u64>) -> Result<Vec<u64>> {
    let sk = if path.extension().is_some_and(|e| e == "gsk") {
        read_sketch(path)?
    } else {
        sketch_file(path, ix.ksize, ix.scaled, min_abund > 1, max_reads)?
            .with_context(|| format!("{}: no sequence records", path.display()))?
    };
    if sk.ksize != ix.ksize || sk.seed != ix.seed {
        bail!("{}: k/seed ({}/{}) do not match the index ({}/{})", path.display(), sk.ksize, sk.seed, ix.ksize, ix.seed);
    }
    if sk.scaled > ix.scaled {
        bail!("{}: query scaled {} is coarser than the index scaled {}", path.display(), sk.scaled, ix.scaled);
    }
    let mut pairs: Vec<(u64, u32)> = match &sk.abunds {
        Some(a) => sk.hashes.iter().copied().zip(a.iter().copied()).collect(),
        None => sk.hashes.iter().map(|&h| (h, 1)).collect(),
    };
    pairs.retain(|&(_, a)| a >= min_abund);
    let h: Vec<u64> = pairs.into_iter().map(|(h, _)| h).collect();
    Ok(if sk.scaled < ix.scaled { downsample(&h, ix.scaled) } else { h })
}
