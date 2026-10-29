use anyhow::{bail, Context, Result};
use clap::{CommandFactory, Parser, Subcommand};
pub(crate) use gmg_core::classify::{status_of, withheld_json, withheld_text};
use gmg_core::{
    downsample, guess, FracSketch, GenomeMeta, IndexBuilder, RefIndex, SketchFile, Thresholds,
    DEFAULT_SEED, RANKS,
};
use needletail::parse_fastx_file;
mod errors;
mod run;

use rayon::prelude::*;
use std::fs::File;
use std::io::{BufReader, BufWriter, Write};
use std::path::{Path, PathBuf};

#[derive(Parser)]
#[command(
    name = "gmg",
    version = concat!(env!("CARGO_PKG_VERSION"), " (", env!("GMG_GIT"), ")"),
    about = "Guess My Genome: guess the fungal taxonomy of genomes or reads from sketches",
    after_long_help = "\
QUICK START
  gmg run my_genomes/ --index bfd-s10000.gmgi -o results/      classify a folder, write results.tsv + report.html
  gmg run reads/ --index bfd-s10000.gmgi -o results/ --report-format all   FASTQ too; also report.md and report.pdf
  gmg guess --index bfd-s10000.gmgi genome.fa.gz               one genome, JSON on stdout

EXIT CODES
  0  the command finished (inputs that could not be read are rows with status=error, not failures)
  1  the command could not run (index missing or unreadable, no input found, output not writable)
  2  usage error (unknown option, missing argument)

Per-command examples: gmg <command> --help.  Tutorial: docs/tutorial.md.  Shell completions: gmg completions bash|zsh|fish."
)]
struct Cli {
    /// worker threads (default: all cores)
    #[arg(long, global = true)]
    threads: Option<usize>,
    /// print no progress or summary lines (errors still go to stderr)
    #[arg(short, long, global = true, conflicts_with = "verbose")]
    quiet: bool,
    /// also print, per input, why the first uncalled rank was withheld
    #[arg(short, long, global = true)]
    verbose: bool,
    #[command(subcommand)]
    cmd: Cmd,
}

#[derive(Subcommand)]
enum Cmd {
    /// Sketch FASTA/FASTQ files (optionally gzipped) into .gsk files
    #[command(
        after_long_help = "EXAMPLES\n  gmg sketch -s 1000 -o sketches/ genomes/*.fa.gz        one .gsk per genome\n  gmg sketch --abund --max-reads 2000000 -s 10000 -o sk/ reads_R1.fastq.gz"
    )]
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
    #[command(
        after_long_help = "EXAMPLES\n  gmg build --meta meta.tsv --sketches sketches/ -o bfd.gmgi\n  gmg build --meta meta.tsv --sketches sketches/ --scaled 10000 -o bfd-s10000.gmgi --taxonomy-version \"NCBI 2026-10\"\n\nmeta.tsv columns: id taxid phylum class order family genus species [cluster]"
    )]
    Build {
        /// TSV with a header: id taxid phylum class order family genus species, optional cluster, <rank>_taxid x6, label (scripts/normalize_taxonomy.py writes all of them)
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
        /// label recorded in the index manifest, e.g. "NCBI taxonomy 2026-10"
        #[arg(long, default_value = "")]
        taxonomy_version: String,
        /// free-text note recorded in the index manifest
        #[arg(long, default_value = "")]
        note: String,
        /// comma-separated phyla of non-fungal decoy genomes in the reference (e.g. "Streptophyta,Chordata"); inputs that
        /// match them get the `not_fungal` QC flag
        #[arg(long, default_value = "")]
        non_fungal_phyla: String,
    },
    /// Guess the taxonomy of one or more queries (.gsk, FASTA or FASTQ)
    #[command(
        after_long_help = "EXAMPLES\n  gmg guess --index bfd.gmgi genome.fa.gz\n  gmg guess --index bfd.gmgi --thresholds novel query.fa        stricter set for a possibly new genus\n  gmg guess --index bfd.gmgi --format table a.fa b.fa            one readable line per query (or --format tsv)"
    )]
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
        /// `known` (default), `permissive`, `novel`, or a JSON file (min_shared[6], min_ani[6], min_consistency[6], margin, top_n, tentative_genus_min_ani)
        #[arg(long)]
        thresholds: Option<String>,
        /// json (default: one JSON object per query), tsv (header + one row per query) or table (one readable line per query)
        #[arg(long, default_value = "json", value_parser = ["json", "tsv", "table"])]
        format: String,
        queries: Vec<PathBuf>,
    },
    /// Classify every genome / read set in a folder: writes results.tsv, hits.tsv and run.json
    #[command(
        after_long_help = "EXAMPLES\n  gmg run genomes/ --index bfd.gmgi -o out/\n  gmg run reads/ --index bfd.gmgi -o out/ --max-reads 500000       R1/R2 and lanes of a sample are merged\n  gmg run --sample-sheet samples.tsv --index bfd.gmgi -o out/ --report-format all\n  gmg run genomes/ --index bfd.gmgi -o out/ --resume               continue after an interruption\n\nWrites results.tsv (one row per input), hits.tsv, results.jsonl, run.json and report.html (docs/results-schema.md)."
    )]
    Run {
        /// directories (searched recursively) and/or files: FASTA (.fa .fna .fasta .fas .ffn) or FASTQ (.fq .fastq), optionally .gz
        #[arg(required_unless_present_any = ["list", "sample_sheet"])]
        inputs: Vec<PathBuf>,
        #[arg(long)]
        index: PathBuf,
        /// output directory
        #[arg(short, long)]
        out: PathBuf,
        /// `known` (default), `permissive` (older, lower genus cutoff), `novel` (queries that may belong to a genus absent from the reference), or a JSON file
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
        /// text file with one input path per line (in addition to the positional inputs)
        #[arg(long)]
        list: Option<PathBuf>,
        /// record the SHA-256 of every input file in run.json (slow for large inputs)
        #[arg(long)]
        checksums: bool,
        /// TSV with columns id, path (relative to the sheet) and optional expected; rows with the same id are one sample (R1/R2, lanes)
        #[arg(long)]
        sample_sheet: Option<PathBuf>,
        /// do not merge read files of one sample (R1/R2, lanes) into a single input
        #[arg(long)]
        no_group: bool,
        /// reference hits per input in hits.tsv
        #[arg(long, default_value_t = 5)]
        top: usize,
        /// reuse finished results from a previous run in the same output directory
        #[arg(long)]
        resume: bool,
        /// do not write report.html
        #[arg(long)]
        no_report: bool,
        /// detail cards in the report (inputs needing attention first)
        #[arg(long, default_value_t = 400)]
        report_cards: usize,
        /// report format: html (default), md (report.md + figures/*.svg), pdf (report.pdf) or all
        #[arg(long, default_value = "html", value_parser = ["html", "md", "pdf", "all"])]
        report_format: String,
        /// hide these reference genome ids (one per line, first column) from every query
        #[arg(long, value_name = "FILE")]
        exclude_ids: Option<PathBuf>,
        /// validation of queries that are in the reference: hide each query's own genome (sample sheet column `ref`,
        /// default the sample id), or all genomes of its species or genus (from `ref`, else the `expected` column)
        #[arg(long, value_name = "LEVEL", value_parser = ["self", "species", "genus"])]
        exclude_level: Option<String>,
    },
    /// Write report.html from a finished `gmg run` directory
    #[command(
        after_long_help = "EXAMPLES\n  gmg report out/                   rewrite out/report.html\n  gmg report out/ --format all      html, Markdown (+ figures/) and PDF"
    )]
    Report {
        /// a directory written by `gmg run`
        dir: PathBuf,
        #[arg(long, default_value_t = 400)]
        max_cards: usize,
        /// html (default), md (report.md + figures/*.svg), pdf (report.pdf) or all
        #[arg(long, default_value = "html", value_parser = ["html", "md", "pdf", "all"])]
        format: String,
    },
    /// Print a shell completion script (bash, zsh, fish, elvish, powershell)
    Completions {
        #[arg(value_enum)]
        shell: clap_complete::Shell,
    },
    /// Print the gmg(1) manual page (roff) to stdout
    Man,
    /// Summarise an index
    #[command(
        after_long_help = "EXAMPLES\n  gmg info bfd.gmgi        sketch parameters, genome count, taxa per rank and the build manifest"
    )]
    Info { index: PathBuf },
    /// Greedy ANI clustering of the genomes in an index (centroids picked in priority order)
    #[command(
        after_long_help = "EXAMPLES\n  gmg cluster --index bfd.gmgi --sketches sketches/ --ani 0.95 -o clusters.tsv"
    )]
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
    #[command(
        after_long_help = "EXAMPLES\n  gmg convert bfd.gmgi -o bfd-packed.gmgi      the packed format: about 3x smaller than raw, loads fast, memory stays near the file size"
    )]
    Convert {
        index: PathBuf,
        #[arg(short, long)]
        out: PathBuf,
        /// v3 (default): packed, loads fast and stays small in memory; v1: delta-coded, expanded at load
        #[arg(long, default_value = "v3", value_parser = ["v1", "v3"])]
        format: String,
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

pub(crate) type Stats = gmg_core::SeqStats;

/// Sketch one file and collect basic stats. Returns None if it has no sequence records.
pub(crate) fn sketch_with_stats(
    path: &Path,
    k: u32,
    scaled: u32,
    abund: bool,
    max_reads: Option<u64>,
) -> Result<Option<(SketchFile, Stats)>> {
    sketch_paths(&[path.to_path_buf()], k, scaled, abund, max_reads)
}

/// Sketch several files into ONE sketch (paired reads, lanes). `max_reads` applies to each file.
pub(crate) fn sketch_paths(
    paths: &[PathBuf],
    k: u32,
    scaled: u32,
    abund: bool,
    max_reads: Option<u64>,
) -> Result<Option<(SketchFile, Stats)>> {
    let path = &paths[0];
    let mut sk = FracSketch::new(k as usize, scaled, DEFAULT_SEED);
    let (mut n, mut total, mut non_acgt) = (0u64, 0u64, 0u64);
    let mut lens: Vec<u64> = Vec::new();
    for path in paths {
        let mut r =
            parse_fastx_file(path).with_context(|| format!("opening {}", path.display()))?;
        let mut file_n = 0u64;
        while let Some(rec) = r.next() {
            let rec = rec.with_context(|| format!("parsing {}", path.display()))?;
            let seq = rec.seq();
            sk.add_sequence(&seq);
            n += 1;
            file_n += 1;
            total += seq.len() as u64;
            non_acgt += seq
                .iter()
                .filter(|b| !matches!(b.to_ascii_uppercase(), b'A' | b'C' | b'G' | b'T'))
                .count() as u64;
            if !abund {
                lens.push(seq.len() as u64);
            }
            if max_reads.is_some_and(|m| file_n >= m) {
                break;
            }
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

pub(crate) fn sha256_file(path: &Path) -> Result<String> {
    use sha2::{Digest, Sha256};
    let mut f = File::open(path).with_context(|| format!("opening {}", path.display()))?;
    let mut h = Sha256::new();
    std::io::copy(&mut f, &mut h)?;
    Ok(h.finalize().iter().map(|b| format!("{b:02x}")).collect())
}

pub(crate) fn load_index(path: &Path) -> Result<RefIndex> {
    let bytes = std::fs::read(path).with_context(|| format!("opening {}", path.display()))?;
    RefIndex::read_owned(bytes).with_context(|| format!("reading index {}", path.display()))
}

/// Metadata table with a header. Required columns: id taxid phylum class order family genus species. Optional: cluster,
/// phylum_taxid .. species_taxid (NCBI taxid of each rank) and label (the original label). Columns are found by name.
fn read_meta(path: &Path) -> Result<Vec<GenomeMeta>> {
    let text = std::fs::read_to_string(path)?;
    let mut lines = text
        .lines()
        .enumerate()
        .filter(|(_, l)| !l.trim().is_empty());
    let Some((_, head)) = lines.next() else {
        bail!("{}: empty metadata table", path.display());
    };
    let names: Vec<&str> = head.split('\t').map(str::trim).collect();
    let col = |n: &str| names.iter().position(|h| h.eq_ignore_ascii_case(n));
    let req = |n: &str| {
        col(n).with_context(|| format!("{}: header lacks the column '{n}' (needs: id taxid phylum class order family genus species)", path.display()))
    };
    let (c_id, c_taxid) = (req("id")?, req("taxid")?);
    let mut c_rank = [0usize; 6];
    for (r, name) in RANKS.iter().enumerate() {
        c_rank[r] = req(name)?;
    }
    let c_cluster = col("cluster");
    let c_rank_taxid: Vec<Option<usize>> =
        RANKS.iter().map(|r| col(&format!("{r}_taxid"))).collect();
    let c_label = col("label");
    let mut out = Vec::new();
    for (i, line) in lines {
        let f: Vec<&str> = line.split('\t').collect();
        if f.len() < names.len().min(8) {
            bail!(
                "{}:{}: expected {} tab-separated columns, got {}",
                path.display(),
                i + 1,
                names.len(),
                f.len()
            );
        }
        let get = |c: usize| f.get(c).copied().unwrap_or("").trim();
        let mut lineage: [String; 6] = Default::default();
        let mut rank_taxids = [0u32; 6];
        for r in 0..6 {
            lineage[r] = get(c_rank[r]).to_string();
            rank_taxids[r] = c_rank_taxid[r]
                .map(|c| get(c).parse().unwrap_or(0))
                .unwrap_or(0);
        }
        out.push(GenomeMeta {
            id: get(c_id).to_string(),
            taxid: get(c_taxid).parse().unwrap_or(0),
            lineage,
            n_hashes: 0,
            cluster: c_cluster.map(|c| get(c).to_string()).unwrap_or_default(),
            rank_taxids,
            label: c_label.map(|c| get(c).to_string()).unwrap_or_default(),
        });
    }
    Ok(out)
}

/// `known` (default; the balanced set), `permissive` (the earlier, lower genus cutoff), `novel` (stricter: queries that may belong
/// to a genus absent from the reference) or a JSON file.
pub(crate) fn load_thresholds(spec: Option<&str>) -> Result<(Thresholds, String)> {
    match spec {
        None | Some("known") => Ok((Thresholds::default(), "known".to_string())),
        Some("permissive") => Ok((
            Thresholds::from_json_str(include_str!(
                "../../../data/thresholds/bfd_s10000.permissive.json"
            ))
            .map_err(|e| anyhow::anyhow!(e))?,
            "permissive".to_string(),
        )),
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
    let (quiet, verbose) = (cli.quiet, cli.verbose);
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
            taxonomy_version,
            note,
            non_fungal_phyla,
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
            let mut ix = b.build(buckets);
            ix.set_manifest(
                serde_json::json!({
                    "gmg": env!("CARGO_PKG_VERSION"), "commit": env!("GMG_GIT"), "built_utc": run::utc_now(),
                    "format": "GMGIDX01", "k": k, "scaled": target, "seed": DEFAULT_SEED,
                    "genomes": used, "skipped": missing,
                    "meta_sha256": sha256_file(&meta)?, "meta_name": meta.file_name().map(|n| n.to_string_lossy().into_owned()),
                    "taxonomy_version": taxonomy_version, "note": note,
                    "non_fungal_phyla": non_fungal_phyla.split(',').map(str::trim).filter(|p| !p.is_empty()).collect::<Vec<_>>(),
                })
                .to_string(),
            );
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
            format,
            queries,
        } => {
            let ix = load_index(&index)?;
            let (th, _) = load_thresholds(thresholds.as_deref())?;
            if format == "tsv" {
                let mut h = vec!["query".to_string(), "status".into(), "best_rank".into()];
                h.extend(RANKS.iter().map(|r| format!("call_{r}")));
                h.extend(
                    [
                        "tentative_genus",
                        "withheld",
                        "top_hit",
                        "top_hit_species",
                        "top_hit_taxid",
                        "ani",
                    ]
                    .map(String::from),
                );
                println!("{}", h.join("\t"));
            }
            for q in &queries {
                let hashes = query_hashes(&ix, q, min_abund, max_reads)?;
                let g = guess(&ix, &hashes, &|_| false, &th);
                let name = q.display().to_string();
                let (best, status) = status_of(&g);
                let w = withheld_json(&g, &th);
                let tophit = g.top.as_ref().map(|t| {
                    let m = &ix.genomes[t.genome as usize];
                    (m.id.clone(), m.lineage[5].clone(), m.taxid, t.ani)
                });
                match format.as_str() {
                    "tsv" => {
                        let mut row = vec![
                            name,
                            status.to_string(),
                            best.map(|b| RANKS[b].to_string()).unwrap_or_default(),
                        ];
                        row.extend(g.calls.iter().map(|c| c.clone().unwrap_or_default()));
                        row.push(
                            g.tentative_genus
                                .as_ref()
                                .map(|t| t.taxon.clone())
                                .unwrap_or_default(),
                        );
                        row.push(withheld_text(&w));
                        match &tophit {
                            Some((id, sp, tid, ani)) => row.extend([
                                id.clone(),
                                sp.clone(),
                                tid.to_string(),
                                format!("{ani:.4}"),
                            ]),
                            None => row.extend(std::iter::repeat_n(String::new(), 4)),
                        }
                        println!("{}", row.join("\t"));
                    }
                    "table" => {
                        let call = best.and_then(|b| g.calls[b].clone());
                        let mut line = match (&call, best) {
                            (Some(c), Some(b)) => format!("{name}  {c}  ({}", RANKS[b]),
                            _ => format!("{name}  no call  (no_match"),
                        };
                        if let Some((_, _, _, ani)) = &tophit {
                            line.push_str(&format!(", ANI {ani:.4}"));
                        }
                        if !w.is_null() && best != Some(5) {
                            line.push_str(&format!(
                                "; {} withheld: {}",
                                w["rank"].as_str().unwrap_or(""),
                                w["detail"].as_str().unwrap_or("")
                            ));
                        }
                        if let Some(t) = &g.tentative_genus {
                            line.push_str(&format!("; tentative genus {}", t.taxon));
                        }
                        line.push(')');
                        println!("{line}");
                    }
                    _ => {
                        let mut j = gmg_core::guess_json(&ix, &g, &name, top);
                        j["status"] = status.into();
                        j["best_rank"] = best.map(|b| RANKS[b]).into();
                        j["withheld"] = w;
                        println!("{j}");
                    }
                }
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
            let ix = load_index(&index)?.into_memory()?;
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
        Cmd::Convert { index, out, format } => {
            let mut ix = load_index(&index)?;
            if let Some(m) = &ix.manifest {
                if let Ok(mut v) = serde_json::from_str::<serde_json::Value>(m) {
                    v["format"] = if format == "v1" {
                        "GMGIDX02"
                    } else {
                        "GMGIDX03"
                    }
                    .into();
                    v["converted_utc"] = run::utc_now().into();
                    ix.set_manifest(v.to_string());
                }
            }
            let mut w = BufWriter::new(File::create(&out)?);
            if format == "v1" {
                ix.write_v1(&mut w)?;
            } else {
                ix.write_v3(&mut w)?;
            }
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
            no_report,
            report_cards,
            report_format,
            exclude_ids,
            exclude_level,
            list,
            checksums,
            sample_sheet,
            no_group,
        } => {
            run::run(run::RunArgs {
                quiet,
                verbose,
                sample_sheet,
                checksums,
                report_format,
                exclude_ids,
                exclude_level: exclude_level.map(|l| l.parse().expect("validated by clap")),
                list,
                group: !no_group,
                inputs,
                index,
                out,
                thresholds,
                kind,
                min_abund,
                max_reads,
                top,
                resume,
                report: !no_report,
                max_cards: report_cards,
            })?;
        }
        Cmd::Report {
            dir,
            max_cards,
            format,
        } => {
            run::write_report(&dir, max_cards, &format)?;
            eprintln!("wrote the {format} report in {}", dir.display());
        }
        Cmd::Completions { shell } => {
            clap_complete::generate(shell, &mut Cli::command(), "gmg", &mut std::io::stdout());
        }
        Cmd::Man => {
            clap_mangen::Man::new(Cli::command()).render(&mut std::io::stdout())?;
        }
        Cmd::Info { index } => {
            let ix = load_index(&index)?;
            let size = std::fs::metadata(&index)?.len();
            println!("file\t{}\t{} bytes", index.display(), size);
            println!(
                "layout\t{}",
                if ix.is_packed() {
                    "packed (format v3)"
                } else {
                    "expanded at load (format v0/v1)"
                }
            );
            println!("k\t{}\nscaled\t{}\nseed\t{}", ix.ksize, ix.scaled, ix.seed);
            println!("genomes\t{}", ix.genomes.len());
            match &ix.manifest {
                Some(m) => match serde_json::from_str::<serde_json::Value>(m) {
                    Ok(serde_json::Value::Object(o)) => {
                        for (k, v) in o {
                            println!(
                                "manifest.{k}\t{}",
                                v.as_str()
                                    .map(str::to_string)
                                    .unwrap_or_else(|| v.to_string())
                            );
                        }
                    }
                    _ => println!("manifest\t{m}"),
                },
                None => println!("manifest\t(none: built before manifests existed)"),
            }
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
            let ix = load_index(&index)?.into_memory()?;
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
    gmg_core::classify::query_hashes_from(ix, sk, min_abund).map_err(anyhow::Error::msg)
}

#[cfg(test)]
mod tests {
    use super::*;

    /// Split on spaces, keeping "double quoted" words together.
    fn shell_words(s: &str) -> Vec<String> {
        let (mut out, mut cur, mut q) = (Vec::new(), String::new(), false);
        for c in s.chars() {
            match c {
                '"' => q = !q,
                ' ' if !q => {
                    if !cur.is_empty() {
                        out.push(std::mem::take(&mut cur));
                    }
                }
                _ => cur.push(c),
            }
        }
        if !cur.is_empty() {
            out.push(cur);
        }
        out
    }

    /// Every `gmg ...` example in the help texts must parse with the real option definitions, so the documentation cannot drift.
    #[test]
    fn help_examples_parse() {
        let cmd = Cli::command();
        let mut texts = vec![cmd
            .get_after_long_help()
            .map(|s| s.to_string())
            .unwrap_or_default()];
        for sub in cmd.get_subcommands() {
            texts.push(
                sub.get_after_long_help()
                    .map(|s| s.to_string())
                    .unwrap_or_default(),
            );
        }
        let mut n = 0;
        for t in texts {
            for line in t.lines() {
                let l = line.trim();
                if !l.starts_with("gmg ")
                    || l.starts_with("gmg <")
                    || l.contains("completions bash")
                {
                    continue;
                }
                // drop the trailing description (two or more spaces) and expand the one shell glob
                let cmdline = l
                    .split("  ")
                    .next()
                    .unwrap()
                    .replace("genomes/*.fa.gz", "genomes/a.fa.gz");
                let args = shell_words(&cmdline);
                if let Err(e) = Cli::try_parse_from(&args) {
                    panic!("help example does not parse: {cmdline}\\n{e}");
                }
                n += 1;
            }
        }
        assert!(n >= 15, "found only {n} examples");
    }
}
