//! `gmg run`: classify every genome / read set in a folder and write `results.tsv`, `hits.tsv`, `run.json`.
//!
//! Each finished input is appended to `results.jsonl` immediately, so a killed run can be resumed (`--resume`);
//! the TSV files and the manifest are generated from that file at the end. Schema version: see `SCHEMA_VERSION`.

use crate::{load_index, load_thresholds, sketch_paths, strip_ext};
use anyhow::{bail, Context, Result};
use gmg_core::classify::{
    classify, error_record, exclusion_mask, hits_tsv, relationship_label, relationships,
    results_tsv, Context as RefContext, Item, READS_SPECIES_MIN_ANI, SCHEMA_VERSION,
};
use gmg_core::{RefIndex, Thresholds};
use rayon::prelude::*;
use serde_json::{json, Value};
use std::collections::HashMap;
use std::collections::HashSet;
use std::fs::{self, File, OpenOptions};
use std::io::{BufWriter, Write};
use std::path::{Path, PathBuf};
use std::sync::atomic::{AtomicUsize, Ordering};
use std::sync::Mutex;
use std::time::{Instant, SystemTime, UNIX_EPOCH};

pub struct RunArgs {
    pub inputs: Vec<PathBuf>,
    /// file with one input path per line
    pub list: Option<PathBuf>,
    /// merge read files of one sample (R1/R2, lanes)
    pub group: bool,
    /// record the SHA-256 of every input file in run.json (slow for large inputs)
    pub checksums: bool,
    pub quiet: bool,
    pub verbose: bool,
    /// TSV naming inputs (id, path, expected)
    pub sample_sheet: Option<PathBuf>,
    pub index: PathBuf,
    pub out: PathBuf,
    pub thresholds: Option<String>,
    pub kind: String,
    pub min_abund: Option<u32>,
    pub max_reads: Option<u64>,
    pub top: usize,
    pub resume: bool,
    /// write report.html (default); `--no-report` turns it off
    pub report: bool,
    pub max_cards: usize,
    /// html (default), md or all
    pub report_format: String,
    /// reference genome ids to hide from every query (one per line)
    pub exclude_ids: Option<PathBuf>,
    /// hide each query's own genome / species / genus from the reference (sample sheet `ref` and `expected` columns)
    pub exclude_level: Option<gmg_core::classify::ExcludeLevel>,
}

struct Input {
    /// one file, or the files of one read sample (R1/R2, lanes) merged into one sketch; `path` is the first
    path: PathBuf,
    files: Vec<PathBuf>,
    id: String,
    kind: &'static str, // "assembly" | "reads"
    /// from a sample sheet: the label the user expects (carried to results.tsv for validation)
    expected: Option<String>,
    /// from a sample sheet: the reference genome id of this query, for `--exclude-level`
    ref_id: Option<String>,
    /// named by a sample sheet: never merged automatically
    explicit: bool,
}

const FASTA_EXT: [&str; 5] = ["fa", "fna", "fasta", "fas", "ffn"];
const FASTQ_EXT: [&str; 2] = ["fq", "fastq"];

/// File extension without a trailing `.gz`, lowercased.
fn base_ext(p: &Path) -> Option<String> {
    let name = p.file_name()?.to_string_lossy().to_lowercase();
    let name = name.strip_suffix(".gz").unwrap_or(&name).to_string();
    name.rsplit_once('.').map(|(_, e)| e.to_string())
}

fn kind_of(p: &Path, forced: &str) -> Option<&'static str> {
    let ext = base_ext(p)?;
    let auto = if FASTA_EXT.contains(&ext.as_str()) {
        Some("assembly")
    } else if FASTQ_EXT.contains(&ext.as_str()) {
        Some("reads")
    } else {
        None
    };
    match forced {
        "assembly" => auto.map(|_| "assembly"),
        "reads" => auto.map(|_| "reads"),
        _ => auto,
    }
}

fn walk(dir: &Path, out: &mut Vec<PathBuf>) -> Result<()> {
    let mut entries: Vec<PathBuf> = fs::read_dir(dir)
        .with_context(|| format!("reading directory {}", dir.display()))?
        .filter_map(|e| e.ok().map(|e| e.path()))
        .collect();
    entries.sort();
    for p in entries {
        if p.is_dir() {
            walk(&p, out)?;
        } else {
            out.push(p);
        }
    }
    Ok(())
}

/// Sample key of a read file: the stem without the read-direction (`_R1`, `_R2`, `_R1_001`, `.1`/`_1`/`_2`) and lane (`_L001`) parts.
fn sample_key(stem: &str) -> String {
    let mut parts: Vec<&str> = stem.split('_').collect();
    // drop lane tokens such as L001
    parts.retain(|t| {
        !(t.len() == 4 && t.starts_with('L') && t[1..].chars().all(|c| c.is_ascii_digit()))
    });
    // a trailing 3-digit chunk after R1/R2 (Illumina `_R1_001`)
    if parts.len() >= 2 {
        let last = parts[parts.len() - 1];
        let prev = parts[parts.len() - 2];
        if last.len() == 3
            && last.chars().all(|c| c.is_ascii_digit())
            && matches!(prev, "R1" | "R2")
        {
            parts.truncate(parts.len() - 2);
            return parts.join("_");
        }
    }
    if let Some(last) = parts.last().copied() {
        if matches!(last, "R1" | "R2" | "1" | "2") && parts.len() > 1 {
            parts.pop();
        } else if let Some(b) = last.strip_suffix(".1").or_else(|| last.strip_suffix(".2")) {
            let n = parts.len();
            parts[n - 1] = b;
        }
    }
    parts.join("_")
}

/// Merge read files of one sample into a single input (same directory, same key, kind reads). A lone file keeps its own id.
fn group_reads(inputs: Vec<Input>) -> Vec<Input> {
    let mut order: Vec<(PathBuf, String)> = Vec::new();
    let mut groups: HashMap<(PathBuf, String), Vec<Input>> = HashMap::new();
    let mut out: Vec<Input> = Vec::new();
    for i in inputs {
        if i.kind != "reads" || i.explicit {
            out.push(i);
            continue;
        }
        let key = (
            i.path.parent().map(|p| p.to_path_buf()).unwrap_or_default(),
            sample_key(&i.id),
        );
        if !groups.contains_key(&key) {
            order.push(key.clone());
        }
        groups.entry(key).or_default().push(i);
    }
    for key in order {
        let mut g = groups.remove(&key).unwrap();
        if g.len() == 1 {
            out.push(g.pop().unwrap());
        } else {
            let files: Vec<PathBuf> = g.iter().map(|i| i.path.clone()).collect();
            out.push(Input {
                path: files[0].clone(),
                files,
                id: key.1.clone(),
                kind: "reads",
                expected: None,
                ref_id: None,
                explicit: false,
            });
        }
    }
    out.sort_by(|a, b| a.id.cmp(&b.id));
    out
}

/// Sample sheet: TSV with a header and columns `id`, `path` (relative to the sheet) and optionally `expected`.
/// Rows with the same id are the files (R1/R2, lanes) of one sample.
fn read_sheet(sheet: &Path, forced: &str) -> Result<Vec<Input>> {
    let text =
        fs::read_to_string(sheet).with_context(|| "reading the --sample-sheet file".to_string())?;
    let base = sheet.parent().map(Path::to_path_buf).unwrap_or_default();
    let mut lines = text
        .lines()
        .filter(|l| !l.trim().is_empty() && !l.starts_with('#'));
    let header: Vec<&str> = lines
        .next()
        .unwrap_or("")
        .split('\t')
        .map(str::trim)
        .collect();
    let col = |n: &str| header.iter().position(|h| h.eq_ignore_ascii_case(n));
    let (Some(ci), Some(cp)) = (col("id"), col("path")) else {
        bail!("sample sheet needs a header with columns 'id' and 'path' (and optionally 'expected'), separated by tabs");
    };
    let ce = col("expected");
    let cr = col("ref");
    let mut out: Vec<Input> = Vec::new();
    for (n, line) in lines.enumerate() {
        let f: Vec<&str> = line.split('\t').map(str::trim).collect();
        let (Some(id), Some(path)) = (
            f.get(ci).filter(|x| !x.is_empty()),
            f.get(cp).filter(|x| !x.is_empty()),
        ) else {
            bail!("sample sheet row {}: id and path are required", n + 2);
        };
        let path = {
            let p = PathBuf::from(path);
            if p.is_absolute() {
                p
            } else {
                base.join(p)
            }
        };
        if !path.is_file() {
            bail!(
                "sample sheet row {}: file not found: {}",
                n + 2,
                path.display()
            );
        }
        let Some(kind) = kind_of(&path, forced) else {
            bail!(
                "sample sheet row {}: {} does not look like FASTA/FASTQ",
                n + 2,
                path.display()
            );
        };
        let expected = ce
            .and_then(|c| f.get(c))
            .filter(|x| !x.is_empty())
            .map(|x| x.to_string());
        let ref_id = cr
            .and_then(|c| f.get(c))
            .filter(|x| !x.is_empty())
            .map(|x| x.to_string());
        match out.iter_mut().find(|i| i.id == *id) {
            Some(i) => {
                if i.kind != kind {
                    bail!("sample sheet: '{id}' mixes assemblies and reads");
                }
                i.files.push(path);
            }
            None => out.push(Input {
                files: vec![path.clone()],
                path,
                id: id.to_string(),
                kind,
                expected,
                ref_id,
                explicit: true,
            }),
        }
    }
    Ok(out)
}

fn discover(args: &RunArgs) -> Result<Vec<Input>> {
    let mut files = Vec::new();
    let mut given: Vec<PathBuf> = args.inputs.clone();
    if let Some(l) = &args.list {
        let text = fs::read_to_string(l).with_context(|| "reading the --list file".to_string())?;
        given.extend(
            text.lines()
                .map(str::trim)
                .filter(|l| !l.is_empty() && !l.starts_with('#'))
                .map(PathBuf::from),
        );
    }
    for p in &given {
        if p.is_dir() {
            walk(p, &mut files)?;
        } else if p.is_file() {
            files.push(p.clone());
        } else {
            bail!("input not found: {}", p.display());
        }
    }
    let mut inputs: Vec<Input> = Vec::new();
    let mut seen: HashMap<String, PathBuf> = HashMap::new();
    for p in files {
        // directories are filtered by extension; explicitly named files must also look like sequence files
        let Some(kind) = kind_of(&p, &args.kind) else {
            continue;
        };
        let id = strip_ext(
            &p.file_name()
                .map(|s| s.to_string_lossy().into_owned())
                .unwrap_or_default(),
        );
        if let Some(prev) = seen.insert(id.clone(), p.clone()) {
            bail!(
                "two inputs would get the id '{id}': {} and {}; rename one",
                prev.display(),
                p.display()
            );
        }
        inputs.push(Input {
            files: vec![p.clone()],
            path: p,
            id,
            kind,
            expected: None,
            ref_id: None,
            explicit: false,
        });
    }
    if args.group {
        inputs = group_reads(inputs);
        let merged: Vec<&Input> = inputs.iter().filter(|i| i.files.len() > 1).collect();
        if !merged.is_empty() && !args.quiet {
            eprintln!(
                "merged {} read files into {} samples (R1/R2 and lanes; --no-group keeps them separate)",
                merged.iter().map(|i| i.files.len()).sum::<usize>(),
                merged.len()
            );
        }
    }
    if let Some(sheet) = &args.sample_sheet {
        for e in read_sheet(sheet, &args.kind)? {
            if let Some(prev) = inputs.iter().find(|i| i.id == e.id) {
                bail!(
                    "sample sheet id '{}' is also used by {}; ids must be unique",
                    e.id,
                    prev.path.display()
                );
            }
            inputs.push(e);
        }
    }
    if inputs.is_empty() {
        bail!("no FASTA/FASTQ files found in the given inputs (extensions: fa fna fasta fas ffn fq fastq, optionally .gz)");
    }
    Ok(inputs)
}

pub(crate) fn utc_now() -> String {
    let secs = SystemTime::now()
        .duration_since(UNIX_EPOCH)
        .map(|d| d.as_secs() as i64)
        .unwrap_or(0);
    let (days, rem) = (secs.div_euclid(86_400), secs.rem_euclid(86_400));
    // civil-from-days (Howard Hinnant)
    let z = days + 719_468;
    let era = z.div_euclid(146_097);
    let doe = z.rem_euclid(146_097);
    let yoe = (doe - doe / 1460 + doe / 36_524 - doe / 146_096) / 365;
    let (y, doy) = (yoe + era * 400, doe - (365 * yoe + yoe / 4 - yoe / 100));
    let mp = (5 * doy + 2) / 153;
    let d = doy - (153 * mp + 2) / 5 + 1;
    let m = if mp < 10 { mp + 3 } else { mp - 9 };
    let y = if m <= 2 { y + 1 } else { y };
    format!(
        "{y:04}-{m:02}-{d:02}T{:02}:{:02}:{:02}Z",
        rem / 3600,
        rem % 3600 / 60,
        rem % 60
    )
}

fn thresholds_json(th: &Thresholds) -> Value {
    json!({"min_shared": th.min_shared, "min_ani": th.min_ani, "min_consistency": th.min_consistency,
           "margin": th.margin, "top_n": th.top_n,
           "tentative_genus_min_ani": th.tentative_genus_min_ani})
}

fn input_size(i: &Input) -> u64 {
    i.files
        .iter()
        .map(|f| fs::metadata(f).map(|m| m.len()).unwrap_or(0))
        .sum()
}

fn process(
    ix: &RefIndex,
    th: &Thresholds,
    inp: &Input,
    a: &RunArgs,
    keep_hashes: bool,
    ctx: &RefContext,
    excl: &HashSet<String>,
) -> (Value, Option<Vec<u64>>) {
    let size = input_size(inp);
    let path = inp.path.display().to_string();
    let attempt = || -> Result<(Value, Option<Vec<u64>>)> {
        let reads = inp.kind == "reads";
        let min_abund = a.min_abund.unwrap_or(if reads { 2 } else { 1 });
        let (sk, st) = sketch_paths(&inp.files, ix.ksize, ix.scaled, reads, a.max_reads)?
            .with_context(|| "no sequence records".to_string())?;
        let item = Item {
            id: &inp.id,
            path: &path,
            size_bytes: size,
            kind: inp.kind,
        };
        let hidden = (!excl.is_empty() || a.exclude_level.is_some()).then(|| {
            exclusion_mask(
                ix,
                excl,
                a.exclude_level,
                inp.ref_id.as_deref().or(Some(inp.id.as_str())),
                inp.expected.as_deref(),
            )
        });
        let skip = |g: u32| hidden.as_ref().is_some_and(|m| m[g as usize]);
        let (mut rec, hashes) = classify(ix, th, ctx, &item, sk, &st, min_abund, a.top, &skip)
            .map_err(anyhow::Error::msg)?;
        if let Some(m) = &hidden {
            rec["excluded"] = json!({"level": a.exclude_level.map(|l| format!("{l:?}")),
                                     "n_genomes": m.iter().filter(|x| **x).count()});
        }
        if let Some(e) = &inp.expected {
            rec["expected"] = json!(e);
        }
        if inp.files.len() > 1 {
            rec["files"] = json!(inp
                .files
                .iter()
                .filter_map(|f| f.file_name().map(|n| n.to_string_lossy().into_owned()))
                .collect::<Vec<_>>());
        }
        Ok((
            rec,
            if keep_hashes && !reads {
                Some(hashes)
            } else {
                None
            },
        ))
    };
    // a bug or an unforeseen file must cost one row, never the whole run
    let result = std::panic::catch_unwind(std::panic::AssertUnwindSafe(attempt))
        .unwrap_or_else(|_| Err(anyhow::anyhow!("internal error while processing this file")));
    match result {
        Ok((v, h)) => (v, h),
        Err(e) => {
            let ie = crate::errors::classify(&inp.path, size, &e);
            let item = Item {
                id: &inp.id,
                path: &path,
                size_bytes: size,
                kind: inp.kind,
            };
            (error_record(&item, ie.code, &ie.message, ie.action), None)
        }
    }
}

fn write_tsv(path: &Path, header: &[String], rows: &[Vec<String>]) -> Result<()> {
    let mut w = BufWriter::new(File::create(path)?);
    writeln!(w, "{}", header.join("\t"))?;
    for r in rows {
        writeln!(w, "{}", r.join("\t"))?;
    }
    w.flush()?;
    Ok(())
}

/// `75` -> `1m15s`, `4000` -> `1h06m`.
fn fmt_secs(s: f64) -> String {
    let s = s.round() as u64;
    match s {
        0..=59 => format!("{s}s"),
        60..=3599 => format!("{}m{:02}s", s / 60, s % 60),
        _ => format!("{}h{:02}m", s / 3600, s % 3600 / 60),
    }
}

pub fn run(args: RunArgs) -> Result<()> {
    let t0 = Instant::now();
    let started = utc_now();
    let inputs = discover(&args)?;
    let ix = load_index(&args.index)?;
    let (th, th_name) = load_thresholds(args.thresholds.as_deref())?;
    fs::create_dir_all(&args.out)?;
    let jsonl = args.out.join("results.jsonl");

    // resume: keep finished, non-error results whose file (path, size) is unchanged
    let mut done: HashMap<String, Value> = HashMap::new();
    if args.resume && jsonl.exists() {
        for line in fs::read_to_string(&jsonl)?.lines() {
            if let Ok(v) = serde_json::from_str::<Value>(line) {
                if v["schema_version"] == SCHEMA_VERSION && v["status"] != "error" {
                    if let Some(id) = v["id"].as_str() {
                        done.insert(id.to_string(), v);
                    }
                }
            }
        }
    }
    let todo: Vec<&Input> = inputs
        .iter()
        .filter(|i| {
            !done.get(&i.id).is_some_and(|v| {
                v["path"] == i.path.display().to_string() && v["size_bytes"] == input_size(i)
            })
        })
        .collect();
    done.retain(|id, _| inputs.iter().any(|i| &i.id == id));
    let kept: Vec<String> = inputs
        .iter()
        .filter(|i| !todo.iter().any(|t| t.id == i.id))
        .map(|i| i.id.clone())
        .collect();
    done.retain(|id, _| kept.contains(id));
    if !args.quiet {
        eprintln!(
        "gmg run: {} inputs ({} to process, {} reused from a previous run); index {} genomes, k={}, scaled={}; thresholds: {}",
        inputs.len(), todo.len(), done.len(), ix.genomes.len(), ix.ksize, ix.scaled, th_name
    );
    }

    // rewrite results.jsonl with the reused results, then append as inputs finish
    {
        let mut w = BufWriter::new(File::create(&jsonl)?);
        for i in &inputs {
            if let Some(v) = done.get(&i.id) {
                writeln!(w, "{v}")?;
            }
        }
    }
    let writer = Mutex::new(BufWriter::new(
        OpenOptions::new().append(true).open(&jsonl)?,
    ));
    let ctx = RefContext::new(&ix);
    let excl: HashSet<String> = match &args.exclude_ids {
        Some(f) => fs::read_to_string(f)
            .with_context(|| "reading the --exclude-ids file".to_string())?
            .lines()
            .map(str::trim)
            .filter(|l| !l.is_empty() && !l.starts_with('#'))
            .map(|l| l.split('\t').next().unwrap_or(l).to_string())
            .collect(),
        None => HashSet::new(),
    };
    let counter = AtomicUsize::new(0);
    let t_work = Instant::now();
    let total = todo.len();
    let n_asm = inputs.iter().filter(|i| i.kind == "assembly").count();
    let want_rel = (2..=60).contains(&n_asm);
    let fresh: Vec<(String, Value, Option<Vec<u64>>)> = todo
        .par_iter()
        .map(|inp| {
            let (v, h) = process(&ix, &th, inp, &args, want_rel, &ctx, &excl);
            {
                let mut w = writer.lock().unwrap();
                let _ = writeln!(w, "{v}");
                let _ = w.flush();
            }
            let n = counter.fetch_add(1, Ordering::SeqCst) + 1;
            let el = t_work.elapsed().as_secs_f64();
            if !args.quiet {
                eprintln!(
                    "[{n}/{total}] {}  {}  {}  ({}, ETA {})",
                    inp.id,
                    v["status"].as_str().unwrap_or(""),
                    v["best_rank"].as_str().unwrap_or("-"),
                    fmt_secs(el),
                    fmt_secs(el / n as f64 * (total - n) as f64)
                );
                if args.verbose {
                    let w = crate::withheld_text(&v["guess"]["withheld"]);
                    if !w.is_empty() {
                        eprintln!("      withheld {w}");
                    }
                }
            }
            (inp.id.clone(), v, h)
        })
        .collect();
    let mut hash_by_id: HashMap<String, Vec<u64>> = HashMap::new();
    let mut by_id = done;
    for (id, v, h) in fresh {
        if let Some(h) = h {
            hash_by_id.insert(id.clone(), h);
        }
        by_id.insert(id, v);
    }
    let results: Vec<Value> = inputs.iter().filter_map(|i| by_id.remove(&i.id)).collect();

    let (h, rows) = results_tsv(&results);
    write_tsv(&args.out.join("results.tsv"), &h, &rows)?;
    let (h, rows) = hits_tsv(&results);
    write_tsv(&args.out.join("hits.tsv"), &h, &rows)?;

    // pairwise ANI between assemblies (only when every assembly was sketched in this run and there are 2..=60 of them)
    let rel_path = args.out.join("relationships.json");
    let asm: Vec<&Value> = results
        .iter()
        .filter(|r| r["kind"] == "assembly" && r["status"] != "error")
        .collect();
    if want_rel
        && asm.len() >= 2
        && asm
            .iter()
            .all(|r| hash_by_id.contains_key(r["id"].as_str().unwrap_or("")))
    {
        let items: Vec<(String, String, Vec<u64>)> = asm
            .iter()
            .map(|r| {
                let id = r["id"].as_str().unwrap_or("").to_string();
                let label = relationship_label(r);
                let h = hash_by_id[&id].clone();
                (id, label, h)
            })
            .collect();
        fs::write(
            &rel_path,
            serde_json::to_string(&relationships(ix.ksize, &items))? + "\n",
        )?;
    } else if !args.resume || !want_rel {
        let _ = fs::remove_file(&rel_path);
    }

    let mut counts: HashMap<String, usize> = HashMap::new();
    for r in &results {
        *counts
            .entry(r["status"].as_str().unwrap_or("").to_string())
            .or_default() += 1;
    }
    let manifest = json!({
        "schema_version": SCHEMA_VERSION, "tool": "gmg", "version": env!("CARGO_PKG_VERSION"), "commit": env!("GMG_GIT"),
        "command": std::env::args().collect::<Vec<_>>(), "started_utc": started, "finished_utc": utc_now(),
        "elapsed_seconds": (t0.elapsed().as_secs_f64() * 10.0).round() / 10.0,
        "index": {"path": args.index.display().to_string(), "bytes": fs::metadata(&args.index).map(|m| m.len()).unwrap_or(0),
                  "k": ix.ksize, "scaled": ix.scaled, "seed": ix.seed, "n_genomes": ix.genomes.len(),
                  "manifest": ix.manifest.as_deref().and_then(|m| serde_json::from_str::<Value>(m).ok())},
        "exclusion": {"ids_file": args.exclude_ids.as_ref().map(|p| p.display().to_string()), "n_ids": excl.len(),
                      "level": args.exclude_level.map(|l| format!("{l:?}"))},
        "thresholds": {"name": th_name, "values": thresholds_json(&th), "reads_species_min_ani": READS_SPECIES_MIN_ANI},
        "status_counts": counts,
        "inputs": results.iter().map(|r| {
            let mut v = json!({"id": r["id"], "path": r["path"], "size_bytes": r["size_bytes"],
                                                "kind": r["kind"], "status": r["status"]});
            if args.checksums {
                if let Some(i) = inputs.iter().find(|i| r["id"] == i.id.as_str()) {
                    let sums: Vec<String> = i.files.iter().filter_map(|f| crate::sha256_file(f).ok()).collect();
                    v["sha256"] = if sums.len() == 1 { json!(sums[0]) } else { json!(sums) };
                }
            }
            v
        }).collect::<Vec<_>>(),
    });
    fs::write(
        args.out.join("run.json"),
        serde_json::to_string_pretty(&manifest)? + "\n",
    )?;
    if args.report {
        write_report(&args.out, args.max_cards, &args.report_format)?;
    }
    let mut cs: Vec<_> = counts.iter().collect();
    cs.sort();
    if !args.quiet {
        eprintln!(
            "done in {:.1}s: {}  ->  {}/{{results.tsv,hits.tsv,run.json,results.jsonl{}}}",
            t0.elapsed().as_secs_f64(),
            cs.iter()
                .map(|(k, v)| format!("{k}={v}"))
                .collect::<Vec<_>>()
                .join(" "),
            args.out.display(),
            if args.report { ",report.html" } else { "" }
        );
    }
    Ok(())
}

/// Render report.html from the files of a finished run in `dir` (results.jsonl, run.json, optional relationships.json).
pub fn write_report(dir: &Path, max_cards: usize, format: &str) -> Result<()> {
    let results = gmg_report::parse_results_jsonl(
        &fs::read_to_string(dir.join("results.jsonl")).with_context(|| {
            format!(
                "reading {}/results.jsonl (is this a gmg run directory?)",
                dir.display()
            )
        })?,
    );
    let run: Value = serde_json::from_str(
        &fs::read_to_string(dir.join("run.json"))
            .with_context(|| format!("reading {}/run.json", dir.display()))?,
    )?;
    let rel: Option<Value> = fs::read_to_string(dir.join("relationships.json"))
        .ok()
        .and_then(|t| serde_json::from_str(&t).ok());
    let opts = gmg_report::Options {
        max_cards,
        ..Default::default()
    };
    if matches!(format, "html" | "all") {
        fs::write(
            dir.join("report.html"),
            gmg_report::render_html(&results, &run, rel.as_ref(), &opts),
        )?;
    }
    if matches!(format, "md" | "pdf" | "all") {
        let m = gmg_report::render_markdown(&results, &run, rel.as_ref(), &opts);
        if matches!(format, "md" | "all") {
            fs::create_dir_all(dir.join("figures"))?;
            for (name, svg) in &m.figures {
                fs::write(dir.join("figures").join(name), svg)?;
            }
            fs::write(dir.join("report.md"), &m.md)?;
        }
        if matches!(format, "pdf" | "all") {
            let pdf =
                gmg_pdf::render_pdf(&m, &opts.title).map_err(|e| anyhow::anyhow!("PDF: {e}"))?;
            fs::write(dir.join("report.pdf"), pdf)?;
        }
    }
    Ok(())
}
