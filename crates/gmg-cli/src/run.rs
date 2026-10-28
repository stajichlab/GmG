//! `gmg run`: classify every genome / read set in a folder and write `results.tsv`, `hits.tsv`, `run.json`.
//!
//! Each finished input is appended to `results.jsonl` immediately, so a killed run can be resumed (`--resume`);
//! the TSV files and the manifest are generated from that file at the end. Schema version: see `SCHEMA_VERSION`.

use crate::{load_index, load_thresholds, query_hashes_from, sketch_with_stats, strip_ext, Stats};
use anyhow::{bail, Context, Result};
use gmg_core::{guess, guess_json, RefIndex, Thresholds, RANKS};
use rayon::prelude::*;
use serde_json::{json, Value};
use std::collections::HashMap;
use std::fs::{self, File, OpenOptions};
use std::io::{BufWriter, Write};
use std::path::{Path, PathBuf};
use std::sync::atomic::{AtomicUsize, Ordering};
use std::sync::Mutex;
use std::time::{Instant, SystemTime, UNIX_EPOCH};

/// Bump when a column of results.tsv / hits.tsv or a key of results.jsonl / run.json changes meaning.
pub const SCHEMA_VERSION: u32 = 1;

pub struct RunArgs {
    pub inputs: Vec<PathBuf>,
    pub index: PathBuf,
    pub out: PathBuf,
    pub thresholds: Option<String>,
    pub kind: String,
    pub min_abund: Option<u32>,
    pub max_reads: Option<u64>,
    pub top: usize,
    pub resume: bool,
}

struct Input {
    path: PathBuf,
    id: String,
    kind: &'static str, // "assembly" | "reads"
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

fn discover(args: &RunArgs) -> Result<Vec<Input>> {
    let mut files = Vec::new();
    for p in &args.inputs {
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
        inputs.push(Input { path: p, id, kind });
    }
    if inputs.is_empty() {
        bail!("no FASTA/FASTQ files found in the given inputs (extensions: fa fna fasta fas ffn fq fastq, optionally .gz)");
    }
    Ok(inputs)
}

fn utc_now() -> String {
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
           "margin": th.margin, "top_n": th.top_n})
}

/// Provisional QC flags (issue #29 will refine them); they explain results, they never change a call.
fn qc_flags(kind: &str, st: &Stats, n_hashes: usize, top: Option<&Value>) -> Vec<&'static str> {
    let mut f = Vec::new();
    if kind == "assembly" {
        if st.total_len < 5_000_000 {
            f.push("short_genome");
        }
        if st.total_len > 150_000_000 {
            f.push("long_genome");
        }
        if st.n50 < 10_000 {
            f.push("fragmented");
        }
        if st.total_len > 0 && st.non_acgt as f64 / st.total_len as f64 > 0.05 {
            f.push("high_N");
        }
    }
    if n_hashes < 100 {
        f.push("few_hashes");
    }
    if let Some(t) = top {
        let ani = t["ani"].as_f64().unwrap_or(0.0);
        let cr = t["containment_ref"].as_f64().unwrap_or(0.0);
        if kind == "assembly" && ani >= 0.9999 && cr >= 0.999 {
            f.push("identical_to_reference");
        }
    }
    f
}

fn process(ix: &RefIndex, th: &Thresholds, inp: &Input, a: &RunArgs) -> Value {
    let size = fs::metadata(&inp.path).map(|m| m.len()).unwrap_or(0);
    let path = inp.path.display().to_string();
    let attempt = || -> Result<Value> {
        let reads = inp.kind == "reads";
        let min_abund = a.min_abund.unwrap_or(if reads { 2 } else { 1 });
        let (sk, st) = sketch_with_stats(&inp.path, ix.ksize, ix.scaled, reads, a.max_reads)?
            .with_context(|| "no sequence records".to_string())?;
        let hashes = query_hashes_from(ix, sk, min_abund)?;
        let g = guess(ix, &hashes, &|_| false, th);
        let gj = guess_json(ix, &g, &path, a.top);
        let best = g.calls.iter().rposition(|c| c.is_some());
        let status = match best {
            None => "no_match",
            Some(i) if i >= 4 => "ok",
            Some(_) => "low_confidence",
        };
        let flags = qc_flags(
            inp.kind,
            &st,
            hashes.len(),
            gj.get("top").filter(|t| !t.is_null()),
        );
        Ok(json!({
            "schema_version": SCHEMA_VERSION, "id": inp.id, "path": path, "size_bytes": size, "kind": inp.kind,
            "status": status, "best_rank": best.map(|i| RANKS[i]), "message": null,
            "qc": {"n_records": st.n_records, "total_bases": st.total_len, "n50": st.n50,
                   "non_acgt_fraction": if st.total_len > 0 { st.non_acgt as f64 / st.total_len as f64 } else { 0.0 },
                   "n_hashes": hashes.len(), "min_abund": min_abund},
            "qc_flags": flags, "guess": gj,
        }))
    };
    attempt().unwrap_or_else(|e| {
        json!({"schema_version": SCHEMA_VERSION, "id": inp.id, "path": path, "size_bytes": size, "kind": inp.kind,
               "status": "error", "best_rank": null, "message": format!("{e:#}"), "qc": null, "qc_flags": [], "guess": null})
    })
}

fn tsv_cell(v: &Value) -> String {
    let s = match v {
        Value::Null => String::new(),
        Value::String(s) => s.clone(),
        other => other.to_string(),
    };
    s.replace(['\t', '\n', '\r'], " ")
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

fn results_tsv(results: &[Value]) -> (Vec<String>, Vec<Vec<String>>) {
    let mut header: Vec<String> = ["query_id", "path", "kind", "status", "best_rank"]
        .iter()
        .map(|s| s.to_string())
        .collect();
    for r in RANKS {
        header.push(format!("call_{r}"));
    }
    for r in RANKS {
        header.push(format!("conf_{r}"));
    }
    for c in [
        "top_hit",
        "top_hit_species",
        "top_hit_taxid",
        "ani",
        "containment_query",
        "containment_ref",
        "shared",
        "cluster",
        "n_hashes",
        "total_bases",
        "n50",
        "qc_flags",
        "message",
        "schema_version",
    ] {
        header.push(c.to_string());
    }
    let mut rows = Vec::new();
    for r in results {
        let g = &r["guess"];
        let top = &g["top"];
        let top_hit = &g["hits"][0];
        let mut row = vec![
            tsv_cell(&r["id"]),
            tsv_cell(&r["path"]),
            tsv_cell(&r["kind"]),
            tsv_cell(&r["status"]),
            tsv_cell(&r["best_rank"]),
        ];
        for k in RANKS {
            row.push(tsv_cell(&g["calls"][k]["taxon"]));
        }
        for k in RANKS {
            row.push(if g.is_null() {
                String::new()
            } else {
                format!(
                    "{:.3}",
                    g["calls"][k]["consistency"].as_f64().unwrap_or(0.0)
                )
            });
        }
        let fmt = |v: &Value| v.as_f64().map(|x| format!("{x:.4}")).unwrap_or_default();
        row.extend([
            tsv_cell(&top["id"]),
            tsv_cell(&top_hit["lineage"][5]),
            tsv_cell(&top_hit["taxid"]),
            fmt(&top["ani"]),
            fmt(&top["containment_query"]),
            fmt(&top["containment_ref"]),
            tsv_cell(&top["shared"]),
            tsv_cell(&top["cluster"]),
            tsv_cell(&r["qc"]["n_hashes"]),
            tsv_cell(&r["qc"]["total_bases"]),
            tsv_cell(&r["qc"]["n50"]),
            r["qc_flags"]
                .as_array()
                .map(|a| {
                    a.iter()
                        .filter_map(|x| x.as_str())
                        .collect::<Vec<_>>()
                        .join(",")
                })
                .unwrap_or_default(),
            tsv_cell(&r["message"]),
            SCHEMA_VERSION.to_string(),
        ]);
        rows.push(row);
    }
    (header, rows)
}

fn hits_tsv(results: &[Value]) -> (Vec<String>, Vec<Vec<String>>) {
    let mut header: Vec<String> = ["query_id", "hit_rank", "hit_id", "taxid", "shared"]
        .iter()
        .map(|s| s.to_string())
        .collect();
    header.extend(RANKS.iter().map(|r| r.to_string()));
    header.extend(["cluster".to_string(), "schema_version".to_string()]);
    let mut rows = Vec::new();
    for r in results {
        if let Some(hits) = r["guess"]["hits"].as_array() {
            for (i, h) in hits.iter().enumerate() {
                let mut row = vec![
                    tsv_cell(&r["id"]),
                    (i + 1).to_string(),
                    tsv_cell(&h["id"]),
                    tsv_cell(&h["taxid"]),
                    tsv_cell(&h["shared"]),
                ];
                for k in 0..6 {
                    row.push(tsv_cell(&h["lineage"][k]));
                }
                row.push(tsv_cell(&h["cluster"]));
                row.push(SCHEMA_VERSION.to_string());
                rows.push(row);
            }
        }
    }
    (header, rows)
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
                v["path"] == i.path.display().to_string()
                    && v["size_bytes"] == fs::metadata(&i.path).map(|m| m.len()).unwrap_or(0)
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
    eprintln!(
        "gmg run: {} inputs ({} to process, {} reused from a previous run); index {} genomes, k={}, scaled={}; thresholds: {}",
        inputs.len(), todo.len(), done.len(), ix.genomes.len(), ix.ksize, ix.scaled, th_name
    );

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
    let counter = AtomicUsize::new(0);
    let total = todo.len();
    let fresh: Vec<(String, Value)> = todo
        .par_iter()
        .map(|inp| {
            let v = process(&ix, &th, inp, &args);
            {
                let mut w = writer.lock().unwrap();
                let _ = writeln!(w, "{v}");
                let _ = w.flush();
            }
            let n = counter.fetch_add(1, Ordering::SeqCst) + 1;
            eprintln!(
                "[{n}/{total}] {}  {}  {}",
                inp.id,
                v["status"].as_str().unwrap_or(""),
                v["best_rank"].as_str().unwrap_or("-")
            );
            (inp.id.clone(), v)
        })
        .collect();
    let mut by_id = done;
    by_id.extend(fresh);
    let results: Vec<Value> = inputs.iter().filter_map(|i| by_id.remove(&i.id)).collect();

    let (h, rows) = results_tsv(&results);
    write_tsv(&args.out.join("results.tsv"), &h, &rows)?;
    let (h, rows) = hits_tsv(&results);
    write_tsv(&args.out.join("hits.tsv"), &h, &rows)?;

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
                  "k": ix.ksize, "scaled": ix.scaled, "seed": ix.seed, "n_genomes": ix.genomes.len()},
        "thresholds": {"name": th_name, "values": thresholds_json(&th)},
        "status_counts": counts,
        "inputs": results.iter().map(|r| json!({"id": r["id"], "path": r["path"], "size_bytes": r["size_bytes"],
                                                "kind": r["kind"], "status": r["status"]})).collect::<Vec<_>>(),
    });
    fs::write(
        args.out.join("run.json"),
        serde_json::to_string_pretty(&manifest)? + "\n",
    )?;
    let mut cs: Vec<_> = counts.iter().collect();
    cs.sort();
    eprintln!(
        "done in {:.1}s: {}  ->  {}/{{results.tsv,hits.tsv,run.json,results.jsonl}}",
        t0.elapsed().as_secs_f64(),
        cs.iter()
            .map(|(k, v)| format!("{k}={v}"))
            .collect::<Vec<_>>()
            .join(" "),
        args.out.display()
    );
    Ok(())
}
