//! `gmg run`: classify every genome / read set in a folder and write `results.tsv`, `hits.tsv`, `run.json`.
//!
//! Each finished input is appended to `results.jsonl` immediately, so a killed run can be resumed (`--resume`);
//! the TSV files and the manifest are generated from that file at the end. Schema version: see `SCHEMA_VERSION`.

use crate::{load_index, load_thresholds, query_hashes_from, sketch_with_stats, strip_ext, Stats};
use anyhow::{bail, Context, Result};
use gmg_core::{containment, guess, guess_json, RefIndex, Thresholds, RANKS};
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
pub const SCHEMA_VERSION: u32 = 2; // 2: tentative_genus column; balanced default thresholds

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
    /// write report.html (default); `--no-report` turns it off
    pub report: bool,
    pub max_cards: usize,
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
           "margin": th.margin, "top_n": th.top_n,
           "tentative_genus_min_ani": th.tentative_genus_min_ani})
}

/// Provisional QC flags (issue #29 will refine them); they explain results, they never change a call.
/// Provisional `mixed_lineage` rule (#64): the second phylum owns at least this many exclusive hashes and share of the first's.
const MIXED_MIN_HASHES: u32 = 50;
const MIXED_MIN_RATIO: f64 = 0.10;

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

fn process(
    ix: &RefIndex,
    th: &Thresholds,
    inp: &Input,
    a: &RunArgs,
    keep_hashes: bool,
) -> (Value, Option<Vec<u64>>) {
    let size = fs::metadata(&inp.path).map(|m| m.len()).unwrap_or(0);
    let path = inp.path.display().to_string();
    let attempt = || -> Result<(Value, Option<Vec<u64>>)> {
        let reads = inp.kind == "reads";
        let min_abund = a.min_abund.unwrap_or(if reads { 2 } else { 1 });
        let (sk, st) = sketch_with_stats(&inp.path, ix.ksize, ix.scaled, reads, a.max_reads)?
            .with_context(|| "no sequence records".to_string())?;
        let hist = if reads {
            abundance_hist(&sk, min_abund)
        } else {
            None
        };
        let hashes = query_hashes_from(ix, sk, min_abund)?;
        let g = guess(ix, &hashes, &|_| false, th);
        let gj = guess_json(ix, &g, &path, a.top);
        let best = g.calls.iter().rposition(|c| c.is_some());
        let status = match best {
            None => "no_match",
            Some(i) if i >= 4 => "ok",
            Some(_) => "low_confidence",
        };
        let split = ix.phylum_exclusive(&hashes);
        let mut flags = qc_flags(
            inp.kind,
            &st,
            hashes.len(),
            gj.get("top").filter(|t| !t.is_null()),
        );
        if split.len() >= 2
            && split[1].1 >= MIXED_MIN_HASHES
            && split[1].1 as f64 >= MIXED_MIN_RATIO * split[0].1 as f64
        {
            flags.push("mixed_lineage");
        }
        let mut qc = json!({"n_records": st.n_records, "total_bases": st.total_len, "n50": st.n50,
                   "non_acgt_fraction": if st.total_len > 0 { st.non_acgt as f64 / st.total_len as f64 } else { 0.0 },
                   "n_hashes": hashes.len(), "min_abund": min_abund});
        if let Some(h) = hist {
            qc["abundance_hist"] = h;
        }
        qc["phylum_split"] = json!(split
            .iter()
            .take(3)
            .map(|(p, c)| json!({"phylum": p, "hashes": c}))
            .collect::<Vec<_>>());
        let rec = json!({
            "schema_version": SCHEMA_VERSION, "id": inp.id, "path": path, "size_bytes": size, "kind": inp.kind,
            "status": status, "best_rank": best.map(|i| RANKS[i]), "message": null,
            "qc": qc, "qc_flags": flags, "guess": gj,
        });
        Ok((
            rec,
            if keep_hashes && !reads {
                Some(hashes)
            } else {
                None
            },
        ))
    };
    match attempt() {
        Ok((v, h)) => (v, h),
        Err(e) => (
            json!({"schema_version": SCHEMA_VERSION, "id": inp.id, "path": path, "size_bytes": size, "kind": inp.kind,
               "status": "error", "best_rank": null, "message": format!("{e:#}"), "qc": null, "qc_flags": [], "guess": null}),
            None,
        ),
    }
}

/// How many times each sketch hash of a read set was seen (1..=12, then 13+), the abundance cutoff, and the median
/// abundance of the hashes kept (an estimate of k-mer coverage). Only for sketches that carry abundances.
fn abundance_hist(sk: &gmg_core::SketchFile, min_abund: u32) -> Option<Value> {
    let ab = sk.abunds.as_ref()?;
    let mut counts = [0u64; 12];
    let mut over = 0u64;
    let mut kept: Vec<u32> = Vec::new();
    for &a in ab {
        match a {
            1..=12 => counts[a as usize - 1] += 1,
            _ => over += 1,
        }
        if a >= min_abund {
            kept.push(a);
        }
    }
    kept.sort_unstable();
    let median = kept.get(kept.len() / 2).copied().unwrap_or(0);
    Some(json!({"counts": counts, "over12": over, "min_abund": min_abund, "median_kept": median}))
}

/// Pairwise containment ANI (geometric mean of the two directions) between assemblies, for the report heatmap.
fn relationships(k: u32, items: &[(String, String, Vec<u64>)]) -> Value {
    // order by taxon label so related genomes sit next to each other (species blocks along the diagonal)
    let mut items: Vec<&(String, String, Vec<u64>)> = items.iter().collect();
    items.sort_by(|a, b| a.1.cmp(&b.1).then(a.0.cmp(&b.0)));
    let n = items.len();
    let mut m = vec![vec![1.0f64; n]; n];
    for i in 0..n {
        for j in (i + 1)..n {
            let (a, b) = (&items[i].2, &items[j].2);
            let v = (containment(a, b) * containment(b, a))
                .sqrt()
                .powf(1.0 / k as f64);
            let v = (v * 10_000.0).round() / 10_000.0;
            m[i][j] = v;
            m[j][i] = v;
        }
    }
    json!({"ids": items.iter().map(|t| t.0.clone()).collect::<Vec<_>>(), "labels": items.iter().map(|t| t.1.clone()).collect::<Vec<_>>(), "ani": m})
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
    header.push("tentative_genus".to_string());
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
        row.push(tsv_cell(&g["tentative_genus"]["taxon"]));
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
    let n_asm = inputs.iter().filter(|i| i.kind == "assembly").count();
    let want_rel = (2..=60).contains(&n_asm);
    let fresh: Vec<(String, Value, Option<Vec<u64>>)> = todo
        .par_iter()
        .map(|inp| {
            let (v, h) = process(&ix, &th, inp, &args, want_rel);
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
                let g = &r["guess"]["calls"];
                let label = g["species"]["taxon"]
                    .as_str()
                    .or(g["genus"]["taxon"].as_str())
                    .unwrap_or(&id)
                    .to_string();
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
    if args.report {
        write_report(&args.out, args.max_cards)?;
    }
    let mut cs: Vec<_> = counts.iter().collect();
    cs.sort();
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
    Ok(())
}

/// Render report.html from the files of a finished run in `dir` (results.jsonl, run.json, optional relationships.json).
pub fn write_report(dir: &Path, max_cards: usize) -> Result<()> {
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
    fs::write(
        dir.join("report.html"),
        gmg_report::render_html(&results, &run, rel.as_ref(), &opts),
    )?;
    Ok(())
}
