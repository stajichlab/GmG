//! Building the result record of one input (and the TSV views of many) from a sketch: the single implementation behind
//! `gmg run` and the browser (gmg-wasm), so both produce the same `results.jsonl` / `results.tsv` for the same input.
//! Feature `json`.

use crate::{downsample, guess, guess_json, RefIndex, SketchFile, Thresholds, RANKS};
use serde_json::{json, Value};
use std::collections::HashMap;

/// Bump when a column of results.tsv / hits.tsv or a key of results.jsonl / run.json changes meaning.
pub const SCHEMA_VERSION: u32 = 6; // 6: taxid_<rank> columns, top_hit_label. 5: withheld column. 4: expected, expected_match columns (sample sheet). 2: tentative_genus column; balanced default thresholds. 3: error_code column, plain error messages

/// Species ANI cutoff for read sets (#32): repeated sequencing errors lower the containment ANI of a read set, so the
/// assembly cutoff (0.9938) leaves 8-15% of correct species unnamed. With *C. immitis* left out of the reference its reads
/// match *C. posadasii* at ANI 0.976-0.980 and are named at 0.99 in 0.6% of sets; at 0.985 also 0.6%, at 0.98 10%.
/// Only lowers a cutoff that is above this value and is still reachable (<= 1.0, so `novel` stays off).
pub const READS_SPECIES_MIN_ANI: f64 = 0.99;

pub use crate::fastx::SeqStats;

/// Expected genome size range (Mb) of the called phylum, from the reference (`RefIndex::size_range_mb`), by phylum name.
pub type SizeRanges = HashMap<String, Option<(f64, f64, String)>>;

/// Per-reference facts used when classifying: expected genome size ranges, and the phyla that are not fungi
/// (decoy genomes added to the reference, listed in the manifest as `non_fungal_phyla`).
pub struct Context {
    pub ranges: SizeRanges,
    pub non_fungal: Vec<String>,
}

impl Context {
    pub fn new(ix: &RefIndex) -> Context {
        let non_fungal = ix
            .manifest
            .as_deref()
            .and_then(|m| serde_json::from_str::<Value>(m).ok())
            .and_then(|v| {
                v["non_fungal_phyla"].as_array().map(|a| {
                    a.iter()
                        .filter_map(|x| x.as_str().map(String::from))
                        .collect()
                })
            })
            .unwrap_or_default();
        Context {
            ranges: size_ranges(ix),
            non_fungal,
        }
    }
}

/// What identifies an input in its record.
pub struct Item<'a> {
    pub id: &'a str,
    /// as shown in results (the CLI passes the file path; the browser the file name)
    pub path: &'a str,
    pub size_bytes: u64,
    /// `assembly` or `reads`
    pub kind: &'a str,
}

/// Best-hit ANI to a non-fungal decoy at or above which an input is flagged `not_fungal` (the largest ANI of a real fungus to a
/// decoy in the 1,257-genome post-freeze set plus IndigoAg was 0.881; non-fungal relatives of a decoy are well above 0.9).
pub const NOT_FUNGAL_MIN_ANI: f64 = 0.90;

/// Classify one sketched input. Returns the record (the line of `results.jsonl`) and the query hashes it was classified
/// with, or a message when the input cannot be classified (too short, k/seed mismatch).
#[allow(clippy::too_many_arguments)]
pub fn classify(
    ix: &RefIndex,
    th: &Thresholds,
    ctx: &Context,
    item: &Item,
    sk: SketchFile,
    st: &SeqStats,
    min_abund: u32,
    top: usize,
    skip: &dyn Fn(u32) -> bool,
) -> Result<(Value, Vec<u64>), String> {
    let reads = item.kind == "reads";
    if st.total_len == 0 {
        return Err("no sequence records".into());
    }
    if st.total_len < ix.ksize as u64 {
        return Err("input is too short for one k-mer window".into());
    }
    let hist = if reads {
        abundance_hist(&sk, min_abund)
    } else {
        None
    };
    let hashes = query_hashes_from(ix, sk, min_abund)?;
    let reads_th;
    let th = if reads && th.min_ani[5] > READS_SPECIES_MIN_ANI && th.min_ani[5] <= 1.0 {
        let mut t = th.clone();
        t.min_ani[5] = READS_SPECIES_MIN_ANI;
        reads_th = t;
        &reads_th
    } else {
        th
    };
    let is_decoy = |gi: u32| {
        let ph = &ix.genomes[gi as usize].lineage[0];
        ctx.non_fungal.iter().any(|n| n == ph)
    };
    let mut g = guess(ix, &hashes, skip, th);
    // the input is nearest to non-fungal decoys: its phylum call, or failing that its best hit, is one of the listed phyla.
    // A real relative (ANI >= NOT_FUNGAL_MIN_ANI) says "not a fungus". A distant best hit says nothing either way, because a
    // divergent fungus and a divergent non-fungus both land on whichever decoy is nearest by chance (ANI 0.75-0.88): the call
    // is then made again without the decoys, so a divergent fungus gets the call it would get from an index without them
    // (instead of e.g. a Chordata phylum), and the flag `distant_non_fungal_hit` records what happened.
    let near_decoy = !ctx.non_fungal.is_empty()
        && (g.calls[0]
            .as_deref()
            .is_some_and(|p| ctx.non_fungal.iter().any(|n| n == p))
            || g.top.as_ref().is_some_and(|t| is_decoy(t.genome)));
    let top_ani = g.top.as_ref().map_or(0.0, |t| t.ani);
    let not_fungal = near_decoy && top_ani >= NOT_FUNGAL_MIN_ANI;
    let distant_decoy = near_decoy && !not_fungal;
    if distant_decoy {
        g = guess(ix, &hashes, &|gi| skip(gi) || is_decoy(gi), th);
    }
    let mut gj = guess_json(ix, &g, item.path, top);
    gj["withheld"] = withheld_json(&g, th);
    let (best, status) = status_of(&g);
    let split = ix.phylum_exclusive(&hashes);
    let est_mb = hashes.len() as f64 * ix.scaled as f64 / 1e6;
    let range = if reads {
        None
    } else {
        g.calls[0]
            .as_deref()
            .and_then(|ph| ctx.ranges.get(ph))
            .filter(|r| r.is_some())
            .or_else(|| ctx.ranges.get(""))
            .and_then(|r| r.as_ref())
    };
    let mut flags = qc_flags(
        item.kind,
        st,
        hashes.len(),
        gj.get("top").filter(|t| !t.is_null()),
        range.map(|r| (est_mb, r)),
    );
    if not_fungal {
        flags.push("not_fungal");
    } else if distant_decoy {
        flags.push("distant_non_fungal_hit");
    }
    let mut qc = json!({"n_records": st.n_records, "total_bases": st.total_len, "n50": st.n50,
               "non_acgt_fraction": if st.total_len > 0 { st.non_acgt as f64 / st.total_len as f64 } else { 0.0 },
               "n_hashes": hashes.len(), "min_abund": min_abund});
    if let Some(h) = hist {
        qc["abundance_hist"] = h;
    }
    if let Some((lo, hi, basis)) = range {
        qc["size_est_mb"] = json!((est_mb * 100.0).round() / 100.0);
        qc["size_range_mb"] = json!([(lo * 100.0).round() / 100.0, (hi * 100.0).round() / 100.0]);
        qc["size_range_basis"] = json!(basis);
    }
    qc["phylum_split"] = json!(split
        .iter()
        .take(3)
        .map(|(p, c)| json!({"phylum": p, "hashes": c}))
        .collect::<Vec<_>>());
    let rec = json!({
        "schema_version": SCHEMA_VERSION, "id": item.id, "path": item.path, "size_bytes": item.size_bytes, "kind": item.kind,
        "status": status, "best_rank": best.map(|i| RANKS[i]), "message": null,
        "qc": qc, "qc_flags": flags, "guess": gj,
    });
    Ok((rec, hashes))
}

/// The record of an input that could not be processed.
pub fn error_record(item: &Item, code: &str, message: &str, action: &str) -> Value {
    json!({"schema_version": SCHEMA_VERSION, "id": item.id, "path": item.path, "size_bytes": item.size_bytes, "kind": item.kind,
           "status": "error", "best_rank": null, "error_code": code, "message": message,
           "error_action": action, "qc": null, "qc_flags": [], "guess": null})
}

/// Why a file held in memory cannot be read as FASTA/FASTQ (the browser has no file system to inspect; gzip is
/// decompressed before this). Mirrors the codes of the CLI (`docs/results-schema.md`).
pub fn bytes_error(bytes: &[u8], reason: &str) -> (&'static str, &'static str, &'static str) {
    if bytes.is_empty() {
        return (
            "empty_file",
            "The file is empty (0 bytes).",
            "Re-download or re-copy it.",
        );
    }
    match bytes.iter().find(|b| !b.is_ascii_whitespace()) {
        None => (
            "no_sequence_records",
            "The file has no sequence records (only blank content).",
            "Check that it is the right file.",
        ),
        Some(b'>') | Some(b'@') => {
            if reason.contains("too short") {
                (
                    "too_short",
                    "The input has fewer bases than one k-mer window.",
                    "Provide a real genome or read set.",
                )
            } else {
                (
                    "no_sequence_records",
                    "The file has a header but no sequence.",
                    "Check that the download or conversion finished.",
                )
            }
        }
        Some(_) => (
            "not_sequence_file",
            "The file does not look like FASTA or FASTQ.",
            "GmG reads FASTA (starting with '>') or FASTQ (starting with '@'), optionally gzipped.",
        ),
    }
}

pub fn size_ranges(ix: &RefIndex) -> SizeRanges {
    let mut m: SizeRanges = HashMap::new();
    m.insert(String::new(), ix.size_range_mb(None));
    for ph in ix.taxa[0].iter().skip(1) {
        m.insert(ph.clone(), ix.size_range_mb(Some(ph)));
    }
    m
}

/// QC flags explain a result, they never change a call (docs/qc.md has the rationale for every threshold).
pub fn qc_flags(
    kind: &str,
    st: &SeqStats,
    n_hashes: usize,
    top: Option<&Value>,
    size: Option<(f64, &(f64, f64, String))>,
) -> Vec<&'static str> {
    let mut f = Vec::new();
    if kind == "assembly" {
        // size: the sketch-based estimate against the 1st-99th percentile of the called phylum in the reference
        // (0.5th-99.5th of all genomes when no phylum is called or the phylum is small)
        if let Some((est, (lo, hi, _))) = size {
            if est < *lo {
                f.push("short_genome");
            }
            if est > *hi {
                f.push("long_genome");
            }
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

/// How many times each sketch hash of a read set was seen (1..=12, then 13+), the abundance cutoff, and the median
/// abundance of the hashes kept (an estimate of k-mer coverage). Only for sketches that carry abundances.
pub fn abundance_hist(sk: &SketchFile, min_abund: u32) -> Option<Value> {
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

pub fn tsv_cell(v: &Value) -> String {
    let s = match v {
        Value::Null => String::new(),
        Value::String(s) => s.clone(),
        other => other.to_string(),
    };
    s.replace(['\t', '\n', '\r'], " ")
}

pub fn results_tsv(results: &[Value]) -> (Vec<String>, Vec<Vec<String>>) {
    let mut header: Vec<String> = ["query_id", "path", "kind", "status", "best_rank"]
        .iter()
        .map(|s| s.to_string())
        .collect();
    for r in RANKS {
        header.push(format!("call_{r}"));
    }
    for r in RANKS {
        header.push(format!("taxid_{r}"));
    }
    header.push("tentative_genus".to_string());
    for r in RANKS {
        header.push(format!("conf_{r}"));
    }
    for c in [
        "top_hit",
        "top_hit_species",
        "top_hit_label",
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
        "error_code",
        "withheld",
        "expected",
        "expected_match",
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
            row.push(tsv_cell(&g["calls"][k]["taxid"]));
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
            tsv_cell(&top_hit["label"]),
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
            tsv_cell(&r["error_code"]),
            withheld_text(&r["guess"]["withheld"]),
            tsv_cell(&r["expected"]),
            expected_match(r),
            SCHEMA_VERSION.to_string(),
        ]);
        rows.push(row);
    }
    (header, rows)
}

/// `yes` when the expected label is the call at some rank (case-insensitive), `no` when it is not, blank when none was given.
pub fn expected_match(r: &Value) -> String {
    let Some(e) = r["expected"].as_str().filter(|e| !e.is_empty()) else {
        return String::new();
    };
    let hit = RANKS.iter().any(|k| {
        r["guess"]["calls"][k]["taxon"]
            .as_str()
            .is_some_and(|t| t.eq_ignore_ascii_case(e))
    });
    if hit { "yes" } else { "no" }.to_string()
}

pub fn hits_tsv(results: &[Value]) -> (Vec<String>, Vec<Vec<String>>) {
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

/// Why the first uncalled rank was withheld, as JSON (`null` when a species was called).
pub fn withheld_json(g: &crate::Guess, th: &Thresholds) -> serde_json::Value {
    match g.withheld(th) {
        Some(w) => serde_json::json!({"rank": RANKS[w.rank], "code": w.code, "detail": w.detail}),
        None => serde_json::Value::Null,
    }
}

/// `genus: ani_below_cutoff (ANI 0.8870 < 0.9000)`, or an empty string.
pub fn withheld_text(v: &serde_json::Value) -> String {
    match (v["rank"].as_str(), v["code"].as_str(), v["detail"].as_str()) {
        (Some(r), Some(c), Some(d)) => format!("{r}: {c} ({d})"),
        _ => String::new(),
    }
}

/// `no_match` / `low_confidence` / `ok` from the deepest rank called.
pub fn status_of(g: &crate::Guess) -> (Option<usize>, &'static str) {
    let best = g.calls.iter().rposition(|c| c.is_some());
    let status = match best {
        None => "no_match",
        Some(i) if i >= 4 => "ok",
        Some(_) => "low_confidence",
    };
    (best, status)
}

/// Check a sketch against the index (k, seed, scaled), drop hashes seen fewer than `min_abund` times, downsample to the index.
pub fn query_hashes_from(
    ix: &RefIndex,
    sk: SketchFile,
    min_abund: u32,
) -> Result<Vec<u64>, String> {
    if sk.ksize != ix.ksize || sk.seed != ix.seed {
        return Err(format!(
            "{}: k/seed ({}/{}) do not match the index ({}/{})",
            sk.name, sk.ksize, sk.seed, ix.ksize, ix.seed
        ));
    }
    if sk.scaled > ix.scaled {
        return Err(format!(
            "{}: query scaled {} is coarser than the index scaled {}",
            sk.name, sk.scaled, ix.scaled
        ));
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

/// Label of an assembly in the relationship heatmap: its species call, else its genus call, else its id.
pub fn relationship_label(r: &Value) -> String {
    let g = &r["guess"]["calls"];
    g["species"]["taxon"]
        .as_str()
        .or(g["genus"]["taxon"].as_str())
        .or(r["id"].as_str())
        .unwrap_or("")
        .to_string()
}

/// Pairwise ANI between assemblies (id, label, hashes): the `relationships.json` of a run, ordered by label so related
/// genomes sit next to each other.
pub fn relationships(k: u32, items: &[(String, String, Vec<u64>)]) -> Value {
    // order by taxon label so related genomes sit next to each other (species blocks along the diagonal)
    let mut items: Vec<&(String, String, Vec<u64>)> = items.iter().collect();
    items.sort_by(|a, b| a.1.cmp(&b.1).then(a.0.cmp(&b.0)));
    let n = items.len();
    let mut m = vec![vec![1.0f64; n]; n];
    for i in 0..n {
        for j in (i + 1)..n {
            let (a, b) = (&items[i].2, &items[j].2);
            let v = (crate::containment(a, b) * crate::containment(b, a))
                .sqrt()
                .powf(1.0 / k as f64);
            let v = (v * 10_000.0).round() / 10_000.0;
            m[i][j] = v;
            m[j][i] = v;
        }
    }
    json!({"ids": items.iter().map(|t| t.0.clone()).collect::<Vec<_>>(), "labels": items.iter().map(|t| t.1.clone()).collect::<Vec<_>>(), "ani": m})
}

/// How much of the reference a query may not see (validation of queries that are themselves in the reference).
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum ExcludeLevel {
    /// only the query's own reference genome
    SelfOnly,
    /// every reference genome of the query's species
    Species,
    /// every reference genome of the query's genus
    Genus,
}

impl std::str::FromStr for ExcludeLevel {
    type Err = String;
    fn from_str(s: &str) -> Result<Self, String> {
        match s {
            "self" => Ok(Self::SelfOnly),
            "species" => Ok(Self::Species),
            "genus" => Ok(Self::Genus),
            _ => Err(format!(
                "unknown exclusion level '{s}' (self, species or genus)"
            )),
        }
    }
}

/// Reference genomes to hide from one query: the ids in `ids`, plus (when `level` is given) the query's own genome
/// (`ref_id`), or every genome sharing its species or genus. The species/genus come from the reference genome `ref_id`
/// when it is in the index, otherwise from `expected` ("Genus species" or a genus). Returns a flag per genome.
pub fn exclusion_mask(
    ix: &RefIndex,
    ids: &std::collections::HashSet<String>,
    level: Option<ExcludeLevel>,
    ref_id: Option<&str>,
    expected: Option<&str>,
) -> Vec<bool> {
    let mut mask: Vec<bool> = ix.genomes.iter().map(|g| ids.contains(&g.id)).collect();
    let Some(level) = level else { return mask };
    let own = ref_id.and_then(|r| ix.genomes.iter().position(|g| g.id == r));
    if let Some(i) = own {
        mask[i] = true;
    }
    if level == ExcludeLevel::SelfOnly {
        return mask;
    }
    let rank = if level == ExcludeLevel::Species { 5 } else { 4 };
    let target: Option<String> = match own {
        Some(i) => Some(ix.genomes[i].lineage[rank].clone()).filter(|t| !t.is_empty()),
        None => expected.map(|e| {
            let e = e.trim();
            if rank == 4 {
                e.split_whitespace().next().unwrap_or(e).to_string()
            } else {
                e.to_string()
            }
        }),
    };
    if let Some(t) = target {
        for (i, g) in ix.genomes.iter().enumerate() {
            if g.lineage[rank] == t {
                mask[i] = true;
            }
        }
    }
    mask
}
