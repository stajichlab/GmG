//! JSON rendering of a guess (feature `json`); shared by the CLI and the wasm wrapper so both print the same thing.

use crate::{Guess, RefIndex, Thresholds, RANKS};
use serde_json::{json, Map, Value};

pub fn guess_json(index: &RefIndex, g: &Guess, query: &str, top: usize) -> Value {
    let hits: Vec<Value> = g
        .hits
        .iter()
        .take(top)
        .map(|&(gi, sh)| {
            let m = &index.genomes[gi as usize];
            json!({"id": m.id, "taxid": m.taxid, "shared": sh, "lineage": m.lineage, "cluster": m.cluster})
        })
        .collect();
    let calls: Map<String, Value> = RANKS
        .iter()
        .zip(&g.calls)
        .zip(&g.votes)
        .map(|((r, c), v)| {
            (
                r.to_string(),
                json!({"taxon": c, "vote": v.taxon, "consistency": v.consistency}),
            )
        })
        .collect();
    let top_json = g.top.as_ref().map(|t| {
        json!({"id": index.genomes[t.genome as usize].id, "cluster": index.genomes[t.genome as usize].cluster,
               "shared": t.shared, "ani": t.ani, "containment_query": t.containment_query,
               "containment_ref": t.containment_ref})
    });
    let tentative = g
        .tentative_genus
        .as_ref()
        .map(|t| json!({"taxon": t.taxon, "ani": t.ani, "consistency": t.consistency}));
    json!({"query": query, "query_hashes": g.query_hashes, "top": top_json, "calls": calls, "tentative_genus": tentative, "hits": hits})
}

impl Thresholds {
    /// Parse thresholds from JSON with optional keys `min_shared[6]`, `min_ani[6]`, `min_consistency[6]`, `margin`, `top_n`
    /// (as written by `scripts/analyze_bench.py`); missing keys keep the defaults.
    pub fn from_json_str(text: &str) -> Result<Thresholds, String> {
        let v: Value = serde_json::from_str(text).map_err(|e| e.to_string())?;
        let mut th = Thresholds::default();
        let arr = |k: &str| v.get(k).and_then(|a| a.as_array()).filter(|a| a.len() == 6);
        if let Some(a) = arr("min_shared") {
            for (i, x) in a.iter().enumerate() {
                th.min_shared[i] = x.as_u64().ok_or("min_shared")? as u32;
            }
        }
        if let Some(a) = arr("min_ani") {
            for (i, x) in a.iter().enumerate() {
                th.min_ani[i] = x.as_f64().ok_or("min_ani")?;
            }
        }
        if let Some(a) = arr("min_consistency") {
            for (i, x) in a.iter().enumerate() {
                th.min_consistency[i] = x.as_f64().ok_or("min_consistency")?;
            }
        }
        if let Some(m) = v.get("margin").and_then(|m| m.as_f64()) {
            th.margin = m;
        }
        if let Some(n) = v.get("top_n").and_then(|n| n.as_u64()) {
            th.top_n = n as usize;
        }
        if let Some(t) = v.get("tentative_genus_min_ani").and_then(|t| t.as_f64()) {
            th.tentative_genus_min_ani = t;
        }
        Ok(th)
    }
}
