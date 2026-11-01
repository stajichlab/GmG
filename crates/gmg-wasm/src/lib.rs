//! JavaScript API around `gmg-core`. Build with `pixi run wasm`.
//!
//! ```js
//! const { Index } = require("./pkg-node/gmg_wasm.js");
//! const ix = new Index(new Uint8Array(fs.readFileSync("bfd.gmgi")));   // v0 or v1 index
//! const result = JSON.parse(ix.guess_sequence(fastaBytes, 1, 0n, 5));   // same JSON as `gmg guess`
//! ```

use gmg_core::classify::{
    bytes_error, classify, error_record, hits_tsv, relationship_label, relationships, results_tsv,
    Context, Item, READS_SPECIES_MIN_ANI, SCHEMA_VERSION,
};
use gmg_core::{guess, guess_json, sketch_fastx, sketch_fastx_stats, RefIndex, Thresholds};
use std::cell::RefCell;
use std::collections::HashMap;
use wasm_bindgen::prelude::*;

#[wasm_bindgen]
pub struct Index {
    ix: RefIndex,
    th: Thresholds,
    th_name: String,
    ctx: Context,
    /// hashes of the assemblies classified in this session (for the relationship heatmap)
    session: RefCell<HashMap<String, Vec<u64>>>,
}

#[wasm_bindgen]
impl Index {
    /// Load an index file (as written by `gmg build` / `gmg convert`).
    #[wasm_bindgen(constructor)]
    pub fn new(bytes: Vec<u8>) -> Result<Index, JsError> {
        // owned: a packed (v3) index keeps this buffer as its data, so the file is held once in memory
        let ix = RefIndex::read_owned(bytes).map_err(|e| JsError::new(&e.to_string()))?;
        let ctx = Context::new(&ix);
        Ok(Index {
            ix,
            th: Thresholds::default(),
            th_name: "known".to_string(),
            ctx,
            session: RefCell::new(HashMap::new()),
        })
    }

    /// The reference's manifest (JSON, or an empty string for indexes built before manifests).
    pub fn manifest(&self) -> String {
        self.ix.manifest.clone().unwrap_or_default()
    }

    /// Choose a built-in threshold set: `known` (default), `permissive` or `novel`.
    pub fn use_thresholds(&mut self, name: &str) -> Result<(), JsError> {
        let text = match name {
            "known" => None,
            "permissive" => Some(include_str!(
                "../../../data/thresholds/bfd_s10000.permissive.json"
            )),
            "novel" => Some(include_str!(
                "../../../data/thresholds/bfd_s10000.novel.json"
            )),
            other => return Err(JsError::new(&format!("unknown threshold set '{other}'"))),
        };
        self.th = match text {
            None => Thresholds::default(),
            Some(t) => Thresholds::from_json_str(t).map_err(|e| JsError::new(&e))?,
        };
        self.th_name = name.to_string();
        Ok(())
    }

    /// Classify one file held in memory and return its record, the same line `gmg run` writes to `results.jsonl`
    /// (also for failures: status `error` with a code and message). `kind`: `assembly`, `reads` or `auto` (FASTQ when the
    /// first byte is `@`). `min_abund` 0 = default (1 assemblies, 2 reads); `max_records` 0 = all. gzip must already be
    /// decompressed.
    #[allow(clippy::too_many_arguments)]
    pub fn classify(
        &self,
        id: &str,
        name: &str,
        bytes: &[u8],
        kind: &str,
        min_abund: u32,
        max_records: u64,
        top: usize,
    ) -> String {
        let kind = match kind {
            "reads" => "reads",
            "assembly" => "assembly",
            _ => {
                if bytes.first() == Some(&b'@') {
                    "reads"
                } else {
                    "assembly"
                }
            }
        };
        let item = Item {
            id,
            path: name,
            size_bytes: bytes.len() as u64,
            kind,
        };
        let reads = kind == "reads";
        let min_abund = if min_abund == 0 {
            if reads {
                2
            } else {
                1
            }
        } else {
            min_abund
        };
        let fail = |reason: &str| {
            let (code, msg, action) = bytes_error(bytes, reason);
            error_record(&item, code, msg, action).to_string()
        };
        let Some((sk, st)) =
            sketch_fastx_stats(bytes, self.ix.ksize, self.ix.scaled, reads, max_records)
        else {
            return fail("");
        };
        match classify(
            &self.ix,
            &self.th,
            &self.ctx,
            &item,
            sk,
            &st,
            min_abund,
            top,
            &|_| false,
        ) {
            Ok((rec, hashes)) => {
                if !reads {
                    self.session.borrow_mut().insert(id.to_string(), hashes);
                }
                rec.to_string()
            }
            Err(msg) => fail(&msg),
        }
    }

    /// Pairwise ANI between the assemblies among `results_jsonl` classified in this session, as `relationships.json`
    /// text for `render_report`; empty when there are fewer than 2 or more than 60.
    pub fn relationships(&self, results_jsonl: &str) -> String {
        let session = self.session.borrow();
        let items: Vec<(String, String, Vec<u64>)> = gmg_report::parse_results_jsonl(results_jsonl)
            .iter()
            .filter(|r| r["kind"] == "assembly" && r["status"] != "error")
            .filter_map(|r| {
                let id = r["id"].as_str()?.to_string();
                Some((id.clone(), relationship_label(r), session.get(&id)?.clone()))
            })
            .collect();
        if !(2..=60).contains(&items.len()) {
            return String::new();
        }
        relationships(self.ix.ksize, &items).to_string()
    }

    /// Forget the hashes kept for the heatmap (call when starting a new batch).
    pub fn clear_session(&self) {
        self.session.borrow_mut().clear();
    }

    /// The `run.json` the report needs (index, thresholds, timing), for a run that started at `started_utc`.
    pub fn run_json(&self, started_utc: &str, elapsed_seconds: f64) -> String {
        serde_json::json!({
            "schema_version": SCHEMA_VERSION, "tool": "gmg-wasm", "version": env!("CARGO_PKG_VERSION"),
            "commit": "web", "started_utc": started_utc, "elapsed_seconds": elapsed_seconds,
            "index": {"n_genomes": self.ix.genomes.len(), "k": self.ix.ksize, "scaled": self.ix.scaled, "seed": self.ix.seed,
                      "manifest": self.ix.manifest.as_deref().and_then(|m| serde_json::from_str::<serde_json::Value>(m).ok())},
            "thresholds": {"name": self.th_name,
                           "values": {"min_shared": self.th.min_shared, "min_ani": self.th.min_ani,
                                      "min_consistency": self.th.min_consistency, "margin": self.th.margin,
                                      "top_n": self.th.top_n, "tentative_genus_min_ani": self.th.tentative_genus_min_ani},
                           "reads_species_min_ani": READS_SPECIES_MIN_ANI},
        })
        .to_string()
    }

    pub fn ksize(&self) -> u32 {
        self.ix.ksize
    }
    pub fn scaled(&self) -> u32 {
        self.ix.scaled
    }
    pub fn n_genomes(&self) -> usize {
        self.ix.genomes.len()
    }

    /// Replace the call thresholds with a JSON string (same keys as `gmg guess --thresholds`).
    pub fn set_thresholds(&mut self, json: &str) -> Result<(), JsError> {
        self.th = Thresholds::from_json_str(json).map_err(|e| JsError::new(&e))?;
        Ok(())
    }

    /// Sketch a FASTA or FASTQ file held in memory at this index's k/scaled. `min_abund` 1 for assemblies, 2 for reads;
    /// `max_records` 0 = all (reads: stop after that many).
    pub fn sketch(&self, bytes: &[u8], min_abund: u32, max_records: u64) -> Vec<u64> {
        sketch_fastx(
            bytes,
            self.ix.ksize as usize,
            self.ix.scaled,
            self.ix.seed,
            min_abund,
            max_records,
        )
    }

    /// Guess from sorted unique hashes (e.g. from `sketch`); returns the JSON of `gmg guess`.
    pub fn guess_hashes(&self, hashes: &[u64], top: usize) -> String {
        let g = guess(&self.ix, hashes, &|_| false, &self.th);
        guess_json(&self.ix, &g, "", top).to_string()
    }

    /// Sketch then guess.
    pub fn guess_sequence(
        &self,
        bytes: &[u8],
        min_abund: u32,
        max_records: u64,
        top: usize,
    ) -> String {
        self.guess_hashes(&self.sketch(bytes, min_abund, max_records), top)
    }
}

/// The GmG version this wasm package was built from (shown in the page footer).
#[wasm_bindgen]
pub fn version() -> String {
    env!("CARGO_PKG_VERSION").to_string()
}

/// `results.tsv` text from `results.jsonl` text (identical to the file `gmg run` writes).
#[wasm_bindgen]
pub fn make_results_tsv(results_jsonl: &str) -> String {
    tsv_text(results_tsv(&gmg_report::parse_results_jsonl(results_jsonl)))
}

/// `hits.tsv` text from `results.jsonl` text.
#[wasm_bindgen]
pub fn make_hits_tsv(results_jsonl: &str) -> String {
    tsv_text(hits_tsv(&gmg_report::parse_results_jsonl(results_jsonl)))
}

fn tsv_text((header, rows): (Vec<String>, Vec<Vec<String>>)) -> String {
    let mut out = header.join("\t");
    out.push('\n');
    for r in rows {
        out.push_str(&r.join("\t"));
        out.push('\n');
    }
    out
}

/// Render the standalone HTML report from the files of a `gmg run`: `results.jsonl` text, `run.json` text and
/// (optionally, empty string if absent) `relationships.json` text. Identical to the report the CLI writes, minus the
/// "Download results.tsv" link.
#[wasm_bindgen]
pub fn render_report(
    results_jsonl: &str,
    run_json: &str,
    relationships_json: &str,
) -> Result<String, JsError> {
    let results = gmg_report::parse_results_jsonl(results_jsonl);
    let run: serde_json::Value =
        serde_json::from_str(run_json).map_err(|e| JsError::new(&format!("run.json: {e}")))?;
    let rel: Option<serde_json::Value> = if relationships_json.trim().is_empty() {
        None
    } else {
        serde_json::from_str(relationships_json).ok()
    };
    let opts = gmg_report::Options {
        tsv_link: None,
        ..Default::default()
    };
    Ok(gmg_report::render_html(&results, &run, rel.as_ref(), &opts))
}
