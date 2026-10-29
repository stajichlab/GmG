//! JavaScript API around `gmg-core`. Build with `pixi run wasm`.
//!
//! ```js
//! const { Index } = require("./pkg-node/gmg_wasm.js");
//! const ix = new Index(new Uint8Array(fs.readFileSync("bfd.gmgi")));   // v0 or v1 index
//! const result = JSON.parse(ix.guess_sequence(fastaBytes, 1, 0n, 5));   // same JSON as `gmg guess`
//! ```

use gmg_core::{guess, guess_json, sketch_fastx, RefIndex, Thresholds};
use wasm_bindgen::prelude::*;

#[wasm_bindgen]
pub struct Index {
    ix: RefIndex,
    th: Thresholds,
}

#[wasm_bindgen]
impl Index {
    /// Load an index file (as written by `gmg build` / `gmg convert`).
    #[wasm_bindgen(constructor)]
    pub fn new(bytes: Vec<u8>) -> Result<Index, JsError> {
        // owned: a packed (v3) index keeps this buffer as its data, so the file is held once in memory
        let ix = RefIndex::read_owned(bytes).map_err(|e| JsError::new(&e.to_string()))?;
        Ok(Index {
            ix,
            th: Thresholds::default(),
        })
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
