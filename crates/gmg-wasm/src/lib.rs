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
    pub fn new(bytes: &[u8]) -> Result<Index, JsError> {
        let ix = RefIndex::read(&mut &bytes[..]).map_err(|e| JsError::new(&e.to_string()))?;
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
