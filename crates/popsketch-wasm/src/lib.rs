//! JavaScript API around `popsketch` (single thread, no files). Build with `pixi run popsketch-wasm`.
//!
//! ```js
//! const { Session } = require("./pkg-node/popsketch_wasm.js");
//! const s = new Session(21, 200, 0n, false);                      // k, scaled, max bases per sample (0 = all), R1 only
//! const withHalves = JSON.parse(s.halves_for(JSON.stringify(ids), 20));   // ids that need split halves (noise floor)
//! s.add_sample(id, "reads", r1Bytes, r2BytesOrEmpty, withHalves.includes(id));
//! const out = JSON.parse(s.analyze(labelsTextOrUndefined, JSON.stringify({ min_cov: 4 }), newickTextOrUndefined, pcsTextOrUndefined));
//! // out.summary (the popsketch.json document), out.html, out.clusters, out.pca, out.newick, out.qc
//! ```
//! Reading is the caller's job: the page slices each (decompressed) FASTQ and hands over the bytes, or hands over sketches that
//! `popsketch sketch` made earlier (`add_sketch`), so the reads are not needed at all.

use gmg_core::RefIndex;
use popsketch::floor::pick_for_halves;
use popsketch::pipeline::{analyze, summary_json, Analysis, FlagMode, Options, SampleRow, Status};
use popsketch::sketch_io::{sketch_bytes, sketch_from_bytes, SampleSketch, SketchParams};
use serde_json::{json, Value};
use std::path::PathBuf;
use wasm_bindgen::prelude::*;

#[wasm_bindgen]
pub struct Session {
    params: SketchParams,
    sketches: Vec<SampleSketch>,
    errors: Vec<SampleRow>,
    warnings: Vec<String>,
    index: Option<RefIndex>,
}

#[wasm_bindgen]
impl Session {
    #[wasm_bindgen(constructor)]
    pub fn new(ksize: u32, scaled: u32, max_bases: u64, r1_only: bool) -> Session {
        Session {
            params: SketchParams {
                ksize,
                scaled,
                max_bases,
                r1_only,
                ..SketchParams::default()
            },
            sketches: Vec::new(),
            errors: Vec::new(),
            warnings: Vec::new(),
            index: None,
        }
    }

    /// Load a reference index (a `.gmgi` file) so `analyze` also screens each sample against it: the organism and
    /// mixture calls of the command line's `--index`. The index must use the same k and a scale at least as coarse.
    pub fn set_index(&mut self, bytes: Vec<u8>) -> Result<(), JsError> {
        self.index = Some(RefIndex::read_owned(bytes).map_err(|e| JsError::new(&e.to_string()))?);
        Ok(())
    }

    /// The ids (JSON array in, JSON array out) that should be sketched with split halves, the same pseudo-random subset the CLI picks.
    pub fn halves_for(&self, ids_json: &str, n: usize) -> Result<String, JsError> {
        let ids: Vec<String> =
            serde_json::from_str(ids_json).map_err(|e| JsError::new(&e.to_string()))?;
        let mut v: Vec<String> = pick_for_halves(&ids, n).into_iter().collect();
        v.sort();
        Ok(json!(v).to_string())
    }

    /// Sketch one sample from the bytes of its files (`kind` is `reads` or `assembly`; `r2` may be empty). A sample that
    /// cannot be read is recorded as an error row and reported by `analyze`, not thrown.
    pub fn add_sample(&mut self, id: &str, kind: &str, r1: &[u8], r2: &[u8], want_halves: bool) {
        let mut parts: Vec<&[u8]> = vec![r1];
        if !r2.is_empty() {
            parts.push(r2);
        }
        match sketch_bytes(id, kind, &parts, &self.params, want_halves) {
            Ok(s) => self.sketches.push(s),
            Err(e) => {
                self.warnings.push(format!("{id}: {e:#}"));
                self.errors.push(SampleRow {
                    id: id.to_string(),
                    kind: kind.to_string(),
                    files: Vec::new(),
                    status: Status::Error,
                    qc: None,
                    index: None,
                    message: Some(format!("{e:#}")),
                });
            }
        }
    }

    /// Add a sample that was already sketched, from the contents of its files as `popsketch sketch` writes them: `<id>.json`,
    /// `<id>.gsk` and the two half sketches (`<id>.half_a.gsk`, `<id>.half_b.gsk`; pass empty arrays for a sample without
    /// them). The sketch must have the session's k and scale, or this throws and the sample is not added.
    pub fn add_sketch(
        &mut self,
        meta_json: &[u8],
        whole: &[u8],
        half_a: &[u8],
        half_b: &[u8],
    ) -> Result<(), JsError> {
        let halves = (!half_a.is_empty() && !half_b.is_empty()).then_some((half_a, half_b));
        let s = sketch_from_bytes(meta_json, whole, halves, &self.params)
            .map_err(|e| JsError::new(&format!("{e:#}")))?;
        self.sketches.push(s);
        Ok(())
    }

    pub fn n_samples(&self) -> usize {
        self.sketches.len()
    }

    /// Audit, distances, clusters, hash PCA, tree and label analysis. `options_json` may set `min_cov`, `size_z`,
    /// `size_floor`, `shoulder_z`, `min_floors`, `k`, `components`, `floor_samples`, `tree`, `no_pca`, `flagged`, `strata`, `permutations`, `collapse` (none, mean, max or median), `collapse_within` (label columns, comma-separated, or all)
    /// (`exclude`, `mark`, `decontaminate`), `root`, `strata` (a label column name), `permutations`. Returns one JSON document.
    pub fn analyze(
        &mut self,
        labels: Option<String>,
        options_json: &str,
        compare_tree: Option<String>,
        compare_pcs: Option<String>,
    ) -> Result<String, JsError> {
        let mut o = Options {
            ksize: self.params.ksize,
            scaled: self.params.scaled,
            max_bases: self.params.max_bases,
            r1_only: self.params.r1_only,
            quiet: true,
            ..Options::default()
        };
        apply(&mut o, options_json).map_err(|e| JsError::new(&e))?;
        if labels.is_some() {
            o.labels = Some(PathBuf::from("labels"));
        }
        if self.index.is_some() {
            o.index = Some(PathBuf::from("index"));
        }
        // the report only needs to know that a comparison was asked for
        if compare_tree.is_some() {
            o.compare_tree = Some(PathBuf::from("compare_tree"));
        }
        if compare_pcs.is_some() {
            o.compare_pcs = Some(PathBuf::from("compare_pcs"));
        }
        let res = analyze(
            &o,
            Analysis {
                sketches: std::mem::take(&mut self.sketches),
                errors: std::mem::take(&mut self.errors),
                warnings: std::mem::take(&mut self.warnings),
                index: self.index.as_ref(),
                labels_text: labels,
                compare_tree: compare_tree.map(Ok),
                compare_pcs: compare_pcs.map(Ok),
            },
        )
        .map_err(|e| JsError::new(&format!("{e:#}")))?;
        let ids: Vec<&str> = res
            .included
            .iter()
            .map(|&r| res.rows[r].id.as_str())
            .collect();
        let clusters = res.clusters.as_ref().map(|c| {
            ids.iter()
                .zip(c)
                .map(|(i, k)| json!([i, k]))
                .collect::<Vec<_>>()
        });
        let pca = res.ordination.as_ref().map(|x| {
            ids.iter()
                .zip(&x.coords)
                .map(|(i, v)| json!({"id": i, "pcs": v}))
                .collect::<Vec<_>>()
        });
        let qc: Vec<Value> = res
            .rows
            .iter()
            .map(|r| {
                json!({"id": r.id, "status": r.status.as_str(),
                       "organism_call": r.index.as_ref().and_then(|c| c.species.clone().or(c.genus.clone())),
                       "second_taxon": r.index.as_ref().and_then(|c| c.mixture.as_ref().map(|m| m.0.clone())),
                       "second_fraction": r.index.as_ref().and_then(|c| c.mixture.as_ref().map(|m| m.1)),
                       "flags": r.flags().iter().map(|f| json!({"code": f.code, "detail": f.detail})).collect::<Vec<_>>()})
            })
            .collect();
        let (md, figures) = popsketch::report::markdown(&o, &res);
        let markdown = json!({"md": md, "figures": figures});
        Ok(json!({
            "summary": summary_json(&o, &res),
            "html": popsketch::report::html(&o, &res),
            "markdown": markdown,
            "clusters": clusters,
            "pca": pca,
            "newick": res.newick,
            "qc": qc,
            "warnings": res.warnings,
        })
        .to_string())
    }
}

fn apply(o: &mut Options, json_text: &str) -> Result<(), String> {
    if json_text.trim().is_empty() {
        return Ok(());
    }
    let v: Value = serde_json::from_str(json_text).map_err(|e| format!("options: {e}"))?;
    let obj = v.as_object().ok_or("options must be a JSON object")?;
    for (k, val) in obj {
        let f = || val.as_f64().ok_or(format!("option {k} must be a number"));
        let u = || {
            val.as_u64()
                .ok_or(format!("option {k} must be a whole number"))
        };
        let b = || {
            val.as_bool()
                .ok_or(format!("option {k} must be true or false"))
        };
        match k.as_str() {
            "min_cov" => o.min_cov = f()?,
            "size_z" => o.size_z = f()?,
            "size_floor" => o.size_floor = f()?,
            "shoulder_z" => o.shoulder_z = f()?,
            "min_floors" => o.min_floors = f()?,
            "k" => o.k = Some(u()? as usize),
            "components" => o.components = u()? as usize,
            "floor_samples" => o.floor_samples = u()? as usize,
            "tree" => o.tree = b()?,
            "no_pca" => o.no_pca = b()?,
            "root" => o.root = val.as_str().map(str::to_string),
            "strata" => o.strata = val.as_str().map(str::to_string),
            "permutations" => o.permutations = u()? as usize,
            "collapse" => {
                o.collapse = match val.as_str() {
                    None | Some("") | Some("none") => None,
                    Some(r) => Some(
                        popsketch::collapse::Rule::parse(r)
                            .ok_or("collapse must be none, mean, max or median")?,
                    ),
                }
            }
            "collapse_within" => {
                o.collapse_within = val
                    .as_str()
                    .unwrap_or("")
                    .split(',')
                    .map(|s| s.trim().to_string())
                    .filter(|s| !s.is_empty())
                    .collect()
            }
            "flagged" => {
                o.flagged = match val.as_str() {
                    Some("exclude") => FlagMode::Exclude,
                    Some("mark") => FlagMode::Mark,
                    Some("decontaminate") => FlagMode::Decontaminate,
                    _ => return Err("flagged must be exclude, mark or decontaminate".into()),
                }
            }
            other => return Err(format!("unknown option {other}")),
        }
    }
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn collapse_options_are_parsed_and_bad_ones_refused() {
        let mut o = Options::default();
        apply(
            &mut o,
            r#"{"collapse":"max","collapse_within":" study, country ,"}"#,
        )
        .unwrap();
        assert_eq!(o.collapse, Some(popsketch::collapse::Rule::Max));
        assert_eq!(
            o.collapse_within,
            vec!["study".to_string(), "country".to_string()]
        );
        // an empty choice (the page's "no") turns it off, and the within list is cleared by an empty string
        apply(&mut o, r#"{"collapse":"","collapse_within":""}"#).unwrap();
        assert_eq!(o.collapse, None);
        assert!(o.collapse_within.is_empty());
        let err = apply(&mut o, r#"{"collapse":"fast"}"#).unwrap_err();
        assert!(err.contains("collapse must be"), "{err}");
    }
}
