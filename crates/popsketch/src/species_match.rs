//! Do the reads of a run belong to the species it is filed under? (PopSeq stajichlab/PopSeq#52, found with SRR30534096: filed
//! *Cryptococcus neoformans*, reads of *Aspergillus fumigatus*.)
//!
//! One verdict per sample from the reference index's own call (the same classifier as `gmg guess`), compared with the filed name:
//!
//! * `match`: the called species is the filed species;
//! * `genus_match`: the called genus is the filed genus but the species is not called or differs (kept, reported);
//! * `mismatch`: the index makes a genus call that differs from the filed genus, and the best hit holds at least `min_share` of
//!   the sample's hashes (a weak call is not enough to refuse a run);
//! * `no_fungal_call_suspect`: no confident call, and the index has references for the filed genus, so it should have found them;
//! * `no_fungal_call`: no confident call, and the index has no references for the filed genus: it cannot judge.
//!
//! `unchecked` (no sketch yet) is written by the consumer, not here. `purity` is the share of the sample's hashes found in the best
//! reference of the filed genus (not used for the status). The default `min_share` is a placeholder until the post-hoc pass over the
//! sketched runs has been read (the order of work in #52).

use serde_json::Value;
use std::collections::HashSet;

#[derive(Clone, Debug)]
pub struct MatchParams {
    /// the best hit must hold at least this share of the sample's hashes for a genus call to count as `mismatch`
    pub min_share: f64,
}

impl Default for MatchParams {
    fn default() -> Self {
        MatchParams { min_share: 0.02 }
    }
}

#[derive(Clone, Debug, PartialEq)]
pub struct MatchRow {
    pub id: String,
    pub filed: String,
    pub status: &'static str,
    pub detected_genus: String,
    pub detected_species: String,
    pub best_shared: u64,
    pub best_share: f64,
    pub purity: Option<f64>,
    pub reason: String,
}

pub const TSV_HEADER: &str =
    "id\tfiled_name\tstatus\tdetected_genus\tdetected_species\tbest_shared_hashes\tbest_share\tpurity\treason";

impl MatchRow {
    pub fn tsv(&self) -> String {
        let f = |x: Option<f64>| x.map_or(String::new(), |v| format!("{v:.4}"));
        let c = |x: &str| x.replace(['\t', '\n', '\r'], " ");
        format!(
            "{}\t{}\t{}\t{}\t{}\t{}\t{:.4}\t{}\t{}",
            c(&self.id),
            c(&self.filed),
            self.status,
            c(&self.detected_genus),
            c(&self.detected_species),
            self.best_shared,
            self.best_share,
            f(self.purity),
            c(&self.reason)
        )
    }
}

/// First word of a species name: the genus.
pub fn genus_of(name: &str) -> &str {
    name.split_whitespace().next().unwrap_or("")
}

/// The verdict for one sample. `rec` is the classifier's record (`guess.calls`, `guess.hits` with a `lineage` and `shared` each),
/// `n_query` the number of the sample's hashes that were classified, `index_genera` the genera the index holds references for.
pub fn judge(
    id: &str,
    filed: &str,
    rec: &Value,
    n_query: usize,
    index_genera: &HashSet<String>,
    p: &MatchParams,
) -> MatchRow {
    let filed_genus = genus_of(filed);
    let call = |rank: &str| {
        rec["guess"]["calls"][rank]["taxon"]
            .as_str()
            .unwrap_or("")
            .to_string()
    };
    let (genus, species) = (call("genus"), call("species"));
    let hits = rec["guess"]["hits"].as_array().cloned().unwrap_or_default();
    let shared = |h: &Value| h["shared"].as_u64().unwrap_or(0);
    let lineage_genus = |h: &Value| h["lineage"][4].as_str().unwrap_or("").to_string();
    let best_shared = hits.first().map_or(0, shared);
    let share = |n: u64| {
        if n_query == 0 {
            0.0
        } else {
            n as f64 / n_query as f64
        }
    };
    let purity = (!filed_genus.is_empty()).then(|| {
        share(
            hits.iter()
                .filter(|h| lineage_genus(h) == filed_genus)
                .map(shared)
                .max()
                .unwrap_or(0),
        )
    });
    let (status, reason) = if genus.is_empty() {
        if index_genera.contains(filed_genus) {
            (
                "no_fungal_call_suspect",
                format!("no confident call, but the index has references for {filed_genus}"),
            )
        } else {
            (
                "no_fungal_call",
                format!("no confident call, and the index has no references for {filed_genus}: cannot judge"),
            )
        }
    } else if genus == filed_genus {
        if species == filed {
            ("match", String::new())
        } else if species.is_empty() {
            ("genus_match", format!("genus {genus}; species not called"))
        } else {
            (
                "genus_match",
                format!("genus {genus}; species called {species}, filed {filed}"),
            )
        }
    } else if share(best_shared) >= p.min_share {
        (
            "mismatch",
            format!(
                "reads are {genus}{}, filed {filed}; best hit holds {:.1}% of the hashes",
                if species.is_empty() {
                    String::new()
                } else {
                    format!(" ({species})")
                },
                100.0 * share(best_shared)
            ),
        )
    } else {
        (
            "no_fungal_call",
            format!(
                "genus call {genus} is weak (best hit {:.2}% of the hashes, below {:.2}%): not enough to refuse the run",
                100.0 * share(best_shared),
                100.0 * p.min_share
            ),
        )
    };
    MatchRow {
        id: id.into(),
        filed: filed.into(),
        status,
        detected_genus: genus,
        detected_species: species,
        best_shared,
        best_share: share(best_shared),
        purity,
        reason,
    }
}

/// Verdicts for every sample of a sketch folder whose id is in `filed` (id -> filed name), against the index in `index_path`.
/// Samples of the folder with no filed name are skipped (and counted in the second value). Reads only the sketches' own
/// cached files, so a run that is not sketched yet is simply absent (`unchecked` for the consumer).
#[cfg(feature = "native")]
pub fn run_folder(
    index_path: &std::path::Path,
    sketches: &std::path::Path,
    filed: &std::collections::HashMap<String, String>,
    p: &MatchParams,
) -> anyhow::Result<(Vec<MatchRow>, usize)> {
    use crate::sketch_io::cache_read_checked;
    use anyhow::Context;
    use gmg_core::classify::Context as RefContext;
    use gmg_core::{max_hash_for_scaled, RefIndex, Thresholds};
    let bytes = std::fs::read(index_path).with_context(|| index_path.display().to_string())?;
    let ix = RefIndex::read_owned(bytes)
        .with_context(|| format!("reading index {}", index_path.display()))?;
    let genera: HashSet<String> = ix.taxa[4]
        .iter()
        .filter(|g| !g.is_empty())
        .cloned()
        .collect();
    let th = Thresholds::default();
    let ctx = RefContext::new(&ix);
    let mh = max_hash_for_scaled(ix.scaled);
    let mut ids: Vec<&String> = filed.keys().collect();
    ids.sort();
    let (mut rows, mut missing) = (Vec::new(), 0usize);
    for id in ids {
        let Some(s) = cache_read_checked(sketches, id, false)? else {
            missing += 1;
            continue;
        };
        let reads = s.kind == "reads";
        let n_query = s
            .spectrum
            .hashes
            .iter()
            .zip(&s.spectrum.abunds)
            .filter(|(&h, &a)| h <= mh && a >= if reads { 2 } else { 1 })
            .count();
        let rec =
            crate::pipeline::screen_record(&ix, &th, &ctx, &s, 20).unwrap_or(serde_json::json!({}));
        rows.push(judge(id, &filed[id], &rec, n_query, &genera, p));
    }
    Ok((rows, missing))
}

#[cfg(test)]
mod tests {
    use super::*;
    use serde_json::json;

    fn rec(genus: &str, species: &str, hits: &[(&str, u64)]) -> Value {
        let hits: Vec<Value> = hits
            .iter()
            .map(|(g, n)| json!({"shared": n, "lineage": ["Asc", "Cl", "Or", "Fam", g, "sp"]}))
            .collect();
        json!({"guess": {"calls": {"genus": {"taxon": (!genus.is_empty()).then_some(genus)},
            "species": {"taxon": (!species.is_empty()).then_some(species)}}, "hits": hits}})
    }
    fn genera() -> HashSet<String> {
        ["Cryptococcus", "Aspergillus"]
            .iter()
            .map(|s| s.to_string())
            .collect()
    }
    fn judge1(filed: &str, r: &Value) -> MatchRow {
        judge("x", filed, r, 30_000, &genera(), &MatchParams::default())
    }

    /// SRR30534096: filed Cryptococcus neoformans, reads Aspergillus fumigatus with 25.8k of about 30k hashes shared.
    #[test]
    fn another_genus_with_strong_evidence_is_a_mismatch() {
        let r = judge1(
            "Cryptococcus neoformans",
            &rec(
                "Aspergillus",
                "Aspergillus fumigatus",
                &[("Aspergillus", 25_821)],
            ),
        );
        assert_eq!(r.status, "mismatch");
        assert!(
            r.reason.contains("Aspergillus") && r.reason.contains("filed Cryptococcus neoformans")
        );
        assert_eq!(r.best_shared, 25_821);
        assert_eq!(r.purity, Some(0.0));
    }

    #[test]
    fn species_agreement_is_a_match_and_other_species_of_the_genus_a_genus_match() {
        let same = judge1(
            "Cryptococcus neoformans",
            &rec(
                "Cryptococcus",
                "Cryptococcus neoformans",
                &[("Cryptococcus", 20_000)],
            ),
        );
        assert_eq!(same.status, "match");
        assert!((same.purity.unwrap() - 20_000.0 / 30_000.0).abs() < 1e-9);
        let other = judge1(
            "Cryptococcus neoformans",
            &rec(
                "Cryptococcus",
                "Cryptococcus gattii",
                &[("Cryptococcus", 20_000)],
            ),
        );
        assert_eq!(other.status, "genus_match");
        let uncalled = judge1(
            "Cryptococcus neoformans",
            &rec("Cryptococcus", "", &[("Cryptococcus", 20_000)]),
        );
        assert_eq!(uncalled.status, "genus_match");
    }

    #[test]
    fn a_weak_other_genus_call_does_not_refuse_the_run() {
        // 100 of 30,000 hashes (0.3%) is below the 2% placeholder
        let r = judge1(
            "Cryptococcus neoformans",
            &rec("Aspergillus", "", &[("Aspergillus", 100)]),
        );
        assert_eq!(r.status, "no_fungal_call");
        assert!(r.reason.contains("weak"));
        // and the cutoff is a parameter
        let strict = judge(
            "x",
            "Cryptococcus neoformans",
            &rec("Aspergillus", "", &[("Aspergillus", 100)]),
            30_000,
            &genera(),
            &MatchParams { min_share: 0.001 },
        );
        assert_eq!(strict.status, "mismatch");
    }

    #[test]
    fn no_call_is_suspect_only_when_the_index_has_the_filed_genus() {
        let none = rec("", "", &[("Aspergillus", 4)]);
        assert_eq!(
            judge1("Zymoseptoria tritici", &none).status,
            "no_fungal_call"
        );
        assert_eq!(
            judge1("Cryptococcus neoformans", &none).status,
            "no_fungal_call_suspect"
        );
    }

    #[test]
    fn an_empty_record_or_no_hashes_does_not_divide_by_zero() {
        let r = judge(
            "x",
            "Cryptococcus neoformans",
            &json!({}),
            0,
            &genera(),
            &MatchParams::default(),
        );
        assert_eq!(r.status, "no_fungal_call_suspect");
        assert_eq!(r.best_share, 0.0);
    }

    #[test]
    fn tsv_row_has_as_many_cells_as_the_header() {
        let r = judge1(
            "Cryptococcus neoformans",
            &rec("Aspergillus", "x\ty", &[("Aspergillus", 25_000)]),
        );
        assert_eq!(r.tsv().split('\t').count(), TSV_HEADER.split('\t').count());
    }
}
