//! Typed view of one `results.jsonl` record (written by `gmg run`).

use serde_json::Value;

pub const RANKS: [&str; 6] = ["phylum", "class", "order", "family", "genus", "species"];

#[derive(Clone, Debug, Default)]
pub struct Call {
    pub taxon: Option<String>,
    pub vote: Option<String>,
    pub consistency: f64,
}

#[derive(Clone, Debug, Default)]
pub struct Top {
    pub id: String,
    pub cluster: String,
    pub shared: u64,
    pub ani: f64,
    pub containment_query: f64,
    pub containment_ref: f64,
}

#[derive(Clone, Debug, Default)]
pub struct Hit {
    pub id: String,
    pub taxid: u64,
    pub shared: u64,
    pub lineage: [String; 6],
    pub cluster: String,
}

/// How many times each sketch hash of a read set was seen (1..=12, then 13+).
#[derive(Clone, Debug, Default)]
pub struct Hist {
    pub counts: Vec<u64>,
    pub over12: u64,
    pub min_abund: u32,
    pub median_kept: u32,
    /// estimated genome coverage of the sketched reads (`qc.coverage_est`) and how it was found (`peak` | `poisson`)
    pub coverage_est: Option<f64>,
    pub coverage_basis: String,
}

#[derive(Clone, Debug, Default)]
pub struct Rec {
    pub id: String,
    /// file name only (never a directory: reports must not leak paths)
    pub file: String,
    pub kind: String,
    pub status: String,
    pub best_rank: Option<usize>,
    pub message: String,
    pub error_code: String,
    pub error_action: String,
    pub total_bases: f64,
    pub n50: f64,
    pub n_hashes: u64,
    /// expected genome size (low, high in Mb, basis) of the called phylum, from the reference
    pub size_range: Option<(f64, f64, String)>,
    pub flags: Vec<String>,
    pub hist: Option<Hist>,
    pub calls: [Call; 6],
    /// (taxon, ani, consistency) of a genus the evidence points to but does not support as a call
    pub tentative_genus: Option<(String, f64, f64)>,
    pub top: Option<Top>,
    pub hits: Vec<Hit>,
}

fn s(v: &Value) -> String {
    v.as_str().unwrap_or("").to_string()
}

/// The last path component, for display.
pub fn file_name(path: &str) -> String {
    path.rsplit(['/', '\\']).next().unwrap_or(path).to_string()
}

impl Rec {
    pub fn from_value(v: &Value) -> Rec {
        let g = &v["guess"];
        let mut calls: [Call; 6] = Default::default();
        for (i, k) in RANKS.iter().enumerate() {
            let c = &g["calls"][*k];
            calls[i] = Call {
                taxon: c["taxon"].as_str().map(String::from),
                vote: c["vote"].as_str().map(String::from),
                consistency: c["consistency"].as_f64().unwrap_or(0.0),
            };
        }
        let tentative_genus = g["tentative_genus"].as_object().and_then(|_| {
            Some((
                g["tentative_genus"]["taxon"].as_str()?.to_string(),
                g["tentative_genus"]["ani"].as_f64().unwrap_or(0.0),
                g["tentative_genus"]["consistency"].as_f64().unwrap_or(0.0),
            ))
        });
        let top = g["top"].as_object().map(|_| Top {
            id: s(&g["top"]["id"]),
            cluster: s(&g["top"]["cluster"]),
            shared: g["top"]["shared"].as_u64().unwrap_or(0),
            ani: g["top"]["ani"].as_f64().unwrap_or(0.0),
            containment_query: g["top"]["containment_query"].as_f64().unwrap_or(0.0),
            containment_ref: g["top"]["containment_ref"].as_f64().unwrap_or(0.0),
        });
        let hits = g["hits"]
            .as_array()
            .map(|a| {
                a.iter()
                    .map(|h| {
                        let mut lineage: [String; 6] = Default::default();
                        for (i, l) in lineage.iter_mut().enumerate() {
                            *l = s(&h["lineage"][i]);
                        }
                        Hit {
                            id: s(&h["id"]),
                            taxid: h["taxid"].as_u64().unwrap_or(0),
                            shared: h["shared"].as_u64().unwrap_or(0),
                            lineage,
                            cluster: s(&h["cluster"]),
                        }
                    })
                    .collect()
            })
            .unwrap_or_default();
        let hist = v["qc"]["abundance_hist"].as_object().map(|_| Hist {
            counts: v["qc"]["abundance_hist"]["counts"]
                .as_array()
                .map(|a| a.iter().map(|x| x.as_u64().unwrap_or(0)).collect())
                .unwrap_or_default(),
            over12: v["qc"]["abundance_hist"]["over12"].as_u64().unwrap_or(0),
            min_abund: v["qc"]["abundance_hist"]["min_abund"].as_u64().unwrap_or(1) as u32,
            median_kept: v["qc"]["abundance_hist"]["median_kept"]
                .as_u64()
                .unwrap_or(0) as u32,
            coverage_est: v["qc"]["coverage_est"].as_f64(),
            coverage_basis: s(&v["qc"]["coverage_basis"]),
        });
        Rec {
            id: s(&v["id"]),
            file: file_name(&s(&v["path"])),
            kind: s(&v["kind"]),
            status: s(&v["status"]),
            best_rank: v["best_rank"]
                .as_str()
                .and_then(|b| RANKS.iter().position(|r| *r == b)),
            message: file_names_only(&s(&v["message"])),
            error_code: s(&v["error_code"]),
            error_action: s(&v["error_action"]),
            total_bases: v["qc"]["total_bases"].as_f64().unwrap_or(0.0),
            n50: v["qc"]["n50"].as_f64().unwrap_or(0.0),
            n_hashes: v["qc"]["n_hashes"].as_u64().unwrap_or(0),
            size_range: match (
                v["qc"]["size_range_mb"][0].as_f64(),
                v["qc"]["size_range_mb"][1].as_f64(),
            ) {
                (Some(lo), Some(hi)) => Some((lo, hi, s(&v["qc"]["size_range_basis"]))),
                _ => None,
            },
            flags: v["qc_flags"]
                .as_array()
                .map(|a| {
                    a.iter()
                        .filter_map(|x| x.as_str().map(String::from))
                        .collect()
                })
                .unwrap_or_default(),
            hist,
            calls,
            tentative_genus,
            top,
            hits,
        }
    }

    pub fn call(&self, rank: usize) -> &str {
        self.calls[rank].taxon.as_deref().unwrap_or("")
    }

    /// Deepest taxon called (species, else genus, ...), or "".
    pub fn call_name(&self) -> &str {
        (0..6)
            .rev()
            .map(|i| self.call(i))
            .find(|c| !c.is_empty())
            .unwrap_or("")
    }

    pub fn phylum(&self) -> &str {
        self.call(0)
    }

    pub fn ani(&self) -> Option<f64> {
        self.top.as_ref().map(|t| t.ani)
    }
}

/// Replace any absolute path in an error message by its file name (messages come from parsers and contain paths).
pub fn file_names_only(msg: &str) -> String {
    msg.split(' ')
        .map(|w| {
            let trimmed = w.trim_end_matches([':', ',', ';']);
            if trimmed.contains('/')
                && (trimmed.starts_with('/')
                    || trimmed.starts_with("./")
                    || trimmed.starts_with("../")
                    || trimmed.contains("//"))
            {
                let tail = &w[trimmed.len()..];
                format!("{}{}", file_name(trimmed), tail)
            } else {
                w.to_string()
            }
        })
        .collect::<Vec<_>>()
        .join(" ")
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn strips_paths_from_messages() {
        let m = file_names_only("opening /Users/x/data/genomes/empty.fa: Failed to read (line 0)");
        assert_eq!(m, "opening empty.fa: Failed to read (line 0)");
        assert!(!file_names_only("no path here").contains('/'));
    }

    #[test]
    fn parses_a_record() {
        let v: Value = serde_json::json!({
            "id": "g1", "path": "/a/b/g1.fa", "kind": "assembly", "status": "ok", "best_rank": "genus",
            "qc": {"total_bases": 1e6, "n50": 5000, "n_hashes": 123}, "qc_flags": ["fragmented"],
            "guess": {"calls": {"phylum": {"taxon": "Ascomycota", "vote": "Ascomycota", "consistency": 1.0},
                                  "genus": {"taxon": "Aspergillus", "vote": "Aspergillus", "consistency": 0.9}},
                      "top": {"id": "r1", "cluster": "c1", "shared": 10, "ani": 0.99, "containment_query": 0.8, "containment_ref": 0.7},
                      "hits": [{"id": "r1", "taxid": 5, "shared": 10, "lineage": ["P","C","O","F","G","S"], "cluster": "c1"}]}
        });
        let r = Rec::from_value(&v);
        assert_eq!(r.file, "g1.fa");
        assert_eq!(r.best_rank, Some(4));
        assert_eq!(r.call_name(), "Aspergillus");
        assert_eq!(r.phylum(), "Ascomycota");
        assert_eq!(r.hits[0].lineage[5], "S");
        assert_eq!(r.flags, vec!["fragmented"]);
    }
}
