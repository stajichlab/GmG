//! Turn raw search hits into per-rank taxon calls with confidence features.

use crate::index::{Hit, RefIndex, RANKS};

/// Per-rank acceptance thresholds. The defaults are placeholders; calibrate them with the
/// leave-one-out benchmark (`gmg bench`) and load the result with `Thresholds::from_json_like`.
#[derive(Clone, Debug)]
pub struct Thresholds {
    pub min_shared: [u32; 6],
    pub min_ani: [f64; 6],
    pub min_consistency: [f64; 6],
    /// Hits within `margin` x the top shared count vote on the call (at most `top_n` of them).
    pub margin: f64,
    pub top_n: usize,
    /// A genus whose evidence is below `min_ani[4]` but at least this ANI (and passes the other gates) is reported as
    /// `tentative_genus`, never as a call: in that band a named genus is right only ~80% of the time even when the genus is in the
    /// reference, and wrong whenever it is not. Set above 1.0 to switch the field off.
    pub tentative_genus_min_ani: f64,
}

impl Default for Thresholds {
    /// The `balanced` set: calibrated on BFD (23,669 genomes, scaled=10000) to ~95% precision per rank for the common case that
    /// the reference holds close relatives (species/genus from leave-one-out, family and above from leave-species-out), with the genus
    /// cutoff raised from ANI 0.7765 to 0.90 so that a genus absent from the reference is named for 18% of queries instead of 58%
    /// (`reports/validation/novelty_curve.md`, issue #57); the 0.7765-0.90 band is kept as `tentative_genus`.
    /// `data/thresholds/bfd_s10000.known.json`; the earlier permissive genus cutoff is `bfd_s10000.permissive.json`;
    /// for genomes that may belong to a genus absent from the reference use `bfd_s10000.novel.json`.
    fn default() -> Self {
        Self {
            min_shared: [5, 5, 5, 5, 5, 5],
            min_ani: [0.7109, 0.7109, 0.7109, 0.7715, 0.90, 0.9938],
            min_consistency: [0.5, 0.5, 0.6, 0.9, 0.5, 0.9],
            margin: 0.9,
            top_n: 10,
            tentative_genus_min_ani: 0.7765,
        }
    }
}

#[derive(Clone, Debug)]
pub struct RankVote {
    pub rank: &'static str,
    /// taxon with the most vote weight (None if no voting hit has this rank)
    pub taxon: Option<String>,
    /// fraction of the (known-rank) vote weight that agrees with `taxon`
    pub consistency: f64,
}

#[derive(Clone, Debug)]
pub struct TopHit {
    pub genome: u32,
    pub shared: u32,
    pub containment_query: f64,
    pub containment_ref: f64,
    /// containment ANI estimate from the query's containment in the hit, containment_query^(1/k).
    /// (Not max(cq, cr): a small reference fully contained in a large query would then score ~1.)
    pub ani: f64,
}

/// A genus the evidence points to but does not support as a call (see `Thresholds::tentative_genus_min_ani`).
#[derive(Clone, Debug)]
pub struct Tentative {
    pub taxon: String,
    pub ani: f64,
    pub consistency: f64,
}

#[derive(Clone, Debug)]
pub struct Guess {
    pub query_hashes: usize,
    pub top: Option<TopHit>,
    /// raw per-rank vote (before thresholds)
    pub votes: Vec<RankVote>,
    /// thresholded, hierarchical calls: Some(taxon) only if this and every higher rank passed
    pub calls: Vec<Option<String>>,
    pub hits: Vec<(u32, u32)>,
    /// set only when the genus was NOT called because its ANI is below `min_ani[4]`, but the vote, the higher ranks and
    /// `tentative_genus_min_ani` all hold
    pub tentative_genus: Option<Tentative>,
}

/// Why a rank was not called: the first rank whose vote exists but failed its threshold (or no hit at all).
#[derive(Clone, Debug, PartialEq)]
pub struct Withheld {
    /// index into `RANKS`
    pub rank: usize,
    /// `no_hits`, `too_few_shared_hashes`, `ani_below_cutoff` or `low_agreement`
    pub code: &'static str,
    /// plain text with the numbers, e.g. `ANI 0.8870 < 0.9000`
    pub detail: String,
}

impl Guess {
    /// The reason the first uncalled rank was withheld; None when a species was called (or no rank has a taxon to call).
    pub fn withheld(&self, th: &Thresholds) -> Option<Withheld> {
        let Some(top) = &self.top else {
            return Some(Withheld {
                rank: 0,
                code: "no_hits",
                detail: "no reference genome shares a sketch hash with the input".into(),
            });
        };
        for (r, v) in self.votes.iter().enumerate() {
            if v.taxon.is_none() || self.calls[r].is_some() {
                continue;
            }
            return Some(if top.shared < th.min_shared[r] {
                Withheld {
                    rank: r,
                    code: "too_few_shared_hashes",
                    detail: format!("{} shared hashes < {}", top.shared, th.min_shared[r]),
                }
            } else if top.ani < th.min_ani[r] {
                Withheld {
                    rank: r,
                    code: "ani_below_cutoff",
                    detail: format!("ANI {:.4} < {:.4}", top.ani, th.min_ani[r]),
                }
            } else {
                Withheld {
                    rank: r,
                    code: "low_agreement",
                    detail: format!(
                        "the near-best hits agree on {:.2} < {:.2}",
                        v.consistency, th.min_consistency[r]
                    ),
                }
            });
        }
        None
    }
}

pub fn guess(
    index: &RefIndex,
    query: &[u64],
    skip: &dyn Fn(u32) -> bool,
    th: &Thresholds,
) -> Guess {
    let hits = index.search(query, skip);
    guess_from_hits(index, query.len(), &hits, th)
}

pub fn guess_from_hits(index: &RefIndex, qlen: usize, hits: &[Hit], th: &Thresholds) -> Guess {
    let empty_votes = || {
        RANKS
            .iter()
            .map(|r| RankVote {
                rank: r,
                taxon: None,
                consistency: 0.0,
            })
            .collect::<Vec<_>>()
    };
    let Some(top) = hits.first() else {
        return Guess {
            query_hashes: qlen,
            top: None,
            votes: empty_votes(),
            calls: vec![None; 6],
            hits: vec![],
            tentative_genus: None,
        };
    };
    let g = &index.genomes[top.genome as usize];
    let cq = top.shared as f64 / qlen.max(1) as f64;
    let cr = top.shared as f64 / g.n_hashes.max(1) as f64;
    let ani = cq.powf(1.0 / index.ksize as f64);
    let topinfo = TopHit {
        genome: top.genome,
        shared: top.shared,
        containment_query: cq,
        containment_ref: cr,
        ani,
    };

    let cutoff = (top.shared as f64 * th.margin).ceil() as u32;
    let voters: Vec<&Hit> = hits
        .iter()
        .take(th.top_n)
        .filter(|h| h.shared >= cutoff)
        .collect();
    let mut votes = Vec::with_capacity(6);
    let mut calls = vec![None; 6];
    let mut open = true;
    let mut tentative = false;
    for r in 0..6 {
        // weight per interned taxon id (0 = unknown, ignored)
        let mut w: Vec<(u32, f64)> = Vec::new();
        let mut total = 0.0;
        for h in &voters {
            let id = index.lineage_ids[h.genome as usize][r];
            if id == 0 {
                continue;
            }
            let wt = h.shared as f64;
            total += wt;
            match w.iter_mut().find(|e| e.0 == id) {
                Some(e) => e.1 += wt,
                None => w.push((id, wt)),
            }
        }
        let best = w
            .iter()
            .cloned()
            .max_by(|a, b| a.1.partial_cmp(&b.1).unwrap().then(b.0.cmp(&a.0)));
        let (taxon, cons) = match best {
            Some((id, wt)) => (Some(index.taxa[r][id as usize].clone()), wt / total),
            None => (None, 0.0),
        };
        // A rank with no taxon on any voting hit (blank in the metadata) is skipped, not failed,
        // so it does not block calls at deeper ranks. Only a failed threshold closes the gate.
        if taxon.is_some() {
            let pass = top.shared >= th.min_shared[r]
                && ani >= th.min_ani[r]
                && cons >= th.min_consistency[r];
            if r == 4 {
                tentative = open
                    && !pass
                    && top.shared >= th.min_shared[4]
                    && cons >= th.min_consistency[4]
                    && ani >= th.tentative_genus_min_ani;
            }
            if open && pass {
                calls[r] = taxon.clone();
            } else {
                open = false;
            }
        }
        votes.push(RankVote {
            rank: RANKS[r],
            taxon,
            consistency: cons,
        });
    }
    let tentative_genus = if tentative {
        votes[4].taxon.clone().map(|taxon| Tentative {
            taxon,
            ani,
            consistency: votes[4].consistency,
        })
    } else {
        None
    };
    Guess {
        query_hashes: qlen,
        top: Some(topinfo),
        votes,
        calls,
        tentative_genus,
        hits: hits
            .iter()
            .take(th.top_n)
            .map(|h| (h.genome, h.shared))
            .collect(),
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::index::{GenomeMeta, IndexBuilder};

    fn m(id: &str, genus: &str, species: &str) -> GenomeMeta {
        let mut l: [String; 6] = Default::default();
        l[0] = "Ascomycota".into();
        l[4] = genus.into();
        l[5] = species.into();
        GenomeMeta {
            id: id.into(),
            taxid: 1,
            lineage: l,
            n_hashes: 0,
            cluster: String::new(),
            ..Default::default()
        }
    }

    #[test]
    fn hierarchical_calls() {
        let mut b = IndexBuilder::new(21, 1, 42);
        b.add(
            m("a1", "Aspergillus", "Aspergillus fumigatus"),
            (0..100).collect(),
        );
        b.add(
            m("a2", "Aspergillus", "Aspergillus flavus"),
            (50..150).collect(),
        );
        b.add(
            m("f1", "Fusarium", "Fusarium oxysporum"),
            (1000..1100).collect(),
        );
        let ix = b.build(4);
        let th = Thresholds {
            min_shared: [1; 6],
            ..Default::default()
        };
        // identical to a1
        let g = guess(&ix, &(0..100).collect::<Vec<_>>(), &|_| false, &th);
        assert_eq!(g.calls[5].as_deref(), Some("Aspergillus fumigatus"));
        assert_eq!(g.calls[0].as_deref(), Some("Ascomycota"));
        // leave a1 out: best remaining is a2 (shares 50..100) -> genus ok, species is a2's
        let g = guess(&ix, &(0..100).collect::<Vec<_>>(), &|i| i == 0, &th);
        assert_eq!(g.calls[4].as_deref(), Some("Aspergillus"));
        // nothing shared -> no calls
        let g = guess(&ix, &[5000, 5001], &|_| false, &th);
        assert!(g.calls.iter().all(|c| c.is_none()));
    }

    #[test]
    fn withheld_names_the_first_failed_rank_and_why() {
        let mut b = IndexBuilder::new(21, 1, 42);
        b.add(
            m("a1", "Aspergillus", "Aspergillus fumigatus"),
            (0..100).collect(),
        );
        let ix = b.build(1);
        let th = Thresholds {
            min_shared: [1; 6],
            ..Default::default()
        };
        // identical: everything called, nothing withheld
        let g = guess(&ix, &(0..100).collect::<Vec<_>>(), &|_| false, &th);
        assert_eq!(g.withheld(&th), None);
        // 60% of the hashes: containment 0.6 -> ANI 0.976, below the species cutoff
        let g = guess(
            &ix,
            &(0..60).chain(500..540).collect::<Vec<_>>(),
            &|_| false,
            &th,
        );
        let w = g.withheld(&th).unwrap();
        assert_eq!((w.rank, w.code), (5, "ani_below_cutoff"), "{w:?}");
        // nothing shared
        let g = guess(&ix, &[7000, 7001], &|_| false, &th);
        assert_eq!(g.withheld(&th).unwrap().code, "no_hits");
    }

    #[test]
    fn tentative_genus_in_the_band_between_the_cutoffs() {
        let mut b = IndexBuilder::new(21, 1, 42);
        b.add(
            m("a1", "Aspergillus", "Aspergillus fumigatus"),
            (0..1000).collect(),
        );
        let ix = b.build(2);
        let th = Thresholds::default();
        // 33 of 1000 query hashes found in the reference: containment 0.033, ANI = 0.033^(1/21) ~ 0.85
        let q: Vec<u64> = (0..33).chain(5000..5967).collect();
        let g = guess(&ix, &q, &|_| false, &th);
        let ani = g.top.as_ref().unwrap().ani;
        assert!(ani > 0.8 && ani < 0.9, "ani {ani}");
        assert!(
            g.calls[4].is_none(),
            "genus must not be called below the 0.90 cutoff"
        );
        assert_eq!(g.calls[0].as_deref(), Some("Ascomycota"));
        let t = g.tentative_genus.expect("tentative genus");
        assert_eq!(t.taxon, "Aspergillus");
        assert!((t.ani - ani).abs() < 1e-12);
        // identical query: genus is a real call, so no tentative field
        let g = guess(&ix, &(0..1000).collect::<Vec<_>>(), &|_| false, &th);
        assert_eq!(g.calls[4].as_deref(), Some("Aspergillus"));
        assert!(g.tentative_genus.is_none());
        let off = Thresholds {
            tentative_genus_min_ani: 2.0,
            ..Thresholds::default()
        };
        let q: Vec<u64> = (0..33).chain(5000..5967).collect();
        assert!(guess(&ix, &q, &|_| false, &off).tentative_genus.is_none());
        // a failed higher rank closes the gate: no tentative genus either
        let strict = Thresholds {
            min_ani: [0.99, 0.99, 0.99, 0.99, 0.90, 0.9938],
            ..Thresholds::default()
        };
        assert!(guess(&ix, &q, &|_| false, &strict)
            .tentative_genus
            .is_none());
    }

    #[test]
    fn min_shared_stops_deeper_ranks() {
        let mut b = IndexBuilder::new(21, 1, 42);
        b.add(
            m("a1", "Aspergillus", "Aspergillus fumigatus"),
            (0..100).collect(),
        );
        let ix = b.build(1);
        let th = Thresholds {
            min_shared: [1, 1, 1, 1, 1, 500],
            ..Default::default()
        };
        let g = guess(&ix, &(0..100).collect::<Vec<_>>(), &|_| false, &th);
        assert_eq!(g.calls[4].as_deref(), Some("Aspergillus"));
        assert!(g.calls[5].is_none());
    }
}
