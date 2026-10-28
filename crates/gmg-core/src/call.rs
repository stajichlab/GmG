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
}

impl Default for Thresholds {
    fn default() -> Self {
        Self {
            min_shared: [5, 5, 5, 8, 15, 30],
            min_ani: [0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
            min_consistency: [0.6, 0.6, 0.6, 0.6, 0.6, 0.6],
            margin: 0.9,
            top_n: 10,
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
    /// containment ANI estimate, max(containment)^(1/k)
    pub ani: f64,
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
}

pub fn guess(index: &RefIndex, query: &[u64], skip: &dyn Fn(u32) -> bool, th: &Thresholds) -> Guess {
    let hits = index.search(query, skip);
    guess_from_hits(index, query.len(), &hits, th)
}

pub fn guess_from_hits(index: &RefIndex, qlen: usize, hits: &[Hit], th: &Thresholds) -> Guess {
    let empty_votes = || {
        RANKS.iter().map(|r| RankVote { rank: r, taxon: None, consistency: 0.0 }).collect::<Vec<_>>()
    };
    let Some(top) = hits.first() else {
        return Guess { query_hashes: qlen, top: None, votes: empty_votes(), calls: vec![None; 6], hits: vec![] };
    };
    let g = &index.genomes[top.genome as usize];
    let cq = top.shared as f64 / qlen.max(1) as f64;
    let cr = top.shared as f64 / g.n_hashes.max(1) as f64;
    let ani = cq.max(cr).powf(1.0 / index.ksize as f64);
    let topinfo = TopHit { genome: top.genome, shared: top.shared, containment_query: cq, containment_ref: cr, ani };

    let cutoff = (top.shared as f64 * th.margin).ceil() as u32;
    let voters: Vec<&Hit> = hits.iter().take(th.top_n).filter(|h| h.shared >= cutoff).collect();
    let mut votes = Vec::with_capacity(6);
    let mut calls = vec![None; 6];
    let mut open = true;
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
        let best = w.iter().cloned().max_by(|a, b| a.1.partial_cmp(&b.1).unwrap().then(b.0.cmp(&a.0)));
        let (taxon, cons) = match best {
            Some((id, wt)) => (Some(index.taxa[r][id as usize].clone()), wt / total),
            None => (None, 0.0),
        };
        // A rank with no taxon on any voting hit (blank in the metadata) is skipped, not failed,
        // so it does not block calls at deeper ranks. Only a failed threshold closes the gate.
        if taxon.is_some() {
            let pass = top.shared >= th.min_shared[r] && ani >= th.min_ani[r] && cons >= th.min_consistency[r];
            if open && pass {
                calls[r] = taxon.clone();
            } else {
                open = false;
            }
        }
        votes.push(RankVote { rank: RANKS[r], taxon, consistency: cons });
    }
    Guess {
        query_hashes: qlen,
        top: Some(topinfo),
        votes,
        calls,
        hits: hits.iter().take(th.top_n).map(|h| (h.genome, h.shared)).collect(),
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
        GenomeMeta { id: id.into(), taxid: 1, lineage: l, n_hashes: 0 }
    }

    #[test]
    fn hierarchical_calls() {
        let mut b = IndexBuilder::new(21, 1, 42);
        b.add(m("a1", "Aspergillus", "Aspergillus fumigatus"), (0..100).collect());
        b.add(m("a2", "Aspergillus", "Aspergillus flavus"), (50..150).collect());
        b.add(m("f1", "Fusarium", "Fusarium oxysporum"), (1000..1100).collect());
        let ix = b.build(4);
        let th = Thresholds { min_shared: [1; 6], ..Default::default() };
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
    fn min_shared_stops_deeper_ranks() {
        let mut b = IndexBuilder::new(21, 1, 42);
        b.add(m("a1", "Aspergillus", "Aspergillus fumigatus"), (0..100).collect());
        let ix = b.build(1);
        let th = Thresholds { min_shared: [1, 1, 1, 1, 1, 500], ..Default::default() };
        let g = guess(&ix, &(0..100).collect::<Vec<_>>(), &|_| false, &th);
        assert_eq!(g.calls[4].as_deref(), Some("Aspergillus"));
        assert!(g.calls[5].is_none());
    }
}
