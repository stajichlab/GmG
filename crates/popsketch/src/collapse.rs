//! Collapse near-identical samples to one representative each, so that clonal or outbreak sets do not swamp the ordination,
//! the tree and the label tests, and so that the number of rows stays small enough for the fast distance method.
//!
//! Two samples count as indistinguishable when their distance (1 - ANI) is at most the technical noise of the comparison,
//! which is read from the split-half distance of each sample (the distance between two sketches of alternate reads of that
//! sample, see `floor`). Leader clustering: samples are visited from the least noisy (best covered) to the noisiest, ties by
//! input order; an unassigned sample becomes a representative and takes every unassigned sample within the threshold of it.
//! No chaining: every member is within the threshold of its own representative. The result depends on the order of
//! visiting, so it is one valid clustering, not the smallest; the order is deterministic.

/// How the threshold of a pair is formed from the two samples' split-half distances `ni` and `nj`.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Rule {
    /// `(ni + nj) / 2`: the typical noise of the comparison (the default)
    Mean,
    /// `max(ni, nj)`: looser, the noisier sample sets the limit
    Max,
    /// the species-wide median split-half distance for every pair
    Median,
}

impl Rule {
    pub fn parse(s: &str) -> Option<Rule> {
        match s {
            "mean" => Some(Rule::Mean),
            "max" => Some(Rule::Max),
            "median" => Some(Rule::Median),
            _ => None,
        }
    }
    pub fn name(self) -> &'static str {
        match self {
            Rule::Mean => "mean",
            Rule::Max => "max",
            Rule::Median => "median",
        }
    }
}

#[derive(Clone, Debug)]
pub struct Collapse {
    /// for every sample, the index of its representative (a representative maps to itself)
    pub rep_of: Vec<usize>,
    /// representatives, in increasing index order
    pub reps: Vec<usize>,
    /// number of samples each representative stands for (itself included), aligned with `reps`
    pub sizes: Vec<usize>,
}

/// A deterministic pseudo-random key for sample `i` under `seed` (SplitMix64 finaliser).
fn shuffle_key(seed: u64, i: usize) -> u64 {
    let mut z = seed ^ (i as u64).wrapping_mul(0x9E37_79B9_7F4A_7C15);
    z = (z ^ (z >> 30)).wrapping_mul(0xBF58_476D_1CE4_E5B9);
    z = (z ^ (z >> 27)).wrapping_mul(0x94D0_49BB_1331_11EB);
    z ^ (z >> 31)
}

/// `d1` is the n x n row-major matrix of 1 - ANI; `noise[i]` is the split-half distance of sample i; `median` is the
/// species-wide value used by `Rule::Median`. `shuffle_seed` replaces the visiting order (least noisy first) by a
/// pseudo-random one, to test how much the result depends on the order.
pub fn collapse(
    d1: &[f64],
    n: usize,
    noise: &[f64],
    rule: Rule,
    median: f64,
    shuffle_seed: Option<u64>,
) -> Collapse {
    assert_eq!(d1.len(), n * n);
    assert_eq!(noise.len(), n);
    let mut order: Vec<usize> = (0..n).collect();
    match shuffle_seed {
        None => order.sort_by(|&a, &b| noise[a].total_cmp(&noise[b]).then(a.cmp(&b))),
        Some(s) => order.sort_by_key(|&i| (shuffle_key(s, i), i)),
    }
    let mut rep_of = vec![usize::MAX; n];
    let mut chosen: Vec<usize> = Vec::new();
    for &i in &order {
        if rep_of[i] != usize::MAX {
            continue;
        }
        rep_of[i] = i;
        chosen.push(i);
        for j in 0..n {
            if rep_of[j] != usize::MAX {
                continue;
            }
            let tau = match rule {
                Rule::Mean => (noise[i] + noise[j]) / 2.0,
                Rule::Max => noise[i].max(noise[j]),
                Rule::Median => median,
            };
            if d1[i * n + j] <= tau {
                rep_of[j] = i;
            }
        }
    }
    chosen.sort_unstable();
    let sizes = chosen
        .iter()
        .map(|&r| rep_of.iter().filter(|&&x| x == r).count())
        .collect();
    Collapse {
        rep_of,
        reps: chosen,
        sizes,
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    /// samples 0,1 are near-identical, 2 is far from everyone, 3 is a near copy of 2
    fn fixture() -> (Vec<f64>, usize) {
        let n = 4;
        let mut d = vec![0.2; n * n];
        for i in 0..n {
            d[i * n + i] = 0.0;
        }
        let mut set = |a: usize, b: usize, v: f64| {
            d[a * n + b] = v;
            d[b * n + a] = v;
        };
        set(0, 1, 0.0005);
        set(2, 3, 0.0008);
        (d, n)
    }

    #[test]
    fn near_identical_samples_share_a_representative_and_distinct_ones_do_not() {
        let (d, n) = fixture();
        let c = collapse(&d, n, &[0.001; 4], Rule::Mean, 0.001, None);
        assert_eq!(c.reps, vec![0, 2]);
        assert_eq!(c.sizes, vec![2, 2]);
        assert_eq!(c.rep_of, vec![0, 0, 2, 2]);
    }

    #[test]
    fn the_least_noisy_sample_is_the_representative() {
        let (d, n) = fixture();
        // sample 1 is the best covered of the pair (0, 1); sample 3 of the pair (2, 3)
        let c = collapse(
            &d,
            n,
            &[0.002, 0.0005, 0.002, 0.0006],
            Rule::Max,
            0.001,
            None,
        );
        assert_eq!(c.reps, vec![1, 3]);
        assert_eq!(c.rep_of, vec![1, 1, 3, 3]);
    }

    #[test]
    fn tighter_noise_keeps_more_representatives_and_no_chaining() {
        let (d, n) = fixture();
        // noise below both pair distances: nothing collapses
        let c = collapse(&d, n, &[0.0001; 4], Rule::Mean, 0.0001, None);
        assert_eq!(c.reps, vec![0, 1, 2, 3]);
        // a chain a-b-c with ab, bc within the threshold but ac not: c is not pulled in through b
        let m = 3;
        let mut e = vec![0.0; m * m];
        for (a, b, v) in [(0, 1, 0.001), (1, 2, 0.001), (0, 2, 0.002)] {
            e[a * m + b] = v;
            e[b * m + a] = v;
        }
        let c = collapse(&e, m, &[0.0, 0.001, 0.002], Rule::Max, 0.001, None);
        // order 0, 1, 2: 0 takes 1 (0.001 <= max(0, 0.001)); 2 is not within max(0, 0.002) = 0.002? it is (0.002 <= 0.002)
        assert_eq!(c.rep_of, vec![0, 0, 0]);
        let c = collapse(&e, m, &[0.0, 0.001, 0.0015], Rule::Max, 0.001, None);
        // now 0 takes 1 but not 2 (0.002 > 0.0015): 2 stays alone although 1-2 is within threshold
        assert_eq!(c.reps, vec![0, 2]);
    }

    #[test]
    fn the_median_rule_uses_one_threshold_for_every_pair() {
        let (d, n) = fixture();
        let c = collapse(&d, n, &[0.0; 4], Rule::Median, 0.001, None);
        assert_eq!(c.reps, vec![0, 2]);
        let c = collapse(&d, n, &[0.0; 4], Rule::Median, 0.0001, None);
        assert_eq!(c.reps.len(), 4);
    }

    #[test]
    fn deterministic_with_ties_and_sizes_add_up() {
        let (d, n) = fixture();
        let a = collapse(&d, n, &[0.001; 4], Rule::Mean, 0.001, None);
        let b = collapse(&d, n, &[0.001; 4], Rule::Mean, 0.001, None);
        assert_eq!(a.rep_of, b.rep_of);
        assert_eq!(a.sizes.iter().sum::<usize>(), n);
        assert!(Rule::parse("mean") == Some(Rule::Mean) && Rule::parse("x").is_none());
    }

    #[test]
    fn a_shuffled_order_is_deterministic_valid_and_can_differ_from_the_noise_order() {
        // 12 samples in a line 0.0004 apart: a leader takes its neighbours within 0.001, so the result depends on the order
        let n = 12;
        let mut d = vec![0.0; n * n];
        for a in 0..n {
            for b in 0..n {
                d[a * n + b] = (a as f64 - b as f64).abs() * 0.0004;
            }
        }
        let noise = vec![0.001; n];
        let base = collapse(&d, n, &noise, Rule::Mean, 0.001, None);
        let mut counts = std::collections::BTreeSet::new();
        for seed in 0..20u64 {
            let a = collapse(&d, n, &noise, Rule::Mean, 0.001, Some(seed));
            let b = collapse(&d, n, &noise, Rule::Mean, 0.001, Some(seed));
            assert_eq!(a.rep_of, b.rep_of, "same seed, same answer");
            assert_eq!(a.sizes.iter().sum::<usize>(), n);
            // every member is within the threshold of its representative
            for (j, &r) in a.rep_of.iter().enumerate() {
                assert!(d[r * n + j] <= 0.001 + 1e-12);
            }
            counts.insert(a.reps.len());
        }
        assert!(
            counts.len() > 1,
            "the order matters on a chain: {counts:?} vs {}",
            base.reps.len()
        );
    }
}
