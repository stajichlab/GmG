//! UPGMA on 1 - ANI, and the largest-gap cut.
//!
//! The cut is made where the merge heights jump most. It is called *clear* only if that jump is at least `min_floors` times
//! the technical noise floor (the split-half distance, in 1 - ANI units); otherwise the top candidate cuts are reported
//! instead of a forced answer. A cohort whose last merge is already below the floor is one population (k = 1).
//!
//! The ratio gap1 / gap2 is reported with its probability under a structureless population (a Kingman coalescent, for which
//! the gaps are independent waiting times), but it is not the gate: that ratio is at least 2 in about 60% of structureless
//! populations, so it cannot tell structure from the natural asymmetry of a coalescent tree (analysis/popgen_cocci/26_gap_null.py).

#[derive(Clone, Debug, PartialEq)]
pub struct Merge {
    /// cluster ids: leaves are 0..n, merge i creates id n + i
    pub a: usize,
    pub b: usize,
    /// average distance between the two clusters
    pub height: f64,
    pub size: usize,
}

/// Average-linkage clustering of a symmetric n x n distance matrix (row-major).
pub fn upgma(d: &[f64], n: usize) -> Vec<Merge> {
    let mut dist = d.to_vec();
    let mut active = vec![true; n];
    let mut id: Vec<usize> = (0..n).collect();
    let mut size = vec![1usize; n];
    let scan = |dist: &[f64], active: &[bool], i: usize| -> (usize, f64) {
        let mut best = (usize::MAX, f64::INFINITY);
        for k in 0..n {
            if k != i && active[k] && dist[i * n + k] < best.1 {
                best = (k, dist[i * n + k]);
            }
        }
        best
    };
    let mut nn: Vec<(usize, f64)> = (0..n).map(|i| scan(&dist, &active, i)).collect();
    let mut merges = Vec::with_capacity(n.saturating_sub(1));
    for step in 0..n.saturating_sub(1) {
        let (mut i, mut best) = (usize::MAX, f64::INFINITY);
        for k in 0..n {
            if active[k] && nn[k].1 < best {
                (i, best) = (k, nn[k].1);
            }
        }
        let j = nn[i].0;
        let (i, j) = (i.min(j), i.max(j));
        let (si, sj) = (size[i] as f64, size[j] as f64);
        merges.push(Merge {
            a: id[i].min(id[j]),
            b: id[i].max(id[j]),
            height: best,
            size: size[i] + size[j],
        });
        for k in 0..n {
            if active[k] && k != i && k != j {
                let v = (si * dist[i * n + k] + sj * dist[j * n + k]) / (si + sj);
                dist[i * n + k] = v;
                dist[k * n + i] = v;
            }
        }
        active[j] = false;
        size[i] += size[j];
        id[i] = n + step;
        for k in 0..n {
            if !active[k] || k == i {
                continue;
            }
            if nn[k].0 == i || nn[k].0 == j {
                nn[k] = scan(&dist, &active, k);
            } else if dist[k * n + i] < nn[k].1 {
                nn[k] = (i, dist[k * n + i]);
            }
        }
        nn[i] = scan(&dist, &active, i);
    }
    merges
}

/// Cluster labels (0..k, in order of each cluster's first member) after cutting the tree into `k` clusters.
pub fn cut_k(merges: &[Merge], n: usize, k: usize) -> Vec<usize> {
    let k = k.clamp(1, n.max(1));
    let mut parent: Vec<usize> = (0..n + merges.len()).collect();
    fn find(p: &mut [usize], x: usize) -> usize {
        let mut r = x;
        while p[r] != r {
            r = p[r];
        }
        let mut c = x;
        while p[c] != r {
            let nx = p[c];
            p[c] = r;
            c = nx;
        }
        r
    }
    for (m, mg) in merges.iter().enumerate().take(n - k) {
        let (ra, rb) = (find(&mut parent, mg.a), find(&mut parent, mg.b));
        parent[ra] = n + m;
        parent[rb] = n + m;
    }
    let mut label_of: Vec<Option<usize>> = vec![None; n + merges.len()];
    let mut next = 0;
    (0..n)
        .map(|i| {
            let r = find(&mut parent, i);
            *label_of[r].get_or_insert_with(|| {
                next += 1;
                next - 1
            })
        })
        .collect()
}

#[derive(Clone, Debug, PartialEq)]
pub struct Candidate {
    /// number of clusters this cut gives
    pub k: usize,
    /// jump in merge height at the cut
    pub gap: f64,
}

#[derive(Clone, Debug, PartialEq)]
pub enum Verdict {
    /// everything merges below the noise floor
    Single,
    /// the largest gap beats the second by at least the ratio
    Clear(usize),
    /// no single cut stands out
    NoClearCut,
    /// too few samples to say
    TooFew,
}

#[derive(Clone, Debug)]
pub struct GapCall {
    pub verdict: Verdict,
    /// gap1 / gap2 (infinity when the second gap is zero)
    pub ratio: f64,
    /// probability of a ratio at least this large in a structureless (coalescent) population of this size
    pub null_p: f64,
    /// the largest gap in units of the noise floor
    pub gap_over_floor: f64,
    /// the cuts with the largest gaps, best first (at most 3)
    pub candidates: Vec<Candidate>,
    pub top_height: f64,
}

/// Probability that gap1 / gap2 is at least `observed` when there is no structure: the gap at the cut that leaves k clusters
/// is the coalescent waiting time T_k ~ Exp(k (k - 1) / 2); take the two largest of T_2 ..= T_kmax. Monte Carlo with a fixed
/// seed, so the answer is reproducible.
pub fn null_ratio_p(n: usize, kmax: usize, observed: f64) -> f64 {
    let kmax = kmax.min(n);
    if kmax < 3 || !observed.is_finite() {
        return if observed.is_finite() { 1.0 } else { 0.0 };
    }
    let mut state: u64 = 0x1234_5678_9ABC_DEF1;
    let mut next = || {
        state = state.wrapping_add(0x9E37_79B9_7F4A_7C15);
        let mut z = state;
        z = (z ^ (z >> 30)).wrapping_mul(0xBF58_476D_1CE4_E5B9);
        z = (z ^ (z >> 27)).wrapping_mul(0x94D0_49BB_1331_11EB);
        ((z ^ (z >> 31)) >> 11) as f64 / (1u64 << 53) as f64
    };
    let reps = 20_000;
    let mut hits = 0usize;
    for _ in 0..reps {
        let (mut first, mut second) = (0.0f64, 0.0f64);
        for k in 2..=kmax {
            let rate = (k * (k - 1)) as f64 / 2.0;
            let t = -(1.0 - next()).ln() / rate;
            if t > first {
                second = first;
                first = t;
            } else if t > second {
                second = t;
            }
        }
        if second > 0.0 && first / second >= observed {
            hits += 1;
        }
    }
    hits as f64 / reps as f64
}

/// Apply the largest-gap rule. `kmax` bounds the number of clusters considered; `floor` is the noise floor in distance units
/// and `min_floors` how many floors the jump must be to be called clear.
pub fn gap_rule(merges: &[Merge], n: usize, min_floors: f64, floor: f64, kmax: usize) -> GapCall {
    let none = |v: Verdict, top: f64, cands: Vec<Candidate>| GapCall {
        verdict: v,
        ratio: f64::NAN,
        null_p: f64::NAN,
        gap_over_floor: f64::NAN,
        candidates: cands,
        top_height: top,
    };
    if n < 4 || merges.len() + 1 != n {
        return none(Verdict::TooFew, f64::NAN, vec![]);
    }
    let h: Vec<f64> = merges.iter().map(|m| m.height).collect();
    let top = *h.last().unwrap();
    let mut cands: Vec<Candidate> = (0..h.len() - 1)
        .filter_map(|i| {
            let k = n - 1 - i;
            (k >= 2 && k <= kmax).then(|| Candidate {
                k,
                gap: h[i + 1] - h[i],
            })
        })
        .collect();
    cands.sort_by(|a, b| b.gap.partial_cmp(&a.gap).unwrap().then(a.k.cmp(&b.k)));
    cands.truncate(3);
    if top <= floor {
        return none(Verdict::Single, top, cands);
    }
    let Some(best) = cands.first().cloned() else {
        return none(Verdict::TooFew, top, cands);
    };
    let ratio = match cands.get(1) {
        Some(b) if b.gap > 0.0 => best.gap / b.gap,
        _ => f64::INFINITY,
    };
    let gap_over_floor = if floor > 0.0 {
        best.gap / floor
    } else {
        f64::INFINITY
    };
    let verdict = if gap_over_floor >= min_floors {
        Verdict::Clear(best.k)
    } else {
        Verdict::NoClearCut
    };
    GapCall {
        verdict,
        ratio,
        null_p: null_ratio_p(n, kmax, ratio),
        gap_over_floor,
        candidates: cands,
        top_height: top,
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn matrix(n: usize, f: impl Fn(usize, usize) -> f64) -> Vec<f64> {
        let mut d = vec![0.0; n * n];
        for i in 0..n {
            for j in 0..n {
                if i != j {
                    d[i * n + j] = f(i.min(j), i.max(j));
                }
            }
        }
        d
    }

    #[test]
    fn upgma_known_heights() {
        // (0,1) at 1, (2,3) at 2, then the clusters are 6 apart
        let d = matrix(4, |i, j| match (i, j) {
            (0, 1) => 1.0,
            (2, 3) => 2.0,
            _ => 6.0,
        });
        let m = upgma(&d, 4);
        let h: Vec<f64> = m.iter().map(|x| x.height).collect();
        assert_eq!(h, vec![1.0, 2.0, 6.0]);
        assert_eq!(m[2].size, 4);
        assert_eq!((m[0].a, m[0].b), (0, 1));
        assert_eq!((m[2].a, m[2].b), (4, 5));
    }

    #[test]
    fn upgma_averages_unequal_clusters() {
        // 0,1 close; 2 is 4 from 0 and 6 from 1: the merge with {0,1} is at (4 + 6) / 2 = 5
        let d = matrix(3, |i, j| match (i, j) {
            (0, 1) => 1.0,
            (0, 2) => 4.0,
            _ => 6.0,
        });
        assert_eq!(upgma(&d, 3)[1].height, 5.0);
    }

    #[test]
    fn cut_gives_the_requested_number_of_clusters() {
        let d = matrix(6, |i, j| if i / 2 == j / 2 { 0.1 } else { 5.0 });
        let m = upgma(&d, 6);
        assert_eq!(cut_k(&m, 6, 3), vec![0, 0, 1, 1, 2, 2]);
        assert_eq!(cut_k(&m, 6, 1), vec![0; 6]);
        assert_eq!(cut_k(&m, 6, 6), vec![0, 1, 2, 3, 4, 5]);
    }

    fn three_groups(within: f64, between: f64) -> Vec<f64> {
        matrix(30, |i, j| {
            let jitter = ((i * 7 + j * 13) % 10) as f64 * within * 0.05;
            if i / 10 == j / 10 {
                within + jitter
            } else {
                between + jitter * 3.0
            }
        })
    }

    #[test]
    fn a_jump_far_above_the_floor_is_cut_at_three() {
        let m = upgma(&three_groups(0.001, 0.02), 30);
        let g = gap_rule(&m, 30, 3.0, 0.0005, 30);
        assert_eq!(g.verdict, Verdict::Clear(3), "{g:?}");
        assert_eq!(g.candidates[0].k, 3);
        assert!(g.gap_over_floor > 20.0);
        let labels = cut_k(&m, 30, 3);
        assert!(labels
            .iter()
            .enumerate()
            .all(|(i, &l)| l == labels[i / 10 * 10]));
    }

    #[test]
    fn everything_below_the_floor_is_one_population() {
        let m = upgma(&three_groups(0.0002, 0.0004), 30);
        let g = gap_rule(&m, 30, 3.0, 0.001, 30);
        assert_eq!(g.verdict, Verdict::Single, "{g:?}");
    }

    #[test]
    fn jumps_below_the_floor_give_candidates_not_a_cut() {
        // the Bd situation: the biggest jump is real but smaller than the technical noise
        let m = upgma(&three_groups(0.0003, 0.0009), 30);
        let g = gap_rule(&m, 30, 3.0, 0.002, 30);
        assert_eq!(g.verdict, Verdict::Single, "{g:?}"); // top merge is itself under the floor
        let m = upgma(&three_groups(0.0006, 0.0030), 30);
        let g = gap_rule(&m, 30, 3.0, 0.002, 30);
        assert_eq!(g.verdict, Verdict::NoClearCut, "{g:?}");
        assert_eq!(g.candidates.len(), 3);
        assert!(g.gap_over_floor < 3.0);
    }

    #[test]
    fn a_uniform_continuum_still_has_a_gap_larger_than_the_floor() {
        // average linkage on equally spaced points makes the top merge a large jump: with a tiny floor it is "clear". The
        // ratio is not the gate precisely because structureless data produce large ratios too
        let m = upgma(&matrix(30, |i, j| (j - i) as f64 * 0.01), 30);
        let g = gap_rule(&m, 30, 3.0, 0.0001, 30);
        assert!(matches!(g.verdict, Verdict::Clear(_)), "{g:?}");
        assert!(
            g.null_p > 0.2,
            "this ratio is unremarkable for a structureless population: {}",
            g.null_p
        );
    }

    #[test]
    fn null_probability_is_monotone_and_matches_the_simulation() {
        // from analysis/popgen_cocci/26_gap_null.py: P(ratio >= 2) about 0.62, >= 5 about 0.26, >= 10 about 0.10
        let (p2, p5, p10) = (
            null_ratio_p(100, 30, 2.0),
            null_ratio_p(100, 30, 5.0),
            null_ratio_p(100, 30, 10.0),
        );
        assert!(p2 > p5 && p5 > p10);
        assert!(
            (p2 - 0.62).abs() < 0.03 && (p5 - 0.26).abs() < 0.03 && (p10 - 0.10).abs() < 0.03,
            "{p2} {p5} {p10}"
        );
        assert_eq!(null_ratio_p(100, 30, 5.0), p5); // reproducible
        assert_eq!(null_ratio_p(100, 30, f64::INFINITY), 0.0);
    }

    #[test]
    fn too_few_samples_are_not_called() {
        let m = upgma(&matrix(3, |_, _| 1.0), 3);
        assert_eq!(gap_rule(&m, 3, 3.0, 0.0, 30).verdict, Verdict::TooFew);
    }
}
