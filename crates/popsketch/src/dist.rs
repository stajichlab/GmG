//! Pairwise containment ANI and mash distance between sorted hash sets.

use crate::par::*;

#[derive(Clone, Debug)]
pub struct DistMat {
    pub n: usize,
    /// containment ANI estimate, C^(1/k) with C = shared / smaller set; diagonal 1
    pub ani: Vec<f64>,
    /// mash distance from the Jaccard index; diagonal 0
    pub mash: Vec<f64>,
}

impl DistMat {
    pub fn ani(&self, i: usize, j: usize) -> f64 {
        self.ani[i * self.n + j]
    }
    pub fn mash(&self, i: usize, j: usize) -> f64 {
        self.mash[i * self.n + j]
    }
    /// 1 - ANI as a dense row-major matrix (the distance used for clustering and PCoA)
    pub fn one_minus_ani(&self) -> Vec<f64> {
        self.ani.iter().map(|a| 1.0 - a).collect()
    }
    pub fn subset(&self, keep: &[usize]) -> DistMat {
        let m = keep.len();
        let mut ani = Vec::with_capacity(m * m);
        let mut mash = Vec::with_capacity(m * m);
        for &i in keep {
            for &j in keep {
                ani.push(self.ani(i, j));
                mash.push(self.mash(i, j));
            }
        }
        DistMat { n: m, ani, mash }
    }
}

/// Number of values two sorted, de-duplicated slices share.
pub fn intersect_count(a: &[u64], b: &[u64]) -> usize {
    let (mut i, mut j, mut n) = (0, 0, 0);
    while i < a.len() && j < b.len() {
        match a[i].cmp(&b[j]) {
            std::cmp::Ordering::Less => i += 1,
            std::cmp::Ordering::Greater => j += 1,
            std::cmp::Ordering::Equal => {
                n += 1;
                i += 1;
                j += 1;
            }
        }
    }
    n
}

/// Containment ANI and mash distance of one pair from their set sizes and intersection.
pub fn pair_stats(na: usize, nb: usize, inter: usize, k: u32) -> (f64, f64) {
    let c = if na.min(nb) == 0 {
        0.0
    } else {
        inter as f64 / na.min(nb) as f64
    };
    let ani = c.powf(1.0 / k as f64);
    let union = na + nb - inter;
    let j = if union == 0 {
        0.0
    } else {
        inter as f64 / union as f64
    };
    let mash = if j > 0.0 {
        -(2.0 * j / (1.0 + j)).ln() / k as f64
    } else {
        1.0
    };
    (ani, mash)
}

fn assemble(n: usize, rows: Vec<Vec<(f64, f64)>>) -> DistMat {
    let mut ani = vec![0.0; n * n];
    let mut mash = vec![0.0; n * n];
    for (i, row) in rows.iter().enumerate() {
        ani[i * n + i] = 1.0;
        for (off, &(a, m)) in row.iter().enumerate() {
            let j = i + 1 + off;
            ani[i * n + j] = a;
            ani[j * n + i] = a;
            mash[i * n + j] = m;
            mash[j * n + i] = m;
        }
    }
    DistMat { n, ani, mash }
}

/// All pairs by merging the two sorted hash lists: the cost of a pair is the size of the sets.
pub fn compute_merge(sets: &[Vec<u64>], k: u32) -> DistMat {
    let n = sets.len();
    let rows: Vec<Vec<(f64, f64)>> = (0..n)
        .into_par_iter()
        .map(|i| {
            ((i + 1)..n)
                .map(|j| {
                    let inter = intersect_count(&sets[i], &sets[j]);
                    pair_stats(sets[i].len(), sets[j].len(), inter, k)
                })
                .collect()
        })
        .collect();
    assemble(n, rows)
}

/// Memory the bitset method may use for the presence matrix.
/// Default cap on the presence matrix (`n` samples x variable hashes / 8 bytes); above it the pairwise merge is used.
pub const BITSET_MAX_BYTES: usize = 4 << 30;
/// Below this many samples the merge is fast enough and simpler.
const BITSET_MIN_SAMPLES: usize = 64;

/// How `compute` got its answer, so a slow run can be explained: the method, the number of variable hashes (present in
/// at least two but not all samples) and the size the presence matrix has or would have had.
#[derive(Clone, Debug, PartialEq)]
pub struct Plan {
    /// `bitset` or `merge`
    pub method: &'static str,
    pub variable_hashes: usize,
    pub bitset_bytes: usize,
    pub limit_bytes: usize,
}

impl Plan {
    pub fn describe(&self) -> String {
        format!(
            "distance method {}: {} variable hashes, presence matrix {:.2} GiB (limit {:.2} GiB)",
            self.method,
            self.variable_hashes,
            self.bitset_bytes as f64 / (1u64 << 30) as f64,
            self.limit_bytes as f64 / (1u64 << 30) as f64
        )
    }
}

/// All pairs, exactly, from bitsets. In population data most hashes are shared by everyone, so merging costs the whole set size
/// for every pair. Here the hashes present in all samples are counted once (`core`), hashes seen in only one sample cannot be
/// shared, and the rest (the variable hashes) become one bit per sample: the intersection of a pair is `core` plus the popcount
/// of the AND of their bits. Returns `Err(plan)` when the presence matrix would be larger than `limit` bytes.
pub fn compute_bitset_limit(
    sets: &[Vec<u64>],
    k: u32,
    limit: usize,
) -> Result<(DistMat, Plan), Plan> {
    let n = sets.len();
    if n == 0 {
        let plan = Plan {
            method: "bitset",
            variable_hashes: 0,
            bitset_bytes: 0,
            limit_bytes: limit,
        };
        return Ok((
            DistMat {
                n: 0,
                ani: vec![],
                mash: vec![],
            },
            plan,
        ));
    }
    let mut freq: std::collections::HashMap<u64, u32> = std::collections::HashMap::new();
    for s in sets {
        for &h in s {
            *freq.entry(h).or_insert(0) += 1;
        }
    }
    let core = freq.values().filter(|&&f| f as usize == n).count();
    let mut cols: Vec<u64> = freq
        .iter()
        .filter(|(_, &f)| f >= 2 && (f as usize) < n)
        .map(|(&h, _)| h)
        .collect();
    cols.sort_unstable();
    let words = cols.len().div_ceil(64);
    let bytes = n.saturating_mul(words).saturating_mul(8);
    let plan = Plan {
        method: "bitset",
        variable_hashes: cols.len(),
        bitset_bytes: bytes,
        limit_bytes: limit,
    };
    if bytes > limit {
        return Err(Plan {
            method: "merge",
            ..plan
        });
    }
    drop(freq);
    let bits: Vec<Vec<u64>> = sets
        .par_iter()
        .map(|s| {
            let mut b = vec![0u64; words];
            for h in s {
                if let Ok(c) = cols.binary_search(h) {
                    b[c / 64] |= 1 << (c % 64);
                }
            }
            b
        })
        .collect();
    let rows: Vec<Vec<(f64, f64)>> = (0..n)
        .into_par_iter()
        .map(|i| {
            ((i + 1)..n)
                .map(|j| {
                    let shared: u32 = bits[i]
                        .iter()
                        .zip(&bits[j])
                        .map(|(x, y)| (x & y).count_ones())
                        .sum();
                    pair_stats(sets[i].len(), sets[j].len(), core + shared as usize, k)
                })
                .collect()
        })
        .collect();
    Ok((assemble(n, rows), plan))
}

/// Bitset with the default limit (None when the presence matrix would be larger than `BITSET_MAX_BYTES`).
pub fn compute_bitset(sets: &[Vec<u64>], k: u32) -> Option<DistMat> {
    compute_bitset_limit(sets, k, BITSET_MAX_BYTES)
        .ok()
        .map(|(d, _)| d)
}

/// Exact pairwise containment ANI and mash distance of all samples (bitsets for hundreds of samples, merging for few),
/// with the method that was used.
pub fn compute_planned(sets: &[Vec<u64>], k: u32, limit: usize) -> (DistMat, Plan) {
    if sets.len() >= BITSET_MIN_SAMPLES {
        match compute_bitset_limit(sets, k, limit) {
            Ok(r) => return r,
            Err(plan) => return (compute_merge(sets, k), plan),
        }
    }
    let plan = Plan {
        method: "merge",
        variable_hashes: 0,
        bitset_bytes: 0,
        limit_bytes: limit,
    };
    (compute_merge(sets, k), plan)
}

/// `compute_planned` with the default limit, without the plan.
pub fn compute(sets: &[Vec<u64>], k: u32) -> DistMat {
    compute_planned(sets, k, BITSET_MAX_BYTES).0
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn intersection_of_sorted_slices() {
        assert_eq!(intersect_count(&[1, 3, 5, 7], &[2, 3, 7, 9]), 2);
        assert_eq!(intersect_count(&[], &[1]), 0);
    }

    #[test]
    fn identical_sets_have_ani_one_and_distance_zero() {
        let a: Vec<u64> = (0..100).collect();
        let (ani, mash) = pair_stats(100, 100, 100, 21);
        assert!((ani - 1.0).abs() < 1e-12 && mash.abs() < 1e-12);
        let m = compute(&[a.clone(), a], 21);
        assert!((m.ani(0, 1) - 1.0).abs() < 1e-12);
    }

    #[test]
    fn containment_uses_the_smaller_set() {
        // a small set fully inside a big one is ANI 1, but its Jaccard (and mash distance) is not
        let big: Vec<u64> = (0..1000).collect();
        let small: Vec<u64> = (0..100).collect();
        let m = compute(&[big, small], 21);
        assert!((m.ani(0, 1) - 1.0).abs() < 1e-12);
        assert!(m.mash(0, 1) > 0.05);
    }

    #[test]
    fn ani_known_value() {
        // 80% of the smaller set shared at k=21: ANI = 0.8^(1/21)
        let (ani, _) = pair_stats(100, 100, 80, 21);
        assert!((ani - 0.8f64.powf(1.0 / 21.0)).abs() < 1e-12);
        assert_eq!(pair_stats(0, 5, 0, 21).0, 0.0);
    }

    #[test]
    fn subset_keeps_pairs() {
        let sets = vec![
            (0..50).collect::<Vec<u64>>(),
            (25..75).collect(),
            (0..10).collect(),
        ];
        let m = compute(&sets, 21);
        let s = m.subset(&[2, 0]);
        assert_eq!(s.n, 2);
        assert_eq!(s.ani(0, 1), m.ani(2, 0));
    }

    /// `n` samples: a core present in all, a pool of variable hashes each sample carries with probability `p`, and private hashes
    fn population(
        n: usize,
        core: u64,
        pool: u64,
        p: f64,
        private: u64,
        seed: u64,
    ) -> Vec<Vec<u64>> {
        let mut x = seed;
        let mut next = move || {
            x ^= x << 13;
            x ^= x >> 7;
            x ^= x << 17;
            (x >> 11) as f64 / (1u64 << 53) as f64
        };
        (0..n as u64)
            .map(|i| {
                let mut h: Vec<u64> = (0..core).collect();
                h.extend((0..pool).filter(|_| next() < p).map(|v| 1_000_000 + v));
                h.extend((0..private).map(|v| 100_000_000 + i * 10_000 + v));
                h
            })
            .collect()
    }

    #[test]
    fn the_bitset_method_gives_exactly_the_merge_result() {
        for (n, p) in [(2usize, 0.5), (3, 0.5), (70, 0.4), (90, 0.9)] {
            let sets = population(n, 400, 600, p, 30, 7 + n as u64);
            let (a, b) = (compute_merge(&sets, 21), compute_bitset(&sets, 21).unwrap());
            assert_eq!(a.n, b.n);
            for (x, y) in a.ani.iter().zip(&b.ani).chain(a.mash.iter().zip(&b.mash)) {
                assert!((x - y).abs() < 1e-12, "n={n}: {x} vs {y}");
            }
        }
    }

    #[test]
    fn compute_uses_bitsets_for_many_samples_and_merging_for_few_and_both_agree() {
        let few = population(10, 200, 300, 0.5, 10, 3);
        assert_eq!(compute(&few, 21).ani, compute_merge(&few, 21).ani);
        let many = population(80, 200, 300, 0.5, 10, 4);
        assert_eq!(
            compute(&many, 21).ani,
            compute_bitset(&many, 21).unwrap().ani
        );
        assert!(compute_bitset(&[], 21).unwrap().n == 0);
    }

    /// cargo test -p popsketch --release -- --ignored --nocapture: the two methods on a population-sized problem
    #[test]
    #[ignore]
    fn timing_on_a_population_sized_set() {
        for n in [300usize, 1500] {
            let sets = population(n, 120_000, 200_000, 0.3, 500, 11);
            let t = std::time::Instant::now();
            let b = compute_bitset(&sets, 21).unwrap();
            let tb = t.elapsed();
            let merge = if n == 300 {
                let t = std::time::Instant::now();
                let m = compute_merge(&sets, 21);
                assert_eq!(m.ani, b.ani);
                Some(t.elapsed())
            } else {
                None
            };
            eprintln!(
                "n={n}: bitset {tb:.1?}{}",
                merge.map_or(String::new(), |m| format!(", merge {m:.1?}"))
            );
        }
    }

    #[test]
    fn plan_reports_the_method_and_a_lower_limit_forces_the_merge_with_the_same_answer() {
        let many: Vec<Vec<u64>> = (0..70u64)
            .map(|i| {
                (0..400u64)
                    .filter(|h| h % 7 != i % 7)
                    .chain(1000 + i * 1000..1000 + i * 1000 + 50)
                    .collect()
            })
            .collect();
        let (d1, p1) = compute_planned(&many, 21, BITSET_MAX_BYTES);
        assert_eq!(p1.method, "bitset");
        assert!(p1.variable_hashes > 0 && p1.bitset_bytes > 0);
        let (d2, p2) = compute_planned(&many, 21, 8);
        assert_eq!(p2.method, "merge");
        assert_eq!(p2.variable_hashes, p1.variable_hashes);
        assert_eq!(d1.ani, d2.ani);
        let few = &many[..10];
        assert_eq!(compute_planned(few, 21, BITSET_MAX_BYTES).1.method, "merge");
    }
}
