//! Pairwise containment ANI and mash distance between sorted hash sets.

use rayon::prelude::*;

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

pub fn compute(sets: &[Vec<u64>], k: u32) -> DistMat {
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
    let mut ani = vec![0.0; n * n];
    let mut mash = vec![0.0; n * n];
    for i in 0..n {
        ani[i * n + i] = 1.0;
        for (off, &(a, m)) in rows[i].iter().enumerate() {
            let j = i + 1 + off;
            ani[i * n + j] = a;
            ani[j * n + i] = a;
            mash[i * n + j] = m;
            mash[j * n + i] = m;
        }
    }
    DistMat { n, ani, mash }
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
}
