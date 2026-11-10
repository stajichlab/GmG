//! The split-half noise floor: the distance (1 - ANI) between two sketches of alternate reads of the same sample is the
//! technical noise of a comparison. The floor is the 90th percentile over a subset of samples (the median is too low: a few
//! low-coverage halves have floors ten times the median, and the median floor left 75% of the pairs inside Bd lineages 2/3/4
//! above it). Each half has half the coverage, so the floor is pessimistic. Without enough reads to split, a fixed floor
//! measured on Coccidioides and Bd is used.

use crate::dist::{intersect_count, pair_stats};
use crate::qc::kept_hashes;
use crate::sketch_io::SampleSketch;
use crate::stats::{median, quantile};
use gmg_core::coverage::auto_min_abund;
use std::collections::HashSet;

/// Used when fewer than `MIN_SAMPLES` samples have halves. Between the median (0.0002-0.0005) and the 90th percentile
/// (0.0016-0.0020) of the split-half distances measured on the two datasets; this is the conservative side of the 90th.
pub const FALLBACK_FLOOR: f64 = 0.001;
pub const MIN_SAMPLES: usize = 8;

#[derive(Clone, Debug)]
pub struct Floor {
    pub value: f64,
    pub n: usize,
    pub median: f64,
    pub p90: f64,
    /// `split_half` or `fallback`
    pub source: &'static str,
    /// each sample's own split-half distance (1 - ANI between its two halves), in input order
    pub per_sample: Vec<(String, f64)>,
}

fn fnv1a(s: &str) -> u64 {
    let mut h: u64 = 0xcbf2_9ce4_8422_2325;
    for b in s.bytes() {
        h ^= b as u64;
        h = h.wrapping_mul(0x0100_0000_01b3);
    }
    h
}

/// The ids that get half-sketches: a fixed pseudo-random subset of `n` (so a rerun picks the same samples).
pub fn pick_for_halves(ids: &[String], n: usize) -> HashSet<String> {
    let mut v: Vec<&String> = ids.iter().collect();
    v.sort_by_key(|s| (fnv1a(s), (*s).clone()));
    v.into_iter().take(n).cloned().collect()
}

/// Estimate the floor from the samples that carry halves.
pub fn estimate(samples: &[&SampleSketch], k: u32) -> Floor {
    let mut d: Vec<f64> = Vec::new();
    let mut per_sample: Vec<(String, f64)> = Vec::new();
    for s in samples {
        let Some((a, b)) = &s.halves else { continue };
        let (ca, cb) = (auto_min_abund(&a.abunds), auto_min_abund(&b.abunds));
        let keep = |sp: &crate::sketch_io::Spectrum, c: u32| -> Vec<u64> {
            sp.hashes
                .iter()
                .zip(&sp.abunds)
                .filter(|(_, &x)| x >= c)
                .map(|(&h, _)| h)
                .collect()
        };
        let (ha, hb) = (keep(a, ca), keep(b, cb));
        if ha.is_empty() || hb.is_empty() {
            continue;
        }
        let inter = intersect_count(&ha, &hb);
        let (ani, _) = pair_stats(ha.len(), hb.len(), inter, k);
        d.push((1.0 - ani).max(0.0));
        per_sample.push((s.id.clone(), (1.0 - ani).max(0.0)));
        let _ = kept_hashes;
    }
    if d.len() >= MIN_SAMPLES {
        let p90 = quantile(&d, 0.9);
        Floor {
            value: p90,
            n: d.len(),
            median: median(&d),
            p90,
            source: "split_half",
            per_sample,
        }
    } else {
        Floor {
            value: FALLBACK_FLOOR,
            n: d.len(),
            median: median(&d),
            p90: quantile(&d, 0.9),
            source: "fallback",
            per_sample,
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::sketch_io::{Counts, Spectrum};

    fn with_halves(id: &str, shared: usize, private: usize) -> SampleSketch {
        // two halves sharing `shared` hashes at high abundance and `private` hashes each (also high abundance)
        let mk = |off: u64| {
            let mut hashes: Vec<u64> = (0..shared as u64).collect();
            hashes.extend((0..private as u64).map(|i| 1_000_000 * off + i));
            hashes.sort_unstable();
            let abunds = vec![20; hashes.len()];
            Spectrum { hashes, abunds }
        };
        SampleSketch {
            id: id.into(),
            kind: "reads".into(),
            files: vec![],
            spectrum: Spectrum::default(),
            counts: Counts::default(),
            halves: Some((mk(1), mk(2))),
        }
    }

    #[test]
    fn floor_is_the_90th_percentile_of_half_distances() {
        let s: Vec<SampleSketch> = (0..10)
            .map(|i| with_halves(&format!("s{i}"), 20_000, 20 + 40 * i))
            .collect();
        let r: Vec<&SampleSketch> = s.iter().collect();
        let f = estimate(&r, 21);
        assert_eq!(f.source, "split_half");
        assert_eq!(f.n, 10);
        assert!(f.value > f.median && (f.value - f.p90).abs() < 1e-12);
        // 20,380 vs 20,020... containment of the larger-private half is shared/min: dist = 1 - (shared/min)^(1/21) is small
        assert!(f.value > 1e-5 && f.value < 0.01, "{}", f.value);
        // each sample's own distance is kept, in input order, and grows with the private hashes of its halves
        assert_eq!(f.per_sample.len(), 10);
        assert_eq!(f.per_sample[0].0, "s0");
        assert!(f.per_sample.windows(2).all(|w| w[0].1 < w[1].1));
        let d: Vec<f64> = f.per_sample.iter().map(|x| x.1).collect();
        assert!((median(&d) - f.median).abs() < 1e-12);
    }

    #[test]
    fn too_few_samples_use_the_fallback() {
        let s: Vec<SampleSketch> = (0..3)
            .map(|i| with_halves(&format!("s{i}"), 20_000, 50))
            .collect();
        let r: Vec<&SampleSketch> = s.iter().collect();
        let f = estimate(&r, 21);
        assert_eq!((f.source, f.value), ("fallback", FALLBACK_FLOOR));
        assert_eq!(estimate(&[], 21).source, "fallback");
    }

    #[test]
    fn the_subset_is_stable_and_the_requested_size() {
        let ids: Vec<String> = (0..100).map(|i| format!("id{i}")).collect();
        let a = pick_for_halves(&ids, 20);
        let mut rev = ids.clone();
        rev.reverse();
        assert_eq!(a, pick_for_halves(&rev, 20));
        assert_eq!(a.len(), 20);
        assert_eq!(pick_for_halves(&ids[..5], 20).len(), 5);
    }
}
