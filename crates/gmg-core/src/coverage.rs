//! Coverage estimate and abundance cutoff from the k-mer abundance histogram of a read set (#32).
//!
//! A sketch of reads carries, for every kept hash, how often its k-mer was seen. For a genome sequenced to coverage `c`
//! the abundances of its k-mers are roughly Poisson with mean `lambda = c * (L - k + 1) / L` (reads of length `L`), while
//! sequencing errors add k-mers seen once or twice. The histogram is therefore an error tail that falls from 1, a valley,
//! and (from about 3x) a peak at `lambda`.
//!
//! * [`profile`] finds that structure and returns the k-mer coverage `lambda`:
//!   * **peak**: the histogram has a valley followed by a clear peak; `lambda` is the mode of the peak (counts smoothed over
//!     three bins);
//!   * **poisson**: no peak (coverage below about 3x, or too few hashes): `lambda = 3 * h3 / h2` (the ratio of the counts
//!     of hashes seen 3 and 2 times of a Poisson distribution, which errors hardly touch) when both counts are at least
//!     10, else `lambda = 2 * h2 / h1` when `h2 >= 8`; otherwise there is no estimate.
//! * [`base_coverage`] converts it to genome coverage with the mean read length.
//! * [`auto_min_abund`] is the abundance cutoff from the valley (hashes seen fewer times are errors): the valley position,
//!   at least 2, and 2 when the histogram has no peak.
//!
//! The estimate is the coverage of the reads that were sketched (so it follows `--max-reads`). For RNA-seq it is only the
//! dominant k-mer abundance, not a genome coverage. Heterozygous or polyploid genomes can show a second peak at a multiple;
//! the first (lowest) strong peak after the valley is reported.

/// Longest abundance resolved; larger counts are pooled in the last bin (coverage above this is reported as this).
const MAX_BIN: usize = 4096;

/// A peak must hold at least this many hashes (smoothed) and exceed the valley by this factor.
const PEAK_MIN_HASHES: f64 = 8.0;
const PEAK_OVER_VALLEY: f64 = 1.3;

#[derive(Clone, Debug, PartialEq)]
pub struct Coverage {
    /// k-mer coverage: the mean abundance of the genome's k-mers
    pub kmer: f64,
    /// `peak` or `poisson`
    pub basis: &'static str,
    /// the abundance at the bottom of the valley between the error tail and the peak (only with `peak`)
    pub valley: Option<u32>,
}

fn histogram(abunds: &[u32]) -> Vec<u64> {
    let max = abunds.iter().copied().max().unwrap_or(0) as usize;
    let mut h = vec![0u64; max.min(MAX_BIN) + 2];
    for &a in abunds {
        if a >= 1 {
            h[(a as usize).min(MAX_BIN)] += 1;
        }
    }
    h
}

/// Count at `a` smoothed over `a - 1 ..= a + 1` (bins beyond the end are 0).
fn smooth(h: &[u64], a: usize) -> f64 {
    let g = |i: usize| h.get(i).copied().unwrap_or(0) as f64;
    (g(a.saturating_sub(1)) + g(a) + g(a + 1)) / 3.0
}

/// The valley and peak of the histogram, if it has a clear peak after an error tail.
fn valley_peak(h: &[u64]) -> Option<(usize, usize)> {
    let last = h.len().saturating_sub(2);
    // walk up from 2 while the smoothed count keeps falling; the first rise ends the error tail
    let mut a = 2usize;
    while a < last && smooth(h, a + 1) <= smooth(h, a) {
        a += 1;
    }
    if a >= last {
        return None;
    }
    let valley = a;
    let mut peak = valley + 1;
    for b in valley + 1..=last {
        if smooth(h, b) > smooth(h, peak) {
            peak = b;
        }
    }
    let (sp, sv) = (smooth(h, peak), smooth(h, valley));
    (sp >= PEAK_MIN_HASHES && sp >= PEAK_OVER_VALLEY * sv).then_some((valley, peak))
}

/// The structure of the abundance histogram: k-mer coverage and how it was found. None when there are too few hashes.
pub fn profile(abunds: &[u32]) -> Option<Coverage> {
    if abunds.len() < 30 {
        return None;
    }
    let h = histogram(abunds);
    if let Some((valley, peak)) = valley_peak(&h) {
        return Some(Coverage {
            kmer: peak as f64,
            basis: "peak",
            valley: Some(valley as u32),
        });
    }
    let g = |i: usize| h.get(i).copied().unwrap_or(0) as f64;
    let (h1, h2, h3) = (g(1), g(2), g(3));
    let lambda = if h3 >= 10.0 && h2 >= 10.0 {
        3.0 * h3 / h2
    } else if h2 >= 8.0 && h1 > 0.0 {
        2.0 * h2 / h1
    } else {
        return None;
    };
    Some(Coverage {
        kmer: lambda,
        basis: "poisson",
        valley: None,
    })
}

/// Genome coverage from the k-mer coverage: `lambda * L / (L - k + 1)` for reads of mean length `L`.
pub fn base_coverage(kmer_cov: f64, mean_read_len: f64, k: u32) -> f64 {
    let window = mean_read_len - k as f64 + 1.0;
    if window < 1.0 {
        return kmer_cov;
    }
    kmer_cov * mean_read_len / window
}

/// Abundance cutoff chosen from the histogram: the valley between the error tail and the coverage peak (at least 2), or 2
/// when there is no peak (low coverage: the fixed default of reads).
pub fn auto_min_abund(abunds: &[u32]) -> u32 {
    match profile(abunds) {
        Some(Coverage {
            valley: Some(v), ..
        }) => v.max(2),
        _ => 2,
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    /// Deterministic abundances: `n` k-mers with Poisson(`lambda`) counts (inverse-CDF on a stratified grid, so no randomness),
    /// plus `errors` k-mers seen once and `errors / 4` seen twice.
    fn poisson_hist(n: usize, lambda: f64, errors: usize) -> Vec<u32> {
        let mut out = Vec::new();
        let mut cdf = 0.0;
        let mut p = (-lambda).exp();
        let mut assigned = 0usize;
        let mut a = 0u32;
        while assigned < n && a < 400 {
            cdf += p;
            let upto = ((cdf * n as f64).round() as usize).min(n);
            for _ in assigned..upto {
                if a >= 1 {
                    out.push(a);
                }
            }
            assigned = upto.max(assigned);
            a += 1;
            p *= lambda / a as f64;
        }
        out.extend(std::iter::repeat_n(1, errors));
        out.extend(std::iter::repeat_n(2, errors / 4));
        out
    }

    #[test]
    fn peak_is_found_at_high_coverage() {
        for lambda in [6.0f64, 12.0, 25.0, 60.0] {
            let c = profile(&poisson_hist(4000, lambda, 500)).unwrap();
            assert_eq!(c.basis, "peak", "lambda {lambda}");
            assert!(
                (c.kmer - lambda).abs() <= (lambda * 0.15).max(1.0),
                "lambda {lambda}: got {}",
                c.kmer
            );
            assert!(c.valley.unwrap() < lambda as u32);
        }
    }

    #[test]
    fn poisson_ratio_is_used_at_low_coverage() {
        for lambda in [0.1f64, 0.3, 0.6, 1.0, 1.6, 2.2] {
            let c = profile(&poisson_hist(20000, lambda, 100)).unwrap();
            assert_eq!(c.basis, "poisson", "lambda {lambda}");
            assert!(
                (c.kmer - lambda).abs() <= (lambda * 0.2).max(0.05),
                "lambda {lambda}: got {}",
                c.kmer
            );
        }
    }

    #[test]
    fn too_few_hashes_or_no_pairs_give_no_estimate() {
        assert!(profile(&[1, 1, 2]).is_none());
        assert!(profile(&vec![1; 500]).is_none()); // only singletons
    }

    #[test]
    fn auto_cutoff_follows_the_valley() {
        assert_eq!(auto_min_abund(&poisson_hist(20000, 0.5, 100)), 2); // no peak
        let v = auto_min_abund(&poisson_hist(4000, 20.0, 800));
        assert!((2..=8).contains(&v), "valley {v}");
    }

    #[test]
    fn base_coverage_scales_with_read_length() {
        let c = base_coverage(10.0, 150.0, 21);
        assert!((c - 10.0 * 150.0 / 130.0).abs() < 1e-9);
        assert_eq!(base_coverage(10.0, 10.0, 21), 10.0); // reads shorter than k: unchanged
    }
}
