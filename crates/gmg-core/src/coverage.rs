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
/// A valley is only accepted if the hashes seen at least that often make up at least this share of the hashes seen twice or
/// more. A real coverage peak holds most of the non-error hashes (about 98% in a healthy 40x sample); a "valley" with a
/// sliver above it is noise in the tail of an error-dominated histogram (sub-1x samples, or a plateau at low coverage with no
/// rise after it), where the walk below would otherwise run on to a spurious minimum and leave almost no hashes.
const PEAK_MASS_SHARE: f64 = 0.5;
/// Deep samples (about 150x and more) fail that test although the peak is clear: recurrent sequencing errors seen 2 to 30
/// times form a long tail that holds more hashes than the peak (share 0.37 to 0.46 in *C. neoformans* at 4 Gb, #157). A peak
/// is also accepted when the hashes within 40% of its position make up at least this share of the hashes seen twice or more.
/// False peaks (high-copy bumps, noise) hold 0.10 or less of them in the sets checked.
const PEAK_WINDOW_SHARE: f64 = 0.3;

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
    let above: f64 = h[valley..].iter().map(|&c| c as f64).sum();
    let twice_or_more: f64 = h[2..].iter().map(|&c| c as f64).sum();
    let mut peak = valley + 1;
    for b in valley + 1..=last {
        if smooth(h, b) > smooth(h, peak) {
            peak = b;
        }
    }
    let (lo, hi) = (
        ((peak as f64 * 0.6) as usize).max(valley),
        (peak as f64 * 1.4).ceil() as usize,
    );
    let window: f64 = h[lo..=hi.min(h.len() - 1)].iter().map(|&c| c as f64).sum();
    if above < PEAK_MASS_SHARE * twice_or_more && window < PEAK_WINDOW_SHARE * twice_or_more {
        return None;
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

    /// Deep sample: a Poisson peak at `lambda` over a long error tail (errors seen 1 to 30 times, falling slowly), as in the
    /// real 4 Gb histograms of #157. Without the window rule the tail outweighs the peak and no peak is reported.
    fn deep_hist(n: usize, lambda: f64) -> Vec<u32> {
        let mut v = poisson_hist(n, lambda, 0);
        for a in 1..=30u32 {
            let count = (n as f64 * 3.0 / (a as f64).powf(2.2)) as usize;
            v.extend(std::iter::repeat_n(a, count));
        }
        v
    }

    #[test]
    fn peak_is_found_at_very_high_coverage_with_a_heavy_error_tail() {
        for lambda in [100.0f64, 150.0, 300.0] {
            let c = profile(&deep_hist(60000, lambda)).unwrap();
            assert_eq!(c.basis, "peak", "lambda {lambda}");
            assert!(
                (c.kmer - lambda).abs() <= lambda * 0.1,
                "lambda {lambda}: got {}",
                c.kmer
            );
            assert!(c.valley.unwrap() >= 2 && c.valley.unwrap() < lambda as u32);
        }
    }

    #[test]
    fn a_small_high_copy_bump_over_an_error_tail_is_not_a_peak() {
        // error tail with 400 hashes in a bump at 90: far below the share of the hashes seen twice or more
        let mut v: Vec<u32> = (1..=30u32)
            .flat_map(|a| std::iter::repeat_n(a, (20000.0 / (a as f64).powf(2.0)) as usize))
            .collect();
        v.extend(std::iter::repeat_n(90, 400));
        assert_ne!(profile(&v).map(|c| c.basis), Some("peak"));
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

    /// Expand the first bins of a real histogram (counts of hashes seen 1, 2, 3, ... times) into abundances.
    fn from_bins(bins: &[u32], tail: &[(u32, u32)]) -> Vec<u32> {
        let mut out = Vec::new();
        for (i, &n) in bins.iter().enumerate() {
            out.extend(std::iter::repeat_n(i as u32 + 1, n as usize));
        }
        for &(a, n) in tail {
            out.extend(std::iter::repeat_n(a, n as usize));
        }
        out
    }

    /// Histograms of read sets where the valley search used to land on noise far above the real (missing or very low) peak and
    /// leave 0.1-3 Mb of hashes (scaled 200): a sub-1x sample, a low-coverage plateau with no rise, and a Bd sub-2x sample.
    #[test]
    fn no_spurious_valley_in_error_dominated_histograms() {
        // sub-1x: error tail plus a thin scatter of repeat k-mers out to high abundance
        let sub1x = from_bins(
            &[
                47862, 12802, 3049, 727, 253, 137, 63, 51, 51, 46, 24, 22, 24, 16, 12, 19, 21, 25,
                21, 29, 22, 31, 26, 28, 25, 6, 20, 27, 18, 14, 19, 9, 10, 8, 8, 7, 8, 3, 4, 1,
            ],
            &[(60, 10), (80, 5), (100, 3), (150, 2), (400, 1)],
        );
        assert_eq!(auto_min_abund(&sub1x), 2, "sub-1x");
        assert_eq!(profile(&sub1x).unwrap().basis, "poisson");
        // plateau at 4-7 over a huge error tail, no rise after it
        let plateau = from_bins(
            &[
                1252617, 67902, 28616, 21793, 21254, 20054, 17468, 13978, 10408, 6812, 4323, 2622,
                1631, 1003, 660, 425, 306, 249, 242, 183, 159, 145, 120, 120, 98, 106, 86, 53, 75,
                63, 67, 51, 40, 61, 38, 38, 40, 31, 23, 25,
            ],
            &[(50, 600), (100, 100)],
        );
        assert_eq!(auto_min_abund(&plateau), 2, "plateau");
        // sub-2x with a long thin tail
        let sub2x = from_bins(
            &[
                379049, 43550, 14214, 5400, 2372, 1241, 794, 535, 354, 263, 214, 192, 147, 138,
                106, 107, 77, 80, 69, 70, 34, 31, 32, 24, 34, 29, 20, 28, 28, 28, 31, 32, 32, 24,
                22, 12, 20, 22, 26, 27,
            ],
            &[(60, 400), (80, 700), (100, 2000), (150, 3000)],
        );
        assert_eq!(auto_min_abund(&sub2x), 2, "sub-2x");
    }

    #[test]
    fn a_real_peak_still_gets_its_valley() {
        // shape of a healthy 14x read set (scaled 200): error tail, valley near 10, peak near 38 (first 40 bins) then the rest
        let healthy = from_bins(
            &[
                356742, 9265, 1485, 581, 304, 236, 179, 143, 146, 105, 114, 107, 119, 134, 118,
                122, 144, 155, 162, 211, 219, 291, 334, 403, 474, 619, 769, 979, 1312, 1597, 1976,
                2302, 2742, 3284, 3809, 4185, 4809, 5243, 5699, 6024,
            ],
            &[
                (43, 25616),
                (47, 22450),
                (51, 15436),
                (55, 8424),
                (59, 3911),
                (63, 1477),
                (67, 523),
                (71, 222),
            ],
        );
        let v = auto_min_abund(&healthy);
        assert!((8..=20).contains(&v), "valley {v}");
        assert_eq!(profile(&healthy).unwrap().basis, "peak");
    }

    #[test]
    fn base_coverage_scales_with_read_length() {
        let c = base_coverage(10.0, 150.0, 21);
        assert!((c - 10.0 * 150.0 / 130.0).abs() < 1e-9);
        assert_eq!(base_coverage(10.0, 10.0, 21), 10.0); // reads shorter than k: unchanged
    }
}
