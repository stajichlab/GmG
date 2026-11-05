//! Cohort-relative QC of sketched samples. Nothing here needs a reference: each sample is compared with the rest of the
//! cohort, and every flag carries its reason and numbers.
//!
//! * `low_coverage`: the k-mer coverage from the abundance histogram is below `min_cov`, or there is no coverage peak.
//! * `size_high` / `size_low`: the estimated size (hashes above the abundance cutoff x scaled) is far from the cohort median
//!   (robust z, with a floor so a tight cohort does not flag everything). Heavy contamination (a second genome) shows as
//!   `size_high`; a sample missing half its content as `size_low`.
//! * `abundance_shoulder`: hashes at low-but-real abundance (0.1c-0.6c) beyond what a single genome predicts, standardised
//!   within coverage bins. Catches a contaminant whose own coverage is several-fold; it cannot see lighter contamination.

use crate::sketch_io::SampleSketch;
use crate::stats::{mad_sd, median, robust_z};
use gmg_core::coverage::{auto_min_abund, profile};

#[derive(Clone, Debug, PartialEq)]
pub struct Flag {
    pub code: &'static str,
    pub detail: String,
}

#[derive(Clone, Debug)]
pub struct QcParams {
    /// k-mer coverage below this is `low_coverage`
    pub min_cov: f64,
    /// |robust z| of the size above this is a size flag
    pub size_z: f64,
    /// floor of the size spread, as a fraction of the cohort median size
    pub size_floor: f64,
    /// robust z of the abundance shoulder above this is `abundance_shoulder`
    pub spectrum_z: f64,
    /// no relative flag is raised for a cohort smaller than this
    pub min_cohort: usize,
}

impl Default for QcParams {
    fn default() -> Self {
        QcParams {
            min_cov: 4.0,
            size_z: 5.0,
            size_floor: 0.03,
            spectrum_z: 5.0,
            min_cohort: 8,
        }
    }
}

#[derive(Clone, Debug)]
pub struct SampleQc {
    pub id: String,
    pub kind: String,
    pub reads: u64,
    pub bases: u64,
    pub n_hashes: usize,
    /// abundance cutoff applied (1 for assemblies)
    pub cutoff: u32,
    pub kept: usize,
    pub est_mb: f64,
    /// k-mer coverage and how it was found (`peak`, `poisson`); None for assemblies and when there is no estimate
    pub kmer_cov: Option<f64>,
    pub cov_basis: &'static str,
    pub shoulder: Option<f64>,
    pub shoulder_z: Option<f64>,
    pub size_z: Option<f64>,
    pub flags: Vec<Flag>,
}

impl SampleQc {
    pub fn flagged(&self) -> bool {
        !self.flags.is_empty()
    }
    pub fn has(&self, code: &str) -> bool {
        self.flags.iter().any(|f| f.code == code)
    }
}

/// Poisson probabilities for 0..=amax (lambda is clamped so exp(-lambda) does not underflow).
fn poisson_pmf(lambda: f64, amax: usize) -> Vec<f64> {
    let l = lambda.clamp(1e-9, 600.0);
    let mut p = vec![0.0; amax + 1];
    p[0] = (-l).exp();
    for a in 1..=amax {
        p[a] = p[a - 1] * l / a as f64;
    }
    p
}

/// Shoulder statistic: hashes with abundance in [lo, hi), lo = max(3, 0.1c), hi = ceil(0.6c), as a fraction of the genome
/// mass (hashes at or above hi). Also returns the excess over the Poisson(c) expectation, for tests and the JSON.
pub fn shoulder_features(abunds: &[u32], c: f64) -> Option<(f64, f64)> {
    let hi = (0.6 * c).ceil() as usize;
    let lo = ((0.1 * c).floor() as usize).max(3);
    if hi <= lo + 1 {
        return None;
    }
    const BINS: usize = 4001;
    let mut h = vec![0f64; BINS];
    for &a in abunds {
        h[(a as usize).min(BINS - 1)] += 1.0;
    }
    let mass: f64 = h[hi.min(BINS - 1)..].iter().sum();
    let pmf = poisson_pmf(c, hi + 1);
    let below_hi: f64 = pmf[..hi].iter().sum();
    let tail = 1.0 - below_hi;
    if mass <= 0.0 || tail <= 1e-9 {
        return None;
    }
    let genome = mass / tail;
    let obs: f64 = h[lo..hi].iter().sum();
    let exp = genome * pmf[lo..hi].iter().sum::<f64>();
    Some((obs / mass, (obs - exp) / mass))
}

const COV_BINS: [f64; 7] = [0.0, 6.0, 10.0, 15.0, 25.0, 40.0, f64::INFINITY];

fn cov_bin(c: f64) -> usize {
    COV_BINS
        .windows(2)
        .position(|w| c >= w[0] && c < w[1])
        .unwrap_or(COV_BINS.len() - 2)
}

/// QC every sample against its cohort (reads and assemblies are compared within their own kind).
pub fn assess(samples: &[SampleSketch], scaled: u32, p: &QcParams) -> Vec<SampleQc> {
    let mut out: Vec<SampleQc> = Vec::with_capacity(samples.len());
    let mut shoulder_raw: Vec<Option<f64>> = Vec::new();
    for s in samples {
        let reads = s.kind == "reads";
        let (cutoff, cov) = if reads {
            let prof = profile(&s.spectrum.abunds);
            (
                auto_min_abund(&s.spectrum.abunds),
                prof.map(|c| (c.kmer, c.basis)),
            )
        } else {
            (1, None)
        };
        let kept = s.spectrum.abunds.iter().filter(|&&a| a >= cutoff).count();
        let (kmer_cov, basis) = match cov {
            Some((c, b)) => (Some(c), b),
            None => (None, "none"),
        };
        let shoulder = match kmer_cov {
            Some(c) if reads && c >= 6.0 => shoulder_features(&s.spectrum.abunds, c).map(|x| x.0),
            _ => None,
        };
        shoulder_raw.push(shoulder);
        out.push(SampleQc {
            id: s.id.clone(),
            kind: s.kind.clone(),
            reads: s.counts.reads,
            bases: s.counts.bases,
            n_hashes: s.spectrum.hashes.len(),
            cutoff,
            kept,
            est_mb: kept as f64 * scaled as f64 / 1e6,
            kmer_cov,
            cov_basis: basis,
            shoulder,
            shoulder_z: None,
            size_z: None,
            flags: Vec::new(),
        });
    }
    // coverage
    for q in out.iter_mut().filter(|q| q.kind == "reads") {
        match q.kmer_cov {
            Some(c) if c >= p.min_cov && q.cov_basis == "peak" => {}
            Some(c) if c >= p.min_cov => {} // Poisson basis at or above the minimum: no peak yet, but enough coverage
            Some(c) => q.flags.push(Flag {
                code: "low_coverage",
                detail: format!("k-mer coverage {c:.1} is below {:.1}", p.min_cov),
            }),
            None => q.flags.push(Flag {
                code: "low_coverage",
                detail: "too few repeated k-mers to estimate coverage".into(),
            }),
        }
    }
    // size, within each kind, against the samples that have enough coverage to be sized
    for kind in ["reads", "assembly"] {
        let ok: Vec<usize> = out
            .iter()
            .enumerate()
            .filter(|(_, q)| q.kind == kind && !q.has("low_coverage"))
            .map(|(i, _)| i)
            .collect();
        if ok.len() < p.min_cohort {
            continue;
        }
        let sizes: Vec<f64> = ok.iter().map(|&i| out[i].est_mb).collect();
        let med = median(&sizes);
        let floor = p.size_floor * med;
        for &i in &ok {
            let z = robust_z(out[i].est_mb, &sizes, floor);
            out[i].size_z = Some(z);
            if z.abs() > p.size_z {
                let (code, what) = if z > 0.0 {
                    (
                        "size_high",
                        "larger than the cohort: a second genome or extra content",
                    )
                } else {
                    ("size_low", "smaller than the cohort: content missing")
                };
                // small genomes (a few Mb or less) need more digits than large ones
                let digits = if med < 5.0 { 2 } else { 1 };
                let detail = format!(
                    "estimated size {:.digits$} Mb vs cohort median {:.digits$} Mb (z {z:+.1}); {what}",
                    out[i].est_mb, med
                );
                out[i].flags.push(Flag { code, detail });
            }
        }
        let _ = mad_sd(&sizes);
    }
    // abundance shoulder, standardised within coverage bins
    for b in 0..COV_BINS.len() - 1 {
        let idx: Vec<usize> = out
            .iter()
            .enumerate()
            .filter(|(i, q)| {
                shoulder_raw[*i].is_some()
                    && !q.has("size_low")
                    && cov_bin(q.kmer_cov.unwrap_or(0.0)) == b
            })
            .map(|(i, _)| i)
            .collect();
        if idx.len() < p.min_cohort {
            continue;
        }
        let vals: Vec<f64> = idx.iter().map(|&i| shoulder_raw[i].unwrap()).collect();
        for &i in &idx {
            let z = robust_z(shoulder_raw[i].unwrap(), &vals, 0.001);
            out[i].shoulder_z = Some(z);
            if z > p.spectrum_z {
                let c = out[i].kmer_cov.unwrap_or(0.0);
                out[i].flags.push(Flag {
                    code: "abundance_shoulder",
                    detail: format!(
                        "k-mers at 0.1-0.6x the main coverage ({c:.0}x) exceed the cohort (z {z:+.1}): a second organism at several-fold coverage"
                    ),
                });
            }
        }
    }
    out
}

/// The hashes of a sample kept after its abundance cutoff.
pub fn kept_hashes(s: &SampleSketch, cutoff: u32) -> Vec<u64> {
    s.spectrum
        .hashes
        .iter()
        .zip(&s.spectrum.abunds)
        .filter(|(_, &a)| a >= cutoff)
        .map(|(&h, _)| h)
        .collect()
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::sketch_io::{Counts, Spectrum};

    struct Rng(u64);
    impl Rng {
        fn next(&mut self) -> f64 {
            self.0 = self.0.wrapping_add(0x9E37_79B9_7F4A_7C15);
            let mut z = self.0;
            z = (z ^ (z >> 30)).wrapping_mul(0xBF58_476D_1CE4_E5B9);
            z = (z ^ (z >> 27)).wrapping_mul(0x94D0_49BB_1331_11EB);
            ((z ^ (z >> 31)) >> 11) as f64 / (1u64 << 53) as f64
        }
        fn poisson(&mut self, l: f64) -> u32 {
            if l < 30.0 {
                let (limit, mut k, mut p) = ((-l).exp(), 0u32, 1.0);
                loop {
                    p *= self.next();
                    if p <= limit {
                        return k;
                    }
                    k += 1;
                }
            }
            let (u1, u2) = (self.next().max(1e-12), self.next());
            let z = (-2.0 * u1.ln()).sqrt() * (2.0 * std::f64::consts::PI * u2).cos();
            (l + l.sqrt() * z).round().max(0.0) as u32
        }
    }

    /// `g` genome hashes at Poisson(c), optional contaminant (hashes, coverage), `err` singletons.
    fn fake(
        id: &str,
        g: usize,
        c: f64,
        contam: Option<(usize, f64)>,
        err: usize,
        seed: u64,
    ) -> SampleSketch {
        let mut r = Rng(seed);
        let mut abunds: Vec<u32> = (0..g).map(|_| r.poisson(c)).filter(|&a| a > 0).collect();
        if let Some((n, c2)) = contam {
            abunds.extend((0..n).map(|_| r.poisson(c2)).filter(|&a| a > 0));
        }
        abunds.extend(std::iter::repeat_n(1, err));
        let hashes: Vec<u64> = (0..abunds.len() as u64).collect();
        SampleSketch {
            id: id.into(),
            kind: "reads".into(),
            files: vec![],
            spectrum: Spectrum { hashes, abunds },
            counts: Counts {
                reads: 1,
                bases: 1,
                non_acgt: 0,
            },
            halves: None,
        }
    }

    fn cohort() -> Vec<SampleSketch> {
        (0..24)
            .map(|i| {
                fake(
                    &format!("c{i}"),
                    20_000,
                    18.0 + (i % 6) as f64 * 5.0,
                    None,
                    40_000,
                    100 + i,
                )
            })
            .collect()
    }

    fn flags_of<'a>(q: &'a [SampleQc], id: &str) -> Vec<&'a str> {
        q.iter()
            .find(|x| x.id == id)
            .unwrap()
            .flags
            .iter()
            .map(|f| f.code)
            .collect()
    }

    #[test]
    fn a_clean_cohort_raises_no_flags() {
        let q = assess(&cohort(), 200, &QcParams::default());
        let flagged: Vec<_> = q
            .iter()
            .filter(|x| x.flagged())
            .map(|x| x.id.clone())
            .collect();
        assert!(flagged.is_empty(), "{flagged:?}");
        assert!(
            q.iter().all(|x| (x.est_mb - 4.0).abs() < 0.3),
            "{:?}",
            q[0].est_mb
        );
    }

    #[test]
    fn planted_problems_are_flagged_and_clean_samples_are_not() {
        let mut s = cohort();
        s.push(fake("lowcov", 20_000, 0.8, None, 60_000, 7));
        s.push(fake("half", 10_000, 30.0, None, 40_000, 8)); // half the content
        s.push(fake("mixed", 20_000, 40.0, Some((20_000, 10.0)), 40_000, 9)); // a second genome at 10x
        s.push(fake("heavy", 20_000, 12.0, Some((20_000, 6.0)), 40_000, 10)); // a second genome at half the coverage
        let q = assess(&s, 200, &QcParams::default());
        assert!(flags_of(&q, "lowcov").contains(&"low_coverage"));
        assert!(
            flags_of(&q, "half").contains(&"size_low"),
            "{:?}",
            flags_of(&q, "half")
        );
        // a second genome as large as the host: the size test or the shoulder test must catch it (which one depends on
        // where the coverage peak is found; the shoulder statistic itself is tested directly below)
        assert!(
            flags_of(&q, "mixed")
                .iter()
                .any(|c| matches!(*c, "size_high" | "abundance_shoulder")),
            "{:?}",
            flags_of(&q, "mixed")
        );
        assert!(
            flags_of(&q, "heavy")
                .iter()
                .any(|c| matches!(*c, "size_high" | "abundance_shoulder")),
            "{:?}",
            flags_of(&q, "heavy")
        );
        for i in 0..24 {
            assert!(
                flags_of(&q, &format!("c{i}")).is_empty(),
                "c{i}: {:?}",
                flags_of(&q, &format!("c{i}"))
            );
        }
    }

    #[test]
    fn light_contamination_is_not_assessable() {
        // a contaminant at 1x next to a 40x genome looks like errors: no flag, by design
        let mut s = cohort();
        s.push(fake("light", 20_000, 40.0, Some((20_000, 1.0)), 40_000, 11));
        let q = assess(&s, 200, &QcParams::default());
        assert!(!flags_of(&q, "light").contains(&"abundance_shoulder"));
    }

    #[test]
    fn a_tiny_cohort_gets_no_relative_flags() {
        let s: Vec<SampleSketch> = (0..4)
            .map(|i| fake(&format!("c{i}"), 20_000, 30.0, None, 40_000, i))
            .collect();
        let mut s2 = s;
        s2.push(fake("half", 5_000, 30.0, None, 40_000, 99));
        assert!(assess(&s2, 200, &QcParams::default())
            .iter()
            .all(|x| x.size_z.is_none()));
    }

    #[test]
    fn shoulder_feature_grows_with_a_second_organism() {
        let mut r = Rng(5);
        let clean: Vec<u32> = (0..20_000)
            .map(|_| r.poisson(40.0))
            .filter(|&a| a > 0)
            .collect();
        let mut mixed = clean.clone();
        mixed.extend((0..20_000).map(|_| r.poisson(10.0)).filter(|&a| a > 0));
        let (a, _) = shoulder_features(&clean, 40.0).unwrap();
        let (b, _) = shoulder_features(&mixed, 40.0).unwrap();
        assert!(b > a + 0.3, "{a} vs {b}");
        assert!(shoulder_features(&clean, 3.0).is_none());
    }
}
