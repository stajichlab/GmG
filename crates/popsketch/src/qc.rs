//! Cohort-relative QC of sketched samples. Nothing here needs a reference: each sample is compared with the rest of the
//! cohort, and every flag carries its reason and numbers.
//!
//! * `low_coverage`: the k-mer coverage from the abundance histogram is below `min_cov`, or there is no coverage peak.
//! * `size_high` / `size_low`: the estimated size (hashes above the abundance cutoff x scaled) is far from the cohort median
//!   (robust z, with a floor so a tight cohort does not flag everything). Heavy contamination (a second genome) shows as
//!   `size_high`; a sample missing half its content as `size_low`.
//! * `abundance_shoulder`: hashes at low-but-real abundance (0.1c-0.4c) beyond what a single genome predicts, standardised
//!   within coverage bins. Catches a contaminant whose own coverage is several-fold; it cannot see lighter contamination.
//!   The window stops at 0.4c on purpose: a heterozygous diploid has a peak at about half the main coverage (on 334 Bd read
//!   sets the shoulder sat at 0.57c in flagged and unflagged samples alike), so a contaminant near half the coverage cannot
//!   be told from heterozygosity and is not looked for.

use crate::sketch_io::SampleSketch;
use crate::stats::{median, robust_z};
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
    /// observations that do not exclude the sample (an abundance shoulder in a cohort where it is widespread)
    pub notes: Vec<Flag>,
}

impl SampleQc {
    pub fn flagged(&self) -> bool {
        !self.flags.is_empty()
    }
    pub fn has(&self, code: &str) -> bool {
        self.flags.iter().any(|f| f.code == code)
    }
}

/// (shoulder fraction, heterozygous-window fraction) of an abundance histogram with k-mer coverage `c`.
///
/// Shoulder: hashes with abundance in [lo, hi), lo = max(3, 0.1c), hi = ceil(0.4c) (below the heterozygous peak at about
/// 0.5c), as a fraction of the genome mass (hashes at or above 0.6c). Heterozygous window: the share of the genome mass at 0.4c
/// to 0.65c, where a heterozygous peak sits. Heterozygosity raises both (the second through the lower tail of its peak); a
/// contaminant raises the shoulder without raising the window. Needs c of about 12 or more for the shoulder window to hold more
/// than a couple of bins; None when it does not, or when no hash is at or above 0.6c.
pub fn window_fractions(abunds: &[u32], c: f64) -> Option<(f64, f64)> {
    let mass_from = (0.6 * c).ceil() as usize;
    let hi = (0.4 * c).ceil() as usize;
    let lo = ((0.1 * c).floor() as usize).max(3);
    if hi <= lo + 1 {
        return None;
    }
    const BINS: usize = 4001;
    let mut h = vec![0f64; BINS];
    for &a in abunds {
        h[(a as usize).min(BINS - 1)] += 1.0;
    }
    let mass: f64 = h[mass_from.min(BINS - 1)..].iter().sum();
    if mass <= 0.0 {
        return None;
    }
    let obs: f64 = h[lo..hi].iter().sum();
    let het_lo = (0.4 * c).ceil() as usize;
    let het_hi = (0.65 * c).ceil() as usize;
    let het: f64 = h[het_lo.min(BINS - 1)..het_hi.min(BINS - 1)].iter().sum();
    Some((obs / mass, het / mass))
}

/// Up to this many points the fit takes every pair; above it a random sample of pairs.
const SLOPE_ALL_PAIRS_MAX_N: usize = 90;
const SLOPE_SAMPLED_PAIRS: usize = 4000;

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
    let mut shoulder_raw: Vec<Option<(f64, f64)>> = Vec::new();
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
            Some(c) if reads && c >= 6.0 => window_fractions(&s.spectrum.abunds, c),
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
            shoulder: shoulder.map(|x| x.0),
            shoulder_z: None,
            size_z: None,
            flags: Vec::new(),
            notes: Vec::new(),
        });
    }
    // coverage
    for q in out.iter_mut().filter(|q| q.kind == "reads") {
        match q.kmer_cov {
            // at or above the minimum, from a peak or (no peak yet) from the Poisson basis: enough coverage
            Some(c) if c >= p.min_cov => {}
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
    }
    // abundance shoulder: the part of the shoulder that the heterozygous-window mass does not explain, standardised within
    // coverage bins. A heterozygous (or polyploid) genome raises both; a contaminant raises only the shoulder.
    let cand: Vec<usize> = (0..out.len())
        .filter(|&i| shoulder_raw[i].is_some() && !out[i].has("size_low"))
        .collect();
    if cand.len() >= p.min_cohort {
        let pts: Vec<(f64, f64)> = cand.iter().map(|&i| shoulder_raw[i].unwrap()).collect();
        // robust line s ~ a + b h: median of pairwise slopes (all pairs up to PAIRS_ALL points, otherwise PAIRS_SAMPLED random
        // pairs from a fixed seed), then the median intercept
        let n = pts.len();
        let mut slopes = Vec::new();
        let mut add = |i: usize, j: usize| {
            if i != j && (pts[i].1 - pts[j].1).abs() > 1e-3 {
                slopes.push((pts[i].0 - pts[j].0) / (pts[i].1 - pts[j].1));
            }
        };
        if n <= SLOPE_ALL_PAIRS_MAX_N {
            for i in 0..n {
                for j in i + 1..n {
                    add(i, j);
                }
            }
        } else {
            let mut rng = crate::rng::SplitMix(0x5107_E5E1);
            for _ in 0..SLOPE_SAMPLED_PAIRS {
                add(rng.below(n), rng.below(n));
            }
        }
        let slope = if slopes.is_empty() {
            0.0
        } else {
            median(&slopes)
        };
        let icpt = median(&pts.iter().map(|(s, h)| s - slope * h).collect::<Vec<_>>());
        let resid: Vec<f64> = pts.iter().map(|(s, h)| s - (icpt + slope * h)).collect();
        let mut zs: Vec<Option<f64>> = vec![None; n];
        for b in 0..COV_BINS.len() - 1 {
            let idx: Vec<usize> = (0..n)
                .filter(|&k| cov_bin(out[cand[k]].kmer_cov.unwrap_or(0.0)) == b)
                .collect();
            if idx.len() < p.min_cohort {
                continue;
            }
            let vals: Vec<f64> = idx.iter().map(|&k| resid[k]).collect();
            for &k in &idx {
                zs[k] = Some(robust_z(resid[k], &vals, 0.001));
            }
        }
        let assessed = zs.iter().flatten().count();
        let hits = zs.iter().flatten().filter(|&&z| z > p.spectrum_z).count();
        // many samples with a shoulder means heterozygosity or varying ploidy, not contamination: note it, do not exclude
        let widespread = hits > 3 && hits as f64 > 0.02 * assessed as f64;
        for (k, &i) in cand.iter().enumerate() {
            let Some(z) = zs[k] else { continue };
            out[i].shoulder_z = Some(z);
            if z > p.spectrum_z {
                let c = out[i].kmer_cov.unwrap_or(0.0);
                let base = format!(
                    "k-mers at 0.1-0.4x the main coverage ({c:.0}x) exceed what the heterozygous-window mass explains (z {z:+.1})"
                );
                if widespread {
                    out[i].notes.push(Flag {
                        code: "abundance_shoulder",
                        detail: format!("{base}; {hits} of {assessed} samples ({:.1}%) show this, which points to heterozygosity or varying ploidy in this organism, not contamination, so the sample is not excluded", 100.0 * hits as f64 / assessed as f64),
                    });
                } else {
                    out[i].flags.push(Flag {
                        code: "abundance_shoulder",
                        detail: format!("{base}: a second organism at several-fold coverage"),
                    });
                }
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

    /// like `fake`, with any number of extra components (hashes, coverage): a heterozygous peak, a contaminant, a triploid third
    fn fake_with(
        id: &str,
        g: usize,
        c: f64,
        comps: &[(usize, f64)],
        err: usize,
        seed: u64,
    ) -> SampleSketch {
        let mut r = Rng(seed);
        let mut abunds: Vec<u32> = (0..g).map(|_| r.poisson(c)).filter(|&a| a > 0).collect();
        for &(n, c2) in comps {
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

    #[test]
    fn varying_heterozygosity_is_not_a_shoulder_but_a_contaminant_is() {
        // 36 samples whose heterozygous peak (at half the coverage) varies from none to large; one also carries a contaminant at a quarter
        let mut s: Vec<SampleSketch> = (0..36u64)
            .map(|i| {
                fake_with(
                    &format!("h{i}"),
                    20_000,
                    36.0 + (i % 5) as f64 * 2.0,
                    &[(((i % 6) as usize) * 1500, 20.0)],
                    40_000,
                    300 + i,
                )
            })
            .collect();
        s.push(fake_with(
            "contam",
            20_000,
            40.0,
            &[(3000, 20.0), (8000, 10.0)],
            40_000,
            999,
        ));
        let q = assess(&s, 200, &QcParams::default());
        assert!(
            flags_of(&q, "contam").contains(&"abundance_shoulder"),
            "{:?} z={:?}",
            flags_of(&q, "contam"),
            q.last().unwrap().shoulder_z
        );
        for i in 0..36 {
            let f = flags_of(&q, &format!("h{i}"));
            assert!(
                !f.contains(&"abundance_shoulder"),
                "h{i} (het {}): {f:?}",
                (i % 6) * 1500
            );
        }
        assert!(
            q.iter().all(|x| x.notes.is_empty()),
            "one hit is not widespread"
        );
    }

    #[test]
    fn a_widespread_shoulder_is_a_note_and_not_an_exclusion() {
        // a polyploid-like cohort: a minority of samples carry a component at a third of the coverage, in varying amounts
        let s: Vec<SampleSketch> = (0..36u64)
            .map(|i| {
                let extra = if i % 6 == 0 {
                    vec![(1200 + 200 * (i as usize / 6), 40.0 / 3.0)]
                } else {
                    vec![]
                };
                fake_with(&format!("p{i}"), 20_000, 40.0, &extra, 40_000, 500 + i)
            })
            .collect();
        // the rule is about how many samples hit, not the z value: a lower z puts six of 36 over it
        let q = assess(
            &s,
            200,
            &QcParams {
                spectrum_z: 2.5,
                ..QcParams::default()
            },
        );
        let noted: Vec<&SampleQc> = q.iter().filter(|x| !x.notes.is_empty()).collect();
        assert!(
            !noted.is_empty(),
            "z values: {:?}",
            q.iter()
                .map(|x| x.shoulder_z.map(|z| (z * 10.0).round() / 10.0))
                .collect::<Vec<_>>()
        );
        assert!(noted.iter().all(|x| !x.has("abundance_shoulder")
            && x.notes[0]
                .detail
                .contains("heterozygosity or varying ploidy")));
        assert!(
            noted.iter().all(|x| !x.flagged()),
            "a note must not exclude the sample"
        );
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
    fn a_heterozygous_peak_at_half_coverage_is_not_a_shoulder_but_a_contaminant_at_a_quarter_is() {
        let mut r = Rng(11);
        let host: Vec<u32> = (0..20_000)
            .map(|_| r.poisson(40.0))
            .filter(|&a| a > 0)
            .collect();
        let (clean, _) = window_fractions(&host, 40.0).unwrap();
        // 8,000 heterozygous sites: their k-mers are seen at half the coverage
        let mut het = host.clone();
        het.extend((0..8_000).map(|_| r.poisson(20.0)).filter(|&a| a > 0));
        let (with_het, _) = window_fractions(&het, 40.0).unwrap();
        // a second organism of 8,000 hashes at a quarter of the coverage
        let mut mixed = host.clone();
        mixed.extend((0..8_000).map(|_| r.poisson(10.0)).filter(|&a| a > 0));
        let (with_contaminant, _) = window_fractions(&mixed, 40.0).unwrap();
        // the lower tail of a half-coverage peak reaches into the window, but far less than a real contaminant does
        let (het_effect, contaminant_effect) = (with_het - clean, with_contaminant - clean);
        assert!(
            contaminant_effect > 0.15,
            "contaminant {with_contaminant} vs clean {clean}"
        );
        assert!(
            het_effect < 0.25 * contaminant_effect,
            "het effect {het_effect} vs contaminant effect {contaminant_effect}"
        );
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
        let (a, _) = window_fractions(&clean, 40.0).unwrap();
        let (b, _) = window_fractions(&mixed, 40.0).unwrap();
        assert!(b > a + 0.3, "{a} vs {b}");
        assert!(window_fractions(&clean, 3.0).is_none());
    }

    #[test]
    fn the_shoulder_is_defined_at_very_high_coverage_too() {
        // the Poisson comparison this replaced returned None above about 1,250x (its lambda was clamped to 600)
        let abunds: Vec<u32> = (0..20_000).map(|i| 1_480 + i % 41).collect();
        let (shoulder, het) = window_fractions(&abunds, 1_500.0).expect("defined at 1,500x");
        assert!(shoulder.abs() < 1e-12 && het.abs() < 1e-12);
    }
}
