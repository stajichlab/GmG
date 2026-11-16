//! Check one sketched sample against the genome it should cover: was it sampled evenly?
//!
//! `popseq fetch` can take 25x of a run as contiguous windows. If the run is sorted by genome position, that covers only
//! `25 / depth` of the genome, at a coverage far above 25x (phase 1 pilot: two of ten *C. auris* runs came out as a 2.2 and a
//! 2.3 Mb genome at 106x, against 12.2 Mb at 12 to 20x for the other eight). The audit flags the size afterwards; this check
//! gives the same answer for one sample straight after sketching, with an exit code a caller can branch on, so the run can be
//! resampled (every Nth read) before it enters the set.
//!
//! Two rules, either of which fails the sample (`Undersampled`, exit 75: resampling by stride may help):
//! * the estimated size is below `min_size_ratio` of the genome size (default 0.5): the reads miss part of the genome;
//! * the k-mer coverage is above `max_cov_ratio` times the k-mer coverage the kept reads should give (default 2.5): the reads
//!   pile up on a part of the genome. Needs the expected coverage. The expectation is in bases (kept bases / genome size); it
//!   is converted to k-mers with the mean read length `L`, since a read of `L` bases holds `L - k + 1` k-mers.
//!
//! Two situations in which the rules cannot say whether the sampling was even, so the verdict is `NotAssessed` (exit 65) and
//! resampling would not change it:
//! * the k-mer coverage is below `min_kmer_cov` (default 3) or not found: the abundance cutoff then removes real k-mers (at
//!   a mean of 3 per k-mer 20% of genome k-mers are seen once, under a Poisson model), so the size is an underestimate whatever
//!   the sampling, and a run with mostly unique k-mers is not a deep read set of this genome (another organism, or a mixture);
//! * the estimated size is above `max_size_ratio` times the genome size (default 5): the sample is not this genome, or the
//!   genome size given is wrong.
//!
//! A failing sample is only called `NotAssessed` when it would have failed a rule; a sample that passes stays `Ok`
//! (except one far larger than the genome).
//!
//! The estimators are the audit's: kept hashes (abundance at or above the automatic cutoff) times `scaled`, and the k-mer
//! coverage of the abundance histogram. The thresholds come from one species (n = 10: 0.97 to 1.01 and 0.5 to 0.8 for sound
//! samples, 0.18 and 4.2 for the biased ones) and have since been compared with 29 sound *P. oryzae* runs (size ratio 0.90 to
//! 1.07) and 299 *B. bruxellensis* runs (size ratio 0.5 to 2.35, median 1.09); beyond those they are UNMEASURED.

use crate::sketch_io::SampleSketch;
use gmg_core::coverage::{auto_min_abund, profile};

/// Exit code of `popsketch check` for an undersampled sample (sysexits EX_TEMPFAIL: "try again", here with another sampling).
pub const EXIT_UNDERSAMPLED: i32 = 75;

/// Exit code of `popsketch check` when the rules cannot judge the sampling (too little k-mer coverage, or not this genome):
/// sysexits EX_DATAERR, "the input data is incorrect in some way". Resampling would not change the verdict.
pub const EXIT_NOT_ASSESSED: i32 = 65;

/// Exit code of `popsketch check` for an error (bad argument, no such sketch, unreadable folder).
pub const EXIT_ERROR: i32 = 2;

#[derive(Clone, Debug)]
pub struct CheckParams {
    pub min_size_ratio: f64,
    pub max_cov_ratio: f64,
    /// below this k-mer coverage a failing sample is `NotAssessed`, not `Undersampled`
    pub min_kmer_cov: f64,
    /// above this multiple of the genome size a failing sample is `NotAssessed` and a passing one is too (it is not this genome)
    pub max_size_ratio: f64,
}

impl Default for CheckParams {
    fn default() -> Self {
        CheckParams {
            min_size_ratio: 0.5,
            max_cov_ratio: 2.5,
            min_kmer_cov: 3.0,
            max_size_ratio: 5.0,
        }
    }
}

#[derive(Clone, Debug, PartialEq, Eq)]
pub enum Verdict {
    Ok,
    /// the reads miss part of the genome or pile up on part of it: resample (exit 75)
    Undersampled,
    /// the rules cannot judge (exit 65): too little k-mer coverage, or the sample is not this genome
    NotAssessed,
}

#[derive(Clone, Debug)]
pub struct SizeCheck {
    pub id: String,
    pub est_size_mb: f64,
    pub genome_mb: f64,
    pub size_ratio: f64,
    pub kmer_cov: Option<f64>,
    pub expected_coverage: Option<f64>,
    /// the expectation in k-mers: `expected_coverage` x (L - k + 1) / L for the mean read length L
    pub expected_kmer_coverage: Option<f64>,
    /// mean read length (bases / reads), when the sketch knows it
    pub read_length: Option<f64>,
    pub cov_ratio: Option<f64>,
    pub verdict: Verdict,
    pub reasons: Vec<String>,
}

/// Columns 1 to 9 are stable (callers cut `reasons` as field 9); new columns are added at the end.
pub const TSV_HEADER: &str = "id\test_size_mb\tgenome_mb\tsize_ratio\tkmer_coverage\texpected_coverage\tcoverage_ratio\tverdict\treasons\tread_length\texpected_kmer_coverage";

impl SizeCheck {
    pub fn tsv(&self) -> String {
        let f = |x: Option<f64>| x.map_or(String::new(), |v| format!("{v:.2}"));
        format!(
            "{}\t{:.2}\t{:.2}\t{:.3}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
            self.id,
            self.est_size_mb,
            self.genome_mb,
            self.size_ratio,
            f(self.kmer_cov),
            f(self.expected_coverage),
            f(self.cov_ratio),
            match self.verdict {
                Verdict::Ok => "ok",
                Verdict::Undersampled => "undersampled",
                Verdict::NotAssessed => "not_assessed",
            },
            self.reasons.join("; "),
            f(self.read_length),
            f(self.expected_kmer_coverage)
        )
    }
}

/// Reject inputs that would make a verdict meaningless: a genome size or coverage that is zero, negative, infinite or NaN
/// (a ratio against zero is infinite, and `NaN < 0.5` is false, so such a sample would pass), and thresholds outside their
/// range. The caller maps an error to exit code 2.
pub fn validate(
    genome_mb: f64,
    expected_coverage: Option<f64>,
    p: &CheckParams,
) -> Result<(), String> {
    let pos = |name: &str, x: f64| {
        if x.is_finite() && x > 0.0 {
            Ok(())
        } else {
            Err(format!("{name} must be a positive number, got {x}"))
        }
    };
    pos("--genome-mb", genome_mb)?;
    if let Some(c) = expected_coverage {
        pos("--coverage", c)?;
    }
    pos("--max-cov-ratio", p.max_cov_ratio)?;
    if !(p.min_kmer_cov.is_finite() && p.min_kmer_cov >= 0.0) {
        return Err(format!(
            "--min-kmer-cov must be a number of 0 or more, got {}",
            p.min_kmer_cov
        ));
    }
    if !(p.max_size_ratio.is_finite() && p.max_size_ratio > 1.0) {
        return Err(format!(
            "--max-size-ratio must be above 1, got {}",
            p.max_size_ratio
        ));
    }
    if !(p.min_size_ratio.is_finite() && p.min_size_ratio > 0.0 && p.min_size_ratio <= 1.0) {
        return Err(format!(
            "--min-size-ratio must be in (0, 1], got {}",
            p.min_size_ratio
        ));
    }
    Ok(())
}

/// `expected_coverage`: the x the kept reads should give if sampled evenly (bases kept / genome size); None skips the
/// coverage rule. `ksize` converts it to k-mers with the mean read length of the sketch's counts (no conversion when the
/// sketch has no read counts).
pub fn check(
    s: &SampleSketch,
    scaled: u32,
    ksize: u32,
    genome_mb: f64,
    expected_coverage: Option<f64>,
    p: &CheckParams,
) -> SizeCheck {
    let abunds = &s.spectrum.abunds;
    let is_reads = s.kind == "reads";
    let cutoff = if is_reads { auto_min_abund(abunds) } else { 1 };
    let kept = abunds.iter().filter(|&&a| a >= cutoff).count();
    let est_size_mb = kept as f64 * scaled as f64 / 1e6;
    let kmer_cov = if is_reads {
        profile(abunds).map(|c| c.kmer)
    } else {
        None
    };
    let size_ratio = est_size_mb / genome_mb;
    // a read of L bases holds L - k + 1 k-mers, so k-mer coverage is base coverage times (L - k + 1) / L
    let read_length = (s.counts.reads > 0).then(|| s.counts.bases as f64 / s.counts.reads as f64);
    let kmer_factor = match read_length {
        Some(l) if l >= ksize as f64 => (l - ksize as f64 + 1.0) / l,
        _ => 1.0,
    };
    let expected_kmer_coverage = expected_coverage.map(|e| e * kmer_factor);
    let cov_ratio = match (kmer_cov, expected_kmer_coverage) {
        (Some(k), Some(e)) if e > 0.0 => Some(k / e),
        _ => None,
    };
    let mut reasons = Vec::new();
    if size_ratio < p.min_size_ratio {
        reasons.push(format!(
            "estimated size {est_size_mb:.1} Mb is {:.0}% of the genome ({genome_mb:.1} Mb), below {:.0}%",
            100.0 * size_ratio,
            100.0 * p.min_size_ratio
        ));
    }
    if let (Some(r), Some(k), Some(e)) = (cov_ratio, kmer_cov, expected_kmer_coverage) {
        if r > p.max_cov_ratio {
            reasons.push(format!(
                "k-mer coverage {k:.0}x is {r:.1} times the expected {e:.0}x, above {:.1} times",
                p.max_cov_ratio
            ));
        }
    }
    let mut verdict = if reasons.is_empty() {
        Verdict::Ok
    } else {
        Verdict::Undersampled
    };
    // the situations the two rules cannot judge (see the module doc)
    if size_ratio > p.max_size_ratio {
        verdict = Verdict::NotAssessed;
        reasons.push(format!(
            "estimated size {est_size_mb:.1} Mb is {size_ratio:.1} times the genome ({genome_mb:.1} Mb), above {:.1} times: not this genome, or the genome size is wrong",
            p.max_size_ratio
        ));
    } else if is_reads
        && verdict == Verdict::Undersampled
        && kmer_cov.is_none_or(|k| k < p.min_kmer_cov)
    {
        verdict = Verdict::NotAssessed;
        reasons.push(match kmer_cov {
            Some(k) => format!(
                "k-mer coverage {k:.1}x is below {:.1}x: the size is not reliable and resampling cannot fix that",
                p.min_kmer_cov
            ),
            None => format!(
                "no k-mer coverage peak found (below {:.1}x or a mixture): the size is not reliable and resampling cannot fix that",
                p.min_kmer_cov
            ),
        });
    }
    SizeCheck {
        id: s.id.clone(),
        est_size_mb,
        genome_mb,
        size_ratio,
        kmer_cov,
        expected_coverage,
        expected_kmer_coverage,
        read_length,
        cov_ratio,
        verdict,
        reasons,
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::sketch_io::{Counts, Spectrum};

    /// `n` genome hashes at abundances around `c` (a fixed spread, no randomness), plus `err` singletons.
    fn sample(id: &str, n: usize, c: u32, err: usize, kind: &str) -> SampleSketch {
        let mut abunds: Vec<u32> = (0..n)
            .map(|i| c.saturating_sub(3) + (i as u32 % 7))
            .collect();
        abunds.extend(std::iter::repeat_n(1, err));
        let hashes: Vec<u64> = (0..abunds.len() as u64).collect();
        SampleSketch {
            id: id.into(),
            kind: kind.into(),
            files: vec![],
            spectrum: Spectrum { hashes, abunds },
            counts: Counts::default(),
            halves: None,
        }
    }

    #[test]
    fn an_evenly_sampled_run_passes() {
        // a 12 Mb genome at scaled 200 is 60,000 hashes; 25x of it gives about 20x k-mer coverage after errors
        let s = sample("good", 60_000, 20, 30_000, "reads");
        let r = check(&s, 200, 21, 12.25, Some(25.0), &CheckParams::default());
        assert_eq!(r.verdict, Verdict::Ok, "{:?}", r.reasons);
        assert!((r.est_size_mb - 12.0).abs() < 0.5, "{}", r.est_size_mb);
        assert!(r.size_ratio > 0.9 && r.cov_ratio.unwrap() < 1.0);
    }

    #[test]
    fn contiguous_windows_of_a_position_sorted_run_fail_on_size_and_on_coverage() {
        // the pilot's biased runs: 18% of the genome at about 106x
        let s = sample("sorted", 11_500, 106, 30_000, "reads");
        let r = check(&s, 200, 21, 12.25, Some(25.0), &CheckParams::default());
        assert_eq!(r.verdict, Verdict::Undersampled);
        assert_eq!(r.reasons.len(), 2, "{:?}", r.reasons);
        assert!(r.size_ratio < 0.25 && r.cov_ratio.unwrap() > 4.0);
    }

    #[test]
    fn either_rule_alone_fails_the_sample_and_the_coverage_rule_is_skipped_without_an_expectation()
    {
        // right size, too much coverage: reads pile up (for example a duplicated library)
        let hot = sample("hot", 60_000, 80, 30_000, "reads");
        let r = check(&hot, 200, 21, 12.25, Some(25.0), &CheckParams::default());
        assert_eq!(r.verdict, Verdict::Undersampled);
        assert_eq!(r.reasons.len(), 1);
        // no expected coverage: only the size rule runs
        assert_eq!(
            check(&hot, 200, 21, 12.25, None, &CheckParams::default()).verdict,
            Verdict::Ok
        );
        // right coverage, too small: part of the genome is missing
        let small = sample("small", 20_000, 20, 30_000, "reads");
        let r = check(&small, 200, 21, 12.25, Some(25.0), &CheckParams::default());
        assert_eq!((r.verdict, r.reasons.len()), (Verdict::Undersampled, 1));
    }

    #[test]
    fn thresholds_are_options_and_an_assembly_has_no_coverage_rule() {
        let s = sample("g", 40_000, 20, 30_000, "reads"); // 8 Mb of a 12.25 Mb genome: ratio 0.65
        assert_eq!(
            check(&s, 200, 21, 12.25, Some(25.0), &CheckParams::default()).verdict,
            Verdict::Ok
        );
        let strict = CheckParams {
            min_size_ratio: 0.8,
            ..CheckParams::default()
        };
        assert_eq!(
            check(&s, 200, 21, 12.25, Some(25.0), &strict).verdict,
            Verdict::Undersampled
        );
        let a = sample("asm", 60_000, 1, 0, "assembly");
        let r = check(&a, 200, 21, 12.0, Some(25.0), &CheckParams::default());
        assert!(r.kmer_cov.is_none() && r.cov_ratio.is_none());
        assert_eq!(r.verdict, Verdict::Ok);
    }

    #[test]
    fn a_zero_nan_or_negative_genome_size_or_coverage_is_rejected_not_passed() {
        let p = CheckParams::default();
        assert!(validate(12.0, Some(25.0), &p).is_ok());
        assert!(validate(12.0, None, &p).is_ok());
        for g in [0.0, -1.0, f64::NAN, f64::INFINITY] {
            assert!(validate(g, Some(25.0), &p).is_err(), "genome {g}");
        }
        for c in [0.0, -5.0, f64::NAN, f64::INFINITY] {
            assert!(validate(12.0, Some(c), &p).is_err(), "coverage {c}");
        }
        let bad = |min: f64, max: f64| {
            validate(
                12.0,
                None,
                &CheckParams {
                    min_size_ratio: min,
                    max_cov_ratio: max,
                    ..CheckParams::default()
                },
            )
        };
        assert!(bad(0.0, 2.5).is_err() && bad(1.5, 2.5).is_err() && bad(f64::NAN, 2.5).is_err());
        assert!(bad(0.5, 0.0).is_err() && bad(0.5, f64::NAN).is_err());
        assert!(bad(1.0, 2.5).is_ok());
    }

    /// like `sample`, with read counts so that the mean read length is `len`
    fn with_length(mut s: SampleSketch, len: u64) -> SampleSketch {
        s.counts = Counts {
            reads: 1_000_000,
            bases: 1_000_000 * len,
            non_acgt: 0,
        };
        s
    }

    #[test]
    fn the_expected_coverage_is_converted_from_bases_to_k_mers_with_the_read_length() {
        // 36 bp reads at k = 21 hold 16 k-mers each: 16/36 = 0.44 of the base coverage. A pile-up at 40x k-mer coverage
        // against 25x of bases is 1.6 times in bases (passed) but 3.6 times the 11x k-mers it should give
        let hot = with_length(sample("hot36", 60_000, 40, 30_000, "reads"), 36);
        let r = check(&hot, 200, 21, 12.25, Some(25.0), &CheckParams::default());
        assert_eq!(r.verdict, Verdict::Undersampled, "{:?}", r.reasons);
        assert!((r.expected_kmer_coverage.unwrap() - 25.0 * 16.0 / 36.0).abs() < 1e-9);
        assert!(r.cov_ratio.unwrap() > 3.0, "{:?}", r.cov_ratio);
        // an evenly sampled 36 bp run gives about 11x and passes
        let sound = with_length(sample("sound36", 60_000, 11, 30_000, "reads"), 36);
        let r = check(&sound, 200, 21, 12.25, Some(25.0), &CheckParams::default());
        assert_eq!(r.verdict, Verdict::Ok, "{:?}", r.reasons);
        assert!(r.cov_ratio.unwrap() < 1.5);
        // without read counts nothing is converted (the old behaviour)
        let r = check(&hot, 200, 21, 12.25, Some(25.0), &CheckParams::default());
        assert!(r.read_length.is_none() || r.read_length == Some(36.0));
        let bare = sample("bare", 60_000, 40, 30_000, "reads");
        let r = check(&bare, 200, 21, 12.25, Some(25.0), &CheckParams::default());
        assert_eq!(r.expected_kmer_coverage, Some(25.0));
    }

    #[test]
    fn a_failing_sample_with_too_little_k_mer_coverage_is_not_assessed_not_undersampled() {
        // 8,000 genome hashes at about 2x plus errors: the size is far below the genome, but at 2x the abundance cutoff
        // removes real k-mers, so evenly sampled reads would look the same; resampling cannot help
        let mut low = sample("low", 8_000, 2, 30_000, "reads");
        low.spectrum.abunds = (0..8_000u32)
            .map(|i| 1 + i % 3) // 1, 2, 3: a mean of 2
            .chain(std::iter::repeat_n(1, 30_000))
            .collect();
        let r = check(&low, 200, 21, 12.25, Some(25.0), &CheckParams::default());
        assert!(
            r.kmer_cov.is_none_or(|k| k < 3.0),
            "k-mer coverage {:?}",
            r.kmer_cov
        );
        assert!(r.size_ratio < 0.5);
        assert_eq!(r.verdict, Verdict::NotAssessed, "{:?}", r.reasons);
        assert!(
            r.reasons
                .iter()
                .any(|x| x.contains("resampling cannot fix")),
            "{:?}",
            r.reasons
        );
        // the same sample with a good size is not touched by the new rule
        let mut fine = sample("fine", 60_000, 2, 0, "reads");
        fine.spectrum.abunds = (0..60_000u32).map(|i| 1 + i % 3).collect();
        let r = check(&fine, 200, 21, 12.25, None, &CheckParams::default());
        assert_eq!(r.verdict, Verdict::Ok, "{:?}", r.reasons);
        // and the threshold is a parameter
        let off = CheckParams {
            min_kmer_cov: 0.0,
            ..CheckParams::default()
        };
        assert_eq!(
            check(&low, 200, 21, 12.25, Some(25.0), &off).verdict,
            Verdict::Undersampled
        );
    }

    #[test]
    fn a_sample_far_larger_than_the_genome_is_not_assessed_even_though_no_rule_failed() {
        // 600,000 hashes at scaled 200 is 120 Mb against a 12.25 Mb genome: another organism, or a wrong genome size
        let big = sample("big", 600_000, 20, 0, "reads");
        let r = check(&big, 200, 21, 12.25, Some(25.0), &CheckParams::default());
        assert_eq!(r.verdict, Verdict::NotAssessed);
        assert!(
            r.reasons.iter().any(|x| x.contains("not this genome")),
            "{:?}",
            r.reasons
        );
        // a heterozygous diploid at 2.4 times (the largest of 299 B. bruxellensis runs) passes
        let diploid = sample("dip", 150_000, 20, 0, "reads"); // 30 Mb of 12.25
        assert_eq!(
            check(
                &diploid,
                200,
                21,
                12.25,
                Some(25.0),
                &CheckParams::default()
            )
            .verdict,
            Verdict::Ok
        );
    }

    #[test]
    fn the_tsv_line_has_the_header_fields_and_scaled_is_read_from_params() {
        let s = sample("x", 60_000, 20, 30_000, "reads");
        let line = check(&s, 200, 21, 12.25, Some(25.0), &CheckParams::default()).tsv();
        assert_eq!(
            line.split('\t').count(),
            TSV_HEADER.split('\t').count(),
            "{line}"
        );
        assert_eq!(
            crate::sketch_io::ksize_scaled_from_tag(
                "k=21;scaled=200;max_bases=0;r1_only=false;target_cov=none\n"
            ),
            Some((21, 200))
        );
        assert_eq!(crate::sketch_io::ksize_scaled_from_tag("k=21"), None);
    }
}
