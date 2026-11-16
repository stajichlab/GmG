//! PERMANOVA (does the group label explain distances?) and PERMDISP (do groups differ in spread?) on a distance matrix,
//! with permutations that can stay inside strata (for example study), and warnings for the ways the test misleads.
//!
//! Why both: PERMANOVA's pseudo-F also rises when groups have equal centres but unequal spread, and more so when the
//! group sizes are unequal (a country with 187 of 265 samples and one with 3). PERMDISP tests the spread alone, so a
//! significant PERMANOVA next to a significant PERMDISP is not evidence of different centres.
//!
//! Why strata: samples from one study are often near-clonal or sampled together, so they are not exchangeable across
//! labels. Permuting the label only among the samples of one stratum keeps that structure and answers "does the label
//! explain more than the study does". If the label never varies inside a stratum (country is the same as study), the
//! test has nothing to permute and returns no p-value, with a warning.
//!
//! Distances are the matrix as given (use `DistMat::one_minus_ani`). Missing labels are left out.

use crate::par::*;
use crate::rng::SplitMix;
use std::collections::BTreeMap;

/// Smallest group for which PERMDISP is reported.
pub const PERMDISP_MIN_GROUP: usize = 3;

pub struct Options {
    pub n_perm: usize,
    pub seed: u64,
}

impl Default for Options {
    fn default() -> Self {
        Options {
            n_perm: 999,
            seed: 1,
        }
    }
}

#[derive(Debug, Clone)]
pub struct GroupStat {
    pub label: String,
    pub n: usize,
    /// mean distance to the group centroid (PERMDISP)
    pub mean_dispersion: f64,
}

#[derive(Debug, Clone)]
pub struct TestResult {
    pub n: usize,
    pub n_groups: usize,
    pub f: f64,
    /// share of the total sum of squares explained by the label
    pub r2: f64,
    pub p: Option<f64>,
    pub disp_f: f64,
    pub disp_p: Option<f64>,
    pub groups: Vec<GroupStat>,
    /// share of the samples that sit in a stratum where the label varies (1 without strata)
    pub informative_fraction: f64,
    pub strata_used: bool,
    /// the sums of squares were weighted (collapsed representatives)
    pub weighted: bool,
    pub warnings: Vec<String>,
}

/// Run both tests. `d` is n x n row-major; `labels` and `strata` have length n. Returns `None` when fewer than two
/// groups with at least two samples remain.
///
/// `weights` (length n, positive; for example the number of samples each representative stands for after collapsing
/// near-identical samples) turns the sums of squares into their weighted forms: a pair counts `w_a w_b` times and a
/// group's size is the sum of its weights. With all weights 1 the result is the unweighted test. The degrees of freedom
/// still use the number of rows (not the sum of weights), so collapsed samples are not counted as independent replicates;
/// the PERMDISP centroids are weighted but its ANOVA on the distances to the centroids is not.
pub fn test(
    d: &[f64],
    n: usize,
    labels: &[Option<String>],
    strata: Option<&[Option<String>]>,
    weights: Option<&[f64]>,
    o: &Options,
) -> Option<TestResult> {
    assert_eq!(d.len(), n * n);
    assert_eq!(labels.len(), n);
    if let Some(w) = weights {
        assert_eq!(w.len(), n);
        assert!(
            w.iter().all(|&x| x > 0.0 && x.is_finite()),
            "weights must be positive"
        );
    }
    let mut warnings = Vec::new();

    // keep labelled samples whose group has at least 2 members
    let mut count: BTreeMap<&str, usize> = BTreeMap::new();
    for l in labels.iter().flatten() {
        *count.entry(l.as_str()).or_default() += 1;
    }
    let singles: Vec<&str> = count
        .iter()
        .filter(|(_, &c)| c < 2)
        .map(|(&l, _)| l)
        .collect();
    if !singles.is_empty() {
        warnings.push(format!(
            "{} group(s) with one sample were left out: {}",
            singles.len(),
            singles.join(", ")
        ));
    }
    let names: Vec<&str> = count
        .iter()
        .filter(|(_, &c)| c >= 2)
        .map(|(&l, _)| l)
        .collect();
    if names.len() < 2 {
        return None;
    }
    let idx: Vec<usize> = (0..n)
        .filter(|&i| labels[i].as_deref().is_some_and(|l| names.contains(&l)))
        .collect();
    let m = idx.len();
    let g = names.len();
    let lab: Vec<usize> = idx
        .iter()
        .map(|&i| {
            names
                .iter()
                .position(|&x| Some(x) == labels[i].as_deref())
                .unwrap()
        })
        .collect();
    let w: Vec<f64> = match weights {
        Some(w) => idx.iter().map(|&i| w[i]).collect(),
        None => vec![1.0; m],
    };
    let mut d2 = vec![0.0f64; m * m];
    for a in 0..m {
        for b in 0..m {
            let x = d[idx[a] * n + idx[b]];
            d2[a * m + b] = x * x;
        }
    }

    // strata: samples with no stratum are each their own, so they can swap with nobody
    let strata_ids: Vec<usize> = match strata {
        Some(s) => {
            assert_eq!(s.len(), n);
            let mut map: BTreeMap<String, usize> = BTreeMap::new();
            idx.iter()
                .enumerate()
                .map(|(k, &i)| match &s[i] {
                    Some(v) => {
                        let nx = map.len();
                        *map.entry(v.clone()).or_insert(nx)
                    }
                    None => 1_000_000 + k,
                })
                .collect()
        }
        None => vec![0; m],
    };
    let mut by_stratum: BTreeMap<usize, Vec<usize>> = BTreeMap::new();
    for (k, &s) in strata_ids.iter().enumerate() {
        by_stratum.entry(s).or_default().push(k);
    }
    let informative: usize = by_stratum
        .values()
        .filter(|v| {
            v.iter()
                .map(|&k| lab[k])
                .collect::<std::collections::BTreeSet<_>>()
                .len()
                > 1
        })
        .map(|v| v.len())
        .sum();
    let informative_fraction = informative as f64 / m as f64;
    let strata_used = strata.is_some();
    if strata_used && informative == 0 {
        warnings.push(
            "the label never varies inside a stratum (it is confounded with the study): no p-value"
                .to_string(),
        );
    } else if strata_used && informative_fraction < 0.5 {
        warnings.push(format!(
            "the label varies inside a stratum for only {:.0}% of the samples; the p-value rests on those",
            100.0 * informative_fraction
        ));
    }

    let sizes: Vec<usize> = (0..g)
        .map(|j| lab.iter().filter(|&&l| l == j).count())
        .collect();
    let (smin, smax) = (*sizes.iter().min().unwrap(), *sizes.iter().max().unwrap());
    if smin < 5 {
        warnings.push(format!(
            "smallest group has {smin} samples; the test has little power for it"
        ));
    }

    // PERMANOVA
    let wsum: f64 = w.iter().sum();
    let total: f64 = (0..m)
        .map(|a| (a + 1..m).map(|b| w[a] * w[b] * d2[a * m + b]).sum::<f64>())
        .sum::<f64>()
        / wsum;
    let ss_w = |labs: &[usize]| -> f64 {
        let mut members: Vec<Vec<usize>> = vec![Vec::new(); g];
        for (k, &l) in labs.iter().enumerate() {
            members[l].push(k);
        }
        members
            .iter()
            .map(|mem| {
                let mut s = 0.0;
                for (x, &a) in mem.iter().enumerate() {
                    for &b in &mem[x + 1..] {
                        s += w[a] * w[b] * d2[a * m + b];
                    }
                }
                s / mem
                    .iter()
                    .map(|&a| w[a])
                    .sum::<f64>()
                    .max(f64::MIN_POSITIVE)
            })
            .sum()
    };
    let pseudo_f = |sw: f64| ((total - sw) / (g - 1) as f64) / (sw / (m - g) as f64);
    let sw_obs = ss_w(&lab);
    let f_obs = pseudo_f(sw_obs);
    // The total sum of squares does not change when labels are permuted, so F is a decreasing function of the within-group
    // sum of squares alone: "F at least as large as observed" is "SS_W at most as large as observed". Comparing the sums of
    // squares keeps the test exact when F is infinite (SS_W of 0: every group is made of exact copies) and removes the need
    // to compare NaN or inf. If the total is 0 (all distances 0) there is nothing to test.
    let degenerate_total = !(total.is_finite() && total > 0.0);
    if degenerate_total {
        warnings.push(
            "all distances are 0 or not finite: there is no structure to test, no p-value"
                .to_string(),
        );
    }
    let r2 = if total > 0.0 {
        (total - sw_obs) / total
    } else {
        0.0
    };

    // PERMDISP: distance of each sample to its group centroid, from the squared distances alone
    let centroid_dist = |labs: &[usize]| -> (Vec<f64>, usize) {
        let mut z = vec![0.0f64; m];
        let mut negative = 0usize;
        let mut members: Vec<Vec<usize>> = vec![Vec::new(); g];
        for (k, &l) in labs.iter().enumerate() {
            members[l].push(k);
        }
        for mem in &members {
            let ng: f64 = mem.iter().map(|&a| w[a]).sum();
            let mut within = 0.0;
            for &a in mem {
                for &b in mem {
                    within += w[a] * w[b] * d2[a * m + b];
                }
            }
            for &a in mem {
                let row: f64 = mem.iter().map(|&b| w[b] * d2[a * m + b]).sum::<f64>() / ng;
                let v = row - within / (2.0 * ng * ng);
                if v < 0.0 {
                    negative += 1;
                }
                z[a] = v.max(0.0).sqrt();
            }
        }
        (z, negative)
    };
    let (z, negative) = centroid_dist(&lab);
    if negative > 0 {
        warnings.push(format!(
            "{negative} sample(s) have a negative squared distance to the centroid (the distances are not Euclidean); set to 0"
        ));
    }
    // one-way ANOVA of the distances to the centroids: (between SS, within SS)
    let anova_ss = |zz: &[f64], labs: &[usize]| -> (f64, f64) {
        let mean: f64 = zz.iter().sum::<f64>() / m as f64;
        let mut gm = vec![0.0; g];
        let mut gn = vec![0usize; g];
        for (k, &l) in labs.iter().enumerate() {
            gm[l] += zz[k];
            gn[l] += 1;
        }
        for j in 0..g {
            gm[j] /= gn[j].max(1) as f64;
        }
        let ssb: f64 = (0..g).map(|j| gn[j] as f64 * (gm[j] - mean).powi(2)).sum();
        let ssw: f64 = zz
            .iter()
            .enumerate()
            .map(|(k, &x)| (x - gm[labs[k]]).powi(2))
            .sum();
        (ssb, ssw)
    };
    let anova_f = |ssb: f64, ssw: f64| (ssb / (g - 1) as f64) / (ssw / (m - g) as f64);
    let (ssb_obs, ssw_obs) = anova_ss(&z, &lab);
    let disp_obs = anova_f(ssb_obs, ssw_obs);
    // PERMDISP needs at least PERMDISP_MIN_GROUP samples in every group: in a group of 2 both samples are at the same distance
    // from the centroid, so the within-group SS is 0 and F infinite for any data, and the p-value is the minimum
    let sst_z = ssb_obs + ssw_obs;
    let disp_small = smin < PERMDISP_MIN_GROUP;
    let degenerate_disp = !(sst_z.is_finite() && sst_z > 0.0);
    if disp_small {
        warnings.push(format!(
            "a group has fewer than {PERMDISP_MIN_GROUP} samples: no PERMDISP p-value (the spread of a group of 2 is not defined)"
        ));
    } else if degenerate_disp {
        warnings.push(
            "every sample is at the same distance from its group centroid: no PERMDISP p-value"
                .to_string(),
        );
    }
    let degenerate_disp = degenerate_disp || disp_small;
    let groups: Vec<GroupStat> = (0..g)
        .map(|j| {
            let zs: Vec<f64> = (0..m).filter(|&k| lab[k] == j).map(|k| z[k]).collect();
            GroupStat {
                label: names[j].to_string(),
                n: zs.len(),
                mean_dispersion: zs.iter().sum::<f64>() / zs.len() as f64,
            }
        })
        .collect();

    // permutations of the label vector, within strata
    let (p, disp_p) = if (strata_used && informative == 0) || degenerate_total {
        (None, None)
    } else {
        let tol = 1e-12;
        let one = |t: usize| -> (bool, bool) {
            let mut rng = SplitMix(o.seed ^ (t as u64).wrapping_mul(0xD6E8_FEB8_6659_FD93));
            let mut labs = lab.clone();
            for v in by_stratum.values() {
                for i in (1..v.len()).rev() {
                    let j = (rng.next_u64() % (i as u64 + 1)) as usize;
                    labs.swap(v[i], v[j]);
                }
            }
            (
                ss_w(&labs) <= sw_obs + tol * total,
                !degenerate_disp && {
                    // the distances to the centroids are recomputed for the permuted groups: with the observed ones held
                    // fixed the null is too narrow for small or unequal groups (false-positive rate 0.10 to 0.35 in the
                    // designs measured), recomputing them keeps it at about 0.05 (see the calibration test)
                    let (zp, _) = centroid_dist(&labs);
                    let (b, wv) = anova_ss(&zp, &labs);
                    let fp = anova_f(b, wv);
                    if disp_obs.is_infinite() {
                        fp.is_infinite()
                    } else {
                        fp >= disp_obs * (1.0 - 1e-12)
                    }
                },
            )
        };
        let hits: Vec<(bool, bool)> = (0..o.n_perm).into_par_iter().map(one).collect();
        let a = hits.iter().filter(|h| h.0).count();
        let b = hits.iter().filter(|h| h.1).count();
        (
            Some((1 + a) as f64 / (1 + o.n_perm) as f64),
            (!degenerate_disp).then(|| (1 + b) as f64 / (1 + o.n_perm) as f64),
        )
    };

    if let (Some(pp), Some(pd)) = (p, disp_p) {
        if pp < 0.05 && pd < 0.05 {
            warnings.push(
                "PERMANOVA and PERMDISP are both significant: the groups differ in spread, so the PERMANOVA result is not evidence of different centres"
                    .to_string(),
            );
        } else if pp < 0.05 && smax as f64 >= 3.0 * smin as f64 && pd < 0.2 {
            warnings.push(
                "group sizes are unequal (largest at least 3 x smallest) and the spreads may differ; read the PERMANOVA with care"
                    .to_string(),
            );
        }
    }

    Some(TestResult {
        n: m,
        n_groups: g,
        f: f_obs,
        r2,
        p,
        disp_f: disp_obs,
        disp_p,
        groups,
        informative_fraction,
        strata_used,
        weighted: weights.is_some(),
        warnings,
    })
}

impl TestResult {
    pub fn to_json(&self) -> serde_json::Value {
        serde_json::json!({
            "n": self.n, "n_groups": self.n_groups, "pseudo_f": self.f, "r2": self.r2, "p": self.p,
            "permdisp_f": self.disp_f, "permdisp_p": self.disp_p,
            "strata_used": self.strata_used, "weighted": self.weighted, "informative_fraction": self.informative_fraction,
            "groups": self.groups.iter().map(|g| serde_json::json!({
                "label": g.label, "n": g.n, "mean_dispersion": g.mean_dispersion})).collect::<Vec<_>>(),
            "warnings": self.warnings,
        })
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn rng_f(r: &mut SplitMix) -> f64 {
        r.unit()
    }

    /// Euclidean distances between points in `dim` dimensions; group j has centre `centres[j]` and spread `spreads[j]`.
    fn planted(
        sizes: &[usize],
        centres: &[f64],
        spreads: &[f64],
        dim: usize,
        seed: u64,
    ) -> (Vec<f64>, usize, Vec<Option<String>>) {
        let mut r = SplitMix(seed);
        let mut pts: Vec<Vec<f64>> = Vec::new();
        let mut labs = Vec::new();
        for (j, &s) in sizes.iter().enumerate() {
            for _ in 0..s {
                let mut p = vec![0.0; dim];
                for x in p.iter_mut() {
                    *x = (rng_f(&mut r) - 0.5) * 2.0 * spreads[j];
                }
                p[0] += centres[j];
                pts.push(p);
                labs.push(Some(format!("g{j}")));
            }
        }
        let n = pts.len();
        let mut d = vec![0.0; n * n];
        for a in 0..n {
            for b in 0..n {
                d[a * n + b] = pts[a]
                    .iter()
                    .zip(&pts[b])
                    .map(|(x, y)| (x - y).powi(2))
                    .sum::<f64>()
                    .sqrt();
            }
        }
        (d, n, labs)
    }

    #[test]
    fn matches_a_hand_computed_pseudo_f() {
        // 4 points on a line at 0, 1, 10, 11; groups {0,1} and {10,11}
        let x: [f64; 4] = [0.0, 1.0, 10.0, 11.0];
        let n = 4;
        let mut d = vec![0.0; 16];
        for a in 0..n {
            for b in 0..n {
                d[a * n + b] = (x[a] - x[b]).abs();
            }
        }
        let labs: Vec<Option<String>> = ["a", "a", "b", "b"]
            .iter()
            .map(|s| Some(s.to_string()))
            .collect();
        let r = test(
            &d,
            n,
            &labs,
            None,
            None,
            &Options {
                n_perm: 99,
                seed: 3,
            },
        )
        .unwrap();
        // SS_T = (sum d^2)/4 = (1+100+121+81+100+1)/4 = 101; SS_W = 1/2 + 1/2 = 1; F = (100/1)/(1/2) = 200
        assert!((r.f - 200.0).abs() < 1e-9, "F = {}", r.f);
        assert!((r.r2 - 100.0 / 101.0).abs() < 1e-9);
        // of the 6 ways to give 2 of 4 points the label a, 2 (the observed one and its mirror) separate the pairs: the exact p
        // is 1/3, and 99 permutations estimate it as (1 + hits) / 100
        let p = r.p.unwrap();
        assert!((0.2..=0.5).contains(&p), "p = {p}");
    }

    #[test]
    fn detects_a_location_effect_and_not_a_null() {
        let (d, n, l) = planted(&[30, 30], &[0.0, 3.0], &[1.0, 1.0], 3, 11);
        let r = test(&d, n, &l, None, None, &Options::default()).unwrap();
        assert!(r.p.unwrap() <= 0.01, "p = {:?}", r.p);
        assert!(r.disp_p.unwrap() > 0.01, "equal spreads: {:?}", r.disp_p);
        let (d, n, l) = planted(&[30, 30], &[0.0, 0.0], &[1.0, 1.0], 3, 12);
        let r = test(&d, n, &l, None, None, &Options::default()).unwrap();
        assert!(r.p.unwrap() > 0.01, "null p = {:?}", r.p);
    }

    #[test]
    fn null_p_values_are_calibrated() {
        // 200 null data sets, alpha 0.1: the false-positive rate should be near 0.1 (binomial sd 0.021)
        let mut hits = 0;
        let o = Options {
            n_perm: 199,
            seed: 5,
        };
        for s in 0..200u64 {
            let (d, n, l) = planted(&[12, 12, 12], &[0.0; 3], &[1.0; 3], 4, 1000 + s);
            if test(&d, n, &l, None, None, &o).unwrap().p.unwrap() <= 0.1 {
                hits += 1;
            }
        }
        let rate = hits as f64 / 200.0;
        assert!((0.04..=0.17).contains(&rate), "false-positive rate {rate}");
    }

    #[test]
    fn dispersion_only_difference_is_caught_by_permdisp_and_warned() {
        // same centre, spread 0.3 vs 3.0, very unequal sizes
        let mut found = 0;
        for s in 0..20u64 {
            let (d, n, l) = planted(&[40, 6], &[0.0, 0.0], &[0.3, 3.0], 3, 50 + s);
            let r = test(
                &d,
                n,
                &l,
                None,
                None,
                &Options {
                    n_perm: 499,
                    seed: s,
                },
            )
            .unwrap();
            assert!(r.disp_p.unwrap() < 0.05, "PERMDISP missed: {:?}", r.disp_p);
            if r.p.unwrap() < 0.05 {
                found += 1;
                assert!(
                    r.warnings.iter().any(|w| w.contains("PERMDISP")),
                    "no warning: {:?}",
                    r.warnings
                );
            }
        }
        // the point of the warning: PERMANOVA alone often calls this a location effect
        assert!(
            found > 0,
            "expected PERMANOVA to be misled at least once in 20"
        );
    }

    #[test]
    fn permuting_within_strata_removes_a_study_effect() {
        // 4 studies, each with its own centre; the label tracks the study but has no effect inside a study
        let mut r = SplitMix(7);
        let n = 80;
        let mut pts = Vec::new();
        let mut labs = Vec::new();
        let mut st = Vec::new();
        for i in 0..n {
            let study = i / 20;
            pts.push(vec![study as f64 * 5.0 + rng_f(&mut r), rng_f(&mut r)]);
            // the label follows the study (a in studies 0-1, b in 2-3) except for every 7th sample, which is flipped
            let base = study < 2;
            labs.push(Some(
                if base != (i % 7 == 0) { "a" } else { "b" }.to_string(),
            ));
            st.push(Some(format!("s{study}")));
        }
        let mut d = vec![0.0; n * n];
        for a in 0..n {
            for b in 0..n {
                d[a * n + b] =
                    ((pts[a][0] - pts[b][0]).powi(2) + (pts[a][1] - pts[b][1]).powi(2)).sqrt();
            }
        }
        let o = Options {
            n_perm: 499,
            seed: 1,
        };
        let free = test(&d, n, &labs, None, None, &o).unwrap();
        let strat = test(&d, n, &labs, Some(&st), None, &o).unwrap();
        // the unrestricted test calls the label an effect because it tracks the study; staying inside studies does not
        assert!(free.p.unwrap() <= 0.01, "unrestricted p = {:?}", free.p);
        assert!(strat.p.unwrap() > 0.05, "within-study p = {:?}", strat.p);
        assert!(strat.informative_fraction > 0.99);
    }

    #[test]
    fn a_label_that_is_the_study_gets_no_p_value() {
        let (d, n, l) = planted(&[15, 15], &[0.0, 4.0], &[1.0, 1.0], 3, 21);
        let st: Vec<Option<String>> = l.clone(); // study == label
        let r = test(&d, n, &l, Some(&st), None, &Options::default()).unwrap();
        assert!(r.p.is_none() && r.disp_p.is_none());
        assert!(
            r.warnings.iter().any(|w| w.contains("confounded")),
            "{:?}",
            r.warnings
        );
    }

    #[test]
    fn missing_and_singleton_labels_are_left_out() {
        let (d, n, mut l) = planted(&[10, 10, 1], &[0.0, 3.0, 9.0], &[1.0; 3], 3, 31);
        l[3] = None;
        let r = test(&d, n, &l, None, None, &Options::default()).unwrap();
        assert_eq!(r.n_groups, 2);
        assert_eq!(r.n, 19);
        assert!(r.warnings.iter().any(|w| w.contains("one sample")));
        // fewer than two usable groups
        let only: Vec<Option<String>> = (0..n)
            .map(|i| if i < 2 { Some("x".into()) } else { None })
            .collect();
        assert!(test(&d, n, &only, None, None, &Options::default()).is_none());
    }

    #[test]
    fn same_seed_same_answer() {
        let (d, n, l) = planted(&[20, 20], &[0.0, 1.0], &[1.0, 1.0], 3, 41);
        let o = Options {
            n_perm: 199,
            seed: 9,
        };
        let a = test(&d, n, &l, None, None, &o).unwrap();
        let b = test(&d, n, &l, None, None, &o).unwrap();
        assert_eq!(a.p, b.p);
        assert_eq!(a.disp_p, b.disp_p);
    }

    /// 3 well separated groups; one sample of group 0 is duplicated `k` times (copies at distance 0 from it)
    fn dup_fixture(k: usize) -> (Vec<f64>, usize, Vec<Option<String>>, Vec<f64>) {
        let (d0, n0, l0) = planted(&[6, 6, 6], &[0.0, 1.0, 2.5], &[0.4, 0.4, 0.4], 3, 77);
        let n = n0 + k;
        let src = |i: usize| if i < n0 { i } else { 0 };
        let mut d = vec![0.0; n * n];
        for a in 0..n {
            for b in 0..n {
                let (x, y) = (src(a), src(b));
                // two copies of the same sample are at distance 0
                d[a * n + b] = if x == y { 0.0 } else { d0[x * n0 + y] };
            }
        }
        let mut l = l0.clone();
        l.extend(std::iter::repeat_n(l0[0].clone(), k));
        let mut w = vec![1.0; n0];
        w[0] = 1.0 + k as f64;
        (d, n, l, w)
    }

    #[test]
    fn weights_of_one_are_the_unweighted_test() {
        let (d, n, l) = planted(&[8, 8], &[0.0, 1.0], &[0.5, 0.5], 3, 5);
        let o = Options {
            n_perm: 99,
            seed: 2,
        };
        let a = test(&d, n, &l, None, None, &o).unwrap();
        let ones = vec![1.0; n];
        let b = test(&d, n, &l, None, Some(&ones), &o).unwrap();
        assert!((a.f - b.f).abs() < 1e-9 && (a.r2 - b.r2).abs() < 1e-12);
        assert_eq!(a.p, b.p);
        assert!(!a.weighted && b.weighted);
    }

    #[test]
    fn a_weight_of_k_plus_one_gives_the_same_sums_of_squares_as_k_exact_copies() {
        let (d, n, l, _) = dup_fixture(4);
        let (d0, n0, l0, w) = {
            let (d0, n0, l0) = planted(&[6, 6, 6], &[0.0, 1.0, 2.5], &[0.4, 0.4, 0.4], 3, 77);
            let mut w = vec![1.0; n0];
            w[0] = 5.0;
            (d0, n0, l0, w)
        };
        let o = Options { n_perm: 9, seed: 1 };
        let dup = test(&d, n, &l, None, None, &o).unwrap();
        let wtd = test(&d0, n0, &l0, None, Some(&w), &o).unwrap();
        assert!((dup.r2 - wtd.r2).abs() < 1e-9, "{} vs {}", dup.r2, wtd.r2);
    }

    #[test]
    fn collapsed_copies_do_not_inflate_the_degrees_of_freedom() {
        // the weighted test keeps n rows: its pseudo-F is the duplicated one scaled by the degrees of freedom only
        let (d, n, l, _) = dup_fixture(30);
        let (d0, n0, l0) = planted(&[6, 6, 6], &[0.0, 1.0, 2.5], &[0.4, 0.4, 0.4], 3, 77);
        let mut w = vec![1.0; n0];
        w[0] = 31.0;
        let o = Options { n_perm: 9, seed: 1 };
        let dup = test(&d, n, &l, None, None, &o).unwrap();
        let wtd = test(&d0, n0, &l0, None, Some(&w), &o).unwrap();
        assert_eq!(dup.n, n0 + 30);
        assert_eq!(wtd.n, n0);
        // same R2, but the weighted F is smaller: it does not count the copies as independent samples
        assert!((dup.r2 - wtd.r2).abs() < 1e-9);
        assert!(wtd.f < dup.f, "{} !< {}", wtd.f, dup.f);
    }

    fn labels(spec: &[&str]) -> Vec<Option<String>> {
        spec.iter().map(|s| Some(s.to_string())).collect()
    }

    #[test]
    fn all_zero_distances_give_no_p_value_not_the_smallest_one() {
        let n = 8;
        let d = vec![0.0; n * n];
        let l = labels(&["a", "a", "a", "a", "b", "b", "b", "b"]);
        let r = test(&d, n, &l, None, None, &Options::default()).unwrap();
        assert!(
            r.p.is_none() && r.disp_p.is_none(),
            "p = {:?}, {:?}",
            r.p,
            r.disp_p
        );
        assert!(
            r.warnings.iter().any(|w| w.contains("no p-value")),
            "{:?}",
            r.warnings
        );
    }

    #[test]
    fn groups_of_exact_copies_get_the_exact_p_not_the_smallest() {
        // two groups of 4 identical points, 0.01 apart: F is infinite. Of the C(8,4) = 70 ways to give 4 samples label a,
        // exactly 2 (the observed one and its mirror) separate the copies, so the exact p is 2/70 = 0.0286
        let n = 8;
        let mut d = vec![0.0; n * n];
        for a in 0..n {
            for b in 0..n {
                if (a < 4) != (b < 4) {
                    d[a * n + b] = 0.01;
                }
            }
        }
        let l = labels(&["a", "a", "a", "a", "b", "b", "b", "b"]);
        let r = test(
            &d,
            n,
            &l,
            None,
            None,
            &Options {
                n_perm: 4999,
                seed: 4,
            },
        )
        .unwrap();
        assert!(r.f.is_infinite());
        let p = r.p.unwrap();
        assert!((0.02..=0.04).contains(&p), "p = {p}");
        // every sample sits on its centroid: no spread to test
        assert!(r.disp_p.is_none(), "{:?}", r.disp_p);
    }

    #[test]
    fn a_real_effect_is_still_found_and_a_finite_f_still_matches_the_sum_of_squares_rule() {
        // the p-value is now decided on SS_W, which must give the same answer as F for finite F
        let (d, n, l) = planted(&[15, 15], &[0.0, 3.0], &[1.0, 1.0], 3, 99);
        let r = test(&d, n, &l, None, None, &Options::default()).unwrap();
        assert!(r.f.is_finite() && r.p.unwrap() <= 0.002, "{:?}", r.p);
    }

    /// share of `reps` null data sets (same centre and spread in every group) with PERMDISP p <= 0.05
    fn permdisp_rate(sizes: &[usize], centres: &[f64], spreads: &[f64], reps: u64) -> f64 {
        let o = |s| Options {
            n_perm: 199,
            seed: s,
        };
        let mut hits = 0;
        for s in 0..reps {
            let (d, n, l) = planted(sizes, centres, spreads, 5, 7000 + s);
            if test(&d, n, &l, None, None, &o(s)).unwrap().disp_p.unwrap() <= 0.05 {
                hits += 1;
            }
        }
        hits as f64 / reps as f64
    }

    #[test]
    fn permdisp_false_positive_rate_is_near_5_percent_for_small_and_unequal_groups() {
        // with the distances to the centroids held fixed the rate was 0.105 (5, 60), 0.14 (3, 40) and 0.35 (3, 3, 30, 30)
        // (code review 2026-11-12); binomial sd at 300 sets is 0.013, so 0.02 to 0.10 is about +-3.5 sd around 0.05
        for sizes in [
            vec![5, 60],
            vec![10, 60],
            vec![3, 40],
            vec![3, 3, 30, 30],
            vec![8, 8, 40, 40],
        ] {
            let c = vec![0.0; sizes.len()];
            let sp = vec![1.0; sizes.len()];
            let r = permdisp_rate(&sizes, &c, &sp, 300);
            assert!(
                (0.02..=0.10).contains(&r),
                "{sizes:?}: false-positive rate {r}"
            );
        }
    }

    #[test]
    fn permdisp_is_not_fooled_by_a_location_difference_and_still_sees_a_spread_difference() {
        // groups differ in centre only: PERMDISP must not reject much, even with unequal sizes
        let r = permdisp_rate(&[10, 60], &[0.0, 3.0], &[1.0, 1.0], 300);
        assert!(r <= 0.10, "location only: {r}");
        // same centre, spread 1 against 2: it must find it
        let r = permdisp_rate(&[10, 60], &[0.0, 0.0], &[1.0, 2.0], 100);
        assert!(r >= 0.95, "spread 1 vs 2: {r}");
    }

    #[test]
    fn groups_of_two_get_no_permdisp_p_value_but_permanova_still_runs() {
        // every group has 2 samples: both are the same distance from their centroid, so F was infinite and p minimal always
        let (d, n, l) = planted(&[2, 2, 2, 2, 2, 2, 2, 2], &[0.0; 8], &[1.0; 8], 4, 5);
        let r = test(&d, n, &l, None, None, &Options::default()).unwrap();
        assert!(r.disp_p.is_none(), "{:?}", r.disp_p);
        assert!(r.p.is_some());
        assert!(
            r.warnings.iter().any(|w| w.contains("PERMDISP")),
            "{:?}",
            r.warnings
        );
        // and the two spread warnings do not fire on a missing p-value
        assert!(!r.warnings.iter().any(|w| w.contains("both significant")));
    }
}
