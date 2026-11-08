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
use std::collections::BTreeMap;

pub struct Options {
    pub n_perm: usize,
    pub seed: u64,
}

impl Default for Options {
    fn default() -> Self {
        Options { n_perm: 999, seed: 1 }
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
    pub warnings: Vec<String>,
}

struct SplitMix(u64);
impl SplitMix {
    fn next(&mut self) -> u64 {
        self.0 = self.0.wrapping_add(0x9E37_79B9_7F4A_7C15);
        let mut z = self.0;
        z = (z ^ (z >> 30)).wrapping_mul(0xBF58_476D_1CE4_E5B9);
        z = (z ^ (z >> 27)).wrapping_mul(0x94D0_49BB_1331_11EB);
        z ^ (z >> 31)
    }
}

/// Run both tests. `d` is n x n row-major; `labels` and `strata` have length n. Returns `None` when fewer than two
/// groups with at least two samples remain.
pub fn test(
    d: &[f64],
    n: usize,
    labels: &[Option<String>],
    strata: Option<&[Option<String>]>,
    o: &Options,
) -> Option<TestResult> {
    assert_eq!(d.len(), n * n);
    assert_eq!(labels.len(), n);
    let mut warnings = Vec::new();

    // keep labelled samples whose group has at least 2 members
    let mut count: BTreeMap<&str, usize> = BTreeMap::new();
    for l in labels.iter().flatten() {
        *count.entry(l.as_str()).or_default() += 1;
    }
    let singles: Vec<&str> = count.iter().filter(|(_, &c)| c < 2).map(|(&l, _)| l).collect();
    if !singles.is_empty() {
        warnings.push(format!(
            "{} group(s) with one sample were left out: {}",
            singles.len(),
            singles.join(", ")
        ));
    }
    let names: Vec<&str> = count.iter().filter(|(_, &c)| c >= 2).map(|(&l, _)| l).collect();
    if names.len() < 2 {
        return None;
    }
    let idx: Vec<usize> = (0..n)
        .filter(|&i| labels[i].as_deref().map_or(false, |l| names.contains(&l)))
        .collect();
    let m = idx.len();
    let g = names.len();
    let lab: Vec<usize> = idx
        .iter()
        .map(|&i| names.iter().position(|&x| Some(x) == labels[i].as_deref()).unwrap())
        .collect();
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
        .filter(|v| v.iter().map(|&k| lab[k]).collect::<std::collections::BTreeSet<_>>().len() > 1)
        .map(|v| v.len())
        .sum();
    let informative_fraction = informative as f64 / m as f64;
    let strata_used = strata.is_some();
    if strata_used && informative == 0 {
        warnings.push(
            "the label never varies inside a stratum (it is confounded with the study): no p-value".to_string(),
        );
    } else if strata_used && informative_fraction < 0.5 {
        warnings.push(format!(
            "the label varies inside a stratum for only {:.0}% of the samples; the p-value rests on those",
            100.0 * informative_fraction
        ));
    }

    let sizes: Vec<usize> = (0..g).map(|j| lab.iter().filter(|&&l| l == j).count()).collect();
    let (smin, smax) = (*sizes.iter().min().unwrap(), *sizes.iter().max().unwrap());
    if smin < 5 {
        warnings.push(format!("smallest group has {smin} samples; the test has little power for it"));
    }

    // PERMANOVA
    let total: f64 = (0..m).map(|a| (a + 1..m).map(|b| d2[a * m + b]).sum::<f64>()).sum::<f64>() / m as f64;
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
                        s += d2[a * m + b];
                    }
                }
                s / mem.len().max(1) as f64
            })
            .sum()
    };
    let pseudo_f = |sw: f64| ((total - sw) / (g - 1) as f64) / (sw / (m - g) as f64);
    let sw_obs = ss_w(&lab);
    let f_obs = pseudo_f(sw_obs);
    let r2 = if total > 0.0 { (total - sw_obs) / total } else { 0.0 };

    // PERMDISP: distance of each sample to its group centroid, from the squared distances alone
    let mut z = vec![0.0f64; m];
    let mut negative = 0usize;
    {
        let mut members: Vec<Vec<usize>> = vec![Vec::new(); g];
        for (k, &l) in lab.iter().enumerate() {
            members[l].push(k);
        }
        for mem in &members {
            let ng = mem.len() as f64;
            let mut within = 0.0;
            for &a in mem {
                for &b in mem {
                    within += d2[a * m + b];
                }
            }
            for &a in mem {
                let row: f64 = mem.iter().map(|&b| d2[a * m + b]).sum::<f64>() / ng;
                let v = row - within / (2.0 * ng * ng);
                if v < 0.0 {
                    negative += 1;
                }
                z[a] = v.max(0.0).sqrt();
            }
        }
    }
    if negative > 0 {
        warnings.push(format!(
            "{negative} sample(s) have a negative squared distance to the centroid (the distances are not Euclidean); set to 0"
        ));
    }
    let anova_f = |zz: &[f64], labs: &[usize]| -> f64 {
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
        let ssw: f64 = zz.iter().enumerate().map(|(k, &x)| (x - gm[labs[k]]).powi(2)).sum();
        (ssb / (g - 1) as f64) / (ssw / (m - g) as f64)
    };
    let disp_obs = anova_f(&z, &lab);
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
    let (p, disp_p) = if strata_used && informative == 0 {
        (None, None)
    } else {
        let one = |t: usize| -> (bool, bool) {
            let mut rng = SplitMix(o.seed ^ (t as u64).wrapping_mul(0xD6E8_FEB8_6659_FD93));
            let mut labs = lab.clone();
            for v in by_stratum.values() {
                for i in (1..v.len()).rev() {
                    let j = (rng.next() % (i as u64 + 1)) as usize;
                    labs.swap(v[i], v[j]);
                }
            }
            (
                pseudo_f(ss_w(&labs)) >= f_obs - 1e-12 * f_obs.abs(),
                anova_f(&z, &labs) >= disp_obs - 1e-12 * disp_obs.abs(),
            )
        };
        let hits: Vec<(bool, bool)> = (0..o.n_perm).into_par_iter().map(one).collect();
        let a = hits.iter().filter(|h| h.0).count();
        let b = hits.iter().filter(|h| h.1).count();
        (
            Some((1 + a) as f64 / (1 + o.n_perm) as f64),
            Some((1 + b) as f64 / (1 + o.n_perm) as f64),
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
        warnings,
    })
}

impl TestResult {
    pub fn to_json(&self) -> serde_json::Value {
        serde_json::json!({
            "n": self.n, "n_groups": self.n_groups, "pseudo_f": self.f, "r2": self.r2, "p": self.p,
            "permdisp_f": self.disp_f, "permdisp_p": self.disp_p,
            "strata_used": self.strata_used, "informative_fraction": self.informative_fraction,
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
        (r.next() >> 11) as f64 / (1u64 << 53) as f64
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
                d[a * n + b] = pts[a].iter().zip(&pts[b]).map(|(x, y)| (x - y).powi(2)).sum::<f64>().sqrt();
            }
        }
        (d, n, labs)
    }

    #[test]
    fn matches_a_hand_computed_pseudo_f() {
        // 4 points on a line at 0, 1, 10, 11; groups {0,1} and {10,11}
        let x = [0.0, 1.0, 10.0, 11.0];
        let n = 4;
        let mut d = vec![0.0; 16];
        for a in 0..n {
            for b in 0..n {
                d[a * n + b] = (x[a] - x[b] as f64).abs();
            }
        }
        let labs: Vec<Option<String>> = ["a", "a", "b", "b"].iter().map(|s| Some(s.to_string())).collect();
        let r = test(&d, n, &labs, None, &Options { n_perm: 99, seed: 3 }).unwrap();
        // SS_T = (sum d^2)/4 = (1+100+121+81+100+1)/4 = 101; SS_W = 1/2 + 1/2 = 1; F = (100/1)/(1/2) = 200
        assert!((r.f - 200.0).abs() < 1e-9, "F = {}", r.f);
        assert!((r.r2 - 100.0 / 101.0).abs() < 1e-9);
        // only 3 distinct splits of 4 points into 2 pairs exist, so the smallest possible p is not below 1/3
        assert!(r.p.unwrap() >= 0.2, "p = {:?}", r.p);
    }

    #[test]
    fn detects_a_location_effect_and_not_a_null() {
        let (d, n, l) = planted(&[30, 30], &[0.0, 3.0], &[1.0, 1.0], 3, 11);
        let r = test(&d, n, &l, None, &Options::default()).unwrap();
        assert!(r.p.unwrap() <= 0.01, "p = {:?}", r.p);
        assert!(r.disp_p.unwrap() > 0.01, "equal spreads: {:?}", r.disp_p);
        let (d, n, l) = planted(&[30, 30], &[0.0, 0.0], &[1.0, 1.0], 3, 12);
        let r = test(&d, n, &l, None, &Options::default()).unwrap();
        assert!(r.p.unwrap() > 0.01, "null p = {:?}", r.p);
    }

    #[test]
    fn null_p_values_are_calibrated() {
        // 200 null data sets, alpha 0.1: the false-positive rate should be near 0.1 (binomial sd 0.021)
        let mut hits = 0;
        let o = Options { n_perm: 199, seed: 5 };
        for s in 0..200u64 {
            let (d, n, l) = planted(&[12, 12, 12], &[0.0; 3], &[1.0; 3], 4, 1000 + s);
            if test(&d, n, &l, None, &o).unwrap().p.unwrap() <= 0.1 {
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
            let r = test(&d, n, &l, None, &Options { n_perm: 499, seed: s }).unwrap();
            assert!(r.disp_p.unwrap() < 0.05, "PERMDISP missed: {:?}", r.disp_p);
            if r.p.unwrap() < 0.05 {
                found += 1;
                assert!(r.warnings.iter().any(|w| w.contains("PERMDISP")), "no warning: {:?}", r.warnings);
            }
        }
        // the point of the warning: PERMANOVA alone often calls this a location effect
        assert!(found > 0, "expected PERMANOVA to be misled at least once in 20");
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
            labs.push(Some(if base != (i % 7 == 0) { "a" } else { "b" }.to_string()));
            st.push(Some(format!("s{study}")));
        }
        let mut d = vec![0.0; n * n];
        for a in 0..n {
            for b in 0..n {
                d[a * n + b] = ((pts[a][0] - pts[b][0]).powi(2) + (pts[a][1] - pts[b][1]).powi(2)).sqrt();
            }
        }
        let o = Options { n_perm: 499, seed: 1 };
        let free = test(&d, n, &labs, None, &o).unwrap();
        let strat = test(&d, n, &labs, Some(&st), &o).unwrap();
        // the unrestricted test calls the label an effect because it tracks the study; staying inside studies does not
        assert!(free.p.unwrap() <= 0.01, "unrestricted p = {:?}", free.p);
        assert!(strat.p.unwrap() > 0.05, "within-study p = {:?}", strat.p);
        assert!(strat.informative_fraction > 0.99);
    }

    #[test]
    fn a_label_that_is_the_study_gets_no_p_value() {
        let (d, n, l) = planted(&[15, 15], &[0.0, 4.0], &[1.0, 1.0], 3, 21);
        let st: Vec<Option<String>> = l.clone(); // study == label
        let r = test(&d, n, &l, Some(&st), &Options::default()).unwrap();
        assert!(r.p.is_none() && r.disp_p.is_none());
        assert!(r.warnings.iter().any(|w| w.contains("confounded")), "{:?}", r.warnings);
    }

    #[test]
    fn missing_and_singleton_labels_are_left_out() {
        let (d, n, mut l) = planted(&[10, 10, 1], &[0.0, 3.0, 9.0], &[1.0; 3], 3, 31);
        l[3] = None;
        let r = test(&d, n, &l, None, &Options::default()).unwrap();
        assert_eq!(r.n_groups, 2);
        assert_eq!(r.n, 19);
        assert!(r.warnings.iter().any(|w| w.contains("one sample")));
        // fewer than two usable groups
        let only: Vec<Option<String>> = (0..n).map(|i| if i < 2 { Some("x".into()) } else { None }).collect();
        assert!(test(&d, n, &only, None, &Options::default()).is_none());
    }

    #[test]
    fn same_seed_same_answer() {
        let (d, n, l) = planted(&[20, 20], &[0.0, 1.0], &[1.0, 1.0], 3, 41);
        let o = Options { n_perm: 199, seed: 9 };
        let a = test(&d, n, &l, None, &o).unwrap();
        let b = test(&d, n, &l, None, &o).unwrap();
        assert_eq!(a.p, b.p);
        assert_eq!(a.disp_p, b.disp_p);
    }
}
