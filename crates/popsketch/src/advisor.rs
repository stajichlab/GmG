//! The resolution advisor: what the data can and cannot say, in plain statements, from the noise floor and the distances.

use crate::dist::DistMat;
use crate::floor::Floor;
use crate::stats::{median, quantile};

#[derive(Clone, Debug)]
pub struct Advice {
    /// fraction of all sample pairs closer than the noise floor
    pub frac_pairs_below_floor: f64,
    /// per cluster: (cluster, samples, median within-cluster distance, that in floors)
    pub clusters: Vec<(usize, usize, f64, f64)>,
    pub statements: Vec<String>,
}

/// `clusters`: label per sample of the clustering that is reported (None: no clear cut).
pub fn advise(d: &DistMat, clusters: Option<&[usize]>, floor: &Floor) -> Advice {
    let n = d.n;
    let mut all = Vec::new();
    for i in 0..n {
        for j in (i + 1)..n {
            all.push(1.0 - d.ani(i, j));
        }
    }
    let below = if all.is_empty() {
        0.0
    } else {
        all.iter().filter(|&&x| x < floor.value).count() as f64 / all.len() as f64
    };
    let mut cl = Vec::new();
    if let Some(c) = clusters {
        let k = c.iter().copied().max().map_or(0, |m| m + 1);
        for g in 0..k {
            let m: Vec<usize> = (0..n).filter(|&i| c[i] == g).collect();
            if m.len() < 2 {
                continue;
            }
            let mut v = Vec::new();
            for (x, &i) in m.iter().enumerate() {
                for &j in &m[x + 1..] {
                    v.push(1.0 - d.ani(i, j));
                }
            }
            let med = median(&v);
            cl.push((
                g,
                m.len(),
                med,
                if floor.value > 0.0 {
                    med / floor.value
                } else {
                    f64::INFINITY
                },
            ));
        }
    }
    let mut st = vec![format!(
        "Technical noise floor: distance (1 - ANI) {:.5} ({}{}). Differences smaller than this are not distinguishable from re-sampling the same reads.",
        floor.value,
        if floor.source == "split_half" { "90th percentile of split-half distances of ".to_string() + &floor.n.to_string() + " samples" } else { "fixed fallback: too few read sets to split".to_string() },
        if floor.source == "split_half" { format!(", median {:.5}", floor.median) } else { String::new() }
    )];
    st.push(format!(
        "{:.0}% of all sample pairs are closer than the noise floor.",
        below * 100.0
    ));
    for (g, nn, med, f) in &cl {
        let tone = if *f < 1.0 {
            "within the noise: no structure inside this group is supported"
        } else if *f < 3.0 {
            "close to the noise: fine structure inside this group is weakly supported"
        } else {
            "well above the noise"
        };
        st.push(format!("Cluster {g} ({nn} samples): median within-cluster distance {med:.5} = {f:.1} floors ({tone})."));
    }
    if below > 0.5 {
        st.push("More than half of all pairs are inside the noise: the order of branches within groups is not supported. Use a SNP-based method (mapping or split k-mer alignment) for strain-level relationships.".to_string());
    }
    let _ = quantile;
    Advice {
        frac_pairs_below_floor: below,
        clusters: cl,
        statements: st,
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::dist::compute;

    fn floor(v: f64) -> Floor {
        Floor {
            value: v,
            n: 20,
            median: v / 5.0,
            p90: v,
            source: "split_half",
        }
    }

    #[test]
    fn tight_clones_are_below_the_floor_and_the_advisor_says_so() {
        // 12 samples that share 99.99% of 5000 hashes: distance about 5e-6, far below a floor of 0.001
        let sets: Vec<Vec<u64>> = (0..12u64)
            .map(|i| (0..5000u64).chain(10_000 + i..10_001 + i).collect())
            .collect();
        let d = compute(&sets, 21);
        let a = advise(&d, Some(&[0; 12]), &floor(0.001));
        assert!(a.frac_pairs_below_floor > 0.99);
        assert!(a.statements.iter().any(|s| s.contains("not supported")));
        assert!(a.clusters[0].3 < 1.0);
    }

    #[test]
    fn separated_groups_are_each_within_the_noise_and_far_apart() {
        let sets: Vec<Vec<u64>> = (0..12u64)
            .map(|i| {
                let g = i / 6;
                (g * 100_000..g * 100_000 + 4000)
                    .chain(5_000_000 + i * 10..5_000_000 + i * 10 + 400)
                    .collect()
            })
            .collect();
        let d = compute(&sets, 21);
        let c: Vec<usize> = (0..12).map(|i| i / 6).collect();
        let a = advise(&d, Some(&c), &floor(0.001));
        // 30 of the 66 pairs are within a group (below the floor), 36 are between groups (far above it)
        assert!(
            (a.frac_pairs_below_floor - 30.0 / 66.0).abs() < 1e-9,
            "{}",
            a.frac_pairs_below_floor
        );
        assert!(!a.statements.iter().any(|s| s.contains("More than half")));
        assert_eq!(a.clusters.len(), 2);
        assert!(a.clusters.iter().all(|x| x.3 < 1.0));
        assert!(a
            .statements
            .iter()
            .any(|s| s.contains("no structure inside this group")));
        // between-group distance is far above the floor
        assert!(1.0 - d.ani(0, 11) > 10.0 * 0.001);
    }
}
