//! DAPC: discriminant analysis of principal components. Linear discriminant analysis on the leading PCs of the hash-presence PCA,
//! with the groups of the primary label column. The discriminant axes are the directions that best separate the labelled groups
//! relative to the spread inside them; every sample (labelled or not) gets a score on each axis, and each labelled sample is
//! assigned to the group whose centroid is nearest in discriminant space. Samples assigned to a group other than their own label
//! are a second, independent hint (next to the ANI-based conflict check) that a label may be wrong.

use crate::ordination::jacobi_eigen;

#[derive(Clone, Debug)]
pub struct Dapc {
    /// number of PCs the analysis used
    pub n_pcs: usize,
    /// groups with at least 2 samples, in order; indices below refer to this list
    pub groups: Vec<String>,
    /// discriminant scores, one row per sample, one column per axis (groups - 1 axes at most)
    pub scores: Vec<Vec<f64>>,
    /// eigenvalue (between / within variance ratio) of each axis
    pub eigenvalues: Vec<f64>,
    /// for each labelled sample, the group whose centroid is nearest in discriminant space
    pub assigned: Vec<Option<usize>>,
    /// fraction of labelled samples assigned to their own group
    pub accuracy: f64,
}

/// `pcs`: one row per sample. `group`: index into `names` for each sample (None: unlabelled). Groups with fewer than 2 samples
/// are ignored; returns None with fewer than 2 usable groups or too few samples for the number of PCs.
pub fn dapc(pcs: &[Vec<f64>], group: &[Option<usize>], names: &[String]) -> Option<Dapc> {
    let n = pcs.len();
    let p_all = pcs.first()?.len();
    let mut count = vec![0usize; names.len()];
    for g in group.iter().flatten() {
        count[*g] += 1;
    }
    let usable: Vec<usize> = (0..names.len()).filter(|&g| count[g] >= 2).collect();
    if usable.len() < 2 {
        return None;
    }
    let remap: Vec<Option<usize>> = (0..names.len())
        .map(|g| usable.iter().position(|&u| u == g))
        .collect();
    let lab: Vec<Option<usize>> = group.iter().map(|g| g.and_then(|g| remap[g])).collect();
    let m = lab.iter().flatten().count();
    let ng = usable.len();
    // the within-group scatter needs more samples than PCs
    let p = p_all.min(m.saturating_sub(ng + 1));
    if p < 2 {
        return None;
    }
    let mu_all: Vec<f64> = (0..p)
        .map(|c| {
            lab.iter()
                .enumerate()
                .filter(|(_, l)| l.is_some())
                .map(|(i, _)| pcs[i][c])
                .sum::<f64>()
                / m as f64
        })
        .collect();
    let mut mu_g = vec![vec![0.0; p]; ng];
    let mut n_g = vec![0.0; ng];
    for (i, l) in lab.iter().enumerate() {
        if let Some(g) = l {
            n_g[*g] += 1.0;
            for c in 0..p {
                mu_g[*g][c] += pcs[i][c];
            }
        }
    }
    for g in 0..ng {
        mu_g[g].iter_mut().for_each(|v| *v /= n_g[g]);
    }
    let (mut sw, mut sb) = (vec![0.0; p * p], vec![0.0; p * p]);
    for (i, l) in lab.iter().enumerate() {
        if let Some(g) = l {
            for a in 0..p {
                for b in 0..p {
                    sw[a * p + b] += (pcs[i][a] - mu_g[*g][a]) * (pcs[i][b] - mu_g[*g][b]);
                }
            }
        }
    }
    for g in 0..ng {
        for a in 0..p {
            for b in 0..p {
                sb[a * p + b] += n_g[g] * (mu_g[g][a] - mu_all[a]) * (mu_g[g][b] - mu_all[b]);
            }
        }
    }
    // a tiny ridge keeps the within scatter invertible
    let ridge = 1e-8 * (0..p).map(|a| sw[a * p + a]).sum::<f64>().max(1e-12) / p as f64;
    for a in 0..p {
        sw[a * p + a] += ridge;
    }
    // Sw = L L^T (Cholesky), M = L^-1 Sb L^-T is symmetric; its eigenvectors u give v = L^-T u
    let mut l = vec![0.0; p * p];
    for i in 0..p {
        for j in 0..=i {
            let mut s = sw[i * p + j];
            for k in 0..j {
                s -= l[i * p + k] * l[j * p + k];
            }
            if i == j {
                if s <= 0.0 {
                    return None;
                }
                l[i * p + j] = s.sqrt();
            } else {
                l[i * p + j] = s / l[j * p + j];
            }
        }
    }
    let solve_l = |b: &[f64]| -> Vec<f64> {
        let mut x = vec![0.0; p];
        for i in 0..p {
            let mut s = b[i];
            for k in 0..i {
                s -= l[i * p + k] * x[k];
            }
            x[i] = s / l[i * p + i];
        }
        x
    };
    let solve_lt = |b: &[f64]| -> Vec<f64> {
        let mut x = vec![0.0; p];
        for i in (0..p).rev() {
            let mut s = b[i];
            for k in (i + 1)..p {
                s -= l[k * p + i] * x[k];
            }
            x[i] = s / l[i * p + i];
        }
        x
    };
    // columns of Sb, then M = L^-1 (L^-1 Sb)^T
    let cols: Vec<Vec<f64>> = (0..p)
        .map(|j| solve_l(&(0..p).map(|i| sb[i * p + j]).collect::<Vec<_>>()))
        .collect(); // cols[j] = L^-1 Sb[:, j]
    let mut mm = vec![0.0; p * p];
    for i in 0..p {
        let row_i = solve_l(&(0..p).map(|j| cols[j][i]).collect::<Vec<_>>()); // row i of L^-1 Sb, then L^-1 applied
        for j in 0..p {
            mm[i * p + j] = row_i[j];
        }
    }
    for i in 0..p {
        for j in (i + 1)..p {
            let v = (mm[i * p + j] + mm[j * p + i]) / 2.0;
            mm[i * p + j] = v;
            mm[j * p + i] = v;
        }
    }
    let (evals, evecs) = jacobi_eigen(&mm, p);
    let axes = (ng - 1).min(p);
    let vs: Vec<Vec<f64>> = (0..axes).map(|a| solve_lt(&evecs[a])).collect();
    let scores: Vec<Vec<f64>> = (0..n)
        .map(|i| {
            vs.iter()
                .map(|v| (0..p).map(|c| (pcs[i][c] - mu_all[c]) * v[c]).sum())
                .collect()
        })
        .collect();
    // centroids in discriminant space, then nearest-centroid assignment
    let mut cen = vec![vec![0.0; axes]; ng];
    for (i, lg) in lab.iter().enumerate() {
        if let Some(g) = lg {
            for a in 0..axes {
                cen[*g][a] += scores[i][a] / n_g[*g];
            }
        }
    }
    let assigned: Vec<Option<usize>> = (0..n)
        .map(|i| {
            lab[i].map(|_| {
                (0..ng)
                    .min_by(|&x, &y| {
                        let dx: f64 = (0..axes).map(|a| (scores[i][a] - cen[x][a]).powi(2)).sum();
                        let dy: f64 = (0..axes).map(|a| (scores[i][a] - cen[y][a]).powi(2)).sum();
                        dx.partial_cmp(&dy).unwrap()
                    })
                    .unwrap()
            })
        })
        .collect();
    let right = (0..n)
        .filter(|&i| lab[i].is_some() && assigned[i] == lab[i])
        .count();
    Some(Dapc {
        n_pcs: p,
        groups: usable.iter().map(|&g| names[g].clone()).collect(),
        scores,
        eigenvalues: evals[..axes].to_vec(),
        assigned,
        accuracy: right as f64 / m as f64,
    })
}

#[cfg(test)]
mod tests {
    use super::*;

    struct Rng(u64);
    impl Rng {
        fn gauss(&mut self) -> f64 {
            let mut u = || {
                self.0 ^= self.0 << 13;
                self.0 ^= self.0 >> 7;
                self.0 ^= self.0 << 17;
                ((self.0 >> 11) as f64 / (1u64 << 53) as f64).max(1e-12)
            };
            (-2.0 * u().ln()).sqrt() * (2.0 * std::f64::consts::PI * u()).cos()
        }
    }

    /// three groups of `per` samples in 8 dimensions: means differ along axes 0 and 1 only, unit noise everywhere
    fn data(per: usize) -> (Vec<Vec<f64>>, Vec<Option<usize>>, Vec<String>) {
        let mut r = Rng(99);
        let mean = [[0.0, 0.0], [6.0, 0.0], [3.0, 5.0]];
        let mut x = Vec::new();
        let mut g = Vec::new();
        for (gi, m) in mean.iter().enumerate() {
            for _ in 0..per {
                x.push(
                    (0..8)
                        .map(|c| if c < 2 { m[c] } else { 0.0 } + r.gauss())
                        .collect(),
                );
                g.push(Some(gi));
            }
        }
        (x, g, vec!["A".into(), "B".into(), "C".into()])
    }

    #[test]
    fn separates_groups_and_assigns_every_sample_to_its_own() {
        let (x, g, names) = data(15);
        let d = dapc(&x, &g, &names).unwrap();
        assert_eq!((d.groups.len(), d.scores[0].len(), d.n_pcs), (3, 2, 8));
        assert!(
            d.eigenvalues[0] >= d.eigenvalues[1] && d.eigenvalues[1] > 1.0,
            "{:?}",
            d.eigenvalues
        );
        assert!((d.accuracy - 1.0).abs() < 1e-12, "accuracy {}", d.accuracy);
        // LDA is optimal: the Fisher ratio of the first axis beats that of any single original axis
        let fisher = |vals: Vec<f64>| -> f64 {
            let gm: Vec<f64> = (0..3)
                .map(|k| vals[k * 15..(k + 1) * 15].iter().sum::<f64>() / 15.0)
                .collect();
            let tot = vals.iter().sum::<f64>() / 45.0;
            let between: f64 = gm.iter().map(|m| 15.0 * (m - tot).powi(2)).sum();
            let within: f64 = (0..45).map(|i| (vals[i] - gm[i / 15]).powi(2)).sum();
            between / within
        };
        let ld1 = fisher(d.scores.iter().map(|s| s[0]).collect());
        let best_pc = (0..8)
            .map(|c| fisher(x.iter().map(|r| r[c]).collect()))
            .fold(0.0, f64::max);
        assert!(ld1 >= best_pc - 1e-9, "{ld1} vs {best_pc}");
    }

    #[test]
    fn a_mislabelled_sample_is_assigned_to_the_group_it_resembles() {
        let (x, mut g, names) = data(15);
        g[3] = Some(1); // a sample drawn from group A, labelled B
        let d = dapc(&x, &g, &names).unwrap();
        assert_eq!(d.assigned[3], Some(0));
        let wrong: Vec<usize> = (0..x.len()).filter(|&i| d.assigned[i] != g[i]).collect();
        assert_eq!(wrong, vec![3], "{wrong:?}");
    }

    #[test]
    fn unlabelled_samples_get_scores_and_tiny_groups_are_ignored() {
        let (x, mut g, mut names) = data(15);
        g[0] = None; // unlabelled
        g[1] = Some(3); // a group of one
        names.push("lonely".into());
        let d = dapc(&x, &g, &names).unwrap();
        assert_eq!(
            d.groups,
            vec!["A".to_string(), "B".to_string(), "C".to_string()]
        );
        assert_eq!(d.scores.len(), x.len());
        assert!(d.assigned[0].is_none() && d.assigned[1].is_none());
        assert!(d.scores[0].iter().all(|v| v.is_finite()));
    }

    #[test]
    fn needs_two_groups_and_enough_samples() {
        let (x, g, names) = data(15);
        let one: Vec<Option<usize>> = g.iter().map(|_| Some(0)).collect();
        assert!(dapc(&x, &one, &names).is_none());
        let (x2, g2, n2) = data(2); // 6 samples: no room for 8 PCs (it uses fewer: 6 - 3 - 1 = 2)
        let d = dapc(&x2, &g2, &n2).unwrap();
        assert_eq!(d.n_pcs, 2);
    }
}
