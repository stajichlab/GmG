//! Ordination of the samples: hash-presence PCA (the default) and PCoA of 1 - ANI (an option).
//!
//! Hash PCA treats every hash seen in at least 2 and not in all samples as a binary feature. Because a private or error
//! k-mer adds variance only to its own sample, the feature set is frequency-filtered first. The PCA is done on the centred
//! n x n Gram matrix (n = samples, far fewer than features): with X the 0/1 sample-by-hash matrix and mu the column means,
//! (X - 1 mu^T)(X - 1 mu^T)^T = X X^T - s 1^T - 1 s^T + t 11^T, where s_i = sum of mu over the hashes of sample i and
//! t = sum of mu^2. The eigenvectors of this matrix, scaled by the square root of their eigenvalue, are the PC scores.

use rayon::prelude::*;
use std::collections::HashMap;

#[derive(Clone, Debug)]
pub struct Ordination {
    /// scores, one row per sample, one column per component
    pub coords: Vec<Vec<f64>>,
    pub eigenvalues: Vec<f64>,
    /// eigenvalue / total variance
    pub var_frac: Vec<f64>,
    /// hashes used (hash PCA) or samples (PCoA)
    pub n_features: usize,
    pub method: &'static str,
}

/// Full eigendecomposition of a small symmetric matrix by cyclic Jacobi rotations. Returns eigenvalues in descending order
/// and the matching eigenvectors as rows of the returned matrix.
pub fn jacobi_eigen(a: &[f64], n: usize) -> (Vec<f64>, Vec<Vec<f64>>) {
    let mut m = a.to_vec();
    let mut v = vec![0.0; n * n];
    for i in 0..n {
        v[i * n + i] = 1.0;
    }
    for _sweep in 0..100 {
        let off: f64 = (0..n)
            .flat_map(|i| (0..n).filter(move |&j| j != i).map(move |j| (i, j)))
            .map(|(i, j)| m[i * n + j] * m[i * n + j])
            .sum();
        if off < 1e-22 {
            break;
        }
        for p in 0..n {
            for q in (p + 1)..n {
                let apq = m[p * n + q];
                if apq.abs() < 1e-300 {
                    continue;
                }
                let theta = (m[q * n + q] - m[p * n + p]) / (2.0 * apq);
                let t = theta.signum() / (theta.abs() + (theta * theta + 1.0).sqrt());
                let t = if theta == 0.0 { 1.0 } else { t };
                let c = 1.0 / (t * t + 1.0).sqrt();
                let s = t * c;
                for k in 0..n {
                    let (akp, akq) = (m[k * n + p], m[k * n + q]);
                    m[k * n + p] = c * akp - s * akq;
                    m[k * n + q] = s * akp + c * akq;
                }
                for k in 0..n {
                    let (apk, aqk) = (m[p * n + k], m[q * n + k]);
                    m[p * n + k] = c * apk - s * aqk;
                    m[q * n + k] = s * apk + c * aqk;
                }
                for k in 0..n {
                    let (vkp, vkq) = (v[k * n + p], v[k * n + q]);
                    v[k * n + p] = c * vkp - s * vkq;
                    v[k * n + q] = s * vkp + c * vkq;
                }
            }
        }
    }
    let mut idx: Vec<usize> = (0..n).collect();
    idx.sort_by(|&x, &y| m[y * n + y].partial_cmp(&m[x * n + x]).unwrap());
    let evals = idx.iter().map(|&i| m[i * n + i]).collect();
    let evecs = idx
        .iter()
        .map(|&i| (0..n).map(|k| v[k * n + i]).collect())
        .collect();
    (evals, evecs)
}

struct SplitMix(u64);
impl SplitMix {
    fn next_f64(&mut self) -> f64 {
        self.0 = self.0.wrapping_add(0x9E37_79B9_7F4A_7C15);
        let mut z = self.0;
        z = (z ^ (z >> 30)).wrapping_mul(0xBF58_476D_1CE4_E5B9);
        z = (z ^ (z >> 27)).wrapping_mul(0x94D0_49BB_1331_11EB);
        (((z ^ (z >> 31)) >> 11) as f64 / (1u64 << 53) as f64) - 0.5
    }
}

fn matvec_block(a: &[f64], n: usize, q: &[Vec<f64>]) -> Vec<Vec<f64>> {
    q.par_iter()
        .map(|col| {
            (0..n)
                .map(|i| {
                    a[i * n..(i + 1) * n]
                        .iter()
                        .zip(col)
                        .map(|(x, y)| x * y)
                        .sum()
                })
                .collect()
        })
        .collect()
}

fn orthonormalize(cols: &mut [Vec<f64>]) {
    for j in 0..cols.len() {
        for i in 0..j {
            let dot: f64 = cols[i].iter().zip(&cols[j]).map(|(x, y)| x * y).sum();
            let ci = cols[i].clone();
            for (x, y) in cols[j].iter_mut().zip(&ci) {
                *x -= dot * y;
            }
        }
        let norm = cols[j]
            .iter()
            .map(|x| x * x)
            .sum::<f64>()
            .sqrt()
            .max(1e-300);
        for x in cols[j].iter_mut() {
            *x /= norm;
        }
    }
}

/// The `k` largest eigenvalues (descending) and eigenvectors of a symmetric n x n matrix: block power iteration on A + shift I
/// (the shift makes every eigenvalue non-negative, so the largest algebraic ones are found even when A has large negative
/// eigenvalues, as a PCoA matrix can), then a Rayleigh-Ritz step.
pub fn top_eigs(a: &[f64], n: usize, k: usize) -> (Vec<f64>, Vec<Vec<f64>>) {
    let k = k.min(n);
    if n <= 64 {
        let (e, v) = jacobi_eigen(a, n);
        return (e[..k].to_vec(), v[..k].to_vec());
    }
    let b = (k + 8).min(n);
    let mut rng = SplitMix(0xC0FFEE);
    let mut q: Vec<Vec<f64>> = (0..b)
        .map(|_| (0..n).map(|_| rng.next_f64()).collect())
        .collect();
    orthonormalize(&mut q);
    // spectral radius bound: a few power steps on the first column
    let mut x = q[0].clone();
    let mut rho = 0.0;
    for _ in 0..30 {
        let y = matvec_block(a, n, std::slice::from_ref(&x)).remove(0);
        rho = y.iter().map(|v| v * v).sum::<f64>().sqrt();
        if rho == 0.0 {
            break;
        }
        x = y.iter().map(|v| v / rho).collect();
    }
    let shift = rho * 1.05;
    let mut prev = vec![0.0; k];
    for _ in 0..3000 {
        let mut z = matvec_block(a, n, &q);
        for (zc, qc) in z.iter_mut().zip(&q) {
            for (zv, qv) in zc.iter_mut().zip(qc) {
                *zv += shift * qv;
            }
        }
        orthonormalize(&mut z);
        q = z;
        let aq = matvec_block(a, n, &q);
        let ritz: Vec<f64> = (0..k)
            .map(|j| q[j].iter().zip(&aq[j]).map(|(x, y)| x * y).sum())
            .collect();
        let done = ritz
            .iter()
            .zip(&prev)
            .all(|(r, p)| (r - p).abs() <= 1e-10 * r.abs().max(1e-12));
        prev = ritz;
        if done {
            break;
        }
    }
    // Rayleigh-Ritz
    let aq = matvec_block(a, n, &q);
    let mut t = vec![0.0; b * b];
    for i in 0..b {
        for j in 0..b {
            t[i * b + j] = q[i].iter().zip(&aq[j]).map(|(x, y)| x * y).sum();
        }
    }
    for i in 0..b {
        for j in (i + 1)..b {
            let m = (t[i * b + j] + t[j * b + i]) / 2.0;
            t[i * b + j] = m;
            t[j * b + i] = m;
        }
    }
    let (evals, evecs) = jacobi_eigen(&t, b);
    let vecs: Vec<Vec<f64>> = (0..k)
        .map(|j| {
            (0..n)
                .map(|i| (0..b).map(|c| evecs[j][c] * q[c][i]).sum())
                .collect()
        })
        .collect();
    (evals[..k].to_vec(), vecs)
}

fn finish(
    evals: Vec<f64>,
    evecs: Vec<Vec<f64>>,
    total: f64,
    n_features: usize,
    method: &'static str,
) -> Ordination {
    let n = evecs.first().map_or(0, |v| v.len());
    let mut coords = vec![vec![0.0; evals.len()]; n];
    for (c, (e, v)) in evals.iter().zip(&evecs).enumerate() {
        // sign convention: the component with the largest magnitude is positive
        let big = v
            .iter()
            .cloned()
            .fold(0.0f64, |m, x| if x.abs() > m.abs() { x } else { m });
        let sgn = if big < 0.0 { -1.0 } else { 1.0 };
        let sd = e.max(0.0).sqrt();
        for (i, x) in v.iter().enumerate() {
            coords[i][c] = sgn * x * sd;
        }
    }
    let var_frac = evals
        .iter()
        .map(|e| if total > 0.0 { e.max(0.0) / total } else { 0.0 })
        .collect();
    Ordination {
        coords,
        eigenvalues: evals,
        var_frac,
        n_features,
        method,
    }
}

/// The centred Gram matrix of the frequency-filtered hash-presence matrix, and the number of features.
pub fn gram_matrix(sets: &[Vec<u64>]) -> Option<(Vec<f64>, usize)> {
    let n = sets.len();
    if n < 3 {
        return None;
    }
    let mut freq: HashMap<u64, u32> = HashMap::new();
    for s in sets {
        for &h in s {
            *freq.entry(h).or_insert(0) += 1;
        }
    }
    let mut cols: Vec<(u64, u32)> = freq
        .into_iter()
        .filter(|&(_, f)| f >= 2 && (f as usize) < n)
        .collect();
    cols.sort_unstable();
    let m = cols.len();
    if m == 0 {
        return None;
    }
    let keys: Vec<u64> = cols.iter().map(|c| c.0).collect();
    let mu: Vec<f64> = cols.iter().map(|c| c.1 as f64 / n as f64).collect();
    let lists: Vec<Vec<u32>> = sets
        .par_iter()
        .map(|s| {
            s.iter()
                .filter_map(|h| keys.binary_search(h).ok().map(|i| i as u32))
                .collect()
        })
        .collect();
    let t: f64 = mu.iter().map(|x| x * x).sum();
    let s_i: Vec<f64> = lists
        .iter()
        .map(|l| l.iter().map(|&c| mu[c as usize]).sum())
        .collect();
    // the number of selected hashes two samples share is the popcount of the AND of their bitsets over the selected columns
    let words = m.div_ceil(64);
    let bits: Vec<Vec<u64>> = lists
        .par_iter()
        .map(|l| {
            let mut b = vec![0u64; words];
            for &c in l {
                b[c as usize / 64] |= 1 << (c as usize % 64);
            }
            b
        })
        .collect();
    let rows: Vec<Vec<f64>> = (0..n)
        .into_par_iter()
        .map(|i| {
            (i..n)
                .map(|j| {
                    let c: u32 = bits[i]
                        .iter()
                        .zip(&bits[j])
                        .map(|(x, y)| (x & y).count_ones())
                        .sum();
                    c as f64 - s_i[i] - s_i[j] + t
                })
                .collect()
        })
        .collect();
    let mut g = vec![0.0; n * n];
    for i in 0..n {
        for (off, &v) in rows[i].iter().enumerate() {
            g[i * n + i + off] = v;
            g[(i + off) * n + i] = v;
        }
    }
    Some((g, m))
}

/// Hash-presence PCA of the samples' (abundance-filtered) hash sets.
pub fn hash_pca(sets: &[Vec<u64>], components: usize) -> Option<Ordination> {
    let (g, m) = gram_matrix(sets)?;
    let n = sets.len();
    let total: f64 = (0..n).map(|i| g[i * n + i]).sum();
    let (e, v) = top_eigs(&g, n, components);
    Some(finish(e, v, total, m, "hash_pca"))
}

/// Classical MDS (PCoA) of a symmetric distance matrix.
pub fn pcoa(dist: &[f64], n: usize, components: usize) -> Option<Ordination> {
    if n < 3 {
        return None;
    }
    let d2: Vec<f64> = dist.iter().map(|d| d * d).collect();
    let row: Vec<f64> = (0..n)
        .map(|i| d2[i * n..(i + 1) * n].iter().sum::<f64>() / n as f64)
        .collect();
    let all: f64 = row.iter().sum::<f64>() / n as f64;
    let mut b = vec![0.0; n * n];
    for i in 0..n {
        for j in 0..n {
            b[i * n + j] = -0.5 * (d2[i * n + j] - row[i] - row[j] + all);
        }
    }
    let total: f64 = (0..n).map(|i| b[i * n + i]).sum();
    let (e, v) = top_eigs(&b, n, components);
    Some(finish(e, v, total, n, "pcoa"))
}

#[cfg(test)]
mod tests {
    use super::*;

    fn rand_psd(n: usize, rank: usize, seed: u64) -> Vec<f64> {
        let mut r = SplitMix(seed);
        let f: Vec<Vec<f64>> = (0..n)
            .map(|_| (0..rank).map(|_| r.next_f64() * 2.0).collect())
            .collect();
        let mut a = vec![0.0; n * n];
        for i in 0..n {
            for j in 0..n {
                a[i * n + j] = (0..rank)
                    .map(|c| f[i][c] * f[j][c] * (rank - c) as f64)
                    .sum();
            }
        }
        a
    }

    #[test]
    fn jacobi_known_matrix() {
        // [[2,1],[1,2]] has eigenvalues 3 and 1
        let (e, v) = jacobi_eigen(&[2.0, 1.0, 1.0, 2.0], 2);
        assert!((e[0] - 3.0).abs() < 1e-12 && (e[1] - 1.0).abs() < 1e-12);
        assert!((v[0][0].abs() - v[0][1].abs()).abs() < 1e-9);
    }

    #[test]
    fn power_iteration_matches_jacobi_on_a_large_matrix() {
        let n = 120;
        let a = rand_psd(n, 6, 3);
        let (full, _) = jacobi_eigen(&a, n);
        let (e, v) = top_eigs(&a, n, 4);
        for i in 0..4 {
            assert!(
                (e[i] - full[i]).abs() < 1e-6 * full[0],
                "eig {i}: {} vs {}",
                e[i],
                full[i]
            );
        }
        // eigenvector equation and orthonormality
        for (i, vec) in v.iter().enumerate() {
            let av: Vec<f64> = (0..n)
                .map(|r| (0..n).map(|c| a[r * n + c] * vec[c]).sum())
                .collect();
            let resid: f64 = av
                .iter()
                .zip(vec)
                .map(|(x, y)| (x - e[i] * y).powi(2))
                .sum::<f64>()
                .sqrt();
            assert!(resid < 1e-5 * full[0], "residual {resid}");
            assert!((vec.iter().map(|x| x * x).sum::<f64>() - 1.0).abs() < 1e-9);
        }
    }

    #[test]
    fn matrix_with_large_negative_eigenvalues_still_gives_the_largest_positive_ones() {
        let n = 100;
        let mut a = rand_psd(n, 3, 9);
        // subtract a large rank-one part so the most extreme eigenvalue is negative
        let mut r = SplitMix(77);
        let u: Vec<f64> = (0..n).map(|_| r.next_f64()).collect();
        let big = 40.0 * a.iter().cloned().fold(0.0, f64::max);
        for i in 0..n {
            for j in 0..n {
                a[i * n + j] -= big * u[i] * u[j] / u.iter().map(|x| x * x).sum::<f64>();
            }
        }
        let (full, _) = jacobi_eigen(&a, n);
        assert!(
            full[n - 1].abs() > full[0],
            "setup: the most extreme eigenvalue should be negative"
        );
        let (e, _) = top_eigs(&a, n, 3);
        for i in 0..3 {
            assert!(
                (e[i] - full[i]).abs() < 1e-5 * full[0].abs().max(1.0),
                "eig {i}: {} vs {}",
                e[i],
                full[i]
            );
        }
    }

    /// Three groups that differ in group-specific hashes, with shared hashes and a few private ones.
    fn three_groups() -> (Vec<Vec<u64>>, Vec<usize>) {
        let mut sets = Vec::new();
        let mut lab = Vec::new();
        let mut r = SplitMix(1);
        for g in 0..3u64 {
            for s in 0..8u64 {
                let mut h: Vec<u64> = (0..500).collect(); // shared by all
                h.extend((1000 * (g + 1))..(1000 * (g + 1) + 300)); // this group's own
                for _ in 0..40 {
                    // sample-specific variation inside the group, drawn from a common pool
                    h.push(10_000 + ((r.next_f64() + 0.5) * 200.0) as u64);
                }
                h.push(1_000_000 + g * 100 + s); // private (error-like), must not matter
                h.sort_unstable();
                h.dedup();
                sets.push(h);
                lab.push(g as usize);
            }
        }
        (sets, lab)
    }

    #[test]
    fn gram_matrix_equals_the_explicit_centred_matrix() {
        let (sets, _) = three_groups();
        let n = sets.len();
        let (g, m) = gram_matrix(&sets).unwrap();
        // explicit: frequency-filtered columns, dense X, centre, X X^T
        let mut freq: HashMap<u64, usize> = HashMap::new();
        for s in &sets {
            for &h in s {
                *freq.entry(h).or_default() += 1;
            }
        }
        let mut cols: Vec<u64> = freq
            .iter()
            .filter(|(_, &f)| f >= 2 && f < n)
            .map(|(&h, _)| h)
            .collect();
        cols.sort_unstable();
        assert_eq!(cols.len(), m);
        let x: Vec<Vec<f64>> = sets
            .iter()
            .map(|s| {
                cols.iter()
                    .map(|c| if s.binary_search(c).is_ok() { 1.0 } else { 0.0 })
                    .collect()
            })
            .collect();
        let mu: Vec<f64> = (0..m)
            .map(|c| x.iter().map(|r| r[c]).sum::<f64>() / n as f64)
            .collect();
        for i in 0..n {
            for j in 0..n {
                let d: f64 = (0..m).map(|c| (x[i][c] - mu[c]) * (x[j][c] - mu[c])).sum();
                assert!(
                    (d - g[i * n + j]).abs() < 1e-9,
                    "({i},{j}): {d} vs {}",
                    g[i * n + j]
                );
            }
        }
    }

    #[test]
    fn hash_pca_separates_the_groups_and_ignores_private_hashes() {
        let (sets, lab) = three_groups();
        let o = hash_pca(&sets, 3).unwrap();
        assert_eq!(o.method, "hash_pca");
        assert!(
            o.var_frac[0] > 0.3 && o.var_frac[0] + o.var_frac[1] > 0.85,
            "{:?}",
            o.var_frac
        );
        // within a group the first two PCs hardly vary compared with the gaps between groups
        let centre = |g: usize| -> (f64, f64) {
            let pts: Vec<_> = o
                .coords
                .iter()
                .zip(&lab)
                .filter(|(_, &l)| l == g)
                .map(|(c, _)| (c[0], c[1]))
                .collect();
            (
                pts.iter().map(|p| p.0).sum::<f64>() / pts.len() as f64,
                pts.iter().map(|p| p.1).sum::<f64>() / pts.len() as f64,
            )
        };
        let (c0, c1, c2) = (centre(0), centre(1), centre(2));
        let d = |a: (f64, f64), b: (f64, f64)| ((a.0 - b.0).powi(2) + (a.1 - b.1).powi(2)).sqrt();
        let spread = o
            .coords
            .iter()
            .zip(&lab)
            .map(|(c, &l)| d((c[0], c[1]), centre(l)))
            .fold(0.0, f64::max);
        assert!(
            d(c0, c1) > 3.0 * spread && d(c0, c2) > 3.0 * spread && d(c1, c2) > 3.0 * spread,
            "{spread} {:?}",
            (c0, c1, c2)
        );
    }

    #[test]
    fn pcoa_recovers_a_known_configuration() {
        // four points on a line at 0, 1, 3, 6: distances are exact, PC1 recovers the positions up to sign and shift
        let pos: [f64; 4] = [0.0, 1.0, 3.0, 6.0];
        let n = 4;
        let mut d = vec![0.0; 16];
        for i in 0..n {
            for j in 0..n {
                d[i * n + j] = (pos[i] - pos[j]).abs();
            }
        }
        let o = pcoa(&d, n, 2).unwrap();
        let c: Vec<f64> = o.coords.iter().map(|r| r[0]).collect();
        let mean = c.iter().sum::<f64>() / n as f64;
        let rec: Vec<f64> = c.iter().map(|x| x - mean).collect();
        let pm = pos.iter().sum::<f64>() / n as f64;
        let sign = if rec[3] > 0.0 { 1.0 } else { -1.0 };
        for i in 0..n {
            assert!((sign * rec[i] - (pos[i] - pm)).abs() < 1e-9, "{rec:?}");
        }
        assert!(o.var_frac[0] > 0.999);
    }

    #[test]
    fn too_few_samples_give_no_ordination() {
        assert!(hash_pca(&[vec![1, 2], vec![2, 3]], 2).is_none());
        assert!(pcoa(&[0.0; 4], 2, 2).is_none());
    }
}
