//! Comparison with an external tree (`--compare-tree`, Newick) and external PC scores (`--compare-pcs`, a table), for example a
//! SNP-based tree and PCA of the same samples. Only the samples that are in both are compared.
//!
//! * Tree: Spearman correlation of the tree's patristic distance with 1 - ANI (all pairs, and the closest half of the pairs,
//!   where the sketches are weakest), and for each cluster whether it is a clade of the tree (exact) and the best F1 of any
//!   clade or its complement against the cluster.
//! * PCs: cosines of the principal angles between the top-m subspaces, the Procrustes disparity of the top 3, and the
//!   Spearman correlation of pairwise distances in the top-5 PC space.

use crate::dist::DistMat;
use crate::ordination::{jacobi_eigen, Ordination};
use crate::stats::{median, quantile};
use anyhow::{bail, Result};
use std::collections::HashMap;

/// Replace the characters Newick reserves, as the tree writer does, so ids compare equal.
pub fn clean_name(name: &str) -> String {
    name.chars()
        .map(|c| {
            if matches!(
                c,
                '(' | ')' | ',' | ':' | ';' | '[' | ']' | '\'' | ' ' | '\t'
            ) {
                '_'
            } else {
                c
            }
        })
        .collect()
}

#[derive(Clone, Debug)]
pub struct ParsedTree {
    pub names: Vec<String>,
    /// node arrays: parent (usize::MAX for the root), branch length, tip index
    pub parent: Vec<usize>,
    pub len: Vec<f64>,
    pub tip: Vec<Option<usize>>,
    pub children: Vec<Vec<usize>>,
}

/// Parse a Newick tree: labels (quoted or not), branch lengths, internal labels (support values) are read and ignored.
pub fn parse_newick(text: &str) -> Result<ParsedTree> {
    let b: Vec<u8> = text
        .bytes()
        .filter(|c| !matches!(c, b'\n' | b'\r' | b'\t'))
        .collect();
    let mut t = ParsedTree {
        names: vec![],
        parent: vec![],
        len: vec![],
        tip: vec![],
        children: vec![],
    };
    let mut i = 0usize;
    fn label(b: &[u8], i: &mut usize) -> String {
        let mut s = String::new();
        if *i < b.len() && b[*i] == b'\'' {
            *i += 1;
            while *i < b.len() && b[*i] != b'\'' {
                s.push(b[*i] as char);
                *i += 1;
            }
            *i += 1;
        } else {
            while *i < b.len() && !matches!(b[*i], b'(' | b')' | b',' | b':' | b';') {
                s.push(b[*i] as char);
                *i += 1;
            }
        }
        s.trim().to_string()
    }
    fn number(b: &[u8], i: &mut usize) -> f64 {
        if *i < b.len() && b[*i] == b':' {
            *i += 1;
            let st = *i;
            while *i < b.len() && !matches!(b[*i], b'(' | b')' | b',' | b';') {
                *i += 1;
            }
            std::str::from_utf8(&b[st..*i])
                .ok()
                .and_then(|x| x.trim().parse().ok())
                .unwrap_or(0.0)
        } else {
            0.0
        }
    }
    // iterative parse with an explicit stack of open internal nodes
    let mut stack: Vec<usize> = Vec::new();
    let mut last: Option<usize> = None;
    let new_node = |t: &mut ParsedTree, parent: Option<usize>| -> usize {
        t.parent.push(parent.unwrap_or(usize::MAX));
        t.len.push(0.0);
        t.tip.push(None);
        t.children.push(vec![]);
        let id = t.parent.len() - 1;
        if let Some(p) = parent {
            t.children[p].push(id);
        }
        id
    };
    while i < b.len() {
        match b[i] {
            b'(' => {
                let n = new_node(&mut t, stack.last().copied());
                stack.push(n);
                i += 1;
                let c = new_node(&mut t, Some(n)); // first child, filled by the next token
                last = Some(c);
            }
            b',' => {
                let p = *stack
                    .last()
                    .ok_or_else(|| anyhow::anyhow!("newick: a comma outside parentheses"))?;
                let c = new_node(&mut t, Some(p));
                last = Some(c);
                i += 1;
            }
            b')' => {
                let n = stack
                    .pop()
                    .ok_or_else(|| anyhow::anyhow!("newick: unbalanced parentheses"))?;
                i += 1;
                let _support = label(&b, &mut i);
                t.len[n] = number(&b, &mut i);
                last = Some(n);
            }
            b';' => break,
            _ => {
                let c = last.ok_or_else(|| anyhow::anyhow!("newick: a label with no node"))?;
                let name = label(&b, &mut i);
                t.len[c] = number(&b, &mut i);
                if !name.is_empty() {
                    t.names.push(name);
                    t.tip[c] = Some(t.names.len() - 1);
                }
            }
        }
    }
    if !stack.is_empty() {
        bail!("newick: unbalanced parentheses");
    }
    if t.names.len() < 3 {
        bail!("newick: the tree has fewer than 3 labelled tips");
    }
    Ok(t)
}

impl ParsedTree {
    fn depth(&self, n: usize) -> Vec<(usize, f64)> {
        let mut v = vec![];
        let (mut cur, mut d) = (n, 0.0);
        loop {
            v.push((cur, d));
            if self.parent[cur] == usize::MAX {
                break;
            }
            d += self.len[cur];
            cur = self.parent[cur];
        }
        v
    }
    /// node index of each tip name (cleaned)
    pub fn tip_nodes(&self) -> HashMap<String, usize> {
        (0..self.parent.len())
            .filter_map(|n| self.tip[n].map(|t| (clean_name(&self.names[t]), n)))
            .collect()
    }
    /// patristic distance between two nodes
    pub fn patristic(&self, a: usize, b: usize) -> f64 {
        let da = self.depth(a);
        let db: HashMap<usize, f64> = self.depth(b).into_iter().collect();
        da.iter()
            .filter_map(|(n, x)| db.get(n).map(|y| x + y))
            .fold(f64::INFINITY, f64::min)
    }
}

/// Spearman correlation (average ranks for ties).
pub fn spearman(a: &[f64], b: &[f64]) -> f64 {
    fn ranks(v: &[f64]) -> Vec<f64> {
        let mut idx: Vec<usize> = (0..v.len()).collect();
        idx.sort_by(|&x, &y| v[x].partial_cmp(&v[y]).unwrap());
        let mut r = vec![0.0; v.len()];
        let mut i = 0;
        while i < idx.len() {
            let mut j = i;
            while j + 1 < idx.len() && v[idx[j + 1]] == v[idx[i]] {
                j += 1;
            }
            let avg = (i + j) as f64 / 2.0 + 1.0;
            for &k in &idx[i..=j] {
                r[k] = avg;
            }
            i = j + 1;
        }
        r
    }
    let (ra, rb) = (ranks(a), ranks(b));
    let n = a.len() as f64;
    let (ma, mb) = (ra.iter().sum::<f64>() / n, rb.iter().sum::<f64>() / n);
    let (mut sab, mut saa, mut sbb) = (0.0, 0.0, 0.0);
    for (x, y) in ra.iter().zip(&rb) {
        sab += (x - ma) * (y - mb);
        saa += (x - ma).powi(2);
        sbb += (y - mb).powi(2);
    }
    if saa == 0.0 || sbb == 0.0 {
        f64::NAN
    } else {
        sab / (saa * sbb).sqrt()
    }
}

#[derive(Clone, Debug)]
pub struct ClusterCheck {
    /// cluster number (1-based, as in the output files)
    pub cluster: usize,
    /// samples of the cluster that are in the tree
    pub n: usize,
    pub exact_clade: bool,
    pub best_f1: f64,
}

#[derive(Clone, Debug)]
pub struct TreeComparison {
    pub tree_tips: usize,
    pub shared: usize,
    pub spearman_all: f64,
    pub spearman_closest_half: f64,
    pub clusters: Vec<ClusterCheck>,
}

pub fn compare_tree(
    text: &str,
    ids: &[String],
    d: &DistMat,
    clusters: Option<&[usize]>,
) -> Result<TreeComparison> {
    let t = parse_newick(text)?;
    let nodes = t.tip_nodes();
    let shared: Vec<usize> = (0..ids.len())
        .filter(|&i| nodes.contains_key(&clean_name(&ids[i])))
        .collect();
    if shared.len() < 4 {
        bail!("only {} samples are in both the run and the tree (ids must match; characters Newick reserves are read as underscores)", shared.len());
    }
    let tn: Vec<usize> = shared
        .iter()
        .map(|&i| nodes[&clean_name(&ids[i])])
        .collect();
    let (mut pat, mut sk) = (Vec::new(), Vec::new());
    for x in 0..shared.len() {
        for y in (x + 1)..shared.len() {
            pat.push(t.patristic(tn[x], tn[y]));
            sk.push(1.0 - d.ani(shared[x], shared[y]));
        }
    }
    let cut = median(&pat);
    let close: Vec<usize> = (0..pat.len()).filter(|&i| pat[i] < cut).collect();
    let _ = quantile;
    let sp_close = spearman(
        &close.iter().map(|&i| pat[i]).collect::<Vec<_>>(),
        &close.iter().map(|&i| sk[i]).collect::<Vec<_>>(),
    );
    // clades as bitsets over the shared samples
    let m = shared.len();
    let words = m.div_ceil(64);
    let mut bits: Vec<Vec<u64>> = vec![vec![0; words]; t.parent.len()];
    let pos: HashMap<usize, usize> = tn.iter().enumerate().map(|(k, &n)| (n, k)).collect();
    for (&node, &k) in &pos {
        bits[node][k / 64] |= 1 << (k % 64);
    }
    // children are always created after their parent, so a reverse sweep folds children into parents
    for n in (0..t.parent.len()).rev() {
        let p = t.parent[n];
        if p != usize::MAX {
            let child = bits[n].clone();
            for (w, v) in bits[p].iter_mut().zip(child) {
                *w |= v;
            }
        }
    }
    let mut checks = Vec::new();
    if let Some(cl) = clusters {
        let k = cl.iter().copied().max().map_or(0, |x| x + 1);
        for g in 0..k {
            let mut target = vec![0u64; words];
            let mut size = 0usize;
            for (x, &i) in shared.iter().enumerate() {
                if cl[i] == g {
                    target[x / 64] |= 1 << (x % 64);
                    size += 1;
                }
            }
            if size == 0 || size == m {
                continue;
            }
            let (mut exact, mut best) = (false, 0.0f64);
            for b in &bits {
                let sz: usize = b.iter().map(|w| w.count_ones() as usize).sum();
                if sz == 0 || sz == m {
                    continue;
                }
                for cand in [
                    b.clone(),
                    b.iter()
                        .enumerate()
                        .map(|(w, v)| {
                            !v & if w + 1 == words && m & 63 != 0 {
                                (1u64 << (m % 64)) - 1
                            } else {
                                u64::MAX
                            }
                        })
                        .collect::<Vec<u64>>(),
                ] {
                    let csz: usize = cand.iter().map(|w| w.count_ones() as usize).sum();
                    let tp: usize = cand
                        .iter()
                        .zip(&target)
                        .map(|(a, b)| (a & b).count_ones() as usize)
                        .sum();
                    if tp == 0 {
                        continue;
                    }
                    best = best.max(2.0 * tp as f64 / (csz + size) as f64);
                    exact |= cand == target;
                }
            }
            checks.push(ClusterCheck {
                cluster: g + 1,
                n: size,
                exact_clade: exact,
                best_f1: best,
            });
        }
    }
    Ok(TreeComparison {
        tree_tips: t.names.len(),
        shared: m,
        spearman_all: spearman(&pat, &sk),
        spearman_closest_half: sp_close,
        clusters: checks,
    })
}

#[derive(Clone, Debug)]
pub struct PcsComparison {
    pub shared: usize,
    /// cosines of the principal angles between the top-m subspaces, m = 2, 3, 5 (largest first)
    pub cosines: Vec<(usize, Vec<f64>)>,
    pub procrustes_top3: Option<f64>,
    pub distance_spearman_top5: f64,
}

fn orthonormal_columns(x: &[Vec<f64>], m: usize) -> Vec<Vec<f64>> {
    // columns centred, then Gram-Schmidt; returns m columns of length n
    let n = x.len();
    let mut cols: Vec<Vec<f64>> = (0..m).map(|c| (0..n).map(|i| x[i][c]).collect()).collect();
    for col in cols.iter_mut() {
        let mean = col.iter().sum::<f64>() / n as f64;
        col.iter_mut().for_each(|v| *v -= mean);
    }
    for j in 0..m {
        for i in 0..j {
            let dot: f64 = cols[i].iter().zip(&cols[j]).map(|(a, b)| a * b).sum();
            let ci = cols[i].clone();
            cols[j].iter_mut().zip(&ci).for_each(|(a, b)| *a -= dot * b);
        }
        let norm = cols[j]
            .iter()
            .map(|v| v * v)
            .sum::<f64>()
            .sqrt()
            .max(1e-300);
        cols[j].iter_mut().for_each(|v| *v /= norm);
    }
    cols
}

/// Singular values (largest first) of the m x m matrix a^T b where a, b are lists of m columns.
fn cross_singular_values(a: &[Vec<f64>], b: &[Vec<f64>]) -> Vec<f64> {
    let m = a.len();
    let mut c = vec![0.0; m * m];
    for i in 0..m {
        for j in 0..m {
            c[i * m + j] = a[i].iter().zip(&b[j]).map(|(x, y)| x * y).sum();
        }
    }
    let mut g = vec![0.0; m * m]; // c c^T
    for i in 0..m {
        for j in 0..m {
            g[i * m + j] = (0..m).map(|k| c[i * m + k] * c[j * m + k]).sum();
        }
    }
    jacobi_eigen(&g, m)
        .0
        .into_iter()
        .map(|e| e.max(0.0).sqrt())
        .collect()
}

pub fn compare_pcs(text: &str, ids: &[String], ord: &Ordination) -> Result<PcsComparison> {
    let mut lines = text
        .lines()
        .filter(|l| !l.trim().is_empty() && !l.starts_with('#'));
    let head = lines
        .next()
        .ok_or_else(|| anyhow::anyhow!("the PC table is empty"))?;
    let sep = if head.contains('\t') { '\t' } else { ',' };
    let ncol = head.split(sep).count() - 1;
    if ncol < 2 {
        bail!("the PC table needs an id column and at least 2 numeric columns");
    }
    let mut table: HashMap<String, Vec<f64>> = HashMap::new();
    for (n, l) in lines.enumerate() {
        let f: Vec<&str> = l.split(sep).collect();
        let vals: Vec<f64> = f[1..]
            .iter()
            .map(|x| x.trim().parse::<f64>())
            .collect::<Result<_, _>>()
            .map_err(|_| anyhow::anyhow!("PC table row {}: a non-numeric value", n + 2))?;
        if vals.len() != ncol {
            bail!(
                "PC table row {}: {} values, expected {ncol}",
                n + 2,
                vals.len()
            );
        }
        table.insert(clean_name(f[0].trim()), vals);
    }
    let shared: Vec<usize> = (0..ids.len())
        .filter(|&i| table.contains_key(&clean_name(&ids[i])))
        .collect();
    if shared.len() < 6 {
        bail!(
            "only {} samples are in both the run and the PC table",
            shared.len()
        );
    }
    let ours: Vec<Vec<f64>> = shared.iter().map(|&i| ord.coords[i].clone()).collect();
    let theirs: Vec<Vec<f64>> = shared
        .iter()
        .map(|&i| table[&clean_name(&ids[i])].clone())
        .collect();
    let kmax = ord.eigenvalues.len().min(ncol);
    let mut cosines = Vec::new();
    for m in [2usize, 3, 5] {
        if m <= kmax && m < shared.len() {
            cosines.push((
                m,
                cross_singular_values(
                    &orthonormal_columns(&ours, m),
                    &orthonormal_columns(&theirs, m),
                ),
            ));
        }
    }
    // Procrustes on the top 3 (centred, unit Frobenius norm): disparity = 1 - (sum of singular values)^2
    let procrustes_top3 = (kmax >= 3).then(|| {
        let prep = |x: &[Vec<f64>]| -> Vec<Vec<f64>> {
            let n = x.len() as f64;
            let mut cols: Vec<Vec<f64>> =
                (0..3).map(|c| x.iter().map(|r| r[c]).collect()).collect();
            for col in cols.iter_mut() {
                let mean = col.iter().sum::<f64>() / n;
                col.iter_mut().for_each(|v| *v -= mean);
            }
            let norm = cols
                .iter()
                .flatten()
                .map(|v| v * v)
                .sum::<f64>()
                .sqrt()
                .max(1e-300);
            cols.iter_mut().flatten().for_each(|v| *v /= norm);
            cols
        };
        let sv = cross_singular_values(&prep(&ours), &prep(&theirs));
        1.0 - sv.iter().sum::<f64>().powi(2)
    });
    let k5 = kmax.min(5);
    let (mut da, mut db) = (Vec::new(), Vec::new());
    for x in 0..shared.len() {
        for y in (x + 1)..shared.len() {
            da.push(
                (0..k5)
                    .map(|c| (ours[x][c] - ours[y][c]).powi(2))
                    .sum::<f64>()
                    .sqrt(),
            );
            db.push(
                (0..k5)
                    .map(|c| (theirs[x][c] - theirs[y][c]).powi(2))
                    .sum::<f64>()
                    .sqrt(),
            );
        }
    }
    Ok(PcsComparison {
        shared: shared.len(),
        cosines,
        procrustes_top3,
        distance_spearman_top5: spearman(&da, &db),
    })
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::dist::compute;

    #[test]
    fn newick_with_lengths_support_values_and_quotes() {
        let t = parse_newick("((A:1,B:2)95:3,('C d':4,D:5)80:6,E:7);").unwrap();
        assert_eq!(t.names, vec!["A", "B", "C d", "D", "E"]);
        let nodes = t.tip_nodes();
        assert!((t.patristic(nodes["A"], nodes["B"]) - 3.0).abs() < 1e-12);
        assert!((t.patristic(nodes["A"], nodes["E"]) - 11.0).abs() < 1e-12);
        assert!((t.patristic(nodes["C_d"], nodes["D"]) - 9.0).abs() < 1e-12);
        assert!((t.patristic(nodes["A"], nodes["D"]) - (1.0 + 3.0 + 6.0 + 5.0)).abs() < 1e-12);
    }

    #[test]
    fn newick_errors() {
        assert!(parse_newick("((A,B);").is_err());
        assert!(parse_newick("(A,B);").is_err()); // fewer than 3 tips
        assert!(parse_newick("A,B,C;").is_err());
    }

    #[test]
    fn spearman_known_values() {
        assert!((spearman(&[1.0, 2.0, 3.0, 4.0], &[10.0, 20.0, 35.0, 90.0]) - 1.0).abs() < 1e-12);
        assert!((spearman(&[1.0, 2.0, 3.0, 4.0], &[4.0, 3.0, 2.0, 1.0]) + 1.0).abs() < 1e-12);
        // scipy.stats.spearmanr([1,2,3,4,5],[5,6,7,8,7]) = 0.8207826816681233
        assert!(
            (spearman(&[1.0, 2.0, 3.0, 4.0, 5.0], &[5.0, 6.0, 7.0, 8.0, 7.0])
                - 0.820_782_681_668_123_3)
                .abs()
                < 1e-9
        );
    }

    /// 12 samples in 3 groups whose hash sets differ by group, and a tree that has the groups as clades
    fn data() -> (Vec<String>, DistMat, Vec<usize>, String) {
        let sets: Vec<Vec<u64>> = (0..12u64)
            .map(|i| {
                let g = i / 4;
                (g * 100_000..g * 100_000 + 3000)
                    .chain(5_000_000 + i * 10..5_000_000 + i * 10 + 200)
                    .collect()
            })
            .collect();
        let ids: Vec<String> = (0..12).map(|i| format!("s{i}")).collect();
        let tree = "((s0:0.01,s1:0.01,s2:0.01,s3:0.01):0.02,(s4:0.01,s5:0.01,s6:0.01,s7:0.01):0.02,(s8:0.01,s9:0.01,s10:0.01,s11:0.01):0.02);";
        (
            ids,
            compute(&sets, 21),
            (0..12).map(|i| i / 4).collect(),
            tree.to_string(),
        )
    }

    #[test]
    fn clusters_that_are_clades_are_exact_and_a_wrong_cluster_is_not() {
        let (ids, d, cl, tree) = data();
        let c = compare_tree(&tree, &ids, &d, Some(&cl)).unwrap();
        assert_eq!((c.tree_tips, c.shared), (12, 12));
        assert_eq!(c.clusters.len(), 3);
        assert!(
            c.clusters
                .iter()
                .all(|x| x.exact_clade && (x.best_f1 - 1.0).abs() < 1e-12),
            "{:?}",
            c.clusters
        );
        assert!(c.spearman_all > 0.7, "{}", c.spearman_all);
        // swap two samples between clusters: no longer a clade, and the F1 drops below 1
        let mut bad = cl.clone();
        bad.swap(0, 4);
        let c2 = compare_tree(&tree, &ids, &d, Some(&bad)).unwrap();
        assert!(
            c2.clusters
                .iter()
                .any(|x| !x.exact_clade && x.best_f1 < 1.0),
            "{:?}",
            c2.clusters
        );
    }

    #[test]
    fn samples_missing_from_the_tree_are_skipped_and_too_few_is_an_error() {
        let (ids, d, cl, tree) = data();
        let small = "((s0:1,s1:1,s2:1):1,(s4:1,s5:1,s6:1):1,(s8:1,s9:1,s10:1):1);";
        let c = compare_tree(small, &ids, &d, Some(&cl)).unwrap();
        assert_eq!(c.shared, 9);
        assert!(c.clusters.iter().all(|x| x.exact_clade));
        assert!(compare_tree("(s0:1,s1:1,zz:1);", &ids, &d, None).is_err());
        let _ = tree;
    }

    fn ordination(coords: Vec<Vec<f64>>) -> Ordination {
        let k = coords[0].len();
        Ordination {
            coords,
            eigenvalues: vec![1.0; k],
            var_frac: vec![0.1; k],
            var_basis: "trace",
            negative_share: None,
            n_features: 1,
            method: "hash_pca",
        }
    }

    fn table(rows: &[Vec<f64>]) -> String {
        let mut s = String::from("id\tPC1\tPC2\tPC3\tPC4\tPC5\n");
        for (i, r) in rows.iter().enumerate() {
            s.push_str(&format!(
                "s{i}\t{}\n",
                r.iter()
                    .map(|v| v.to_string())
                    .collect::<Vec<_>>()
                    .join("\t")
            ));
        }
        s
    }

    fn pts() -> Vec<Vec<f64>> {
        let mut r = 7u64;
        let mut f = || {
            r = r
                .wrapping_mul(6364136223846793005)
                .wrapping_add(1442695040888963407);
            ((r >> 33) as f64 / (1u64 << 31) as f64) - 0.5
        };
        (0..14)
            .map(|_| (0..5).map(|c| f() * (5 - c) as f64).collect())
            .collect()
    }

    #[test]
    fn identical_rotated_and_scaled_pcs_match_and_unrelated_ones_do_not() {
        let ids: Vec<String> = (0..14).map(|i| format!("s{i}")).collect();
        let p = pts();
        let ord = ordination(p.clone());
        let same = compare_pcs(&table(&p), &ids, &ord).unwrap();
        assert!(
            same.cosines
                .iter()
                .all(|(_, c)| c.iter().all(|x| (x - 1.0).abs() < 1e-6)),
            "{:?}",
            same.cosines
        );
        assert!(same.procrustes_top3.unwrap().abs() < 1e-6);
        assert!((same.distance_spearman_top5 - 1.0).abs() < 1e-9);
        // a rotation of the first two axes by 40 degrees and a scale: same subspace, same configuration up to rotation
        let (c, s) = (40f64.to_radians().cos(), 40f64.to_radians().sin());
        let rot: Vec<Vec<f64>> = p
            .iter()
            .map(|r| {
                vec![
                    3.0 * (c * r[0] - s * r[1]),
                    3.0 * (s * r[0] + c * r[1]),
                    3.0 * r[2],
                    3.0 * r[3],
                    3.0 * r[4],
                ]
            })
            .collect();
        let r = compare_pcs(&table(&rot), &ids, &ord).unwrap();
        assert!(
            r.cosines
                .iter()
                .all(|(_, c)| c.iter().all(|x| (x - 1.0).abs() < 1e-6)),
            "{:?}",
            r.cosines
        );
        assert!(r.procrustes_top3.unwrap() < 1e-6);
        // unrelated scores
        let mut q = 99u64;
        let other: Vec<Vec<f64>> = (0..14)
            .map(|_| {
                (0..5)
                    .map(|_| {
                        q = q.wrapping_mul(6364136223846793005).wrapping_add(1);
                        ((q >> 33) as f64 / (1u64 << 31) as f64) - 0.5
                    })
                    .collect()
            })
            .collect();
        let u = compare_pcs(&table(&other), &ids, &ord).unwrap();
        assert!(
            u.cosines
                .iter()
                .find(|x| x.0 == 2)
                .unwrap()
                .1
                .iter()
                .any(|x| *x < 0.9),
            "{:?}",
            u.cosines
        );
        assert!(u.procrustes_top3.unwrap() > 0.3);
    }

    #[test]
    fn bad_pc_tables_are_refused() {
        let ids: Vec<String> = (0..14).map(|i| format!("s{i}")).collect();
        let ord = ordination(pts());
        assert!(compare_pcs("id\tPC1\ns0\t1\n", &ids, &ord).is_err());
        assert!(compare_pcs("id\tPC1\tPC2\ns0\t1\tx\n", &ids, &ord).is_err());
        assert!(compare_pcs("id\tPC1\tPC2\nzz\t1\t2\n", &ids, &ord).is_err());
    }
}
