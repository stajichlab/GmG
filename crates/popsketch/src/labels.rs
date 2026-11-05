//! Label file: a table whose first column is the sample id and whose other columns are label sets (species, lineage, source ...).
//! A blank or NA cell means unlabelled. The first label column is the primary one. Each column is compared with the clusters
//! and checked for samples whose ANI says they belong to a different label.

use crate::dist::DistMat;
use crate::stats::{adjusted_rand, median};
use anyhow::{bail, Result};
use std::collections::{BTreeMap, HashMap};

#[derive(Clone, Debug)]
pub struct LabelTable {
    pub columns: Vec<String>,
    pub rows: HashMap<String, Vec<Option<String>>>,
}

fn cell(s: &str) -> Option<String> {
    let t = s.trim();
    if t.is_empty()
        || matches!(
            t.to_ascii_uppercase().as_str(),
            "NA" | "N/A" | "NAN" | "NONE" | "."
        )
    {
        None
    } else {
        Some(t.to_string())
    }
}

/// Parse a tab-separated (or, if the header has no tab, comma-separated) table with a header.
pub fn parse(text: &str) -> Result<LabelTable> {
    let mut lines = text
        .lines()
        .filter(|l| !l.trim().is_empty() && !l.starts_with('#'));
    let Some(head) = lines.next() else {
        bail!("the label file is empty")
    };
    let sep = if head.contains('\t') { '\t' } else { ',' };
    let cols: Vec<&str> = head.split(sep).map(str::trim).collect();
    if cols.len() < 2 {
        bail!("the label file needs a header with the sample id column and at least one label column, separated by tabs");
    }
    let columns: Vec<String> = cols[1..].iter().map(|s| s.to_string()).collect();
    let mut rows = HashMap::new();
    for (n, line) in lines.enumerate() {
        let f: Vec<&str> = line.split(sep).collect();
        let id = f[0].trim();
        if id.is_empty() {
            bail!("label file row {}: empty sample id", n + 2);
        }
        let vals: Vec<Option<String>> = (0..columns.len())
            .map(|c| f.get(c + 1).and_then(|x| cell(x)))
            .collect();
        if rows.insert(id.to_string(), vals).is_some() {
            bail!("label file: sample id '{id}' appears twice");
        }
    }
    Ok(LabelTable { columns, rows })
}

#[derive(Clone, Debug)]
pub struct Conflict {
    pub id: String,
    pub label: String,
    pub own_median_ani: f64,
    pub other_label: String,
    pub other_median_ani: f64,
    pub nearest: String,
    pub nearest_label: Option<String>,
    pub nearest_ani: f64,
}

#[derive(Clone, Debug)]
pub struct ColumnReport {
    pub name: String,
    pub n_labelled: usize,
    /// label and number of included samples carrying it
    pub groups: Vec<(String, usize)>,
    /// (number of clusters, adjusted Rand index of the clusters against this label), labelled samples only
    pub ari: Vec<(usize, f64)>,
    /// label -> cluster -> count, at the clustering used for the contingency table
    pub contingency: BTreeMap<String, BTreeMap<usize, usize>>,
    pub contingency_k: Option<usize>,
    /// fraction of each label's samples that sit in its most common cluster
    pub purity: Vec<(String, f64)>,
    pub conflicts: Vec<Conflict>,
    /// labels with at least 3 samples that do not form a coherent sketch cluster: not checked for conflicts
    pub unchecked: Vec<String>,
}

/// Labels of the included samples for one column.
pub fn column_labels(t: &LabelTable, ids: &[String], col: usize) -> Vec<Option<String>> {
    ids.iter()
        .map(|id| t.rows.get(id).and_then(|v| v[col].clone()))
        .collect()
}

/// `cluster_sets`: clusterings to score, as (k, label per included sample); the first one is used for the contingency table.
pub fn analyze_column(
    name: &str,
    labels: &[Option<String>],
    ids: &[String],
    dist: &DistMat,
    cluster_sets: &[(usize, Vec<usize>)],
    margin: f64,
) -> ColumnReport {
    let n = labels.len();
    let mut counts: BTreeMap<String, usize> = BTreeMap::new();
    for l in labels.iter().flatten() {
        *counts.entry(l.clone()).or_default() += 1;
    }
    let groups: Vec<(String, usize)> = counts.iter().map(|(k, v)| (k.clone(), *v)).collect();
    let n_labelled = counts.values().sum();
    let lab_idx: HashMap<&String, usize> = counts.keys().enumerate().map(|(i, k)| (k, i)).collect();
    let labelled: Vec<usize> = (0..n).filter(|&i| labels[i].is_some()).collect();
    let mut ari = Vec::new();
    for (k, cl) in cluster_sets {
        // an index against one label group is undefined: report none
        if labelled.len() >= 2 && counts.len() >= 2 {
            let a: Vec<usize> = labelled
                .iter()
                .map(|&i| lab_idx[labels[i].as_ref().unwrap()])
                .collect();
            let b: Vec<usize> = labelled.iter().map(|&i| cl[i]).collect();
            ari.push((*k, adjusted_rand(&a, &b)));
        }
    }
    let mut contingency: BTreeMap<String, BTreeMap<usize, usize>> = BTreeMap::new();
    let mut contingency_k = None;
    if let Some((k, cl)) = cluster_sets.first() {
        contingency_k = Some(*k);
        for &i in &labelled {
            *contingency
                .entry(labels[i].clone().unwrap())
                .or_default()
                .entry(cl[i])
                .or_default() += 1;
        }
    }
    let purity = contingency
        .iter()
        .map(|(l, m)| {
            (
                l.clone(),
                *m.values().max().unwrap_or(&0) as f64 / m.values().sum::<usize>().max(1) as f64,
            )
        })
        .collect();
    // conflicts: only for labels with at least 3 samples that form a coherent sketch cluster. A label that does not follow the
    // genetic structure (a sequencing centre, a QC verdict, lineages the sketches cannot tell apart) would flag samples at random.
    let med_to = |i: usize, pred: &dyn Fn(usize) -> bool| -> f64 {
        let v: Vec<f64> = (0..n)
            .filter(|&j| j != i && pred(j))
            .map(|j| dist.ani(i, j))
            .collect();
        median(&v)
    };
    let mut big: Vec<&String> = Vec::new();
    let mut unchecked: Vec<String> = Vec::new();
    for (lab, &c) in counts.iter().filter(|(_, &c)| c >= 3) {
        let members: Vec<usize> = (0..n)
            .filter(|&i| labels[i].as_ref() == Some(lab))
            .collect();
        let agree = members
            .iter()
            .filter(|&&i| {
                let own = med_to(i, &|j| labels[j].as_ref() == Some(lab));
                let other = med_to(i, &|j| labels[j].as_ref() != Some(lab));
                own - other > margin
            })
            .count();
        if agree as f64 >= 0.6 * c as f64 {
            big.push(lab);
        } else {
            unchecked.push(lab.clone());
        }
    }
    let mut conflicts = Vec::new();
    if big.len() >= 2 {
        for &i in &labelled {
            let own = labels[i].as_ref().unwrap();
            if !big.contains(&own) {
                continue;
            }
            let med = |lab: &String| med_to(i, &|j| labels[j].as_ref() == Some(lab));
            let own_m = med(own);
            let (best, best_m) = big
                .iter()
                .filter(|l| **l != own)
                .map(|l| ((*l).clone(), med(l)))
                .fold((String::new(), f64::NEG_INFINITY), |a, b| {
                    if b.1 > a.1 {
                        b
                    } else {
                        a
                    }
                });
            if best_m - own_m > margin {
                let nn = (0..n)
                    .filter(|&j| j != i)
                    .max_by(|&a, &b| dist.ani(i, a).partial_cmp(&dist.ani(i, b)).unwrap())
                    .unwrap();
                conflicts.push(Conflict {
                    id: ids[i].clone(),
                    label: own.clone(),
                    own_median_ani: own_m,
                    other_label: best,
                    other_median_ani: best_m,
                    nearest: ids[nn].clone(),
                    nearest_label: labels[nn].clone(),
                    nearest_ani: dist.ani(i, nn),
                });
            }
        }
    }
    ColumnReport {
        name: name.to_string(),
        n_labelled,
        groups,
        ari,
        contingency,
        contingency_k,
        purity,
        conflicts,
        unchecked,
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::dist::compute;

    #[test]
    fn parses_tabs_commas_blanks_and_na() {
        let t = parse("id\tspecies\tsource\ns1\tA\t\ns2\tNA\tsoil\ns3\tB\tclinical\n").unwrap();
        assert_eq!(t.columns, vec!["species", "source"]);
        assert_eq!(t.rows["s1"], vec![Some("A".into()), None]);
        assert_eq!(t.rows["s2"], vec![None, Some("soil".into())]);
        let c = parse("id,lineage\nx,L1\ny,\n").unwrap();
        assert_eq!(c.rows["y"], vec![None]);
    }

    #[test]
    fn rejects_bad_tables() {
        assert!(parse("").is_err());
        assert!(parse("id\nx\n").is_err());
        assert!(parse("id\tg\na\tx\na\ty\n").is_err());
    }

    fn two_groups() -> (DistMat, Vec<String>) {
        // 6 samples: 0-2 share a big hash core, 3-5 another; sample 5 is really group 0's
        let mut sets: Vec<Vec<u64>> = Vec::new();
        for i in 0..6u64 {
            let g = if i < 3 || i == 5 { 0 } else { 1 };
            let mut h: Vec<u64> = (g * 100_000..g * 100_000 + 2000).collect();
            h.extend(1_000_000 + i * 10..1_000_000 + i * 10 + 5);
            sets.push(h);
        }
        let ids = (0..6).map(|i| format!("s{i}")).collect();
        (compute(&sets, 21), ids)
    }

    #[test]
    fn a_swapped_label_is_a_conflict_and_clean_samples_are_not() {
        let (d, ids) = two_groups();
        let labels: Vec<Option<String>> = ["A", "A", "A", "B", "B", "B"]
            .iter()
            .map(|s| Some(s.to_string()))
            .collect();
        let clusters = vec![(2, vec![0, 0, 0, 1, 1, 0])];
        let r = analyze_column("species", &labels, &ids, &d, &clusters, 0.001);
        assert_eq!(r.conflicts.len(), 1, "{:?}", r.conflicts);
        assert_eq!(r.conflicts[0].id, "s5");
        assert_eq!(r.conflicts[0].other_label, "A");
        assert!(r.conflicts[0].other_median_ani > r.conflicts[0].own_median_ani);
        assert_eq!(r.conflicts[0].nearest_label.as_deref(), Some("A"));
        // clusters agree with the labels except for s5
        assert!(r.ari[0].1 < 1.0 && r.ari[0].1 > 0.3);
        assert_eq!(r.contingency["A"][&0], 3);
        assert!((r.purity.iter().find(|p| p.0 == "B").unwrap().1 - 2.0 / 3.0).abs() < 1e-9);
    }

    #[test]
    fn small_and_missing_labels_are_not_judged() {
        let (d, ids) = two_groups();
        // label B has only two members; one sample is unlabelled
        let labels: Vec<Option<String>> = vec![
            Some("A".into()),
            Some("A".into()),
            Some("A".into()),
            Some("B".into()),
            Some("B".into()),
            None,
        ];
        let r = analyze_column(
            "c",
            &labels,
            &ids,
            &d,
            &[(2, vec![0, 0, 0, 1, 1, 0])],
            0.001,
        );
        assert!(r.conflicts.is_empty());
        assert_eq!(r.n_labelled, 5);
    }

    #[test]
    fn labels_that_do_not_follow_the_structure_are_not_checked() {
        let (d, ids) = two_groups();
        // alternate labels across the two real groups: neither label is a coherent cluster
        let labels: Vec<Option<String>> = ["X", "Y", "X", "Y", "X", "Y"]
            .iter()
            .map(|s| Some(s.to_string()))
            .collect();
        let r = analyze_column(
            "centre",
            &labels,
            &ids,
            &d,
            &[(2, vec![0, 0, 0, 1, 1, 0])],
            0.001,
        );
        assert!(r.conflicts.is_empty(), "{:?}", r.conflicts);
        assert_eq!(r.unchecked, vec!["X".to_string(), "Y".to_string()]);
        // a label that does follow it is checked
        let good: Vec<Option<String>> = ["A", "A", "A", "B", "B", "B"]
            .iter()
            .map(|s| Some(s.to_string()))
            .collect();
        assert!(analyze_column(
            "species",
            &good,
            &ids,
            &d,
            &[(2, vec![0, 0, 0, 1, 1, 0])],
            0.001
        )
        .unchecked
        .is_empty());
    }
}
