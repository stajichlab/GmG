//! The whole run: discover samples, sketch them, audit them against the cohort, and (for the samples that pass) compute
//! distances, an ordination, clusters, a tree and the label comparison. Writes the TSV/JSON outputs and the report.

use crate::advisor::{advise, Advice};
use crate::cluster::{cut_k, gap_rule, upgma, GapCall, Verdict};
use crate::dist::{compute, DistMat};
use crate::floor::{estimate, pick_for_halves, Floor};
use crate::labels::{
    analyze_column, column_labels, parse as parse_labels, ColumnReport, LabelTable,
};
use crate::ordination::{hash_pca, pcoa, Ordination};
use crate::qc::{assess, kept_hashes, Flag, QcParams, SampleQc};
use crate::sketch_io::{
    cache_init, cache_read, cache_ready, cache_write, sketch_sample, SampleSketch, SketchParams,
};
use crate::stats::{mad_sd, median};
use crate::tree::{collapse, neighbor_joining, newick, root_tree};
use anyhow::{bail, Context, Result};
use gmg_core::classify::{classify, Context as RefContext, Item};
use gmg_core::{max_hash_for_scaled, RefIndex, SeqStats, SketchFile, Thresholds, DEFAULT_SEED};
use gmg_input::{discover, DiscoverArgs};
use rayon::prelude::*;
use serde_json::{json, Value};
use std::collections::{HashMap, HashSet};
use std::fmt::Write as _;
use std::fs;
use std::path::PathBuf;

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum FlagMode {
    /// flagged samples are left out of distances, PCA, clusters and tree (default)
    Exclude,
    /// flagged samples stay in and are drawn with their own marker
    Mark,
    /// contaminated samples are cleaned with the cohort vote and put back; other flagged samples stay out
    Decontaminate,
}

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum ReportFormat {
    Html,
    Md,
    All,
}

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum OrdKind {
    HashPca,
    Pcoa,
}

#[derive(Clone, Debug)]
pub struct Options {
    pub inputs: Vec<PathBuf>,
    pub sample_sheet: Option<PathBuf>,
    pub labels: Option<PathBuf>,
    pub out: PathBuf,
    /// `auto`, `assembly` or `reads`
    pub kind: String,
    pub no_group: bool,
    pub ksize: u32,
    pub scaled: u32,
    pub max_bases: u64,
    pub r1_only: bool,
    /// stop reading a sample once its k-mer coverage reaches this (`--max-bases` stays as a ceiling)
    pub target_cov: Option<f64>,
    pub min_cov: f64,
    pub max_samples: usize,
    pub flagged: FlagMode,
    pub index: Option<PathBuf>,
    pub ordination: OrdKind,
    pub components: usize,
    /// PCs used by DAPC on the primary label column (0: no DAPC); needs a label file
    pub dapc_pcs: usize,
    /// force this many clusters
    pub k: Option<usize>,
    /// cut the clustering at this distance (1 - ANI)
    pub cut: Option<f64>,
    /// the largest merge-height jump is called clear when it is at least this many noise floors
    pub min_floors: f64,
    pub tree: bool,
    pub root: Option<String>,
    /// collapse tree branches shorter than this (default: the noise floor)
    pub tree_floor: Option<f64>,
    /// number of samples that get split-half sketches for the noise floor (0: use the fixed fallback)
    pub floor_samples: usize,
    /// Newick tree of (some of) the same samples, for example from SNPs: compared with the sketch distances and clusters
    pub compare_tree: Option<PathBuf>,
    /// PC scores of (some of) the same samples, for example from SNPs: compared with the sketch PCA
    pub compare_pcs: Option<PathBuf>,
    pub resume: bool,
    /// `report.html` (default), `report.md` with `figures/*.svg`, or both
    pub report_format: ReportFormat,
    pub quiet: bool,
    pub threads: Option<usize>,
}

impl Default for Options {
    fn default() -> Self {
        Options {
            inputs: vec![],
            sample_sheet: None,
            labels: None,
            out: PathBuf::from("popsketch_out"),
            kind: "auto".into(),
            no_group: false,
            ksize: 21,
            scaled: 200,
            max_bases: 1_000_000_000,
            r1_only: false,
            target_cov: None,
            min_cov: 4.0,
            max_samples: 1500,
            flagged: FlagMode::Exclude,
            index: None,
            ordination: OrdKind::HashPca,
            components: 5,
            dapc_pcs: 10,
            k: None,
            cut: None,
            min_floors: 3.0,
            tree: true,
            root: None,
            tree_floor: None,
            floor_samples: 20,
            compare_tree: None,
            compare_pcs: None,
            resume: false,
            report_format: ReportFormat::Html,
            quiet: false,
            threads: None,
        }
    }
}

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Status {
    Included,
    /// flagged and left out
    Excluded,
    /// flagged and kept (`--flagged mark`)
    Marked,
    /// contamination removed with the cohort vote, then included
    Decontaminated,
    /// could not be read
    Error,
}

impl Status {
    pub fn as_str(self) -> &'static str {
        match self {
            Status::Included => "included",
            Status::Excluded => "excluded",
            Status::Marked => "marked",
            Status::Decontaminated => "decontaminated",
            Status::Error => "error",
        }
    }
}

#[derive(Clone, Debug, Default)]
pub struct IndexCall {
    pub genus: Option<String>,
    pub species: Option<String>,
    /// second taxon and the fraction of hashes it explains, when the index suspects a mixture
    pub mixture: Option<(String, f64)>,
}

#[derive(Clone, Debug)]
pub struct SampleRow {
    pub id: String,
    pub kind: String,
    pub files: Vec<String>,
    pub status: Status,
    pub qc: Option<SampleQc>,
    pub index: Option<IndexCall>,
    pub message: Option<String>,
}

impl SampleRow {
    pub fn flags(&self) -> &[Flag] {
        self.qc.as_ref().map_or(&[], |q| &q.flags)
    }
}

/// DAPC on the primary label column.
pub struct DapcOut {
    pub column: String,
    pub dapc: crate::dapc::Dapc,
    /// label of each included sample (None: unlabelled)
    pub labels: Vec<Option<String>>,
}

pub struct Results {
    pub rows: Vec<SampleRow>,
    /// indices into `rows` of the samples in the analysis, in matrix order
    pub included: Vec<usize>,
    pub dist: Option<DistMat>,
    pub floor: Floor,
    pub gap: Option<GapCall>,
    pub k_used: Option<usize>,
    /// cluster of each included sample
    pub clusters: Option<Vec<usize>>,
    pub ordination: Option<Ordination>,
    pub newick: Option<String>,
    pub collapsed: usize,
    pub label_table: Option<LabelTable>,
    pub label_reports: Vec<ColumnReport>,
    pub tree_comparison: Option<crate::compare::TreeComparison>,
    pub pcs_comparison: Option<crate::compare::PcsComparison>,
    pub dapc: Option<DapcOut>,
    pub warnings: Vec<String>,
    pub advice: Option<Advice>,
}

fn flag(q: &mut SampleQc, code: &'static str, detail: String) {
    q.flags.push(Flag { code, detail });
}

/// Classify one sample against the index (downsampled to the index's scale); None when it cannot be done.
fn screen(ix: &RefIndex, th: &Thresholds, ctx: &RefContext, s: &SampleSketch) -> Option<IndexCall> {
    let mh = max_hash_for_scaled(ix.scaled);
    let (mut hashes, mut abunds) = (Vec::new(), Vec::new());
    for (&h, &a) in s.spectrum.hashes.iter().zip(&s.spectrum.abunds) {
        if h <= mh {
            hashes.push(h);
            abunds.push(a);
        }
    }
    let sk = SketchFile {
        ksize: ix.ksize,
        scaled: ix.scaled,
        seed: DEFAULT_SEED,
        name: s.id.clone(),
        hashes,
        abunds: Some(abunds),
    };
    let st = SeqStats {
        n_records: s.counts.reads,
        total_len: s.counts.bases,
        non_acgt: s.counts.non_acgt,
        n50: 0,
    };
    let reads = s.kind == "reads";
    let item = Item {
        id: &s.id,
        path: "",
        size_bytes: 0,
        kind: &s.kind,
    };
    let (rec, _) = classify(
        ix,
        th,
        ctx,
        &item,
        sk,
        &st,
        // a fixed cutoff of 2 for reads, as `gmg run` uses: the auto cutoff drops the low-abundance k-mers of a light contaminant
        if reads { 2 } else { 1 },
        1,
        &|_| false,
    )
    .ok()?;
    let taxon = |r: &str| rec["guess"]["calls"][r]["taxon"].as_str().map(String::from);
    let mixture = (rec["qc"]["mixture"]["suspected"] == true).then(|| {
        (
            rec["qc"]["mixture"]["second"]["taxon"]
                .as_str()
                .unwrap_or("?")
                .to_string(),
            rec["qc"]["mixture"]["second"]["explained_fraction"]
                .as_f64()
                .unwrap_or(0.0),
        )
    });
    Some(IndexCall {
        genus: taxon("genus"),
        species: taxon("species"),
        mixture,
    })
}

fn say(o: &Options, msg: &str) {
    if !o.quiet {
        eprintln!("{msg}");
    }
}

pub fn run(o: &Options) -> Result<Results> {
    if let Some(t) = o.threads {
        rayon::ThreadPoolBuilder::new()
            .num_threads(t)
            .build_global()
            .ok();
    }
    let inputs = discover(&DiscoverArgs {
        inputs: &o.inputs,
        list: None,
        group: !o.no_group,
        quiet: o.quiet,
        sample_sheet: o.sample_sheet.as_deref(),
        kind: &o.kind,
    })?;
    if inputs.len() > o.max_samples {
        bail!(
            "{} samples is more than --max-samples {}: exact all-pairs distances grow with the square of the number of samples. \
             Raise --max-samples, use a coarser --scaled, or split the set into groups",
            inputs.len(),
            o.max_samples
        );
    }
    fs::create_dir_all(&o.out).with_context(|| o.out.display().to_string())?;
    let params = SketchParams {
        ksize: o.ksize,
        scaled: o.scaled,
        max_bases: o.max_bases,
        r1_only: o.r1_only,
        target_cov: o.target_cov,
        ..SketchParams::default()
    };
    let read_ids: Vec<String> = inputs
        .iter()
        .filter(|i| i.kind == "reads")
        .map(|i| i.id.clone())
        .collect();
    let halves_for = pick_for_halves(&read_ids, o.floor_samples);
    let cache_dir = o.out.join("sketches");
    let reuse = o.resume && cache_ready(&cache_dir, &params);
    if !reuse {
        cache_init(&cache_dir, &params)?;
    }
    say(
        o,
        &format!(
            "sketching {} samples (k={}, scaled={}, {})",
            inputs.len(),
            o.ksize,
            o.scaled,
            if o.max_bases == 0 {
                "all reads".to_string()
            } else {
                format!("up to {} bases each", o.max_bases)
            }
        ),
    );
    let sketched: Vec<Result<SampleSketch>> = inputs
        .par_iter()
        .map(|inp| {
            let want = halves_for.contains(&inp.id);
            if reuse {
                if let Some(s) = cache_read(&cache_dir, &inp.id, want) {
                    return Ok(s);
                }
            }
            let s = sketch_sample(&inp.id, inp.kind, &inp.files, &params, want)?;
            cache_write(&cache_dir, &s, &params)?;
            Ok(s)
        })
        .collect();
    let mut rows: Vec<SampleRow> = Vec::new();
    let mut sketches: Vec<SampleSketch> = Vec::new();
    let mut warnings: Vec<String> = Vec::new();
    for (inp, r) in inputs.iter().zip(sketched) {
        match r {
            Ok(s) => sketches.push(s),
            Err(e) => {
                warnings.push(format!("{}: {e:#}", inp.id));
                rows.push(SampleRow {
                    id: inp.id.clone(),
                    kind: inp.kind.to_string(),
                    files: inp.files.iter().map(|f| f.display().to_string()).collect(),
                    status: Status::Error,
                    qc: None,
                    index: None,
                    message: Some(format!("{e:#}")),
                });
            }
        }
    }
    if sketches.is_empty() {
        bail!("no sample could be read");
    }
    // ---- audit
    let qp = QcParams {
        min_cov: o.min_cov,
        ..QcParams::default()
    };
    let mut qcs = assess(&sketches, o.scaled, &qp);
    let mut idx_calls: Vec<Option<IndexCall>> = vec![None; sketches.len()];
    if let Some(path) = &o.index {
        let bytes = fs::read(path).with_context(|| format!("opening {}", path.display()))?;
        let ix = RefIndex::read_owned(bytes)
            .with_context(|| format!("reading index {}", path.display()))?;
        if ix.ksize != o.ksize || ix.scaled < o.scaled {
            warnings.push(format!(
                "index screen skipped: the index has k={} scaled={}; samples are sketched at k={} scaled={} (the index must use the same k and a scale at least as coarse)",
                ix.ksize, ix.scaled, o.ksize, o.scaled
            ));
        } else {
            let th = Thresholds::default();
            let ctx = RefContext::new(&ix);
            idx_calls = sketches
                .par_iter()
                .map(|s| screen(&ix, &th, &ctx, s))
                .collect();
            let mut modes: HashMap<&str, usize> = HashMap::new();
            for c in idx_calls.iter().flatten() {
                if let Some(g) = &c.genus {
                    *modes.entry(g.as_str()).or_default() += 1;
                }
            }
            let called: usize = modes.values().sum();
            let top = modes
                .iter()
                .max_by_key(|(g, n)| (**n, std::cmp::Reverse(**g)))
                .map(|(g, n)| (g.to_string(), *n));
            for (q, c) in qcs.iter_mut().zip(&idx_calls) {
                let Some(c) = c else { continue };
                if let Some((m, f)) = &c.mixture {
                    flag(
                        q,
                        "mixture_suspected",
                        format!(
                            "the index finds a second taxon, {m}, explaining {:.1}% of the hashes",
                            f * 100.0
                        ),
                    );
                }
                if let (Some((tg, tn)), Some(g)) = (&top, &c.genus) {
                    if g != tg && tn * 2 >= called {
                        flag(q, "different_organism", format!("the index calls this sample {g}; {tn} of {called} called samples are {tg}"));
                    }
                }
            }
        }
    }
    let noted = qcs
        .iter()
        .filter(|q| q.notes.iter().any(|n| n.code == "abundance_shoulder"))
        .count();
    if noted > 0 {
        warnings.push(format!(
            "{noted} samples have an abundance shoulder that is widespread in this cohort: more likely heterozygosity or varying ploidy than contamination, so they are noted but not excluded"
        ));
    }
    // ---- distances among everything with enough hashes, then the unrelated-to-cohort check
    const MIN_HASHES: usize = 500;
    let mut kept: Vec<Vec<u64>> = sketches
        .iter()
        .zip(&qcs)
        .map(|(s, q)| kept_hashes(s, q.cutoff))
        .collect();
    for (q, k) in qcs.iter_mut().zip(&kept) {
        if k.len() < MIN_HASHES {
            flag(
                q,
                "too_few_hashes",
                format!(
                    "only {} hashes after the abundance cutoff (need {MIN_HASHES})",
                    k.len()
                ),
            );
        }
    }
    let cand: Vec<usize> = (0..sketches.len())
        .filter(|&i| kept[i].len() >= MIN_HASHES)
        .collect();
    let cand_sets: Vec<Vec<u64>> = cand.iter().map(|&i| kept[i].clone()).collect();
    let d_all = compute(&cand_sets, o.ksize);
    let unflagged: Vec<usize> = (0..cand.len())
        .filter(|&c| !qcs[cand[c]].flagged())
        .collect();
    // nearest neighbour among the clean cohort only: a contaminated sample can carry the same foreign hashes as an unrelated one
    let mut clean = vec![false; cand.len()];
    for &c in &unflagged {
        clean[c] = true;
    }
    let nn: Vec<f64> = (0..cand.len())
        .map(|c| {
            (0..cand.len())
                .filter(|&j| j != c && clean[j])
                .map(|j| d_all.ani(c, j))
                .fold(0.0, f64::max)
        })
        .collect();
    if unflagged.len() >= QcParams::default().min_cohort {
        let vals: Vec<f64> = unflagged.iter().map(|&c| nn[c]).collect();
        let (med, spread) = (median(&vals), mad_sd(&vals).max(0.005));
        for c in 0..cand.len() {
            let z = (nn[c] - med) / spread;
            if z < -8.0 {
                let q = &mut qcs[cand[c]];
                flag(q, "unrelated_to_cohort", format!("nearest neighbour ANI {:.4} vs cohort median {:.4} (z {z:+.1}): shares almost nothing with the other samples", nn[c], med));
            }
        }
    }
    // ---- decontamination (opt-in): contamination flags only, by the cohort vote
    let mut statuses: Vec<Status> = vec![Status::Included; sketches.len()];
    let contam = |q: &SampleQc| {
        q.flags.iter().all(|f| {
            matches!(
                f.code,
                "size_high" | "abundance_shoulder" | "mixture_suspected"
            )
        }) && q.flagged()
    };
    if o.flagged == FlagMode::Decontaminate {
        let clean: Vec<usize> = (0..sketches.len())
            .filter(|&i| !qcs[i].flagged() && kept[i].len() >= MIN_HASHES)
            .collect();
        if clean.len() >= 3 {
            let mut freq: HashMap<u64, u32> = HashMap::new();
            for &i in &clean {
                for &h in &kept[i] {
                    *freq.entry(h).or_default() += 1;
                }
            }
            for i in 0..sketches.len() {
                if contam(&qcs[i]) && kept[i].len() >= MIN_HASHES {
                    kept[i].retain(|h| freq.get(h).is_some_and(|&c| c >= 3));
                    statuses[i] = Status::Decontaminated;
                }
            }
        } else {
            warnings.push("--flagged decontaminate needs at least 3 unflagged samples for the vote; flagged samples were excluded".into());
        }
    }
    for i in 0..sketches.len() {
        let q = &qcs[i];
        if q.flagged() && statuses[i] != Status::Decontaminated {
            statuses[i] = if o.flagged == FlagMode::Mark && kept[i].len() >= MIN_HASHES {
                Status::Marked
            } else {
                Status::Excluded
            };
        }
    }
    // a decontaminated sample that ended up too small is excluded
    for i in 0..sketches.len() {
        if statuses[i] == Status::Decontaminated && kept[i].len() < MIN_HASHES {
            statuses[i] = Status::Excluded;
        }
    }
    let incl_s: Vec<usize> = (0..sketches.len())
        .filter(|&i| {
            matches!(
                statuses[i],
                Status::Included | Status::Marked | Status::Decontaminated
            )
        })
        .collect();
    // noise floor from clean samples that carry halves
    let half_refs: Vec<&SampleSketch> = (0..sketches.len())
        .filter(|&i| !qcs[i].flagged() && sketches[i].halves.is_some())
        .map(|i| &sketches[i])
        .collect();
    let floor = estimate(&half_refs, o.ksize);
    say(
        o,
        &format!(
            "{} samples, {} in the analysis; noise floor {:.5} ({})",
            sketches.len(),
            incl_s.len(),
            floor.value,
            floor.source
        ),
    );

    let mut res = Results {
        rows: Vec::new(),
        included: Vec::new(),
        dist: None,
        floor: floor.clone(),
        gap: None,
        k_used: None,
        clusters: None,
        ordination: None,
        newick: None,
        collapsed: 0,
        label_table: None,
        label_reports: Vec::new(),
        tree_comparison: None,
        pcs_comparison: None,
        dapc: None,
        warnings,
        advice: None,
    };
    // rows: errors first were pushed in `rows`; add the sketched ones, then order by id
    let mut all_rows = rows;
    for (i, s) in sketches.iter().enumerate() {
        all_rows.push(SampleRow {
            id: s.id.clone(),
            kind: s.kind.clone(),
            files: s.files.iter().map(|f| f.display().to_string()).collect(),
            status: statuses[i],
            qc: Some(qcs[i].clone()),
            index: idx_calls[i].clone(),
            message: None,
        });
    }
    // map sketch index -> row index after sorting
    let mut order: Vec<usize> = (0..all_rows.len()).collect();
    order.sort_by(|&a, &b| all_rows[a].id.cmp(&all_rows[b].id));
    let n_err = all_rows
        .iter()
        .filter(|r| r.status == Status::Error)
        .count();
    let sketch_to_row: HashMap<usize, usize> = (0..sketches.len())
        .map(|i| (i, order.iter().position(|&x| x == n_err + i).unwrap()))
        .collect();
    res.rows = order.iter().map(|&x| all_rows[x].clone()).collect();
    res.included = {
        let mut v: Vec<usize> = incl_s.iter().map(|i| sketch_to_row[i]).collect();
        v.sort_unstable();
        v
    };
    let inv: HashMap<usize, usize> = (0..sketches.len())
        .map(|i| (sketch_to_row[&i], i))
        .collect();
    let n = res.included.len();
    if n >= 3 {
        let sets: Vec<Vec<u64>> = res.included.iter().map(|r| kept[inv[r]].clone()).collect();
        let decon_changed = res
            .included
            .iter()
            .any(|r| statuses[inv[r]] == Status::Decontaminated);
        let dist = if decon_changed {
            compute(&sets, o.ksize)
        } else {
            let pos: HashMap<usize, usize> =
                cand.iter().enumerate().map(|(c, &s)| (s, c)).collect();
            d_all.subset(
                &res.included
                    .iter()
                    .map(|r| pos[&inv[r]])
                    .collect::<Vec<_>>(),
            )
        };
        // clusters
        let d1 = dist.one_minus_ani();
        let merges = upgma(&d1, n);
        let gap = gap_rule(&merges, n, o.min_floors, floor.value, 30);
        let k_used = if let Some(k) = o.k {
            Some(k.clamp(1, n))
        } else if let Some(c) = o.cut {
            Some(merges.iter().filter(|m| m.height > c).count() + 1)
        } else if let Verdict::Clear(k) = gap.verdict {
            Some(k)
        } else if gap.verdict == Verdict::Single {
            Some(1)
        } else {
            None
        };
        let clusters = k_used.map(|k| cut_k(&merges, n, k));
        // ordination: DAPC needs more PCs than are plotted
        let ncomp = if o.labels.is_some() && o.dapc_pcs > 0 {
            o.components.max(o.dapc_pcs)
        } else {
            o.components
        };
        res.ordination = match o.ordination {
            OrdKind::HashPca => hash_pca(&sets, ncomp),
            OrdKind::Pcoa => pcoa(&d1, n, ncomp),
        };
        // tree
        if o.tree {
            let names: Vec<String> = res
                .included
                .iter()
                .map(|&r| res.rows[r].id.clone())
                .collect();
            let mut t = neighbor_joining(&dist.mash, n, names.clone());
            res.collapsed = collapse(&mut t, o.tree_floor.unwrap_or(floor.value));
            if let Some(r) = &o.root {
                if !names.contains(r) {
                    res.warnings.push(format!(
                        "--root {r}: no such sample in the analysis; the tree is unrooted"
                    ));
                }
            }
            res.newick = Some(newick(&root_tree(&t, o.root.as_deref()), &names));
        }
        // labels
        if let Some(lp) = &o.labels {
            let text =
                fs::read_to_string(lp).with_context(|| format!("reading {}", lp.display()))?;
            let table = parse_labels(&text)?;
            let ids: Vec<String> = res
                .included
                .iter()
                .map(|&r| res.rows[r].id.clone())
                .collect();
            let all_ids: HashSet<&str> = res.rows.iter().map(|r| r.id.as_str()).collect();
            let unknown: Vec<&String> = table
                .rows
                .keys()
                .filter(|k| !all_ids.contains(k.as_str()))
                .collect();
            if !unknown.is_empty() {
                res.warnings.push(format!(
                    "{} ids in the label file are not samples in this run (first: {})",
                    unknown.len(),
                    unknown[0]
                ));
            }
            let margin = floor.value.max(0.001);
            for (c, name) in table.columns.iter().enumerate() {
                let labs = column_labels(&table, &ids, c);
                let distinct: HashSet<&String> = labs.iter().flatten().collect();
                let mut sets: Vec<(usize, Vec<usize>)> = Vec::new();
                if let (Some(k), Some(cl)) = (k_used, &clusters) {
                    sets.push((k, cl.clone()));
                }
                for k in gap
                    .candidates
                    .iter()
                    .map(|c| c.k)
                    .chain(std::iter::once(distinct.len().max(2)))
                {
                    if k <= n && !sets.iter().any(|s| s.0 == k) {
                        sets.push((k, cut_k(&merges, n, k)));
                    }
                }
                res.label_reports
                    .push(analyze_column(name, &labs, &ids, &dist, &sets, margin));
            }
            // DAPC on the primary label column
            if let (Some(ord), true, Some(rep)) =
                (&res.ordination, o.dapc_pcs > 0, res.label_reports.first())
            {
                let labs = column_labels(&table, &ids, 0);
                let names: Vec<String> = rep.groups.iter().map(|g| g.0.clone()).collect();
                let group: Vec<Option<usize>> = labs
                    .iter()
                    .map(|l| l.as_ref().and_then(|l| names.iter().position(|x| x == l)))
                    .collect();
                let p = o.dapc_pcs.min(ord.eigenvalues.len());
                let pcs: Vec<Vec<f64>> = ord.coords.iter().map(|r| r[..p].to_vec()).collect();
                match crate::dapc::dapc(&pcs, &group, &names) {
                    Some(d) => res.dapc = Some(DapcOut { column: rep.name.clone(), dapc: d, labels: labs }),
                    None => res.warnings.push(format!("DAPC on label column {} was skipped: it needs at least 2 label groups with 2 samples each and more samples than PCs", rep.name)),
                }
            }
            res.label_table = Some(table);
        }
        // comparison with an external SNP tree and PCs
        let names: Vec<String> = res
            .included
            .iter()
            .map(|&r| res.rows[r].id.clone())
            .collect();
        if let Some(p) = &o.compare_tree {
            match fs::read_to_string(p)
                .with_context(|| format!("reading {}", p.display()))
                .and_then(|t| crate::compare::compare_tree(&t, &names, &dist, clusters.as_deref()))
            {
                Ok(c) => res.tree_comparison = Some(c),
                Err(e) => res.warnings.push(format!("--compare-tree: {e:#}")),
            }
        }
        if let (Some(p), Some(ord)) = (&o.compare_pcs, &res.ordination) {
            match fs::read_to_string(p)
                .with_context(|| format!("reading {}", p.display()))
                .and_then(|t| crate::compare::compare_pcs(&t, &names, ord))
            {
                Ok(c) => res.pcs_comparison = Some(c),
                Err(e) => res.warnings.push(format!("--compare-pcs: {e:#}")),
            }
        }
        res.advice = Some(advise(&dist, clusters.as_deref(), &floor));
        res.gap = Some(gap);
        res.k_used = k_used;
        res.clusters = clusters;
        res.dist = Some(dist);
    } else {
        res.warnings.push(format!(
            "only {n} samples passed the audit: no distances, ordination or tree"
        ));
    }
    write_outputs(o, &res)?;
    crate::report::write(o, &res)?;
    Ok(res)
}

fn tsv_cell(s: &str) -> String {
    s.replace(['\t', '\n', '\r'], " ")
}

fn write_matrix(path: PathBuf, ids: &[String], m: &[f64], n: usize, prec: usize) -> Result<()> {
    let mut s = String::from("id");
    for id in ids {
        let _ = write!(s, "\t{id}");
    }
    s.push('\n');
    for i in 0..n {
        s.push_str(&ids[i]);
        for j in 0..n {
            let _ = write!(s, "\t{:.prec$}", m[i * n + j]);
        }
        s.push('\n');
    }
    fs::write(&path, s).with_context(|| path.display().to_string())
}

fn write_outputs(o: &Options, r: &Results) -> Result<()> {
    // qc.tsv
    let mut s = String::from("id\tkind\tstatus\tflags\treasons\tnotes\treads\tbases\thashes\tabundance_cutoff\tkept_hashes\test_size_mb\tkmer_coverage\tcoverage_basis\tsize_z\tshoulder_z\torganism_call\tsecond_taxon\tsecond_fraction\tfiles\n");
    let f = |x: Option<f64>, p: usize| x.map_or(String::new(), |v| format!("{v:.p$}"));
    for row in &r.rows {
        let (flags, reasons) = (
            row.flags()
                .iter()
                .map(|x| x.code)
                .collect::<Vec<_>>()
                .join(";"),
            row.flags()
                .iter()
                .map(|x| x.detail.as_str())
                .collect::<Vec<_>>()
                .join("; "),
        );
        let ic = row.index.as_ref();
        match &row.qc {
            Some(q) => {
                let _ = writeln!(
                    s,
                    "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{:.2}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
                    row.id,
                    row.kind,
                    row.status.as_str(),
                    flags,
                    tsv_cell(&reasons),
                    tsv_cell(&q.notes.iter().map(|x| x.detail.as_str()).collect::<Vec<_>>().join("; ")),
                    q.reads,
                    q.bases,
                    q.n_hashes,
                    q.cutoff,
                    q.kept,
                    q.est_mb,
                    f(q.kmer_cov, 1),
                    q.cov_basis,
                    f(q.size_z, 1),
                    f(q.shoulder_z, 1),
                    ic.and_then(|c| c.species.clone().or(c.genus.clone()))
                        .unwrap_or_default(),
                    ic.and_then(|c| c.mixture.as_ref().map(|m| m.0.clone()))
                        .unwrap_or_default(),
                    f(ic.and_then(|c| c.mixture.as_ref().map(|m| m.1)), 3),
                    tsv_cell(&row.files.join(","))
                );
            }
            None => {
                let _ = writeln!(
                    s,
                    "{}\t{}\terror\t\t{}\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t{}",
                    row.id,
                    row.kind,
                    tsv_cell(row.message.as_deref().unwrap_or("")),
                    tsv_cell(&row.files.join(","))
                );
            }
        }
    }
    fs::write(o.out.join("qc.tsv"), s)?;
    let ids: Vec<String> = r.included.iter().map(|&i| r.rows[i].id.clone()).collect();
    if let Some(d) = &r.dist {
        write_matrix(o.out.join("ani.tsv"), &ids, &d.ani, d.n, 5)?;
        write_matrix(o.out.join("mash.tsv"), &ids, &d.mash, d.n, 6)?;
    }
    // pca.tsv, clusters.tsv
    if let Some(ord) = &r.ordination {
        let k = ord.eigenvalues.len();
        let mut s = String::from("id\tstatus");
        for c in 0..k {
            let _ = write!(s, "\tPC{}", c + 1);
        }
        s.push_str("\tcluster\n");
        for (i, &row) in r.included.iter().enumerate() {
            let _ = write!(s, "{}\t{}", r.rows[row].id, r.rows[row].status.as_str());
            for c in 0..k {
                let _ = write!(s, "\t{:.5}", ord.coords[i][c]);
            }
            let cl = r
                .clusters
                .as_ref()
                .map_or("NA".to_string(), |c| (c[i] + 1).to_string());
            let _ = writeln!(s, "\t{cl}");
        }
        fs::write(o.out.join("pca.tsv"), s)?;
        let var: Vec<String> = ord.var_frac.iter().map(|v| format!("{v:.4}")).collect();
        fs::write(
            o.out.join("pca_variance.txt"),
            format!(
                "method\t{}\nvariance_fraction\t{}\nfeatures\t{}\n",
                ord.method,
                var.join(" "),
                ord.n_features
            ),
        )?;
    }
    let mut s = String::from("id\tcluster\n");
    for (i, &row) in r.included.iter().enumerate() {
        let cl = r
            .clusters
            .as_ref()
            .map_or("NA".to_string(), |c| (c[i] + 1).to_string());
        let _ = writeln!(s, "{}\t{cl}", r.rows[row].id);
    }
    fs::write(o.out.join("clusters.tsv"), s)?;
    if let Some(dp) = &r.dapc {
        let axes = dp.dapc.eigenvalues.len();
        let mut s = String::from("id\tlabel\tassigned\tmatches_label");
        for a in 0..axes {
            let _ = write!(s, "\tLD{}", a + 1);
        }
        s.push('\n');
        for (i, &row) in r.included.iter().enumerate() {
            let label = dp.labels[i].clone().unwrap_or_default();
            let (asg, ok) = match dp.dapc.assigned[i] {
                Some(g) => (
                    dp.dapc.groups[g].clone(),
                    if dp.labels[i].as_deref() == Some(dp.dapc.groups[g].as_str()) {
                        "1"
                    } else {
                        "0"
                    },
                ),
                None => (String::new(), ""),
            };
            let _ = write!(s, "{}\t{label}\t{asg}\t{ok}", r.rows[row].id);
            for a in 0..axes {
                let _ = write!(s, "\t{:.4}", dp.dapc.scores[i][a]);
            }
            s.push('\n');
        }
        fs::write(o.out.join("dapc.tsv"), s)?;
    }
    if let Some(nw) = &r.newick {
        fs::write(o.out.join("tree.nwk"), format!("{nw}\n"))?;
    }
    // label_check.tsv
    if !r.label_reports.is_empty() {
        let mut s = String::from("label_column\tid\tlabel\town_median_ani\tbetter_matching_label\tits_median_ani\tnearest_sample\tnearest_label\tnearest_ani\n");
        for rep in &r.label_reports {
            for c in &rep.conflicts {
                let _ = writeln!(
                    s,
                    "{}\t{}\t{}\t{:.4}\t{}\t{:.4}\t{}\t{}\t{:.4}",
                    rep.name,
                    c.id,
                    c.label,
                    c.own_median_ani,
                    c.other_label,
                    c.other_median_ani,
                    c.nearest,
                    c.nearest_label.clone().unwrap_or_default(),
                    c.nearest_ani
                );
            }
        }
        fs::write(o.out.join("label_check.tsv"), s)?;
    }
    fs::write(
        o.out.join("popsketch.json"),
        serde_json::to_string_pretty(&summary_json(o, r))?,
    )?;
    Ok(())
}

pub fn summary_json(o: &Options, r: &Results) -> Value {
    let mut by_status: HashMap<&str, usize> = HashMap::new();
    let mut by_flag: HashMap<&str, Vec<&str>> = HashMap::new();
    let mut by_note: HashMap<&str, Vec<&str>> = HashMap::new();
    for row in &r.rows {
        *by_status.entry(row.status.as_str()).or_default() += 1;
        for fl in row.flags() {
            by_flag.entry(fl.code).or_default().push(&row.id);
        }
        for n in row.qc.iter().flat_map(|q| &q.notes) {
            by_note.entry(n.code).or_default().push(&row.id);
        }
    }
    let gap = r.gap.as_ref().map(|g| {
        json!({"verdict": match &g.verdict { Verdict::Clear(k) => format!("clear:{k}"), Verdict::NoClearCut => "no_clear_cut".into(), Verdict::Single => "single_population".into(), Verdict::TooFew => "too_few".into() },
               "ratio": g.ratio, "null_p": g.null_p, "gap_over_floor": g.gap_over_floor, "top_merge_height": g.top_height,
               "candidates": g.candidates.iter().map(|c| json!({"k": c.k, "gap": c.gap})).collect::<Vec<_>>()})
    });
    let labels: Vec<Value> = r
        .label_reports
        .iter()
        .map(|c| {
            json!({"column": c.name, "labelled": c.n_labelled, "groups": c.groups, "ari_by_k": c.ari, "contingency_k": c.contingency_k,
                   "contingency": c.contingency, "purity": c.purity, "unchecked_labels": c.unchecked,
                   "conflicts": c.conflicts.iter().map(|x| json!({"id": x.id, "label": x.label, "own_median_ani": x.own_median_ani, "better_label": x.other_label, "better_median_ani": x.other_median_ani, "nearest": x.nearest, "nearest_label": x.nearest_label})).collect::<Vec<_>>()})
        })
        .collect();
    json!({
        "tool": "popsketch", "version": env!("CARGO_PKG_VERSION"),
        "options": {"ksize": o.ksize, "scaled": o.scaled, "max_bases": o.max_bases, "r1_only": o.r1_only, "min_cov": o.min_cov,
                    "flagged": format!("{:?}", o.flagged), "ordination": format!("{:?}", o.ordination), "min_floors": o.min_floors,
                    "index": o.index.as_ref().map(|p| p.display().to_string())},
        "samples": {"total": r.rows.len(), "by_status": by_status, "in_analysis": r.included.len()},
        "flags": by_flag,
        "notes": by_note,
        "noise_floor": {"value": r.floor.value, "source": r.floor.source, "n": r.floor.n, "median": r.floor.median, "p90": r.floor.p90},
        "clusters": {"k_used": r.k_used, "gap": gap},
        "ordination": r.ordination.as_ref().map(|x| json!({"method": x.method, "variance_fraction": x.var_frac, "features": x.n_features})),
        "tree": {"written": r.newick.is_some(), "branches_collapsed": r.collapsed},
        "labels": labels,
        "dapc": r.dapc.as_ref().map(|dp| json!({"column": dp.column, "n_pcs": dp.dapc.n_pcs, "axes": dp.dapc.eigenvalues.len(), "eigenvalues": dp.dapc.eigenvalues,
            "groups": dp.dapc.groups, "accuracy": dp.dapc.accuracy,
            "assigned_elsewhere": (0..dp.labels.len()).filter(|&i| dp.dapc.assigned[i].is_some() && dp.labels[i].as_deref() != dp.dapc.assigned[i].map(|g| dp.dapc.groups[g].as_str()))
                .map(|i| json!({"id": r.rows[r.included[i]].id, "label": dp.labels[i], "assigned": dp.dapc.assigned[i].map(|g| dp.dapc.groups[g].clone())})).collect::<Vec<_>>()})),
        "compare": json!({
            "tree": r.tree_comparison.as_ref().map(|c| json!({"tree_tips": c.tree_tips, "shared": c.shared, "spearman_all_pairs": c.spearman_all,
                "spearman_closest_half": c.spearman_closest_half,
                "clusters": c.clusters.iter().map(|x| json!({"cluster": x.cluster, "samples": x.n, "exact_clade": x.exact_clade, "best_f1": x.best_f1})).collect::<Vec<_>>()})),
            "pcs": r.pcs_comparison.as_ref().map(|c| json!({"shared": c.shared, "cosines": c.cosines.iter().map(|(m, v)| json!({"top": m, "cosines": v})).collect::<Vec<_>>(),
                "procrustes_top3": c.procrustes_top3, "distance_spearman_top5": c.distance_spearman_top5})),
        }),
        "advice": r.advice.as_ref().map(|a| json!({"fraction_pairs_below_floor": a.frac_pairs_below_floor, "statements": a.statements})),
        "warnings": r.warnings,
    })
}
