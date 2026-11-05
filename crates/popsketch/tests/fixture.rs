//! End-to-end test on simulated reads with planted truth: three groups of clean samples (a known population) plus planted
//! problems (low coverage, a foreign genome, half the content, a second genome mixed in) and one swapped label. The
//! generator is deterministic (fixed seeds) and needs no external data or tool.

use popsketch::cluster::Verdict;
use popsketch::pipeline::{run, FlagMode, Options, ReportFormat, Results, Status};
use std::collections::HashSet;
use std::fs;
use std::path::{Path, PathBuf};

use popsketch::demo::{generate, mutate, random_genome, reads, Rng};

const GLEN: usize = 120_000;
const COV: f64 = 12.0;

struct Truth {
    dir: PathBuf,
    labels: PathBuf,
    problems: HashSet<String>,
    swapped: String,
    group_of: Vec<(String, usize)>,
}

fn build(name: &str) -> Truth {
    let dir = std::env::temp_dir().join(format!("popsketch_fx_{name}_{}", std::process::id()));
    let _ = fs::remove_dir_all(&dir);
    let d = generate(&dir, GLEN, COV).unwrap();
    Truth {
        dir,
        labels: d.labels,
        problems: d.problems.into_iter().collect(),
        swapped: d.swapped,
        group_of: d.group_of,
    }
}

fn opts(t: &Truth, name: &str) -> Options {
    Options {
        inputs: vec![t.dir.join("reads")],
        labels: Some(t.labels.clone()),
        out: t.dir.join(format!("out_{name}")),
        scaled: 5,
        max_bases: 0,
        quiet: true,
        floor_samples: 20,
        threads: Some(4),
        ..Options::default()
    }
}

fn status<'a>(r: &'a Results, id: &str) -> &'a Status {
    &r.rows
        .iter()
        .find(|x| x.id == id)
        .unwrap_or_else(|| panic!("no sample {id}"))
        .status
}

fn flagged(r: &Results) -> HashSet<String> {
    r.rows
        .iter()
        .filter(|x| x.qc.as_ref().is_some_and(|q| q.flagged()))
        .map(|x| x.id.clone())
        .collect()
}

fn cleanup(t: &Truth) {
    if std::env::var_os("POPSKETCH_KEEP").is_some() {
        eprintln!("kept: {}", t.dir.display());
        return;
    }
    let _ = fs::remove_dir_all(&t.dir);
}

fn exists(p: &Path) -> bool {
    p.exists() && fs::metadata(p).map(|m| m.len() > 0).unwrap_or(false)
}

#[test]
fn audit_finds_exactly_the_planted_problems_and_the_structure_is_recovered() {
    let t = build("main");
    let r = run(&opts(&t, "main")).unwrap();
    // 1. every planted problem and nothing else is flagged
    let f = flagged(&r);
    assert_eq!(f, t.problems, "flagged {:?}", {
        let mut v: Vec<_> = r
            .rows
            .iter()
            .filter(|x| x.qc.as_ref().is_some_and(|q| q.flagged()))
            .map(|x| {
                (
                    x.id.clone(),
                    x.flags().iter().map(|f| f.code).collect::<Vec<_>>(),
                )
            })
            .collect();
        v.sort();
        v
    });
    for id in &t.problems {
        assert_eq!(*status(&r, id), Status::Excluded, "{id}");
    }
    let code = |id: &str, c: &str| {
        r.rows
            .iter()
            .find(|x| x.id == id)
            .unwrap()
            .flags()
            .iter()
            .any(|f| f.code == c)
    };
    assert!(code("lowcov0", "low_coverage") && code("lowcov1", "low_coverage"));
    assert!(code("half", "size_low"));
    assert!(code("mixed", "size_high"));
    assert!(code("foreign", "unrelated_to_cohort"));
    // 2. the 27 clean samples are in the analysis
    assert_eq!(r.included.len(), 27);
    // 3. noise floor measured from split halves
    assert_eq!(r.floor.source, "split_half", "{:?}", r.floor);
    assert!(
        r.floor.n >= 8 && r.floor.value > 0.0 && r.floor.value < 0.01,
        "{:?}",
        r.floor
    );
    // 4. three clear clusters that match the planted groups
    let gap = r.gap.as_ref().unwrap();
    assert_eq!(gap.verdict, Verdict::Clear(3), "{gap:?}");
    let cl = r.clusters.as_ref().unwrap();
    for (i, &row) in r.included.iter().enumerate() {
        let id = &r.rows[row].id;
        let truth = t.group_of.iter().find(|x| &x.0 == id).unwrap().1;
        let mate = r
            .included
            .iter()
            .position(|&rr| t.group_of.iter().find(|x| x.0 == r.rows[rr].id).unwrap().1 == truth)
            .unwrap();
        assert_eq!(cl[i], cl[mate], "{id} is not with its group");
    }
    // 5. the PCA separates the groups (the first two PCs: between-group centroid distance >> within-group spread)
    let ord = r.ordination.as_ref().unwrap();
    let pt = |i: usize| (ord.coords[i][0], ord.coords[i][1]);
    let centroid = |g: usize| {
        let idx: Vec<usize> = (0..r.included.len())
            .filter(|&i| {
                t.group_of
                    .iter()
                    .find(|x| x.0 == r.rows[r.included[i]].id)
                    .unwrap()
                    .1
                    == g
            })
            .collect();
        let n = idx.len() as f64;
        (
            idx.iter().map(|&i| pt(i).0).sum::<f64>() / n,
            idx.iter().map(|&i| pt(i).1).sum::<f64>() / n,
            idx,
        )
    };
    let cs: Vec<_> = (0..3).map(centroid).collect();
    let spread = cs
        .iter()
        .flat_map(|c| {
            c.2.iter()
                .map(|&i| ((pt(i).0 - c.0).powi(2) + (pt(i).1 - c.1).powi(2)).sqrt())
        })
        .fold(0.0, f64::max);
    for a in 0..3 {
        for b in (a + 1)..3 {
            let d = ((cs[a].0 - cs[b].0).powi(2) + (cs[a].1 - cs[b].1).powi(2)).sqrt();
            assert!(
                d > 3.0 * spread,
                "groups {a} and {b}: {d} vs spread {spread}"
            );
        }
    }
    assert!(ord.var_frac[0] + ord.var_frac[1] > 0.5);
    // 6. the swapped label is the only conflict, and the nearest sample is in the group the data say
    let rep = &r.label_reports[0];
    assert_eq!(rep.conflicts.len(), 1, "{:?}", rep.conflicts);
    assert_eq!(rep.conflicts[0].id, t.swapped);
    assert_eq!(rep.conflicts[0].label, "group1");
    assert_eq!(rep.conflicts[0].other_label, "group0");
    // DAPC on the primary label column: two axes for three groups, and only the swapped sample lands in another group
    let dp = r.dapc.as_ref().expect("dapc");
    assert_eq!(dp.column, "group");
    assert_eq!(dp.dapc.eigenvalues.len(), 2);
    let pos = r
        .included
        .iter()
        .position(|&x| r.rows[x].id == t.swapped)
        .unwrap();
    assert_eq!(
        dp.dapc.assigned[pos].map(|g| dp.dapc.groups[g].as_str()),
        Some("group0"),
        "{:?}",
        dp.dapc.groups
    );
    let elsewhere: Vec<usize> = (0..r.included.len())
        .filter(|&i| {
            dp.dapc.assigned[i].map(|g| dp.dapc.groups[g].as_str()) != dp.labels[i].as_deref()
        })
        .collect();
    assert_eq!(elsewhere, vec![pos]);
    assert!(
        (dp.dapc.accuracy - 26.0 / 27.0).abs() < 1e-9,
        "{}",
        dp.dapc.accuracy
    );
    // the second label column has blanks: its groups are judged only from labelled samples
    assert_eq!(r.label_reports.len(), 2);
    assert!(r.label_reports[1].n_labelled < r.included.len());
    // 7. outputs
    let out = t.dir.join("out_main");
    for f in [
        "dapc.tsv",
        "qc.tsv",
        "ani.tsv",
        "mash.tsv",
        "pca.tsv",
        "clusters.tsv",
        "tree.nwk",
        "label_check.tsv",
        "popsketch.json",
        "report.html",
    ] {
        assert!(exists(&out.join(f)), "{f} missing");
    }
    let nwk = fs::read_to_string(out.join("tree.nwk")).unwrap();
    assert_eq!(nwk.matches('(').count(), nwk.matches(')').count());
    for id in &t.problems {
        assert!(!nwk.contains(id.as_str()), "{id} is in the tree");
    }
    assert!(nwk.contains("g1s4"));
    let html = fs::read_to_string(out.join("report.html")).unwrap();
    assert!(
        !html.contains("The set spans distant genomes"),
        "the note appears for a single-species set"
    );
    assert!(
        html.contains("1. Audit")
            && html.contains("DAPC on the labels of group")
            && html.contains("What the audit cannot see")
            && html.contains("group1")
            && html.contains("g0s0")
    );
    // 8. resume reuses the sketches and gives the same answer
    let mut o2 = opts(&t, "main");
    o2.resume = true;
    let r2 = run(&o2).unwrap();
    assert_eq!(flagged(&r2), t.problems);
    assert_eq!(r2.clusters, r.clusters);
    cleanup(&t);
}

#[test]
fn mark_keeps_flagged_samples_in_and_decontaminate_puts_the_mixed_one_back() {
    let t = build("modes");
    let mut o = opts(&t, "mark");
    o.flagged = FlagMode::Mark;
    let r = run(&o).unwrap();
    assert_eq!(*status(&r, "mixed"), Status::Marked);
    assert_eq!(*status(&r, "half"), Status::Marked);
    assert!(r.included.len() > 27);
    // the unrelated genome stays in, so the set spans distant genomes and the report says so
    assert_eq!(*status(&r, "foreign"), Status::Marked);
    let html = fs::read_to_string(t.dir.join("out_mark/report.html")).unwrap();
    assert!(
        html.contains("The set spans distant genomes"),
        "no distant-genomes note"
    );
    let mut o = opts(&t, "decon");
    o.flagged = FlagMode::Decontaminate;
    let r = run(&o).unwrap();
    assert_eq!(
        *status(&r, "mixed"),
        Status::Decontaminated,
        "{:?}",
        r.rows.iter().find(|x| x.id == "mixed").map(|x| x
            .flags()
            .iter()
            .map(|f| f.code)
            .collect::<Vec<_>>())
    );
    assert_eq!(*status(&r, "half"), Status::Excluded); // not a contamination flag
    assert_eq!(*status(&r, "lowcov0"), Status::Excluded);
    assert_eq!(*status(&r, "foreign"), Status::Excluded);
    // the decontaminated sample joins its own group (group 1)
    let cl = r.clusters.as_ref().unwrap();
    let pos = |id: &str| r.included.iter().position(|&x| r.rows[x].id == id).unwrap();
    assert_eq!(cl[pos("mixed")], cl[pos("g1s3")]);
    cleanup(&t);
}

#[test]
fn a_clonal_set_is_one_population_and_too_many_samples_are_refused() {
    let dir = std::env::temp_dir().join(format!("popsketch_fx_clone_{}", std::process::id()));
    let _ = fs::remove_dir_all(&dir);
    fs::create_dir_all(dir.join("reads")).unwrap();
    let base = random_genome(GLEN, 5);
    for s in 0..12 {
        let mut t = String::new();
        reads(
            &mutate(&base, 0.00003, 600 + s),
            0,
            GLEN,
            COV,
            40_000 + s,
            &mut t,
            &format!("c{s}_"),
        );
        fs::write(dir.join("reads").join(format!("clone{s:02}.fastq")), t).unwrap();
    }
    let o = Options {
        inputs: vec![dir.join("reads")],
        out: dir.join("out"),
        scaled: 5,
        max_bases: 0,
        quiet: true,
        threads: Some(4),
        ..Options::default()
    };
    let r = run(&o).unwrap();
    assert_eq!(r.included.len(), 12);
    assert_eq!(
        r.gap.as_ref().unwrap().verdict,
        Verdict::Single,
        "{:?} floor {:?}",
        r.gap,
        r.floor
    );
    assert_eq!(r.k_used, Some(1));
    let html = fs::read_to_string(dir.join("out/report.html")).unwrap();
    assert!(html.contains("one population"));
    let too_many = Options {
        max_samples: 5,
        out: dir.join("out2"),
        ..o.clone()
    };
    let err = run(&too_many).err().unwrap().to_string();
    assert!(err.contains("--max-samples"), "{err}");
    let _ = fs::remove_dir_all(&dir);
}

#[test]
fn comparison_with_an_external_tree_and_pcs_reports_clades_and_matching_planes() {
    let t = build("compare");
    // an external "SNP" tree with the three planted groups as clades, and PC scores that place the groups apart
    let clean: Vec<&(String, usize)> = t.group_of.iter().collect();
    let mut nwk = String::from("(");
    for g in 0..3 {
        if g > 0 {
            nwk.push(',');
        }
        nwk.push('(');
        let members: Vec<String> = clean
            .iter()
            .filter(|x| x.1 == g)
            .map(|x| format!("{}:0.001", x.0))
            .collect();
        nwk.push_str(&members.join(","));
        nwk.push_str("):0.01");
    }
    nwk.push_str(");");
    fs::write(t.dir.join("snp.nwk"), nwk).unwrap();
    let centre = [(-10.0, 0.0), (5.0, 8.0), (5.0, -8.0)];
    let mut pcs = String::from("id\tPC1\tPC2\tPC3\n");
    let mut r = Rng(4242);
    for (id, g) in &t.group_of {
        pcs.push_str(&format!(
            "{id}\t{:.4}\t{:.4}\t{:.4}\n",
            centre[*g].0 + r.f() - 0.5,
            centre[*g].1 + r.f() - 0.5,
            r.f() - 0.5
        ));
    }
    fs::write(t.dir.join("snp_pcs.tsv"), pcs).unwrap();
    let mut o = opts(&t, "compare");
    o.compare_tree = Some(t.dir.join("snp.nwk"));
    o.compare_pcs = Some(t.dir.join("snp_pcs.tsv"));
    o.report_format = ReportFormat::All;
    let r = run(&o).unwrap();
    assert!(r.warnings.is_empty(), "{:?}", r.warnings);
    let tc = r.tree_comparison.as_ref().expect("tree comparison");
    assert_eq!(tc.shared, 27);
    assert_eq!(tc.clusters.len(), 3);
    assert!(
        tc.clusters.iter().all(|c| c.exact_clade),
        "{:?}",
        tc.clusters
    );
    assert!(tc.spearman_all > 0.7, "{}", tc.spearman_all);
    let pc = r.pcs_comparison.as_ref().expect("pcs comparison");
    assert_eq!(pc.shared, 27);
    let top2 = &pc.cosines.iter().find(|c| c.0 == 2).unwrap().1;
    assert!(top2.iter().all(|&c| c > 0.9), "{top2:?}");
    assert!(pc.procrustes_top3.is_some());
    let json: serde_json::Value = serde_json::from_str(
        &fs::read_to_string(t.dir.join("out_compare/popsketch.json")).unwrap(),
    )
    .unwrap();
    assert_eq!(
        json["compare"]["tree"]["clusters"]
            .as_array()
            .unwrap()
            .len(),
        3
    );
    let html = fs::read_to_string(t.dir.join("out_compare/report.html")).unwrap();
    assert!(
        html.contains("Comparison with the supplied SNP results")
            && html.contains("is a clade of the tree")
    );
    // the Markdown report: the same sections, figures as standalone SVG files
    let md = fs::read_to_string(t.dir.join("out_compare/report.md")).unwrap();
    for needle in [
        "# popsketch: sample audit and coarse structure",
        "## 1. Audit",
        "## 2. Coarse structure (secondary)",
        "## 4. Comparison with the supplied SNP results",
        "| flag | samples | which |",
        "![size_vs_coverage](figures/size_vs_coverage.svg)",
        "![pca](figures/pca_1.svg)",
        "![tree](figures/tree.svg)",
        "What the audit cannot see.",
    ] {
        assert!(md.contains(needle), "report.md lacks {needle:?}");
    }
    for f in [
        "size_vs_coverage.svg",
        "pca_1.svg",
        "pca_2.svg",
        "distances.svg",
        "tree.svg",
    ] {
        let svg = fs::read_to_string(t.dir.join("out_compare/figures").join(f))
            .unwrap_or_else(|_| panic!("{f} missing"));
        assert!(
            svg.starts_with("<svg xmlns=\"http://www.w3.org/2000/svg\" width=\""),
            "{f}: {}",
            &svg[..60.min(svg.len())]
        );
        assert!(svg.contains("<style>") && svg.ends_with("</svg>"), "{f}");
    }
    // the PDF report (written by `All`): a real PDF with pages
    let pdf = fs::read(t.dir.join("out_compare/report.pdf")).unwrap();
    assert!(
        pdf.starts_with(b"%PDF-") && pdf.len() > 5_000,
        "report.pdf: {} bytes",
        pdf.len()
    );
    // the page objects are inside compressed streams: check the trailer instead (the pages were checked with pdfinfo and pdftotext)
    assert!(
        pdf.windows(5).rev().take(64).any(|w| w == b"%%EOF"),
        "no %%EOF trailer"
    );
    assert!(String::from_utf8_lossy(&pdf).contains("startxref"));
    // the default writes the HTML only
    assert!(!t.dir.join("out_main").join("report.md").exists());
    // a missing file is a warning, not a failure
    let mut o2 = opts(&t, "compare2");
    o2.compare_tree = Some(t.dir.join("nope.nwk"));
    let r2 = run(&o2).unwrap();
    assert!(r2.warnings.iter().any(|w| w.contains("--compare-tree")));
    cleanup(&t);
}

#[test]
fn labels_that_do_not_follow_the_structure_give_no_conflicts_and_no_dapc() {
    let t = build("noise");
    // alternate two labels across the samples: neither is a coherent cluster of the sketches
    let mut labels = String::from("id\tnoise\n");
    for (i, (id, _)) in t.group_of.iter().enumerate() {
        labels.push_str(&format!("{id}\t{}\n", if i % 2 == 0 { "X" } else { "Y" }));
    }
    let lp = t.dir.join("noise_labels.tsv");
    fs::write(&lp, labels).unwrap();
    let mut o = opts(&t, "noise");
    o.labels = Some(lp);
    let r = run(&o).unwrap();
    let rep = &r.label_reports[0];
    assert!(rep.conflicts.is_empty(), "{:?}", rep.conflicts);
    assert_eq!(rep.unchecked, vec!["X".to_string(), "Y".to_string()]);
    assert!(r.dapc.is_none());
    assert!(
        r.warnings
            .iter()
            .any(|w| w.contains("DAPC") && w.contains("coherent")),
        "{:?}",
        r.warnings
    );
    cleanup(&t);
}

/// A recombining population: every genome is a random mosaic of blocks of three founders (about 1.5% from the root each), with a
/// little private variation. There are no groups, so the tool must not call a clear cut.
#[test]
fn a_recombining_population_is_not_split_into_clusters() {
    let dir = std::env::temp_dir().join(format!("popsketch_fx_recomb_{}", std::process::id()));
    let _ = fs::remove_dir_all(&dir);
    fs::create_dir_all(dir.join("reads")).unwrap();
    let root = random_genome(GLEN, 3);
    let founders: Vec<Vec<u8>> = (0..3).map(|f| mutate(&root, 0.015, 900 + f)).collect();
    let block = 4_000;
    let mut r = Rng(31337);
    for s in 0..24u64 {
        let mut g: Vec<u8> = Vec::with_capacity(GLEN);
        let mut start = 0;
        while start < GLEN {
            let end = (start + block).min(GLEN);
            g.extend_from_slice(&founders[r.below(3)][start..end]);
            start = end;
        }
        let genome = mutate(&g, 0.0015, 4000 + s);
        let mut t = String::new();
        reads(&genome, 0, GLEN, COV, 50_000 + s, &mut t, &format!("m{s}_"));
        fs::write(dir.join("reads").join(format!("mosaic{s:02}.fastq")), t).unwrap();
    }
    let o = Options {
        inputs: vec![dir.join("reads")],
        out: dir.join("out"),
        scaled: 5,
        max_bases: 0,
        quiet: true,
        threads: Some(4),
        ..Options::default()
    };
    let res = run(&o).unwrap();
    assert_eq!(
        res.included.len(),
        24,
        "nothing should be flagged in a clean recombining set"
    );
    let gap = res.gap.as_ref().unwrap();
    assert!(
        !matches!(gap.verdict, Verdict::Clear(_)),
        "a recombining population has no groups, yet: {:?} (jump {:.1} floors, floor {:.5})",
        gap.verdict,
        gap.gap_over_floor,
        res.floor.value
    );
    eprintln!(
        "recombining set: verdict {:?}, jump {:.2} floors, ratio {:.2}, floor {:.5}, top merge {:.5}",
        gap.verdict, gap.gap_over_floor, gap.ratio, res.floor.value, gap.top_height
    );
    let html = fs::read_to_string(dir.join("out/report.html")).unwrap();
    assert!(html.contains("No single clear cut") || html.contains("one population"));
    let _ = fs::remove_dir_all(&dir);
}

#[test]
fn threshold_options_and_no_pca_take_effect() {
    let t = build("thr");
    // a size threshold nothing can reach: the size flags (half, mixed) go away, the other planted problems stay
    let mut o = opts(&t, "thr");
    o.size_z = 1e6;
    o.no_pca = true;
    let r = run(&o).unwrap();
    let f = flagged(&r);
    assert!(!f.contains("half") && !f.contains("mixed"), "{f:?}");
    assert!(f.contains("lowcov0") && f.contains("foreign"), "{f:?}");
    // --no-pca: distances and clusters, no ordination, no DAPC, and the report says so
    assert!(r.dist.is_some() && r.clusters.is_some());
    assert!(r.ordination.is_none() && r.dapc.is_none());
    let html = fs::read_to_string(t.dir.join("out_thr/report.html")).unwrap();
    assert!(html.contains("--no-pca"));
    assert!(!t.dir.join("out_thr/pca.tsv").exists());
    // a lower floor multiple makes the clear cut easier to call; an impossible one makes it impossible
    let mut hard = opts(&t, "thr_hard");
    hard.min_floors = 1e6;
    hard.no_pca = true;
    let r2 = run(&hard).unwrap();
    assert!(
        matches!(r2.gap.as_ref().unwrap().verdict, Verdict::NoClearCut),
        "{:?}",
        r2.gap.as_ref().unwrap().verdict
    );
    assert!(r2.clusters.is_none());
    cleanup(&t);
}
