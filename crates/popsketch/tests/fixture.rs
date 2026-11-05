//! End-to-end test on simulated reads with planted truth: three groups of clean samples (a known population) plus planted
//! problems (low coverage, a foreign genome, half the content, a second genome mixed in) and one swapped label. The
//! generator is deterministic (fixed seeds) and needs no external data or tool.

use popsketch::cluster::Verdict;
use popsketch::pipeline::{run, FlagMode, Options, Results, Status};
use std::collections::HashSet;
use std::fs;
use std::path::{Path, PathBuf};

struct Rng(u64);
impl Rng {
    fn next(&mut self) -> u64 {
        self.0 ^= self.0 << 13;
        self.0 ^= self.0 >> 7;
        self.0 ^= self.0 << 17;
        self.0
    }
    fn below(&mut self, n: usize) -> usize {
        (self.next() % n as u64) as usize
    }
    fn f(&mut self) -> f64 {
        (self.next() >> 11) as f64 / (1u64 << 53) as f64
    }
}

fn random_genome(len: usize, seed: u64) -> Vec<u8> {
    let mut r = Rng(seed);
    (0..len).map(|_| b"ACGT"[r.below(4)]).collect()
}

/// Substitute a fraction `rate` of the positions with a different base.
fn mutate(g: &[u8], rate: f64, seed: u64) -> Vec<u8> {
    let mut r = Rng(seed);
    g.iter()
        .map(|&b| {
            if r.f() < rate {
                loop {
                    let c = b"ACGT"[r.below(4)];
                    if c != b {
                        break c;
                    }
                }
            } else {
                b
            }
        })
        .collect()
}

fn revcomp(s: &[u8]) -> Vec<u8> {
    s.iter()
        .rev()
        .map(|&b| match b {
            b'A' => b'T',
            b'C' => b'G',
            b'G' => b'C',
            _ => b'A',
        })
        .collect()
}

/// FASTQ text of reads at `cov` x from `g[from..to]`, both strands, with substitution errors.
fn reads(g: &[u8], from: usize, to: usize, cov: f64, seed: u64, out: &mut String, tag: &str) {
    let (len, err) = (150usize, 0.003);
    let mut r = Rng(seed);
    let n = (cov * (to - from) as f64 / len as f64) as usize;
    for i in 0..n {
        let start = from + r.below(to - from - len);
        let mut s = g[start..start + len].to_vec();
        for b in s.iter_mut() {
            if r.f() < err {
                *b = b"ACGT"[r.below(4)];
            }
        }
        if r.below(2) == 1 {
            s = revcomp(&s);
        }
        out.push_str(&format!(
            "@{tag}{i}\n{}\n+\n{}\n",
            String::from_utf8_lossy(&s),
            "I".repeat(len)
        ));
    }
}

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
    fs::create_dir_all(dir.join("reads")).unwrap();
    let root = random_genome(GLEN, 1);
    let foreign = random_genome(GLEN, 99);
    let foreign2 = random_genome(GLEN, 77);
    let groups: Vec<Vec<u8>> = (0..3).map(|g| mutate(&root, 0.015, 1000 + g)).collect();
    let mut labels = String::from("id\tgroup\tsource\n");
    let mut problems = HashSet::new();
    let mut group_of = Vec::new();
    let write = |id: &str, text: &str| {
        fs::write(dir.join("reads").join(format!("{id}.fastq")), text).unwrap()
    };
    for (g, gg) in groups.iter().enumerate() {
        for s in 0..9 {
            let id = format!("g{g}s{s}");
            let genome = mutate(gg, 0.0015, 5000 + (g * 100 + s) as u64);
            let mut t = String::new();
            reads(
                &genome,
                0,
                GLEN,
                COV,
                20_000 + (g * 100 + s) as u64,
                &mut t,
                &id,
            );
            write(&id, &t);
            // sample g0s0 is labelled as group 1 (a swap)
            let lab = if id == "g0s0" { 1 } else { g };
            labels.push_str(&format!(
                "{id}\tgroup{lab}\t{}\n",
                if s % 3 == 0 { "NA" } else { "field" }
            ));
            group_of.push((id, g));
        }
    }
    // planted problems
    for (i, grp) in groups.iter().take(2).enumerate() {
        let id = format!("lowcov{i}");
        let mut t = String::new();
        reads(
            &mutate(grp, 0.0015, 7000 + i as u64),
            0,
            GLEN,
            0.6,
            30_000 + i as u64,
            &mut t,
            &id,
        );
        write(&id, &t);
        problems.insert(id.clone());
        labels.push_str(&format!("{id}\tgroup{i}\tNA\n"));
    }
    let mut t = String::new();
    reads(&foreign, 0, GLEN, COV, 31_000, &mut t, "foreign");
    write("foreign", &t);
    problems.insert("foreign".to_string());
    labels.push_str("foreign\t\tNA\n");
    let mut t = String::new();
    let g = mutate(&groups[2], 0.0015, 7100);
    reads(&g, 0, GLEN / 2, COV, 32_000, &mut t, "half"); // only half the genome
    write("half", &t);
    problems.insert("half".to_string());
    labels.push_str("half\tgroup2\tNA\n");
    let mut t = String::new();
    reads(
        &mutate(&groups[1], 0.0015, 7200),
        0,
        GLEN,
        COV,
        33_000,
        &mut t,
        "mixh",
    );
    reads(&foreign2, 0, GLEN, COV * 0.8, 34_000, &mut t, "mixc"); // a second genome at 0.8x the host coverage
    write("mixed", &t);
    problems.insert("mixed".to_string());
    labels.push_str("mixed\tgroup1\tNA\n");
    let lp = dir.join("labels.tsv");
    fs::write(&lp, labels).unwrap();
    Truth {
        dir,
        labels: lp,
        problems,
        swapped: "g0s0".into(),
        group_of,
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
    // the second label column has blanks: its groups are judged only from labelled samples
    assert_eq!(r.label_reports.len(), 2);
    assert!(r.label_reports[1].n_labelled < r.included.len());
    // 7. outputs
    let out = t.dir.join("out_main");
    for f in [
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
        html.contains("1. Audit")
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
