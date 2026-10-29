//! End-to-end test of `gmg sketch/build/run` on tiny synthetic genomes (no data files needed).

use std::fs;
use std::path::Path;
use std::process::{Command, Output};

fn gmg(args: &[&str]) -> Output {
    Command::new(env!("CARGO_BIN_EXE_gmg"))
        .args(args)
        .output()
        .expect("run gmg")
}

fn ok(args: &[&str]) -> Output {
    let o = gmg(args);
    assert!(
        o.status.success(),
        "gmg {args:?} failed: {}",
        String::from_utf8_lossy(&o.stderr)
    );
    o
}

struct Rng(u64);
impl Rng {
    fn next(&mut self) -> u64 {
        self.0 ^= self.0 << 13;
        self.0 ^= self.0 >> 7;
        self.0 ^= self.0 << 17;
        self.0
    }
    fn genome(&mut self, n: usize) -> Vec<u8> {
        (0..n)
            .map(|_| b"ACGT"[(self.next() % 4) as usize])
            .collect()
    }
    fn mutate(&mut self, g: &[u8], per_thousand: u64) -> Vec<u8> {
        g.iter()
            .map(|&b| {
                if self.next() % 1000 < per_thousand {
                    b"ACGT"[(self.next() % 4) as usize]
                } else {
                    b
                }
            })
            .collect()
    }
}

fn fasta(name: &str, seq: &[u8]) -> String {
    let mut s = format!(">{name}\n");
    for c in seq.chunks(70) {
        s.push_str(std::str::from_utf8(c).unwrap());
        s.push('\n');
    }
    s
}

fn read_tsv(p: &Path) -> (Vec<String>, Vec<Vec<String>>) {
    let t = fs::read_to_string(p).unwrap();
    let mut lines = t.lines();
    let header = lines
        .next()
        .unwrap()
        .split('\t')
        .map(String::from)
        .collect();
    (
        header,
        lines
            .map(|l| l.split('\t').map(String::from).collect())
            .collect(),
    )
}

fn col(h: &[String], name: &str) -> usize {
    h.iter()
        .position(|c| c == name)
        .unwrap_or_else(|| panic!("no column {name}"))
}

/// Six random 200 kb genomes in two genera of three phyla, sketched at scaled=10 and indexed; returns (genomes, index path, rng).
fn build_reference(d: &Path) -> (Vec<Vec<u8>>, std::path::PathBuf, Rng) {
    let mut rng = Rng(0x9E3779B97F4A7C15);
    let genomes: Vec<Vec<u8>> = (0..6).map(|_| rng.genome(200_000)).collect();
    let lineage = [
        "Ascomycota\tEurotiomycetes\tEurotiales\tAspergillaceae\tAspergillus\tAspergillus alpha",
        "Ascomycota\tEurotiomycetes\tEurotiales\tAspergillaceae\tAspergillus\tAspergillus beta",
        "Ascomycota\tSordariomycetes\tHypocreales\tNectriaceae\tFusarium\tFusarium gamma",
        "Ascomycota\tSordariomycetes\tHypocreales\tNectriaceae\tFusarium\tFusarium delta",
        "Basidiomycota\tAgaricomycetes\tAgaricales\tAgaricaceae\tAgaricus\tAgaricus eps",
        "Basidiomycota\tAgaricomycetes\tAgaricales\tAgaricaceae\tAgaricus\tAgaricus zeta",
    ];
    let mut meta = String::from("id\ttaxid\tphylum\tclass\torder\tfamily\tgenus\tspecies\n");
    let refs = d.join("refs");
    fs::create_dir(&refs).unwrap();
    for (i, g) in genomes.iter().enumerate() {
        fs::write(refs.join(format!("g{i}.fa")), fasta(&format!("g{i}"), g)).unwrap();
        meta.push_str(&format!("g{i}\t{}\t{}\n", 100 + i, lineage[i]));
    }
    fs::write(d.join("meta.tsv"), meta).unwrap();
    let sk = d.join("sk");
    let mut args = vec!["sketch", "-k", "21", "-s", "10", "-o", sk.to_str().unwrap()];
    let ref_files: Vec<String> = (0..6)
        .map(|i| refs.join(format!("g{i}.fa")).display().to_string())
        .collect();
    args.extend(ref_files.iter().map(|s| s.as_str()));
    ok(&args);
    let ix = d.join("ix.gmgi");
    ok(&[
        "build",
        "--meta",
        d.join("meta.tsv").to_str().unwrap(),
        "--sketches",
        sk.to_str().unwrap(),
        "-o",
        ix.to_str().unwrap(),
    ]);
    (genomes, ix, rng)
}

#[test]
fn run_folder_end_to_end() {
    let dir = tempfile::tempdir().unwrap();
    let d = dir.path();
    let (genomes, ix, mut rng) = build_reference(d);

    // inputs: an exact copy of g2, a 2%-mutated copy of g4 in a sub-folder, an empty FASTA, a text file to ignore
    let inp = d.join("in");
    fs::create_dir_all(inp.join("sub")).unwrap();
    fs::write(inp.join("exact.fa"), fasta("exact", &genomes[2])).unwrap();
    fs::write(
        inp.join("sub/mutant.fasta"),
        fasta("mutant", &rng.mutate(&genomes[4], 20)),
    )
    .unwrap();
    fs::write(inp.join("empty.fa"), "").unwrap();
    fs::write(inp.join("notes.txt"), "not a genome").unwrap();

    let out = d.join("out");
    let o = ok(&[
        "run",
        inp.to_str().unwrap(),
        "--index",
        ix.to_str().unwrap(),
        "-o",
        out.to_str().unwrap(),
    ]);
    assert!(
        String::from_utf8_lossy(&o.stderr).contains("3 inputs"),
        "notes.txt must be ignored"
    );

    let (h, rows) = read_tsv(&out.join("results.tsv"));
    assert_eq!(rows.len(), 3);
    let row = |id: &str| {
        rows.iter()
            .find(|r| r[col(&h, "query_id")] == id)
            .unwrap_or_else(|| panic!("no row {id}"))
    };
    let exact = row("exact");
    assert_eq!(exact[col(&h, "status")], "ok");
    assert_eq!(exact[col(&h, "best_rank")], "species");
    assert_eq!(exact[col(&h, "call_species")], "Fusarium gamma");
    assert_eq!(exact[col(&h, "top_hit")], "g2");
    assert_eq!(exact[col(&h, "top_hit_taxid")], "102");
    // a 2% mutant has ANI ~0.98: below the calibrated species threshold, so the call stops at genus
    let mutant = row("mutant");
    assert_eq!(mutant[col(&h, "best_rank")], "genus");
    assert_eq!(mutant[col(&h, "call_genus")], "Agaricus");
    assert_eq!(mutant[col(&h, "call_species")], "");
    // the empty file is an error row, not a failed run
    let empty = row("empty");
    assert_eq!(empty[col(&h, "status")], "error");
    assert!(!empty[col(&h, "message")].is_empty());
    assert!(h.contains(&"schema_version".to_string()));

    let (hh, hits) = read_tsv(&out.join("hits.tsv"));
    assert_eq!(hh[0], "query_id");
    assert!(
        hits.iter().filter(|r| r[0] == "exact").count() >= 1
            && hits.iter().all(|r| r[0] != "empty")
    );
    let manifest: serde_json::Value =
        serde_json::from_str(&fs::read_to_string(out.join("run.json")).unwrap()).unwrap();
    assert_eq!(manifest["status_counts"]["ok"], 2);
    assert_eq!(manifest["status_counts"]["error"], 1);
    assert_eq!(manifest["index"]["n_genomes"], 6);

    // report.html is written by default, with the heatmap (2 assemblies), and `gmg report` regenerates it
    let report = fs::read_to_string(out.join("report.html")).unwrap();
    assert!(
        report.contains("Classification report") && report.contains("Relationships among inputs")
    );
    assert!(report.contains("id=\"q-empty\""), "error card present");
    assert!(
        !report.contains(d.to_str().unwrap()),
        "no local path in the report"
    );
    assert!(out.join("relationships.json").exists());
    // the Markdown report is produced from the same run on request
    ok(&["report", out.to_str().unwrap(), "--format", "md"]);
    let md = fs::read_to_string(out.join("report.md")).unwrap();
    assert!(md.contains("## All inputs") && md.contains("exact"));
    assert!(out.join("figures/depth.svg").exists());
    assert!(
        !md.contains(d.to_str().unwrap()),
        "no local path in the markdown"
    );
    fs::remove_file(out.join("report.html")).unwrap();
    ok(&["report", out.to_str().unwrap()]);
    assert!(out.join("report.html").exists());

    // --resume reuses the two good results and retries only the failed one
    let o = ok(&[
        "run",
        inp.to_str().unwrap(),
        "--index",
        ix.to_str().unwrap(),
        "-o",
        out.to_str().unwrap(),
        "--resume",
    ]);
    let err = String::from_utf8_lossy(&o.stderr);
    assert!(err.contains("1 to process, 2 reused"), "{err}");

    // `novel` thresholds withhold family/genus/species: the same mutant now stops higher up
    let out2 = d.join("out2");
    ok(&[
        "run",
        inp.to_str().unwrap(),
        "--index",
        ix.to_str().unwrap(),
        "-o",
        out2.to_str().unwrap(),
        "--thresholds",
        "novel",
    ]);
    assert!(out2.join("report.html").exists());
    let out3 = d.join("out3");
    ok(&[
        "run",
        inp.to_str().unwrap(),
        "--index",
        ix.to_str().unwrap(),
        "-o",
        out3.to_str().unwrap(),
        "--no-report",
    ]);
    assert!(!out3.join("report.html").exists());
    let (h2, rows2) = read_tsv(&out2.join("results.tsv"));
    let m2 = rows2
        .iter()
        .find(|r| r[col(&h2, "query_id")] == "mutant")
        .unwrap();
    assert_eq!(m2[col(&h2, "call_genus")], "");
    assert_eq!(m2[col(&h2, "status")], "low_confidence");
}

#[test]
fn divergent_input_gets_a_tentative_genus_not_a_genus_call() {
    let dir = tempfile::tempdir().unwrap();
    let d = dir.path();
    let (genomes, ix, mut rng) = build_reference(d);
    // 15% substitutions: ANI to g4 about 0.85, between the tentative cutoff (0.7765) and the genus cutoff (0.90)
    let inp = d.join("in");
    fs::create_dir_all(&inp).unwrap();
    fs::write(
        inp.join("divergent.fa"),
        fasta("divergent", &rng.mutate(&genomes[4], 150)),
    )
    .unwrap();
    let out = d.join("out");
    ok(&[
        "run",
        inp.to_str().unwrap(),
        "--index",
        ix.to_str().unwrap(),
        "-o",
        out.to_str().unwrap(),
    ]);
    let (h, rows) = read_tsv(&out.join("results.tsv"));
    let r = &rows[0];
    let ani: f64 = r[col(&h, "ani")].parse().unwrap();
    assert!(ani > 0.8 && ani < 0.9, "ani {ani}");
    assert_eq!(
        r[col(&h, "call_genus")],
        "",
        "genus must not be called below ANI 0.90"
    );
    assert_eq!(r[col(&h, "tentative_genus")], "Agaricus");
    assert_eq!(r[col(&h, "status")], "low_confidence");
    let report = fs::read_to_string(out.join("report.html")).unwrap();
    assert!(
        report.contains("Agaricus (tentative)") && report.contains("tentative genus: Agaricus")
    );
    // the permissive set (the earlier 0.7765 genus cutoff) names the genus instead
    let out2 = d.join("out2");
    ok(&[
        "run",
        inp.to_str().unwrap(),
        "--index",
        ix.to_str().unwrap(),
        "-o",
        out2.to_str().unwrap(),
        "--thresholds",
        "permissive",
    ]);
    let (h2, rows2) = read_tsv(&out2.join("results.tsv"));
    assert_eq!(rows2[0][col(&h2, "call_genus")], "Agaricus");
    assert_eq!(rows2[0][col(&h2, "tentative_genus")], "");
}

#[test]
fn failed_inputs_get_plain_error_codes() {
    let dir = tempfile::tempdir().unwrap();
    let d = dir.path();
    let (_, ix, _) = build_reference(d);
    let inp = d.join("in");
    fs::create_dir_all(&inp).unwrap();
    fs::write(inp.join("empty.fa"), b"").unwrap();
    fs::write(inp.join("header_only.fa"), b">only_a_header\n").unwrap();
    fs::write(
        inp.join("notseq.fa"),
        b"this is a text file, not a genome\n",
    )
    .unwrap();
    fs::write(
        inp.join("bad.fa.gz"),
        [
            0x1f, 0x8b, 8, 0, 0, 0, 0, 0, 0, 3, 0xde, 0xad, 0xbe, 0xef, 1, 2, 3, 4,
        ],
    )
    .unwrap();
    fs::write(inp.join("tiny.fa"), b">t\nACGTACGT\n").unwrap();
    let out = d.join("out");
    ok(&[
        "run",
        inp.to_str().unwrap(),
        "--index",
        ix.to_str().unwrap(),
        "-o",
        out.to_str().unwrap(),
    ]);
    let (h, rows) = read_tsv(&out.join("results.tsv"));
    let code = |id: &str| {
        let r = rows.iter().find(|r| r[0] == id).unwrap();
        assert_eq!(r[col(&h, "status")], "error", "{id}");
        let m = &r[col(&h, "message")];
        assert!(
            !m.contains('/') && !m.is_empty(),
            "{id}: plain message, no path: {m}"
        );
        r[col(&h, "error_code")].clone()
    };
    assert_eq!(code("empty"), "empty_file");
    assert_eq!(code("header_only"), "no_sequence_records");
    assert_eq!(code("notseq"), "not_sequence_file");
    assert_eq!(code("bad"), "bad_gzip");
    assert_eq!(code("tiny"), "too_short");
    let report = fs::read_to_string(out.join("report.html")).unwrap();
    assert!(
        report.contains("<code>empty_file</code>") && report.contains("Re-download or re-copy it.")
    );
    assert!(!report.contains(d.to_str().unwrap()));

    #[cfg(unix)]
    {
        use std::os::unix::fs::PermissionsExt;
        let f = inp.join("locked.fa");
        fs::write(&f, b">x\nACGT\n").unwrap();
        fs::set_permissions(&f, fs::Permissions::from_mode(0o000)).unwrap();
        // running as root can still read it; only assert when the OS really denies access
        if fs::File::open(&f).is_err() {
            let out2 = d.join("out2");
            ok(&[
                "run",
                f.to_str().unwrap(),
                "--index",
                ix.to_str().unwrap(),
                "-o",
                out2.to_str().unwrap(),
                "--no-report",
            ]);
            let (h2, r2) = read_tsv(&out2.join("results.tsv"));
            assert_eq!(r2[0][col(&h2, "error_code")], "unreadable");
        }
        fs::set_permissions(&f, fs::Permissions::from_mode(0o644)).unwrap();
    }
}

fn fastq(prefix: &str, g: &[u8], start: usize, n: usize, len: usize, step: usize) -> String {
    let mut s = String::new();
    for i in 0..n {
        let a = (start + i * step) % (g.len() - len);
        s.push_str(&format!("@{prefix}{i}\n"));
        s.push_str(std::str::from_utf8(&g[a..a + len]).unwrap());
        s.push_str(&format!("\n+\n{}\n", "I".repeat(len)));
    }
    s
}

#[test]
fn paired_and_lane_split_reads_equal_the_concatenated_file() {
    let dir = tempfile::tempdir().unwrap();
    let d = dir.path();
    let (genomes, ix, _) = build_reference(d);
    let g = &genomes[2];
    // sample A: R1/R2 x two lanes; sample B: the same reads concatenated in one file
    let parts = [
        ("A_S1_L001_R1_001", fastq("a", g, 0, 3000, 100, 61)),
        ("A_S1_L001_R2_001", fastq("b", g, 30, 3000, 100, 61)),
        ("A_S1_L002_R1_001", fastq("c", g, 7, 3000, 100, 59)),
        ("A_S1_L002_R2_001", fastq("d", g, 40, 3000, 100, 59)),
    ];
    let inp = d.join("in");
    fs::create_dir_all(&inp).unwrap();
    let mut all = String::new();
    for (n, t) in &parts {
        fs::write(inp.join(format!("{n}.fastq")), t).unwrap();
        all.push_str(t);
    }
    fs::write(inp.join("whole.fastq"), &all).unwrap();
    // a lone R1 file keeps its own id
    fs::write(
        inp.join("solo_R1_001.fastq"),
        fastq("e", g, 0, 3000, 100, 61),
    )
    .unwrap();
    let out = d.join("out");
    ok(&[
        "run",
        inp.to_str().unwrap(),
        "--index",
        ix.to_str().unwrap(),
        "-o",
        out.to_str().unwrap(),
        "--no-report",
        "--min-abund",
        "1",
    ]);
    let (h, rows) = read_tsv(&out.join("results.tsv"));
    let get = |id: &str| {
        rows.iter()
            .find(|r| r[0] == id)
            .unwrap_or_else(|| {
                panic!(
                    "no row {id}: {:?}",
                    rows.iter().map(|r| &r[0]).collect::<Vec<_>>()
                )
            })
            .clone()
    };
    let (a, w) = (get("A_S1"), get("whole"));
    for c in [
        "status",
        "call_species",
        "n_hashes",
        "shared",
        "total_bases",
        "top_hit",
    ] {
        assert_eq!(a[col(&h, c)], w[col(&h, c)], "{c}");
    }
    assert_eq!(a[col(&h, "call_species")], "Fusarium gamma");
    assert_eq!(get("solo_R1_001")[col(&h, "status")], "ok");
    assert_eq!(rows.len(), 3, "4 lane/pair files -> 1 sample");
    let jl = fs::read_to_string(out.join("results.jsonl")).unwrap();
    assert!(
        jl.contains("\"files\":[\"A_S1_L001_R1_001.fastq\""),
        "grouped sample lists its files"
    );

    // --no-group keeps the five files separate; --list reads paths from a file
    let list = d.join("list.txt");
    fs::write(
        &list,
        format!("# samples\n{}\n", inp.join("whole.fastq").display()),
    )
    .unwrap();
    let out2 = d.join("out2");
    ok(&[
        "run",
        "--list",
        list.to_str().unwrap(),
        "--index",
        ix.to_str().unwrap(),
        "-o",
        out2.to_str().unwrap(),
        "--no-report",
    ]);
    let (_, r2) = read_tsv(&out2.join("results.tsv"));
    assert_eq!(r2.len(), 1);
    let out3 = d.join("out3");
    ok(&[
        "run",
        inp.to_str().unwrap(),
        "--no-group",
        "--index",
        ix.to_str().unwrap(),
        "-o",
        out3.to_str().unwrap(),
        "--no-report",
    ]);
    assert_eq!(read_tsv(&out3.join("results.tsv")).1.len(), 6);
}

#[test]
fn odd_but_legal_and_hostile_inputs_do_not_crash_a_run() {
    let dir = tempfile::tempdir().unwrap();
    let d = dir.path();
    let (genomes, ix, mut rng) = build_reference(d);
    let inp = d.join("in");
    fs::create_dir_all(&inp).unwrap();
    let g = &genomes[2];
    let lf = fasta("lf", g);
    // same genome with CRLF endings, and with endings mixed line by line
    fs::write(inp.join("crlf.fa"), lf.replace('\n', "\r\n")).unwrap();
    let mixed: String = lf
        .split_inclusive('\n')
        .enumerate()
        .map(|(i, l)| {
            if i % 2 == 0 {
                l.replace('\n', "\r\n")
            } else {
                l.to_string()
            }
        })
        .collect();
    fs::write(inp.join("mixed.fa"), mixed).unwrap();
    fs::write(inp.join("lf.fa"), &lf).unwrap();
    // one 200 kb line (no wrapping) and a long one-line header
    let long: String = format!(
        ">{}\n{}\n",
        "x".repeat(100_000),
        std::str::from_utf8(g).unwrap()
    );
    fs::write(inp.join("longline.fa"), long).unwrap();
    // binary garbage, a text file, a FASTQ cut in the middle of a record, a lone '>'
    let junk: Vec<u8> = (0..5000).map(|_| (rng.next() & 0xff) as u8).collect();
    fs::write(inp.join("binary.fa"), &junk).unwrap();
    let mut gt = vec![b'>'];
    gt.extend(&junk);
    fs::write(inp.join("gt_binary.fa"), gt).unwrap();
    fs::write(
        inp.join("cut.fastq"),
        "@r1\nACGTACGTACGTACGTACGTACGT\n+\nIIIIIIIIIIIIIIIIIIIIIIII\n@r2\nACGTAC",
    )
    .unwrap();
    fs::write(inp.join("gt.fa"), ">").unwrap();
    let out = d.join("out");
    ok(&[
        "run",
        inp.to_str().unwrap(),
        "--index",
        ix.to_str().unwrap(),
        "-o",
        out.to_str().unwrap(),
    ]);
    let (h, rows) = read_tsv(&out.join("results.tsv"));
    assert_eq!(rows.len(), 8, "every input has a row");
    let get = |id: &str| rows.iter().find(|r| r[0] == id).unwrap().clone();
    for id in ["lf", "crlf", "mixed", "longline"] {
        assert_eq!(get(id)[col(&h, "status")], "ok", "{id}");
        assert_eq!(get(id)[col(&h, "call_species")], "Fusarium gamma", "{id}");
        assert_eq!(
            get(id)[col(&h, "n_hashes")],
            get("lf")[col(&h, "n_hashes")],
            "{id}: line endings must not change the sketch"
        );
    }
    for id in ["binary", "gt_binary", "cut", "gt"] {
        let r = get(id);
        assert!(
            matches!(
                r[col(&h, "status")].as_str(),
                "error" | "ok" | "no_match" | "low_confidence"
            ),
            "{id}"
        );
        if r[col(&h, "status")] == "error" {
            assert!(
                !r[col(&h, "error_code")].is_empty(),
                "{id}: has an error code"
            );
        }
    }
    assert_eq!(get("binary")[col(&h, "error_code")], "not_sequence_file");
}

#[test]
fn index_manifest_is_recorded_shown_and_kept_by_convert() {
    let dir = tempfile::tempdir().unwrap();
    let d = dir.path();
    let (genomes, ix, _) = build_reference(d);
    let info = String::from_utf8(ok(&["info", ix.to_str().unwrap()]).stdout).unwrap();
    assert!(
        info.contains("manifest.meta_sha256\t")
            && info.contains("manifest.built_utc\t")
            && info.contains("manifest.format\tGMGIDX01"),
        "{info}"
    );
    let v1 = d.join("ix.v1.gmgi");
    ok(&["convert", ix.to_str().unwrap(), "-o", v1.to_str().unwrap()]);
    let info1 = String::from_utf8(ok(&["info", v1.to_str().unwrap()]).stdout).unwrap();
    assert!(
        info1.contains("manifest.format\tGMGIDX02") && info1.contains("manifest.converted_utc\t"),
        "{info1}"
    );
    let inp = d.join("in");
    fs::create_dir_all(&inp).unwrap();
    fs::write(inp.join("a.fa"), fasta("a", &genomes[0])).unwrap();
    let out = d.join("out");
    ok(&[
        "run",
        inp.to_str().unwrap(),
        "--index",
        v1.to_str().unwrap(),
        "-o",
        out.to_str().unwrap(),
        "--no-report",
    ]);
    let run: serde_json::Value =
        serde_json::from_str(&fs::read_to_string(out.join("run.json")).unwrap()).unwrap();
    assert_eq!(run["index"]["manifest"]["k"], 21);
}

#[test]
fn run_agrees_with_guess_and_can_record_checksums() {
    let dir = tempfile::tempdir().unwrap();
    let d = dir.path();
    let (genomes, ix, mut rng) = build_reference(d);
    let inp = d.join("in");
    fs::create_dir_all(&inp).unwrap();
    fs::write(inp.join("exact.fa"), fasta("exact", &genomes[2])).unwrap();
    fs::write(
        inp.join("mutant.fa"),
        fasta("mutant", &rng.mutate(&genomes[4], 50)),
    )
    .unwrap();
    let out = d.join("out");
    let o = ok(&[
        "run",
        inp.to_str().unwrap(),
        "--index",
        ix.to_str().unwrap(),
        "-o",
        out.to_str().unwrap(),
        "--no-report",
        "--checksums",
    ]);
    assert!(
        String::from_utf8_lossy(&o.stderr).contains("ETA"),
        "progress shows an ETA"
    );
    let results: Vec<serde_json::Value> = fs::read_to_string(out.join("results.jsonl"))
        .unwrap()
        .lines()
        .map(|l| serde_json::from_str(l).unwrap())
        .collect();
    for name in ["exact", "mutant"] {
        let file = inp.join(format!("{name}.fa"));
        let g = ok(&[
            "guess",
            "--index",
            ix.to_str().unwrap(),
            file.to_str().unwrap(),
        ]);
        let guess: serde_json::Value =
            serde_json::from_str(String::from_utf8(g.stdout).unwrap().lines().next().unwrap())
                .unwrap();
        let row = results.iter().find(|r| r["id"] == name).unwrap();
        assert_eq!(row["guess"]["calls"], guess["calls"], "{name}: calls");
        assert_eq!(row["guess"]["top"], guess["top"], "{name}: top hit");
    }
    let run: serde_json::Value =
        serde_json::from_str(&fs::read_to_string(out.join("run.json")).unwrap()).unwrap();
    let sum = run["inputs"][0]["sha256"].as_str().unwrap();
    assert_eq!(sum.len(), 64);
}

#[test]
fn sample_sheet_names_inputs_pairs_files_and_records_expected_labels() {
    let dir = tempfile::tempdir().unwrap();
    let d = dir.path();
    let (genomes, ix, _) = build_reference(d);
    let g = &genomes[2];
    let data = d.join("data");
    fs::create_dir_all(&data).unwrap();
    // names that no automatic rule would pair
    fs::write(data.join("forward.fastq"), fastq("a", g, 0, 3000, 100, 61)).unwrap();
    fs::write(data.join("reverse.fastq"), fastq("b", g, 30, 3000, 100, 61)).unwrap();
    fs::write(data.join("asm.fa"), fasta("asm", &genomes[0])).unwrap();
    let sheet = d.join("sheet.tsv");
    fs::write(
        &sheet,
        "id\tpath\texpected\nstrainA\tdata/forward.fastq\tFusarium gamma\nstrainA\tdata/reverse.fastq\t\nstrainB\tdata/asm.fa\tFusarium\n",
    )
    .unwrap();
    let out = d.join("out");
    ok(&[
        "run",
        "--sample-sheet",
        sheet.to_str().unwrap(),
        "--index",
        ix.to_str().unwrap(),
        "-o",
        out.to_str().unwrap(),
        "--no-report",
        "--min-abund",
        "1",
    ]);
    let (h, rows) = read_tsv(&out.join("results.tsv"));
    assert_eq!(rows.len(), 2);
    let a = rows.iter().find(|r| r[0] == "strainA").unwrap();
    let b = rows.iter().find(|r| r[0] == "strainB").unwrap();
    assert_eq!(a[col(&h, "call_species")], "Fusarium gamma");
    assert_eq!(
        (
            a[col(&h, "expected")].as_str(),
            a[col(&h, "expected_match")].as_str()
        ),
        ("Fusarium gamma", "yes")
    );
    // strainB is Aspergillus alpha; the sheet expected Fusarium, so the mismatch is visible
    assert_eq!(b[col(&h, "expected_match")], "no");
    // a bad sheet is a usage error with a plain message
    let bad = d.join("bad.tsv");
    fs::write(&bad, "name\tfile\nx\ty\n").unwrap();
    let o = gmg(&[
        "run",
        "--sample-sheet",
        bad.to_str().unwrap(),
        "--index",
        ix.to_str().unwrap(),
        "-o",
        d.join("o2").to_str().unwrap(),
    ]);
    assert!(!o.status.success() && String::from_utf8_lossy(&o.stderr).contains("'id' and 'path'"));
}

/// The documented workflow (docs/tutorial.md) on the tiny fixture must keep producing the same results.tsv.
/// Regenerate with `UPDATE_GOLDEN=1 cargo test -p gmg-cli --test run golden`.
#[test]
fn golden_results_tsv_for_the_tutorial_workflow() {
    let dir = tempfile::tempdir().unwrap();
    let d = dir.path();
    let (genomes, ix, mut rng) = build_reference(d);
    let inp = d.join("genomes");
    fs::create_dir_all(&inp).unwrap();
    fs::write(inp.join("exact_fusarium.fa"), fasta("a", &genomes[2])).unwrap();
    fs::write(
        inp.join("mutant_agaricus.fa"),
        fasta("b", &rng.mutate(&genomes[4], 20)),
    )
    .unwrap();
    fs::write(
        inp.join("divergent_agaricus.fa"),
        fasta("c", &rng.mutate(&genomes[5], 150)),
    )
    .unwrap();
    fs::write(inp.join("unrelated.fa"), fasta("d", &rng.genome(150_000))).unwrap();
    fs::write(inp.join("empty.fa"), b"").unwrap();
    let out = d.join("results");
    ok(&[
        "run",
        inp.to_str().unwrap(),
        "--index",
        ix.to_str().unwrap(),
        "-o",
        out.to_str().unwrap(),
        "--no-report",
    ]);
    // drop the path column (machine specific); everything else is deterministic
    let text = fs::read_to_string(out.join("results.tsv")).unwrap();
    let mut norm = String::new();
    for line in text.lines() {
        let f: Vec<&str> = line.split('\t').collect();
        norm.push_str(&[&f[..1], &f[2..]].concat().join("\t"));
        norm.push('\n');
    }
    let golden = Path::new(env!("CARGO_MANIFEST_DIR")).join("tests/golden/results.tsv");
    if std::env::var("UPDATE_GOLDEN").is_ok() {
        fs::write(&golden, &norm).unwrap();
    }
    let want =
        fs::read_to_string(&golden).expect("tests/golden/results.tsv (UPDATE_GOLDEN=1 creates it)");
    assert_eq!(
        norm, want,
        "results.tsv changed; if intended, regenerate with UPDATE_GOLDEN=1"
    );
}

#[test]
fn completions_and_man_page_are_generated() {
    for shell in ["bash", "zsh", "fish"] {
        let o = ok(&["completions", shell]);
        let t = String::from_utf8_lossy(&o.stdout);
        assert!(
            t.contains("gmg") && t.contains("report"),
            "{shell} completion"
        );
    }
    let man = String::from_utf8(ok(&["man"]).stdout).unwrap();
    assert!(
        man.starts_with(".ie") || man.contains(".TH gmg"),
        "roff header"
    );
    assert!(man.contains("gmg-run") || man.contains("run"));
}

fn golden_check(name: &str, got: &str) {
    let golden = Path::new(env!("CARGO_MANIFEST_DIR"))
        .join("tests/golden")
        .join(name);
    if std::env::var("UPDATE_GOLDEN").is_ok() {
        fs::write(&golden, got).unwrap();
    }
    let want = fs::read_to_string(&golden)
        .unwrap_or_else(|_| panic!("tests/golden/{name} (UPDATE_GOLDEN=1 creates it)"));
    assert_eq!(
        got, want,
        "{name} changed; if intended, regenerate with UPDATE_GOLDEN=1"
    );
}

#[test]
fn guess_formats_json_tsv_and_table_are_stable() {
    let dir = tempfile::tempdir().unwrap();
    let d = dir.path();
    let (genomes, ix, mut rng) = build_reference(d);
    let q: Vec<(String, Vec<u8>)> = vec![
        ("exact_fusarium.fa".into(), genomes[2].clone()),
        ("mutant_agaricus.fa".into(), rng.mutate(&genomes[4], 20)),
        ("divergent_agaricus.fa".into(), rng.mutate(&genomes[5], 150)),
        ("unrelated.fa".into(), rng.genome(150_000)),
    ];
    let mut paths = Vec::new();
    for (n, g) in &q {
        let p = d.join(n);
        fs::write(&p, fasta(n, g)).unwrap();
        paths.push(p.to_str().unwrap().to_string());
    }
    let run = |fmt: &str| {
        let mut args = vec!["guess", "--index", ix.to_str().unwrap(), "--format", fmt];
        args.extend(paths.iter().map(|s| s.as_str()));
        String::from_utf8(ok(&args).stdout)
            .unwrap()
            .replace(d.to_str().unwrap(), "<dir>")
    };
    golden_check("guess.tsv", &run("tsv"));
    golden_check("guess.table.txt", &run("table"));
    // json: one object per line with the documented keys
    let json = run("json");
    let first: serde_json::Value = serde_json::from_str(json.lines().next().unwrap()).unwrap();
    for k in [
        "query",
        "status",
        "best_rank",
        "calls",
        "withheld",
        "tentative_genus",
        "top",
        "hits",
    ] {
        assert!(first.get(k).is_some(), "json has {k}");
    }
    assert_eq!(json.lines().count(), 4);
    let third: serde_json::Value = serde_json::from_str(json.lines().nth(2).unwrap()).unwrap();
    assert_eq!(third["withheld"]["rank"], "genus");
    assert_eq!(third["withheld"]["code"], "ani_below_cutoff");
}

#[test]
fn quiet_and_verbose_control_progress_lines() {
    let dir = tempfile::tempdir().unwrap();
    let d = dir.path();
    let (genomes, ix, mut rng) = build_reference(d);
    let inp = d.join("in");
    fs::create_dir_all(&inp).unwrap();
    fs::write(
        inp.join("divergent.fa"),
        fasta("x", &rng.mutate(&genomes[4], 150)),
    )
    .unwrap();
    let go = |extra: &[&str]| {
        let mut a = vec![
            "run",
            inp.to_str().unwrap(),
            "--index",
            ix.to_str().unwrap(),
            "--no-report",
            "-o",
        ];
        let out = d.join(format!("o{}", extra.join("")));
        let outs = out.to_str().unwrap().to_string();
        a.push(&outs);
        a.extend(extra);
        String::from_utf8_lossy(&ok(&a).stderr).to_string()
    };
    assert_eq!(go(&["-q"]).trim(), "", "quiet prints nothing");
    let v = go(&["-v"]);
    assert!(v.contains("withheld genus: ani_below_cutoff"), "{v}");
    assert!(
        !go(&[]).contains("withheld"),
        "default does not print the reason"
    );
}
