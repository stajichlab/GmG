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
    // the same holds with the abundance filter, because abundances are summed over the files of a sample: the default cutoff (2) and
    // `auto` give the same hashes, the same coverage estimate and the same call for the merged sample and the concatenated file (#32)
    for extra in [vec![], vec!["--min-abund", "auto"]] {
        let outc = d.join(format!("out_cut{}", extra.join("")));
        let mut args = vec![
            "run",
            inp.to_str().unwrap(),
            "--index",
            ix.to_str().unwrap(),
            "-o",
            outc.to_str().unwrap(),
            "--no-report",
        ];
        args.extend_from_slice(&extra);
        ok(&args);
        let (hc, rc) = read_tsv(&outc.join("results.tsv"));
        let find = |id: &str| rc.iter().find(|r| r[0] == id).unwrap().clone();
        let (ca, cw) = (find("A_S1"), find("whole"));
        for c in [
            "status",
            "call_species",
            "n_hashes",
            "shared",
            "total_bases",
            "coverage_est",
        ] {
            assert_eq!(ca[col(&hc, c)], cw[col(&hc, c)], "{c} with {extra:?}");
        }
        assert!(!ca[col(&hc, "coverage_est")].is_empty(), "{extra:?}");
    }

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
        info1.contains("manifest.format\tGMGIDX03")
            && info1.contains("manifest.converted_utc\t")
            && info1.contains("packed (format v3)"),
        "{info1}"
    );
    // the older compact format is still available and gives the same answers
    let old = d.join("ix.old.gmgi");
    ok(&[
        "convert",
        ix.to_str().unwrap(),
        "-o",
        old.to_str().unwrap(),
        "--format",
        "v1",
    ]);
    assert!(
        String::from_utf8(ok(&["info", old.to_str().unwrap()]).stdout)
            .unwrap()
            .contains("manifest.format\tGMGIDX02")
    );
    let g = |index: &std::path::Path| {
        let q = d.join("q.fa");
        fs::write(&q, fasta("q", &genomes[1])).unwrap();
        String::from_utf8(
            ok(&[
                "guess",
                "--index",
                index.to_str().unwrap(),
                q.to_str().unwrap(),
            ])
            .stdout,
        )
        .unwrap()
        .replace(d.to_str().unwrap(), "")
    };
    assert_eq!(g(&v1), g(&old));
    assert_eq!(g(&v1), g(&ix));
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

#[test]
fn size_flags_come_from_the_reference_size_range_of_the_called_phylum() {
    let dir = tempfile::tempdir().unwrap();
    let d = dir.path();
    let mut rng = Rng(0xABCDEF1234567);
    // 60 Ascomycota of 150-250 kb, one genus each; scaled 10 -> ~15-25k hashes
    let refs = d.join("refs");
    fs::create_dir(&refs).unwrap();
    let mut meta = String::from("id\ttaxid\tphylum\tclass\torder\tfamily\tgenus\tspecies\n");
    let mut files = Vec::new();
    let mut gs = Vec::new();
    for i in 0..60 {
        let g = rng.genome(150_000 + (i % 10) * 11_000);
        let f = refs.join(format!("g{i}.fa"));
        fs::write(&f, fasta(&format!("g{i}"), &g)).unwrap();
        meta.push_str(&format!(
            "g{i}\t{}\tAscomycota\tC\tO\tF\tGenus{i}\tGenus{i} sp\n",
            1000 + i
        ));
        files.push(f.display().to_string());
        gs.push(g);
    }
    fs::write(d.join("meta.tsv"), meta).unwrap();
    let sk = d.join("sk");
    let mut args = vec!["sketch", "-k", "21", "-s", "10", "-o", sk.to_str().unwrap()];
    args.extend(files.iter().map(|s| s.as_str()));
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
    let inp = d.join("in");
    fs::create_dir_all(&inp).unwrap();
    // a 20 kb piece of g3, the whole of g3, and g3 followed by 700 kb of unrelated sequence
    fs::write(inp.join("piece.fa"), fasta("p", &gs[3][..20_000])).unwrap();
    fs::write(inp.join("whole.fa"), fasta("w", &gs[3])).unwrap();
    let mut big = gs[3].clone();
    big.extend(rng.genome(700_000));
    fs::write(inp.join("big.fa"), fasta("b", &big)).unwrap();
    let out = d.join("out");
    ok(&[
        "run",
        inp.to_str().unwrap(),
        "--index",
        ix.to_str().unwrap(),
        "-o",
        out.to_str().unwrap(),
        "--no-report",
    ]);
    let (h, rows) = read_tsv(&out.join("results.tsv"));
    let flags = |id: &str| rows.iter().find(|r| r[0] == id).unwrap()[col(&h, "qc_flags")].clone();
    assert!(
        flags("piece").contains("short_genome"),
        "20 kb of a 150-250 kb lineage: {}",
        flags("piece")
    );
    assert!(
        !flags("whole").contains("short_genome") && !flags("whole").contains("long_genome"),
        "{}",
        flags("whole")
    );
    assert!(flags("big").contains("long_genome"), "{}", flags("big"));
    // the range and its basis are recorded for the report
    let jl = fs::read_to_string(out.join("results.jsonl")).unwrap();
    let rec: serde_json::Value = jl
        .lines()
        .map(|l| serde_json::from_str::<serde_json::Value>(l).unwrap())
        .find(|v| v["id"] == "piece")
        .unwrap();
    assert!(rec["qc"]["size_range_basis"]
        .as_str()
        .unwrap()
        .starts_with("Ascomycota (n = 60"));
    assert!(rec["qc"]["size_range_mb"][0].as_f64().unwrap() > 0.1);
    // the report says what is typical for the lineage
    ok(&["report", out.to_str().unwrap()]);
    let html = fs::read_to_string(out.join("report.html")).unwrap();
    assert!(
        html.contains("Typical size"),
        "qc panel shows the typical size"
    );
}

#[test]
fn per_rank_taxids_and_original_labels_reach_the_results() {
    let dir = tempfile::tempdir().unwrap();
    let d = dir.path();
    let (genomes, _ix, _) = build_reference(d);
    // the same reference with NCBI taxids per rank and the original label, as scripts/normalize_taxonomy.py writes it
    let lineage = [
        ("Aspergillus alpha", 5052, "Aspergillus", 5052),
        ("Aspergillus beta", 5053, "Aspergillus", 5052),
        ("Fusarium gamma", 5507, "Fusarium", 5506),
        ("Fusarium delta", 5508, "Fusarium", 5506),
        ("Agaricus eps", 5340, "Agaricus", 5341),
        ("Agaricus zeta", 5342, "Agaricus", 5341),
    ];
    let mut meta = String::from("id\ttaxid\tphylum\tclass\torder\tfamily\tgenus\tspecies\tcluster\tphylum_taxid\tclass_taxid\torder_taxid\tfamily_taxid\tgenus_taxid\tspecies_taxid\tlabel\n");
    for (i, (sp, st, genus, gt)) in lineage.iter().enumerate() {
        let (ph, pt) = if i < 4 {
            ("Ascomycota", 4890)
        } else {
            ("Basidiomycota", 5204)
        };
        meta.push_str(&format!("g{i}\t{st}\t{ph}\tC\tO\tF\t{genus}\t{sp}\t\t{pt}\t0\t0\t0\t{gt}\t{st}\toriginal label {i}\n"));
    }
    let meta_path = d.join("meta_ncbi.tsv");
    fs::write(&meta_path, meta).unwrap();
    let ix2 = d.join("ix2.gmgi");
    ok(&[
        "build",
        "--meta",
        meta_path.to_str().unwrap(),
        "--sketches",
        d.join("sk").to_str().unwrap(),
        "-o",
        ix2.to_str().unwrap(),
    ]);
    let inp = d.join("in");
    fs::create_dir_all(&inp).unwrap();
    fs::write(inp.join("q.fa"), fasta("q", &genomes[2])).unwrap();
    let out = d.join("out");
    ok(&[
        "run",
        inp.to_str().unwrap(),
        "--index",
        ix2.to_str().unwrap(),
        "-o",
        out.to_str().unwrap(),
        "--no-report",
    ]);
    let (h, rows) = read_tsv(&out.join("results.tsv"));
    let r = &rows[0];
    assert_eq!(r[col(&h, "call_species")], "Fusarium gamma");
    assert_eq!(r[col(&h, "taxid_phylum")], "4890");
    assert_eq!(r[col(&h, "taxid_genus")], "5506");
    assert_eq!(r[col(&h, "taxid_species")], "5507");
    assert_eq!(
        r[col(&h, "taxid_class")],
        "",
        "a rank without a taxid stays blank"
    );
    assert_eq!(r[col(&h, "top_hit_label")], "original label 2");
    // a header without the required columns is a plain error
    let bad = d.join("bad.tsv");
    fs::write(&bad, "id\ttaxid\nx\t1\n").unwrap();
    let o = gmg(&[
        "build",
        "--meta",
        bad.to_str().unwrap(),
        "--sketches",
        d.join("sk").to_str().unwrap(),
        "-o",
        d.join("b.gmgi").to_str().unwrap(),
    ]);
    assert!(
        !o.status.success()
            && String::from_utf8_lossy(&o.stderr).contains("header lacks the column")
    );
}

#[test]
fn exclusion_hides_the_query_its_species_or_its_genus_from_the_reference() {
    let dir = tempfile::tempdir().unwrap();
    let d = dir.path();
    let (genomes, ix, _) = build_reference(d);
    let data = d.join("data");
    fs::create_dir_all(&data).unwrap();
    // the query is reference genome g2 (Fusarium gamma; its genus mate is g3)
    fs::write(data.join("q.fa"), fasta("q", &genomes[2])).unwrap();
    let sheet = d.join("sheet.tsv");
    fs::write(
        &sheet,
        "id\tpath\tref\texpected\nq\tdata/q.fa\tg2\tFusarium gamma\n",
    )
    .unwrap();
    fs::write(d.join("ids.txt"), "# comment\ng2\n").unwrap();
    let run = |extra: &[&str]| -> serde_json::Value {
        let out = d.join(format!("out{}", extra.join("")));
        let mut a = vec![
            "run",
            "--sample-sheet",
            sheet.to_str().unwrap(),
            "--index",
            ix.to_str().unwrap(),
            "-o",
            out.to_str().unwrap(),
            "--no-report",
        ];
        a.extend_from_slice(extra);
        ok(&a);
        let line = fs::read_to_string(out.join("results.jsonl")).unwrap();
        serde_json::from_str(line.lines().next().unwrap()).unwrap()
    };
    let plain = run(&[]);
    assert_eq!(
        plain["guess"]["calls"]["species"]["taxon"],
        "Fusarium gamma"
    );
    assert!(plain.get("excluded").is_none());
    for (extra, n) in [
        (vec!["--exclude-level", "self"], 1),
        (vec!["--exclude-level", "species"], 1),
        (vec!["--exclude-level", "genus"], 2),
        (
            vec!["--exclude-ids", d.join("ids.txt").to_str().unwrap()],
            1,
        ),
    ] {
        let r = run(&extra);
        assert_eq!(r["excluded"]["n_genomes"], n, "{extra:?}");
        assert_ne!(
            r["guess"]["calls"]["species"]["taxon"], "Fusarium gamma",
            "{extra:?}"
        );
        assert_eq!(r["status"], "no_match", "{extra:?}");
        // the phylum split (hashes found in one phylum only) must not see the hidden genomes either (it used to: the hidden genome
        // then made every leave-out phylum look perfectly resolved)
        if extra[0] == "--exclude-level" && extra[1] == "self" {
            assert!(
                plain["qc"]["phylum_split"]
                    .as_array()
                    .is_some_and(|a| !a.is_empty()),
                "the unhidden run has a split"
            );
            let hidden_hashes: u64 = r["qc"]["phylum_split"]
                .as_array()
                .unwrap()
                .iter()
                .filter(|p| p["phylum"] == "Ascomycota")
                .map(|p| p["hashes"].as_u64().unwrap())
                .sum();
            let plain_hashes: u64 = plain["qc"]["phylum_split"]
                .as_array()
                .unwrap()
                .iter()
                .filter(|p| p["phylum"] == "Ascomycota")
                .map(|p| p["hashes"].as_u64().unwrap())
                .sum();
            assert!(
                hidden_hashes < plain_hashes / 4,
                "Ascomycota-exclusive hashes {hidden_hashes} vs {plain_hashes} unhidden"
            );
        }
    }
}

#[test]
fn stdin_and_max_bases_are_inputs_like_any_other() {
    use std::io::Write;
    use std::process::Stdio;
    let dir = tempfile::tempdir().unwrap();
    let d = dir.path();
    let (genomes, ix, _) = build_reference(d);
    let fa = fasta("q", &genomes[2]);
    fs::write(d.join("q.fa"), &fa).unwrap();
    let guess = |args: &[&str], stdin: Option<&[u8]>| -> serde_json::Value {
        let mut c = Command::new(env!("CARGO_BIN_EXE_gmg"));
        c.args(["guess", "--index", ix.to_str().unwrap(), "--format", "json"])
            .args(args);
        c.stdin(if stdin.is_some() {
            Stdio::piped()
        } else {
            Stdio::null()
        })
        .stdout(Stdio::piped())
        .stderr(Stdio::piped());
        let mut ch = c.spawn().unwrap();
        if let Some(b) = stdin {
            ch.stdin.take().unwrap().write_all(b).unwrap();
        }
        let o = ch.wait_with_output().unwrap();
        assert!(o.status.success(), "{}", String::from_utf8_lossy(&o.stderr));
        serde_json::from_slice(&o.stdout).unwrap()
    };
    let file = guess(&[d.join("q.fa").to_str().unwrap()], None);
    let piped = guess(&["-"], Some(fa.as_bytes()));
    assert_eq!(file["calls"], piped["calls"]);
    assert_eq!(piped["calls"]["species"]["taxon"], "Fusarium gamma");
    // a gzipped stream is detected from its bytes
    let mut gz = flate2::write::GzEncoder::new(Vec::new(), flate2::Compression::fast());
    gz.write_all(fa.as_bytes()).unwrap();
    let zipped = guess(&["-"], Some(&gz.finish().unwrap()));
    assert_eq!(file["calls"], zipped["calls"]);
    // reading stops after the record that reaches the limit: 10 kb records, 20 kb limit = a tenth of the genome
    let many: String = genomes[2]
        .chunks(10_000)
        .enumerate()
        .map(|(i, c)| fasta(&format!("c{i}"), c))
        .collect();
    fs::write(d.join("many.fa"), &many).unwrap();
    let many_p = d.join("many.fa");
    let all = guess(&[many_p.to_str().unwrap()], None);
    let part = guess(&["--max-bases", "20000", many_p.to_str().unwrap()], None);
    assert!(part["query_hashes"].as_u64().unwrap() < all["query_hashes"].as_u64().unwrap() / 4);
}

#[test]
fn decoy_genomes_make_non_fungal_inputs_say_so() {
    let dir = tempfile::tempdir().unwrap();
    let d = dir.path();
    let (genomes, _ix, mut rng) = build_reference(d);
    // two non-fungal decoys added to the same reference, listed in the manifest
    let plant = rng.genome(300_000);
    let worm = rng.genome(250_000);
    let refs = d.join("refs");
    fs::write(refs.join("plant.fa"), fasta("plant", &plant)).unwrap();
    fs::write(refs.join("worm.fa"), fasta("worm", &worm)).unwrap();
    ok(&[
        "sketch",
        "-k",
        "21",
        "-s",
        "10",
        "-o",
        d.join("sk").to_str().unwrap(),
        refs.join("plant.fa").to_str().unwrap(),
        refs.join("worm.fa").to_str().unwrap(),
    ]);
    let mut meta = fs::read_to_string(d.join("meta.tsv")).unwrap();
    meta.push_str("plant\t3702\tStreptophyta\tMagnoliopsida\tBrassicales\tBrassicaceae\tArabidopsis\tArabidopsis thaliana\n");
    meta.push_str("worm\t6239\tNematoda\tChromadorea\tRhabditida\tRhabditidae\tCaenorhabditis\tCaenorhabditis elegans\n");
    let meta_path = d.join("meta_decoys.tsv");
    fs::write(&meta_path, meta).unwrap();
    let ix = d.join("ix_decoys.gmgi");
    ok(&[
        "build",
        "--meta",
        meta_path.to_str().unwrap(),
        "--sketches",
        d.join("sk").to_str().unwrap(),
        "-o",
        ix.to_str().unwrap(),
        "--non-fungal-phyla",
        "Streptophyta, Nematoda",
    ]);
    assert!(
        String::from_utf8(ok(&["info", ix.to_str().unwrap()]).stdout)
            .unwrap()
            .contains("non_fungal_phyla")
    );
    let inp = d.join("in");
    fs::create_dir_all(&inp).unwrap();
    fs::write(
        inp.join("plant_relative.fa"),
        fasta("p", &rng.mutate(&plant, 30)),
    )
    .unwrap(); // ~0.98 ANI: a close relative
    fs::write(inp.join("worm_copy.fa"), fasta("w", &worm)).unwrap();
    fs::write(
        inp.join("plant_distant.fa"),
        fasta("d", &rng.mutate(&plant, 200)),
    )
    .unwrap(); // ~0.85 ANI: a decoy is nearest, but not a relative
    fs::write(inp.join("fungus.fa"), fasta("f", &genomes[2])).unwrap();
    fs::write(inp.join("unrelated.fa"), fasta("u", &rng.genome(150_000))).unwrap();
    let out = d.join("out");
    ok(&[
        "run",
        inp.to_str().unwrap(),
        "--index",
        ix.to_str().unwrap(),
        "-o",
        out.to_str().unwrap(),
        "-q",
    ]);
    let (h, rows) = read_tsv(&out.join("results.tsv"));
    let get = |id: &str| rows.iter().find(|r| r[0] == id).unwrap().clone();
    let flags = |id: &str| get(id)[col(&h, "qc_flags")].clone();
    assert!(
        flags("plant_relative").contains("not_fungal")
            && get("plant_relative")[col(&h, "call_phylum")] == "Streptophyta"
    );
    assert!(
        flags("worm_copy").contains("not_fungal")
            && get("worm_copy")[col(&h, "call_species")] == "Caenorhabditis elegans"
    );
    assert!(
        !flags("plant_distant").contains("not_fungal")
            && flags("plant_distant").contains("distant_non_fungal_hit"),
        "a distant decoy hit is not a verdict: {} (ANI {})",
        flags("plant_distant"),
        get("plant_distant")[col(&h, "ani")]
    );
    assert!(
        get("plant_distant")[col(&h, "call_phylum")] != "Streptophyta",
        "a distant decoy hit does not become a phylum call"
    );
    assert!(
        !flags("fungus").contains("not_fungal"),
        "a fungus is not flagged: {}",
        flags("fungus")
    );
    assert_eq!(get("fungus")[col(&h, "call_species")], "Fusarium gamma");
    assert!(
        !flags("unrelated").contains("not_fungal"),
        "no match is not a decoy match"
    );
    let html = fs::read_to_string(out.join("report.html")).unwrap();
    assert!(html.contains("Probably not a fungus"));
    // an index without the manifest list never flags
    let out2 = d.join("out2");
    ok(&[
        "run",
        inp.to_str().unwrap(),
        "--index",
        d.join("ix.gmgi").to_str().unwrap(),
        "-o",
        out2.to_str().unwrap(),
        "-q",
        "--no-report",
    ]);
    let (h2, r2) = read_tsv(&out2.join("results.tsv"));
    assert!(r2
        .iter()
        .all(|r| !r[col(&h2, "qc_flags")].contains("not_fungal")));
}

#[test]
fn run_reads_standard_input_next_to_files() {
    use std::io::Write;
    use std::process::Stdio;
    let dir = tempfile::tempdir().unwrap();
    let d = dir.path();
    let (genomes, ix, _) = build_reference(d);
    let fa = fasta("q", &genomes[2]);
    fs::create_dir_all(d.join("in")).unwrap();
    fs::write(d.join("in/file.fa"), &fa).unwrap();
    let run = |extra: &[&str], stdin: &[u8]| -> std::process::Output {
        let mut ch = Command::new(env!("CARGO_BIN_EXE_gmg"))
            .args(["run", "--index", ix.to_str().unwrap(), "-q", "--no-report"])
            .args(extra)
            .stdin(Stdio::piped())
            .stdout(Stdio::piped())
            .stderr(Stdio::piped())
            .spawn()
            .unwrap();
        let _ = ch.stdin.take().unwrap().write_all(stdin); // the process may exit first (refused input)
        ch.wait_with_output().unwrap()
    };
    // a file and standard input in one run: two rows, the same call
    let out = d.join("out");
    let o = run(
        &[
            "-o",
            out.to_str().unwrap(),
            d.join("in").to_str().unwrap(),
            "-",
        ],
        fa.as_bytes(),
    );
    assert!(o.status.success(), "{}", String::from_utf8_lossy(&o.stderr));
    let (h, rows) = read_tsv(&out.join("results.tsv"));
    assert_eq!(rows.len(), 2);
    let sp = |id: &str| rows.iter().find(|r| r[0] == id).unwrap()[col(&h, "call_species")].clone();
    assert_eq!(sp("stdin"), "Fusarium gamma");
    assert_eq!(sp("file"), sp("stdin"));
    // FASTQ on standard input needs --kind reads
    let fq = "@r1\nACGTACGTACGTACGTACGTACGTAAA\n+\nIIIIIIIIIIIIIIIIIIIIIIIIIII\n";
    let out2 = d.join("out2");
    let o = run(&["-o", out2.to_str().unwrap(), "-"], fq.as_bytes());
    let (h2, rows2) = read_tsv(&out2.join("results.tsv"));
    assert_eq!(rows2[0][col(&h2, "status")], "error");
    assert!(String::from_utf8_lossy(&o.stderr).len() < 10_000);
    assert!(
        fs::read_to_string(out2.join("results.jsonl"))
            .unwrap()
            .contains("--kind reads"),
        "the error row says how to fix it"
    );
    let out3 = d.join("out3");
    let o = run(
        &["-o", out3.to_str().unwrap(), "--kind", "reads", "-"],
        fq.as_bytes(),
    );
    assert!(o.status.success(), "{}", String::from_utf8_lossy(&o.stderr));
    // standard input twice is refused
    let o = run(&["-o", d.join("out4").to_str().unwrap(), "-", "-"], b"");
    assert!(!o.status.success());
}

#[test]
fn weak_phylum_flags_a_phylum_called_from_a_distant_assembly() {
    let dir = tempfile::tempdir().unwrap();
    let d = dir.path();
    let (genomes, ix, mut rng) = build_reference(d);
    let inp = d.join("in");
    fs::create_dir_all(&inp).unwrap();
    fs::write(inp.join("copy.fa"), fasta("c", &genomes[2])).unwrap();
    fs::write(
        inp.join("distant.fa"),
        fasta("d", &rng.mutate(&genomes[2], 250)),
    )
    .unwrap(); // ~0.81 ANI to its reference genome
    let out = d.join("out");
    ok(&[
        "run",
        inp.to_str().unwrap(),
        "--index",
        ix.to_str().unwrap(),
        "-o",
        out.to_str().unwrap(),
        "-q",
        "--no-report",
    ]);
    let (h, rows) = read_tsv(&out.join("results.tsv"));
    let get = |id: &str| rows.iter().find(|r| r[0] == id).unwrap().clone();
    let (c, dist) = (get("copy"), get("distant"));
    assert!(!c[col(&h, "qc_flags")].contains("weak_phylum"));
    assert_eq!(
        dist[col(&h, "call_phylum")],
        "Ascomycota",
        "the distant input still gets its phylum (flag, not withhold)"
    );
    assert!(
        dist[col(&h, "qc_flags")].contains("weak_phylum"),
        "flags: {} ANI {}",
        dist[col(&h, "qc_flags")],
        dist[col(&h, "ani")]
    );
}

/// Error-bearing single-end reads of length `len` at coverage `cov` from `g` (positions and errors from `rng`), as FASTQ text.
fn fastq_reads(rng: &mut Rng, g: &[u8], cov: f64, len: usize, err_per_thousand: u64) -> String {
    let n = (cov * g.len() as f64 / len as f64).round() as usize;
    let mut out = String::new();
    for i in 0..n {
        let start = (rng.next() as usize) % (g.len() - len);
        let read = rng.mutate(&g[start..start + len], err_per_thousand);
        out.push_str(&format!(
            "@r{i}\n{}\n+\n{}\n",
            String::from_utf8(read).unwrap(),
            "I".repeat(len)
        ));
    }
    out
}

#[test]
fn coverage_is_estimated_from_the_abundance_histogram_and_auto_cutoff_runs() {
    let dir = tempfile::tempdir().unwrap();
    let d = dir.path();
    let (genomes, ix, mut rng) = build_reference(d);
    let g = &genomes[2];
    let inp = d.join("reads");
    fs::create_dir_all(&inp).unwrap();
    let covs = [0.5f64, 2.0, 8.0, 30.0];
    for c in covs {
        // 0.5% substitution errors, 150 bp reads
        fs::write(
            inp.join(format!("c{c}.fastq")),
            fastq_reads(&mut rng, g, c, 150, 5),
        )
        .unwrap();
    }
    let out = d.join("out");
    ok(&[
        "run",
        inp.to_str().unwrap(),
        "--index",
        ix.to_str().unwrap(),
        "-o",
        out.to_str().unwrap(),
        "-q",
        "--no-report",
    ]);
    let (h, rows) = read_tsv(&out.join("results.tsv"));
    let mut seen = vec![];
    for c in covs {
        let r = rows.iter().find(|r| r[0] == format!("c{c}")).unwrap();
        let est: f64 = r[col(&h, "coverage_est")]
            .parse()
            .expect("a coverage estimate");
        seen.push((c, est));
        // the estimate follows the truth within 30% (errors, repeated k-mers and small histograms are noise)
        assert!((est - c).abs() <= 0.3 * c, "true {c}x, estimated {est}x");
    }
    let line = fs::read_to_string(out.join("results.jsonl")).unwrap();
    assert!(
        line.contains("\"coverage_basis\":\"peak\"")
            && line.contains("\"coverage_basis\":\"poisson\""),
        "{seen:?}"
    );
    // --min-abund auto: runs for reads and assemblies, records the cutoff it used, and keeps the call at depth
    let asm = d.join("asm");
    fs::create_dir_all(&asm).unwrap();
    fs::write(asm.join("g.fa"), fasta("g", g)).unwrap();
    for (dir_in, kind) in [(&inp, "reads"), (&asm, "assembly")] {
        let o = d.join(format!("auto_{kind}"));
        ok(&[
            "run",
            dir_in.to_str().unwrap(),
            "--index",
            ix.to_str().unwrap(),
            "-o",
            o.to_str().unwrap(),
            "-q",
            "--no-report",
            "--min-abund",
            "auto",
        ]);
        let (h2, r2) = read_tsv(&o.join("results.tsv"));
        assert!(
            r2.iter().all(|r| r[col(&h2, "status")] != "error"),
            "{kind}"
        );
    }
    let (h3, r3) = read_tsv(&d.join("auto_reads").join("results.tsv"));
    let deep = r3.iter().find(|r| r[0] == "c30").unwrap();
    assert_eq!(deep[col(&h3, "call_species")], "Fusarium gamma");
    let rec: serde_json::Value = serde_json::from_str(
        fs::read_to_string(d.join("auto_reads").join("results.jsonl"))
            .unwrap()
            .lines()
            .find(|l| l.contains("\"c30\""))
            .unwrap(),
    )
    .unwrap();
    let used = rec["qc"]["min_abund"].as_u64().unwrap();
    assert!((2..=10).contains(&used), "auto cutoff at 30x: {used}");
    // guess takes it too
    let o = ok(&[
        "guess",
        "--index",
        ix.to_str().unwrap(),
        "--min-abund",
        "auto",
        "--format",
        "json",
        inp.join("c30.fastq").to_str().unwrap(),
    ]);
    let j: serde_json::Value = serde_json::from_slice(&o.stdout).unwrap();
    assert_eq!(j["calls"]["species"]["taxon"], "Fusarium gamma");
    // a bad value is refused
    assert!(!gmg(&[
        "guess",
        "--index",
        ix.to_str().unwrap(),
        "--min-abund",
        "many",
        inp.join("c30.fastq").to_str().unwrap()
    ])
    .status
    .success());
}

/// Result records of one run on `inputs` (id, FASTA text), with extra `run` arguments.
fn run_records(
    d: &Path,
    ix: &Path,
    tag: &str,
    inputs: &[(&str, String)],
    extra: &[&str],
) -> Vec<serde_json::Value> {
    let data = d.join(format!("data_{tag}"));
    fs::create_dir_all(&data).unwrap();
    for (id, text) in inputs {
        fs::write(data.join(format!("{id}.fa")), text).unwrap();
    }
    let out = d.join(format!("out_{tag}"));
    let mut a = vec![
        "run",
        data.to_str().unwrap(),
        "--index",
        ix.to_str().unwrap(),
        "-o",
        out.to_str().unwrap(),
        "--no-report",
        "-q",
    ];
    a.extend_from_slice(extra);
    ok(&a);
    let mut recs: Vec<serde_json::Value> = fs::read_to_string(out.join("results.jsonl"))
        .unwrap()
        .lines()
        .map(|l| serde_json::from_str(l).unwrap())
        .collect();
    recs.sort_by_key(|r| r["id"].as_str().unwrap().to_string());
    recs
}

fn flags_of(r: &serde_json::Value) -> Vec<String> {
    r["qc_flags"]
        .as_array()
        .unwrap()
        .iter()
        .map(|f| f.as_str().unwrap().to_string())
        .collect()
}

#[test]
fn a_second_genus_in_the_input_is_flagged_without_changing_the_call() {
    let dir = tempfile::tempdir().unwrap();
    let d = dir.path();
    let (genomes, ix, _) = build_reference(d);
    // pure: g0 (Aspergillus alpha). Mixed: g0 plus 50 kb of g4 (Agaricus eps) = 20% of the bases. Tiny: 1 kb of g4 (0.5%).
    // Same genus: g0 plus 50 kb of g1 (Aspergillus beta), which the pass must not call a second genus.
    let mix = |extra: &[u8]| {
        let mut s = genomes[0].clone();
        s.extend_from_slice(extra);
        fasta("q", &s)
    };
    let inputs = [
        ("pure", fasta("q", &genomes[0])),
        ("mixed", mix(&genomes[4][..50_000])),
        ("tiny", mix(&genomes[4][..1_000])),
        ("samegenus", mix(&genomes[1][..50_000])),
    ];
    let recs = run_records(d, &ix, "plain", &inputs, &[]);
    let by = |id: &str| recs.iter().find(|r| r["id"] == id).unwrap();
    for r in &recs {
        assert_eq!(r["schema_version"], 8);
        // the call is the same for every input
        assert_eq!(
            r["guess"]["calls"]["genus"]["taxon"], "Aspergillus",
            "{}",
            r["id"]
        );
    }
    // the species is called for the pure input only (the minority lowers the ANI below the species cutoff), as without the detector
    assert_eq!(
        by("pure")["guess"]["calls"]["species"]["taxon"],
        "Aspergillus alpha"
    );
    let m = &by("mixed")["qc"]["mixture"];
    assert!(flags_of(by("mixed")).contains(&"mixture_suspected".to_string()));
    assert_eq!(m["primary"], "Aspergillus");
    assert_eq!(m["second"]["taxon"], "Agaricus");
    assert_eq!(m["second"]["phylum"], "Basidiomycota");
    let e = m["explained_fraction"].as_f64().unwrap();
    let s = m["second"]["explained_fraction"].as_f64().unwrap();
    assert!((e - 0.8).abs() < 0.03 && (s - 0.2).abs() < 0.03, "{e} {s}");
    for id in ["pure", "samegenus"] {
        assert!(
            !flags_of(by(id)).contains(&"mixture_suspected".to_string()),
            "{id}"
        );
        assert!(by(id)["qc"]["mixture"].get("second").is_none(), "{id}");
    }
    assert!(!flags_of(by("tiny")).contains(&"mixture_suspected".to_string()));

    // the TSV columns
    let (h, rows) = read_tsv(&d.join("out_plain/results.tsv"));
    let row = |id: &str| rows.iter().find(|r| r[col(&h, "query_id")] == id).unwrap();
    let mixed = row("mixed");
    assert_eq!(mixed[col(&h, "second_taxon")], "Agaricus");
    let ts: f64 = mixed[col(&h, "second_explained_fraction")].parse().unwrap();
    assert!((ts - 0.2).abs() < 0.03);
    assert!(mixed[col(&h, "qc_flags")].contains("mixture_suspected"));
    let pure = row("pure");
    assert_eq!(pure[col(&h, "second_taxon")], "");
    assert_eq!(pure[col(&h, "second_explained_fraction")], "");
    let pe: f64 = pure[col(&h, "explained_fraction")].parse().unwrap();
    assert!(pe > 0.99);
    assert_eq!(pure[col(&h, "schema_version")], "8");
}

#[test]
fn the_mixture_pass_respects_the_exclusion_masks() {
    let dir = tempfile::tempdir().unwrap();
    let d = dir.path();
    let (genomes, ix, _) = build_reference(d);
    let mut s = genomes[0].clone();
    s.extend_from_slice(&genomes[4][..50_000]);
    let inputs = [("mixed", fasta("q", &s))];
    let second = |recs: &[serde_json::Value]| recs[0]["qc"]["mixture"]["second"]["taxon"].clone();
    assert_eq!(second(&run_records(d, &ix, "a", &inputs, &[])), "Agaricus");
    // hiding the minority's genus (both Agaricus genomes) leaves nothing to explain it: no second taxon, no flag
    fs::write(d.join("ids.txt"), "g4\ng5\n").unwrap();
    let ids = d.join("ids.txt");
    let recs = run_records(
        d,
        &ix,
        "b",
        &inputs,
        &["--exclude-ids", ids.to_str().unwrap()],
    );
    assert!(second(&recs).is_null(), "{}", recs[0]["qc"]["mixture"]);
    assert!(!flags_of(&recs[0]).contains(&"mixture_suspected".to_string()));
    // hiding a genome that is not the minority's changes nothing
    fs::write(d.join("ids1.txt"), "g5\n").unwrap();
    let ids1 = d.join("ids1.txt");
    let recs = run_records(
        d,
        &ix,
        "c",
        &inputs,
        &["--exclude-ids", ids1.to_str().unwrap()],
    );
    assert_eq!(second(&recs), "Agaricus");
    // hiding the majority's genomes makes the minority the best hit; the genus that is gone is no second taxon
    fs::write(d.join("ids0.txt"), "g0\ng1\n").unwrap();
    let ids0 = d.join("ids0.txt");
    let recs = run_records(
        d,
        &ix,
        "e",
        &inputs,
        &["--exclude-ids", ids0.to_str().unwrap()],
    );
    assert_eq!(recs[0]["qc"]["mixture"]["primary"], "Agaricus");
    assert!(second(&recs).is_null(), "{}", recs[0]["qc"]["mixture"]);
}

#[test]
fn high_error_reads_flag_comes_from_the_kmer_survival_and_never_changes_a_call() {
    let dir = tempfile::tempdir().unwrap();
    let d = dir.path();
    let (genomes, ix, mut rng) = build_reference(d);
    let g = &genomes[2];
    let inp = d.join("reads");
    fs::create_dir_all(&inp).unwrap();
    // 10x reads with 0.5%, 2% and 8% substitution errors (per thousand: 5, 20, 80), and 0.5x reads with 8% (below the coverage the flag is judged from)
    for (name, cov, e) in [
        ("e05", 10.0, 5),
        ("e2", 10.0, 20),
        ("e8", 10.0, 80),
        ("e8_low", 0.5, 80),
    ] {
        fs::write(
            inp.join(format!("{name}.fastq")),
            fastq_reads(&mut rng, g, cov, 150, e),
        )
        .unwrap();
    }
    let out = d.join("out");
    ok(&[
        "run",
        inp.to_str().unwrap(),
        "--index",
        ix.to_str().unwrap(),
        "-o",
        out.to_str().unwrap(),
        "-q",
        "--no-report",
    ]);
    let (h, rows) = read_tsv(&out.join("results.tsv"));
    let flags = |n: &str| rows.iter().find(|r| r[0] == n).unwrap()[col(&h, "qc_flags")].clone();
    assert!(
        !flags("e05").contains("high_error_reads"),
        "{}",
        flags("e05")
    );
    assert!(!flags("e2").contains("high_error_reads"), "{}", flags("e2"));
    assert!(flags("e8").contains("high_error_reads"), "{}", flags("e8"));
    assert!(!flags("e8_low").contains("high_error_reads"));
    // the estimate is in the record: about 8% for the 8% input
    let rec: serde_json::Value = serde_json::from_str(
        fs::read_to_string(out.join("results.jsonl"))
            .unwrap()
            .lines()
            .find(|l| l.contains("\"e8\""))
            .unwrap(),
    )
    .unwrap();
    let est = rec["qc"]["error_est"].as_f64().expect("error_est");
    assert!((0.05..=0.12).contains(&est), "error estimate {est}");
}
