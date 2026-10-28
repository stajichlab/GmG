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

#[test]
fn run_folder_end_to_end() {
    let dir = tempfile::tempdir().unwrap();
    let d = dir.path();
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
    let (h2, rows2) = read_tsv(&out2.join("results.tsv"));
    let m2 = rows2
        .iter()
        .find(|r| r[col(&h2, "query_id")] == "mutant")
        .unwrap();
    assert_eq!(m2[col(&h2, "call_genus")], "");
    assert_eq!(m2[col(&h2, "status")], "low_confidence");
}
