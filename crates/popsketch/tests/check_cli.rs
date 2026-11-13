//! `popsketch check` as a program: the exit code is the contract (0 ok, 75 undersampled, 2 error), so it is tested as one.
//! A sketch is made with `popsketch sketch` from a tiny simulated genome, so no external data is needed.

use std::path::{Path, PathBuf};
use std::process::{Command, Output};

const BIN: &str = env!("CARGO_BIN_EXE_popsketch");

fn run(args: &[&str]) -> Output {
    Command::new(BIN).args(args).output().unwrap()
}

fn code(o: &Output) -> i32 {
    o.status.code().unwrap()
}

/// a folder with one sketched sample `s1` (a 120 kb genome at 12x) and the genome size in Mb
fn folder(name: &str) -> (PathBuf, f64) {
    let dir = std::env::temp_dir().join(format!("popsketch_chk_{name}_{}", std::process::id()));
    let _ = std::fs::remove_dir_all(&dir);
    let d = popsketch::demo::generate(&dir, 120_000, 12.0).unwrap();
    let reads = dir.join("reads");
    let mut files: Vec<_> = std::fs::read_dir(&reads)
        .unwrap()
        .map(|e| e.unwrap().path())
        .filter(|p| p.to_string_lossy().contains(&d.group_of[0].0))
        .collect();
    files.sort();
    let sk = dir.join("sketches");
    let mut args = vec![
        "sketch".to_string(),
        "--id".into(),
        "s1".into(),
        "-o".into(),
        sk.to_string_lossy().into(),
        "--scaled".into(),
        "5".into(),
        "--max-bases".into(),
        "0".into(),
    ];
    args.extend(files.iter().map(|f| f.to_string_lossy().to_string()));
    let o = run(&args.iter().map(String::as_str).collect::<Vec<_>>());
    assert_eq!(code(&o), 0, "{}", String::from_utf8_lossy(&o.stderr));
    (sk, 0.12)
}

fn check(sk: &Path, extra: &[&str]) -> Output {
    let mut a = vec!["check", "--id", "s1", "-s", sk.to_str().unwrap()];
    a.extend_from_slice(extra);
    run(&a)
}

#[test]
fn a_sound_sample_is_exit_0_and_a_too_big_genome_is_exit_75() {
    let (sk, g) = folder("ok");
    assert_eq!(code(&check(&sk, &["--genome-mb", &g.to_string()])), 0);
    // the same sample against a genome 10 times larger: most of it is missing
    assert_eq!(
        code(&check(&sk, &["--genome-mb", &(g * 10.0).to_string()])),
        75
    );
}

#[test]
fn bad_numbers_are_exit_2_and_never_ok() {
    let (sk, _) = folder("bad");
    for a in [
        vec!["--genome-mb", "0"],
        vec!["--genome-mb", "NaN"],
        vec!["--genome-mb", "inf"],
        vec!["--genome-mb", "-3"],
        vec!["--genome-mb", "0.12", "--coverage", "0"],
        vec!["--genome-mb", "0.12", "--coverage", "NaN"],
        vec!["--genome-mb", "0.12", "--min-size-ratio", "2"],
        vec!["--genome-mb", "0.12", "--max-cov-ratio", "0"],
    ] {
        let o = check(&sk, &a);
        assert_eq!(code(&o), 2, "{a:?}: {}", String::from_utf8_lossy(&o.stderr));
        assert!(o.stdout.is_empty(), "{a:?} printed a verdict");
    }
}

#[test]
fn a_missing_sketch_or_folder_is_exit_2() {
    let (sk, g) = folder("missing");
    let g = g.to_string();
    let o = run(&[
        "check",
        "--id",
        "nope",
        "-s",
        sk.to_str().unwrap(),
        "--genome-mb",
        &g,
    ]);
    assert_eq!(code(&o), 2);
    let o = run(&[
        "check",
        "--id",
        "s1",
        "-s",
        "/nonexistent/folder",
        "--genome-mb",
        &g,
    ]);
    assert_eq!(code(&o), 2);
}

#[test]
fn a_mistyped_kind_is_refused_not_treated_as_an_assembly() {
    let (sk, _) = folder("kind");
    let f = sk.parent().unwrap().join("reads");
    let any = std::fs::read_dir(&f)
        .unwrap()
        .next()
        .unwrap()
        .unwrap()
        .path();
    let o = run(&[
        "sketch",
        "--id",
        "x",
        "-o",
        sk.to_str().unwrap(),
        "--kind",
        "read",
        any.to_str().unwrap(),
    ]);
    assert_ne!(code(&o), 0);
    assert!(
        String::from_utf8_lossy(&o.stderr).contains("invalid value"),
        "{}",
        String::from_utf8_lossy(&o.stderr)
    );
}

#[test]
fn a_second_sketch_with_other_parameters_is_refused_and_does_not_change_the_folder() {
    let (sk, g) = folder("params");
    let before = std::fs::read_to_string(sk.join("params.txt")).unwrap();
    let any = std::fs::read_dir(sk.parent().unwrap().join("reads"))
        .unwrap()
        .next()
        .unwrap()
        .unwrap()
        .path();
    let o = run(&[
        "sketch",
        "--id",
        "other",
        "-o",
        sk.to_str().unwrap(),
        "--scaled",
        "50",
        "--max-bases",
        "0",
        any.to_str().unwrap(),
    ]);
    assert_ne!(code(&o), 0);
    assert!(
        String::from_utf8_lossy(&o.stderr).contains("use a new folder"),
        "{}",
        String::from_utf8_lossy(&o.stderr)
    );
    assert_eq!(
        std::fs::read_to_string(sk.join("params.txt")).unwrap(),
        before
    );
    assert!(!sk.join("other.gsk").exists());
    // and the first sample still checks as it did
    assert_eq!(code(&check(&sk, &["--genome-mb", &g.to_string()])), 0);
}
