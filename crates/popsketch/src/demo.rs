//! A small simulated dataset with planted truth, used by the tests and by `popsketch demo`: three groups of clean samples (a
//! population with a known structure) plus planted problems (low coverage, a foreign genome, half the content, a second genome
//! mixed in) and one swapped label. Deterministic (fixed seeds); nothing is read from outside.

use anyhow::{Context, Result};
use std::collections::BTreeSet;
use std::fs;
use std::path::{Path, PathBuf};

pub struct Rng(pub u64);
impl Rng {
    pub fn next_u64(&mut self) -> u64 {
        self.0 ^= self.0 << 13;
        self.0 ^= self.0 >> 7;
        self.0 ^= self.0 << 17;
        self.0
    }
    pub fn below(&mut self, n: usize) -> usize {
        (self.next_u64() % n as u64) as usize
    }
    pub fn f(&mut self) -> f64 {
        (self.next_u64() >> 11) as f64 / (1u64 << 53) as f64
    }
}

pub fn random_genome(len: usize, seed: u64) -> Vec<u8> {
    let mut r = Rng(seed);
    (0..len).map(|_| b"ACGT"[r.below(4)]).collect()
}

/// Substitute a fraction `rate` of the positions with a different base.
pub fn mutate(g: &[u8], rate: f64, seed: u64) -> Vec<u8> {
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
pub fn reads(g: &[u8], from: usize, to: usize, cov: f64, seed: u64, out: &mut String, tag: &str) {
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

/// What was planted.
#[derive(Clone, Debug)]
pub struct Demo {
    pub reads_dir: PathBuf,
    pub labels: PathBuf,
    /// samples with a planted problem: they should be flagged, and nothing else
    pub problems: BTreeSet<String>,
    /// a clean sample whose label in the label table is wrong (it says group1; the reads are group0)
    pub swapped: String,
    /// (sample, true group 0..3) of the clean samples
    pub group_of: Vec<(String, usize)>,
}

/// Write the dataset under `dir` (`dir/reads/*.fastq`, `dir/labels.tsv`): `genome_len` bases per genome and `coverage` for the
/// clean samples. 27 clean samples in 3 groups and 5 samples with planted problems.
pub fn generate(dir: &Path, genome_len: usize, coverage: f64) -> Result<Demo> {
    let reads_dir = dir.join("reads");
    fs::create_dir_all(&reads_dir).with_context(|| reads_dir.display().to_string())?;
    let (glen, cov) = (genome_len, coverage);
    let root = random_genome(glen, 1);
    let foreign = random_genome(glen, 99);
    let foreign2 = random_genome(glen, 77);
    let groups: Vec<Vec<u8>> = (0..3).map(|g| mutate(&root, 0.015, 1000 + g)).collect();
    let mut labels = String::from("id\tgroup\tsource\n");
    let mut problems = BTreeSet::new();
    let mut group_of = Vec::new();
    let write = |id: &str, text: &str| -> Result<()> {
        Ok(fs::write(reads_dir.join(format!("{id}.fastq")), text)?)
    };
    for (g, gg) in groups.iter().enumerate() {
        for s in 0..9 {
            let id = format!("g{g}s{s}");
            let genome = mutate(gg, 0.0015, 5000 + (g * 100 + s) as u64);
            let mut t = String::new();
            reads(
                &genome,
                0,
                glen,
                cov,
                20_000 + (g * 100 + s) as u64,
                &mut t,
                &id,
            );
            write(&id, &t)?;
            // sample g0s0 is labelled as group 1 (a swap)
            let lab = if id == "g0s0" { 1 } else { g };
            labels.push_str(&format!(
                "{id}\tgroup{lab}\t{}\n",
                if s % 3 == 0 { "NA" } else { "field" }
            ));
            group_of.push((id, g));
        }
    }
    for (i, grp) in groups.iter().take(2).enumerate() {
        let id = format!("lowcov{i}");
        let mut t = String::new();
        reads(
            &mutate(grp, 0.0015, 7000 + i as u64),
            0,
            glen,
            0.05 * cov,
            30_000 + i as u64,
            &mut t,
            &id,
        );
        write(&id, &t)?;
        problems.insert(id.clone());
        labels.push_str(&format!("{id}\tgroup{i}\tNA\n"));
    }
    let mut t = String::new();
    reads(&foreign, 0, glen, cov, 31_000, &mut t, "foreign");
    write("foreign", &t)?;
    problems.insert("foreign".to_string());
    labels.push_str("foreign\t\tNA\n");
    let mut t = String::new();
    reads(
        &mutate(&groups[2], 0.0015, 7100),
        0,
        glen / 2,
        cov,
        32_000,
        &mut t,
        "half",
    ); // only half the genome
    write("half", &t)?;
    problems.insert("half".to_string());
    labels.push_str("half\tgroup2\tNA\n");
    let mut t = String::new();
    reads(
        &mutate(&groups[1], 0.0015, 7200),
        0,
        glen,
        cov,
        33_000,
        &mut t,
        "mixh",
    );
    reads(&foreign2, 0, glen, cov * 0.8, 34_000, &mut t, "mixc"); // a second genome at 0.8x the host coverage
    write("mixed", &t)?;
    problems.insert("mixed".to_string());
    labels.push_str("mixed\tgroup1\tNA\n");
    let lp = dir.join("labels.tsv");
    fs::write(&lp, labels)?;
    Ok(Demo {
        reads_dir,
        labels: lp,
        problems,
        swapped: "g0s0".into(),
        group_of,
    })
}
