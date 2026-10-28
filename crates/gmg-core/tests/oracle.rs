//! gmg-core must produce exactly the hashes of the `sourmash` crate (k=21, scaled=1000).
//! Uses the pilot genomes in `data/genomes/`; skipped if they are not present.

use gmg_core::{FracSketch, DEFAULT_SEED};
use needletail::parse_fastx_file;
use sourmash::encodings::HashFunctions;
use sourmash::signature::SigsTrait;
use sourmash::sketch::minhash::KmerMinHash;
use std::path::PathBuf;

fn genome(acc: &str) -> Option<PathBuf> {
    let p = PathBuf::from(env!("CARGO_MANIFEST_DIR")).join(format!("../../data/genomes/{acc}.fna"));
    p.exists().then_some(p)
}

#[test]
fn matches_sourmash_crate() {
    for acc in ["GCA_000150145.1", "GCA_000002655.1"] {
        let Some(path) = genome(acc) else {
            eprintln!("skipping {acc}: data/genomes not present");
            continue;
        };
        let mut ours = FracSketch::new(21, 1000, DEFAULT_SEED);
        let mut theirs =
            KmerMinHash::new(1000, 21, HashFunctions::Murmur64Dna, DEFAULT_SEED, false, 0);
        let mut r = parse_fastx_file(&path).unwrap();
        while let Some(rec) = r.next() {
            let seq = rec.unwrap().seq().to_vec();
            ours.add_sequence(&seq);
            theirs.add_sequence(&seq, true).unwrap();
        }
        let mut want = theirs.mins();
        want.sort_unstable();
        assert_eq!(ours.finish(), want, "hash mismatch for {acc}");
    }
}

/// Deterministic synthetic genome (no data files needed, so this always runs in CI): random ACGT contigs with
/// lowercase stretches, N runs, IUPAC codes and a contig shorter than k, to exercise the cases real assemblies have.
fn synthetic_contigs(seed: u64, n_contigs: usize, len: usize) -> Vec<Vec<u8>> {
    let mut s = seed | 1;
    let mut next = move || {
        s ^= s << 13;
        s ^= s >> 7;
        s ^= s << 17;
        s
    };
    let mut out = Vec::new();
    for c in 0..n_contigs {
        let mut v: Vec<u8> = (0..len).map(|_| b"ACGT"[(next() % 4) as usize]).collect();
        for _ in 0..3 {
            let a = (next() as usize) % (len - 600);
            for b in v[a..a + 200].iter_mut() {
                *b = b.to_ascii_lowercase();
            }
            let a = (next() as usize) % (len - 600);
            for b in v[a..a + 40].iter_mut() {
                *b = b'N';
            }
        }
        v[len / 2] = b'R';
        v[len / 3] = b'Y';
        if c == 0 {
            // a repeat and its reverse complement, to hit the canonical k-mer logic
            let rep: Vec<u8> = v[1000..1300].to_vec();
            let rc: Vec<u8> = rep
                .iter()
                .rev()
                .map(|&b| match b.to_ascii_uppercase() {
                    b'A' => b'T',
                    b'C' => b'G',
                    b'G' => b'C',
                    b'T' => b'A',
                    x => x,
                })
                .collect();
            v[5000..5300].copy_from_slice(&rc);
        }
        out.push(v);
    }
    out.push(b"ACGTACGTAC".to_vec()); // shorter than k=21
    out
}

#[test]
fn synthetic_matches_sourmash_crate() {
    for (k, scaled, seed) in [(21u32, 1000u32, 7u64), (21, 10, 11), (31, 100, 13)] {
        let contigs = synthetic_contigs(seed, 4, 120_000);
        let mut ours = FracSketch::new(k as usize, scaled, DEFAULT_SEED);
        let mut theirs = KmerMinHash::new(
            scaled,
            k,
            HashFunctions::Murmur64Dna,
            DEFAULT_SEED,
            false,
            0,
        );
        for c in &contigs {
            ours.add_sequence(c);
            theirs.add_sequence(c, true).unwrap();
        }
        let mut want = theirs.mins();
        want.sort_unstable();
        let got = ours.finish();
        assert!(!got.is_empty());
        assert_eq!(got, want, "hash mismatch for k={k} scaled={scaled}");
    }
}
