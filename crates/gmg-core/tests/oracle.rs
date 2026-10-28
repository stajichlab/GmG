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
        let mut theirs = KmerMinHash::new(1000, 21, HashFunctions::Murmur64Dna, DEFAULT_SEED, false, 0);
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
