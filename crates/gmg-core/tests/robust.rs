//! Corrupt, truncated and hostile index / sketch files must give an error (or a usable index), never a panic or a huge allocation (#34).
use gmg_core::{guess, GenomeMeta, IndexBuilder, RefIndex, SketchFile, Thresholds};
use std::panic::catch_unwind;

struct Rng(u64);
impl Rng {
    fn next(&mut self) -> u64 {
        self.0 ^= self.0 << 13;
        self.0 ^= self.0 >> 7;
        self.0 ^= self.0 << 17;
        self.0
    }
}

fn meta(id: &str, genus: &str) -> GenomeMeta {
    let mut l: [String; 6] = Default::default();
    l[0] = "Ascomycota".into();
    l[4] = genus.into();
    GenomeMeta {
        id: id.into(),
        taxid: 7,
        lineage: l,
        n_hashes: 0,
        cluster: String::new(),
        ..Default::default()
    }
}

fn index_bytes(v: u8) -> Vec<u8> {
    let mut b = IndexBuilder::new(21, 1, 42);
    b.add(meta("a", "Aspergillus"), (0..300).map(|x| x * 7).collect());
    b.add(meta("b", "Fusarium"), (100..400).map(|x| x * 7).collect());
    b.add(meta("c", "Fusarium"), (350..600).map(|x| x * 7).collect());
    let ix = b.build(4);
    let mut buf = Vec::new();
    match v {
        1 => ix.write_v1(&mut buf).unwrap(),
        3 => ix.write_v3(&mut buf).unwrap(),
        _ => ix.write(&mut buf).unwrap(),
    }
    buf
}

/// Loading must not panic; if it loads, a search and a guess must not panic either.
fn exercise(bytes: &[u8]) {
    let r = catch_unwind(|| {
        let q: Vec<u64> = (0..50).map(|x| x * 7).collect();
        if let Ok(ix) = RefIndex::read(&mut &bytes[..]) {
            let _ = guess(&ix, &q, &|_| false, &Thresholds::default());
        }
        // the owned loader keeps v3 packed and decodes on demand: a corrupt stream must end the walk, not panic
        if let Ok(ix) = RefIndex::read_owned(bytes.to_vec()) {
            let _ = guess(&ix, &q, &|_| false, &Thresholds::default());
        }
    });
    assert!(r.is_ok(), "panicked on {} bytes", bytes.len());
}

#[test]
fn mutated_and_truncated_indexes_never_panic() {
    let mut rng = Rng(0x2545F4914F6CDD1D);
    for v in [0u8, 1, 3] {
        let good = index_bytes(v);
        assert!(RefIndex::read(&mut &good[..]).is_ok());
        // every truncation point
        for n in 0..good.len() {
            exercise(&good[..n]);
        }
        // random byte changes, single and several
        for _ in 0..4000 {
            let mut b = good.clone();
            for _ in 0..1 + rng.next() % 4 {
                let i = (rng.next() as usize) % b.len();
                b[i] = rng.next() as u8;
            }
            exercise(&b);
        }
        // a hostile header: every 8-byte window of the header region set to a huge value
        for i in 0..64.min(good.len().saturating_sub(8)) {
            let mut b = good.clone();
            b[i..i + 8].copy_from_slice(&u64::MAX.to_le_bytes());
            exercise(&b);
        }
    }
}

#[test]
fn hostile_lengths_do_not_allocate() {
    // magic + k + scaled + seed + n_genomes=1 + meta_len = 2^62: must fail fast, not try to allocate
    let mut b = b"GMGIDX01".to_vec();
    b.extend(21u32.to_le_bytes());
    b.extend(1000u32.to_le_bytes());
    b.extend(42u64.to_le_bytes());
    b.extend(1u32.to_le_bytes());
    b.extend((1u64 << 62).to_le_bytes());
    assert!(RefIndex::read(&mut &b[..]).is_err());
    // sketch file with a name length of 4 GB
    let mut s = b"GMGSKT01".to_vec();
    s.extend(21u32.to_le_bytes());
    s.extend(1000u32.to_le_bytes());
    s.extend(42u64.to_le_bytes());
    s.extend(u32::MAX.to_le_bytes());
    assert!(SketchFile::read(&mut &s[..]).is_err());
    // k = 0 and scaled = 0 are rejected
    let mut z = b"GMGIDX01".to_vec();
    z.extend(0u32.to_le_bytes());
    z.extend(0u32.to_le_bytes());
    z.extend([0u8; 40]);
    assert!(RefIndex::read(&mut &z[..]).is_err());
    assert!(RefIndex::read(&mut &b""[..]).is_err());
    assert!(RefIndex::read(&mut &b"not an index at all"[..]).is_err());
}

#[test]
fn mutated_sketch_files_never_panic() {
    let mut rng = Rng(77);
    let sk = SketchFile {
        ksize: 21,
        scaled: 1000,
        seed: 42,
        name: "x".into(),
        hashes: vec![1, 5, 9, 30],
        abunds: Some(vec![3, 1, 2, 9]),
    };
    let mut good = Vec::new();
    sk.write(&mut good).unwrap();
    for n in 0..good.len() {
        assert!(catch_unwind(|| {
            let _ = SketchFile::read(&mut &good[..n]);
        })
        .is_ok());
    }
    for _ in 0..3000 {
        let mut b = good.clone();
        let i = (rng.next() as usize) % b.len();
        b[i] = rng.next() as u8;
        assert!(catch_unwind(|| {
            let _ = SketchFile::read(&mut &b[..]);
        })
        .is_ok());
    }
}

/// FASTA/FASTQ text from anywhere (empty, binary, CRLF, truncated records, one huge line, lone headers, NULs) must sketch or
/// return nothing, never panic.
#[test]
fn hostile_sequence_text_never_panics() {
    let mut rng = Rng(0x9E3779B97F4A7C15);
    let mut cases: Vec<Vec<u8>> = vec![
        b"".to_vec(),
        b">".to_vec(),
        b">\n".to_vec(),
        b"@".to_vec(),
        b"@r\nACGT\n+\n".to_vec(),
        b"@r\nACGT\n+\nIII".to_vec(),
        b">a\r\nACGTACGTACGTACGTACGTACGTAC\r\n>b\r\n\r\n".to_vec(),
        b"ACGTACGTACGTACGTACGTACGTACGTACGT".to_vec(),
        vec![0u8; 4096],
        vec![b'>'; 10_000],
        b">x\n\0\0ACGT\xff\xfe\n".to_vec(),
    ];
    let mut huge = b">h\n".to_vec();
    huge.extend((0..3_000_000).map(|i| b"ACGT"[i % 4]));
    cases.push(huge);
    for len in [1usize, 7, 64, 1000, 20_000] {
        cases.push((0..len).map(|_| rng.next() as u8).collect());
        let mut fa = b">r\n".to_vec();
        fa.extend((0..len).map(|_| b"ACGTN\n\r>@"[(rng.next() % 9) as usize]));
        cases.push(fa);
    }
    for (i, c) in cases.iter().enumerate() {
        let r = catch_unwind(|| {
            let _ = gmg_core::sketch_fastx(c, 21, 10, 42, 1, 0);
            let _ = gmg_core::sketch_fastx_stats(c, 21, 10, false, 0);
            let _ = gmg_core::sketch_fastx_stats(c, 21, 10, true, 5);
        });
        assert!(r.is_ok(), "panicked on case {i} ({} bytes)", c.len());
    }
}
