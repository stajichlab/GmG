//! Minimal FASTA/FASTQ reader for sketching from memory (no I/O, no dependencies).
//!
//! The browser hands over whole files as bytes; the native CLI uses `needletail` (which also reads gzip) instead.
//! FASTA records may span lines; FASTQ records are the usual four lines. Format is detected from the first byte.

use crate::sketch::{FracSketch, DEFAULT_SEED};
use crate::skfile::SketchFile;

/// Per-input counts gathered while reading the sequence (used for QC flags).
#[derive(Clone, Debug, Default)]
pub struct SeqStats {
    pub n_records: u64,
    pub total_len: u64,
    /// bases other than A/C/G/T (case-insensitive)
    pub non_acgt: u64,
    pub n50: u64,
}

/// Sketch a FASTA/FASTQ held in memory and collect the stats `gmg run` records: `abund` tracks abundances (reads),
/// otherwise record lengths are kept for the N50. `max_records` stops early (0 = no limit). None if no record was read.
pub fn sketch_fastx_stats(
    bytes: &[u8],
    ksize: u32,
    scaled: u32,
    abund: bool,
    max_records: u64,
) -> Option<(SketchFile, SeqStats)> {
    let mut sk = FracSketch::new(ksize as usize, scaled, DEFAULT_SEED);
    let (mut n, mut total, mut non_acgt) = (0u64, 0u64, 0u64);
    let mut lens: Vec<u64> = Vec::new();
    let mut each = |seq: &[u8]| -> bool {
        sk.add_sequence(seq);
        n += 1;
        total += seq.len() as u64;
        non_acgt += seq
            .iter()
            .filter(|b| !matches!(b.to_ascii_uppercase(), b'A' | b'C' | b'G' | b'T'))
            .count() as u64;
        if !abund {
            lens.push(seq.len() as u64);
        }
        max_records == 0 || n < max_records
    };
    match bytes.first() {
        Some(b'@') => fastq(bytes, &mut each),
        _ => fasta(bytes, &mut each),
    }
    if n == 0 {
        return None;
    }
    lens.sort_unstable_by(|a, b| b.cmp(a));
    let (mut acc, mut n50) = (0u64, 0u64);
    for l in &lens {
        acc += l;
        if acc * 2 >= total {
            n50 = *l;
            break;
        }
    }
    let (hashes, abunds) = if abund {
        let (h, a) = sk.finish_with_abundance();
        (h, Some(a))
    } else {
        (sk.finish(), None)
    };
    Some((
        SketchFile {
            ksize,
            scaled,
            seed: DEFAULT_SEED,
            name: String::new(),
            hashes,
            abunds,
        },
        SeqStats {
            n_records: n,
            total_len: total,
            non_acgt,
            n50,
        },
    ))
}

/// Sketch every record in `bytes` and return sorted unique hashes, keeping only hashes seen at least `min_abund`
/// times (use 1 for assemblies, 2 for reads). `max_records` stops early (reads subsample); 0 = no limit.
pub fn sketch_fastx(
    bytes: &[u8],
    ksize: usize,
    scaled: u32,
    seed: u64,
    min_abund: u32,
    max_records: u64,
) -> Vec<u64> {
    let mut sk = FracSketch::new(ksize, scaled, seed);
    let mut n = 0u64;
    let mut each = |seq: &[u8]| -> bool {
        sk.add_sequence(seq);
        n += 1;
        max_records == 0 || n < max_records
    };
    match bytes.first() {
        Some(b'@') => fastq(bytes, &mut each),
        _ => fasta(bytes, &mut each),
    }
    if min_abund <= 1 {
        sk.finish()
    } else {
        let (h, a) = sk.finish_with_abundance();
        h.into_iter()
            .zip(a)
            .filter(|&(_, c)| c >= min_abund)
            .map(|(h, _)| h)
            .collect()
    }
}

fn fasta(bytes: &[u8], each: &mut dyn FnMut(&[u8]) -> bool) {
    let mut seq: Vec<u8> = Vec::new();
    let mut in_rec = false;
    for line in bytes.split(|&b| b == b'\n') {
        let line = line.strip_suffix(b"\r").unwrap_or(line);
        if line.first() == Some(&b'>') {
            if in_rec && !each(&seq) {
                return;
            }
            seq.clear();
            in_rec = true;
        } else {
            seq.extend_from_slice(line);
        }
    }
    if in_rec && !seq.is_empty() {
        each(&seq);
    }
}

fn fastq(bytes: &[u8], each: &mut dyn FnMut(&[u8]) -> bool) {
    let mut lines = bytes.split(|&b| b == b'\n');
    while let (Some(_head), Some(seq), Some(_plus), Some(_qual)) =
        (lines.next(), lines.next(), lines.next(), lines.next())
    {
        if !each(seq.strip_suffix(b"\r").unwrap_or(seq)) {
            return;
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::sketch::DEFAULT_SEED;

    fn direct(seqs: &[&[u8]]) -> Vec<u64> {
        let mut s = FracSketch::new(5, 1, DEFAULT_SEED);
        for q in seqs {
            s.add_sequence(q);
        }
        s.finish()
    }

    #[test]
    fn fasta_multiline_and_crlf() {
        let f = b">a desc\r\nACGTACG\r\nTTGCA\r\n>b\nGGGCCCAAATTT\n";
        assert_eq!(
            sketch_fastx(f, 5, 1, DEFAULT_SEED, 1, 0),
            direct(&[b"ACGTACGTTGCA", b"GGGCCCAAATTT"])
        );
    }

    #[test]
    fn fastq_records_and_limit() {
        let q = b"@r1\nACGTACGTTGCA\n+\nIIIIIIIIIIII\n@r2\nGGGCCCAAATTT\n+\nIIIIIIIIIIII\n";
        assert_eq!(
            sketch_fastx(q, 5, 1, DEFAULT_SEED, 1, 0),
            direct(&[b"ACGTACGTTGCA", b"GGGCCCAAATTT"])
        );
        assert_eq!(
            sketch_fastx(q, 5, 1, DEFAULT_SEED, 1, 1),
            direct(&[b"ACGTACGTTGCA"])
        );
    }

    #[test]
    fn min_abundance_drops_rare_hashes() {
        // read A three times, read B once: with min_abund 3 only A's hashes survive (a canonical hash can repeat
        // within one read when a window and its reverse complement both occur, but never 3 times in B here)
        let rec = |name: &str, seq: &str| format!("@{name}\n{seq}\n+\n{}\n", "I".repeat(seq.len()));
        let a = "ACGTACGTTGCA";
        let b = "GGGCCCAAATTT";
        let q = format!(
            "{}{}{}{}",
            rec("a1", a),
            rec("a2", a),
            rec("a3", a),
            rec("b1", b)
        );
        let all = sketch_fastx(q.as_bytes(), 5, 1, DEFAULT_SEED, 1, 0);
        let kept = sketch_fastx(q.as_bytes(), 5, 1, DEFAULT_SEED, 3, 0);
        assert!(kept.len() < all.len());
        assert_eq!(kept, direct(&[a.as_bytes()]));
    }
}
