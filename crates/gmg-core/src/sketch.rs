use murmurhash3::murmurhash3_x64_128;

/// sourmash's default MurmurHash3 seed.
pub const DEFAULT_SEED: u64 = 42;

/// Largest hash kept for a given `scaled` (same arithmetic as sourmash).
pub fn max_hash_for_scaled(scaled: u32) -> u64 {
    match scaled {
        0 => 0,
        1 => u64::MAX,
        _ => (u64::MAX as f64 / scaled as f64) as u64,
    }
}

/// A scaled (FracMinHash) DNA sketch under construction. Feed sequences with
/// [`add_sequence`](Self::add_sequence), then call [`finish`](Self::finish).
pub struct FracSketch {
    ksize: usize,
    seed: u64,
    max_hash: u64,
    hashes: Vec<u64>,
}

impl FracSketch {
    pub fn new(ksize: usize, scaled: u32, seed: u64) -> Self {
        assert!(ksize > 0 && scaled > 0);
        Self { ksize, seed, max_hash: max_hash_for_scaled(scaled), hashes: Vec::new() }
    }

    /// Add one sequence (one FASTA record). Case-insensitive. k-mers containing anything but
    /// ACGT are skipped, like sourmash with `force`.
    pub fn add_sequence(&mut self, seq: &[u8]) {
        let k = self.ksize;
        if seq.len() < k {
            return;
        }
        let fw: Vec<u8> = seq.iter().map(|b| b.to_ascii_uppercase()).collect();
        let rc: Vec<u8> = fw.iter().rev().map(|&b| complement(b)).collect();
        let mut run = 0usize; // consecutive valid bases ending at i
        for (i, &b) in fw.iter().enumerate() {
            if matches!(b, b'A' | b'C' | b'G' | b'T') {
                run += 1;
            } else {
                run = 0;
            }
            if run >= k {
                let start = i + 1 - k;
                let kmer = &fw[start..=i];
                let krc = &rc[fw.len() - 1 - i..fw.len() - start];
                let h = murmurhash3_x64_128(kmer.min(krc), self.seed).0;
                if h <= self.max_hash {
                    self.hashes.push(h);
                }
            }
        }
    }

    /// Sorted, de-duplicated hashes.
    pub fn finish(mut self) -> Vec<u64> {
        self.hashes.sort_unstable();
        self.hashes.dedup();
        self.hashes
    }

    /// Sorted unique hashes with the number of times each was seen (for reads).
    pub fn finish_with_abundance(mut self) -> (Vec<u64>, Vec<u32>) {
        self.hashes.sort_unstable();
        let mut hs: Vec<u64> = Vec::new();
        let mut ab: Vec<u32> = Vec::new();
        for h in self.hashes {
            if hs.last() == Some(&h) {
                *ab.last_mut().unwrap() += 1;
            } else {
                hs.push(h);
                ab.push(1);
            }
        }
        (hs, ab)
    }
}

/// Reduce a sorted sketch made at scale `from` to the coarser scale `to` (to >= from), as
/// sourmash's `downsample` does: keep hashes <= max_hash(to).
pub fn downsample(hashes: &[u64], to_scaled: u32) -> Vec<u64> {
    let mh = max_hash_for_scaled(to_scaled);
    hashes[..hashes.partition_point(|&h| h <= mh)].to_vec()
}

fn complement(b: u8) -> u8 {
    match b {
        b'A' => b'T',
        b'C' => b'G',
        b'G' => b'C',
        b'T' => b'A',
        _ => b'N',
    }
}

fn intersection_size(a: &[u64], b: &[u64]) -> usize {
    let (mut i, mut j, mut n) = (0, 0, 0);
    while i < a.len() && j < b.len() {
        match a[i].cmp(&b[j]) {
            std::cmp::Ordering::Less => i += 1,
            std::cmp::Ordering::Greater => j += 1,
            std::cmp::Ordering::Equal => {
                n += 1;
                i += 1;
                j += 1;
            }
        }
    }
    n
}

/// Fraction of `query` hashes found in `reference` (both sorted, de-duplicated).
pub fn containment(query: &[u64], reference: &[u64]) -> f64 {
    if query.is_empty() {
        return 0.0;
    }
    intersection_size(query, reference) as f64 / query.len() as f64
}

/// Jaccard similarity of two sorted, de-duplicated hash lists.
pub fn jaccard(a: &[u64], b: &[u64]) -> f64 {
    let inter = intersection_size(a, b);
    let union = a.len() + b.len() - inter;
    if union == 0 {
        0.0
    } else {
        inter as f64 / union as f64
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn revcomp_is_canonical() {
        let mut a = FracSketch::new(5, 1, DEFAULT_SEED);
        a.add_sequence(b"ACGTTGCAAT");
        let mut b = FracSketch::new(5, 1, DEFAULT_SEED);
        b.add_sequence(b"ATTGCAACGT"); // reverse complement of the above
        assert_eq!(a.finish(), b.finish());
    }

    #[test]
    fn skips_non_acgt_and_ignores_case() {
        let mut a = FracSketch::new(4, 1, DEFAULT_SEED);
        a.add_sequence(b"acgtNacgt");
        let mut b = FracSketch::new(4, 1, DEFAULT_SEED);
        b.add_sequence(b"ACGT");
        b.add_sequence(b"ACGT");
        assert_eq!(a.finish(), b.finish());
    }

    #[test]
    fn abundance_counts_repeats() {
        let mut a = FracSketch::new(4, 1, DEFAULT_SEED);
        a.add_sequence(b"ACGTACGT");
        a.add_sequence(b"ACGT");
        let plain = {
            let mut b = FracSketch::new(4, 1, DEFAULT_SEED);
            b.add_sequence(b"ACGTACGT");
            b.add_sequence(b"ACGT");
            b.finish()
        };
        let (h, ab) = a.finish_with_abundance();
        assert_eq!(h, plain);
        assert_eq!(ab.iter().sum::<u32>(), 6); // 5 windows in ACGTACGT + 1
    }

    #[test]
    fn downsample_keeps_low_hashes() {
        let mh = max_hash_for_scaled(2);
        let v = vec![1, mh / 2, mh, mh + 1];
        assert_eq!(downsample(&v, 4), vec![1, mh / 2]);
    }

    #[test]
    fn containment_and_jaccard() {
        let a = [1u64, 2, 3, 4];
        let b = [3u64, 4, 5];
        assert_eq!(containment(&a, &b), 0.5);
        assert!((jaccard(&a, &b) - 0.4).abs() < 1e-12);
    }
}
