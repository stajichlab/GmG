//! A small seeded generator (SplitMix64) shared by the permutation tests, the ordination start vectors, the cluster-count null
//! and the pair sampling of the abundance-shoulder fit. Deterministic for a seed, no dependency.

#[derive(Clone, Debug)]
pub struct SplitMix(pub u64);

impl SplitMix {
    pub fn next_u64(&mut self) -> u64 {
        self.0 = self.0.wrapping_add(0x9E37_79B9_7F4A_7C15);
        let mut z = self.0;
        z = (z ^ (z >> 30)).wrapping_mul(0xBF58_476D_1CE4_E5B9);
        z = (z ^ (z >> 27)).wrapping_mul(0x94D0_49BB_1331_11EB);
        z ^ (z >> 31)
    }

    /// Uniform in [0, 1).
    pub fn unit(&mut self) -> f64 {
        (self.next_u64() >> 11) as f64 / (1u64 << 53) as f64
    }

    /// Uniform in 0..n (n > 0); the modulo bias is below 2^-32 for any n used here.
    pub fn below(&mut self, n: usize) -> usize {
        (self.next_u64() % n as u64) as usize
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn the_sequence_is_fixed_for_a_seed_and_unit_stays_in_range() {
        // the first outputs of SplitMix64 from seed 0 (the reference values)
        let mut r = SplitMix(0);
        assert_eq!(r.next_u64(), 0xE220_A839_7B1D_CDAF);
        assert_eq!(r.next_u64(), 0x6E78_9E6A_A1B9_65F4);
        let mut a = SplitMix(7);
        let mut b = SplitMix(7);
        for _ in 0..1000 {
            let (x, y) = (a.unit(), b.unit());
            assert_eq!(x, y);
            assert!((0.0..1.0).contains(&x));
        }
        let mut c = SplitMix(1);
        assert!((0..1000).all(|_| c.below(7) < 7));
    }
}
