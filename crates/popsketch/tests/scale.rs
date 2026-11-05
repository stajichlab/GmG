//! How the stages that grow with the number of samples behave at 1,500 and 3,000 samples (synthetic hash sets with a shared core,
//! a pool of variable hashes and a few private ones). Ignored by default; run it with
//! `cargo test -p popsketch --release --test scale -- --ignored --nocapture`.

use popsketch::cluster::{gap_rule, upgma};
use popsketch::dist::compute;
use popsketch::ordination::{hash_pca, pcoa};
use popsketch::tree::{collapse, neighbor_joining};
use std::time::Instant;

fn population(n: usize, core: u64, pool: u64, p: f64, private: u64) -> Vec<Vec<u64>> {
    let mut x: u64 = 12345;
    let mut next = move || {
        x ^= x << 13;
        x ^= x >> 7;
        x ^= x << 17;
        (x >> 11) as f64 / (1u64 << 53) as f64
    };
    (0..n as u64)
        .map(|i| {
            let mut h: Vec<u64> = (0..core).collect();
            h.extend((0..pool).filter(|_| next() < p).map(|v| 1_000_000 + v));
            h.extend((0..private).map(|v| 100_000_000 + i * 10_000 + v));
            h
        })
        .collect()
}

#[test]
#[ignore]
fn stage_times_at_1500_and_3000_samples() {
    for n in [1500usize, 3000] {
        let sets = population(n, 20_000, 50_000, 0.3, 100);
        let t = Instant::now();
        let d = compute(&sets, 21);
        let t_dist = t.elapsed();
        let d1 = d.one_minus_ani();
        let t = Instant::now();
        let m = upgma(&d1, n);
        let t_upgma = t.elapsed();
        let _ = gap_rule(&m, n, 3.0, 0.001, 30);
        let t = Instant::now();
        let mut tree = neighbor_joining(&d.mash, n, (0..n).map(|i| i.to_string()).collect());
        collapse(&mut tree, 0.001);
        let t_nj = t.elapsed();
        let t = Instant::now();
        let pca = hash_pca(&sets, 5).unwrap();
        let t_pca = t.elapsed();
        let t = Instant::now();
        let pc = pcoa(&d1, n, 5).unwrap();
        let t_pcoa = t.elapsed();
        eprintln!(
            "n={n}: distances {t_dist:.1?}, UPGMA {t_upgma:.1?}, NJ {t_nj:.1?}, hash PCA {t_pca:.1?}, PCoA {t_pcoa:.1?} (PC1 {:.3}, {:.3})",
            pca.var_frac[0], pc.var_frac[0]
        );
    }
}
