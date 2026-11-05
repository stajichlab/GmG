//! Small robust statistics used by the cohort-relative QC.

pub fn median(v: &[f64]) -> f64 {
    if v.is_empty() {
        return f64::NAN;
    }
    let mut s: Vec<f64> = v.iter().copied().filter(|x| x.is_finite()).collect();
    if s.is_empty() {
        return f64::NAN;
    }
    s.sort_by(|a, b| a.partial_cmp(b).unwrap());
    let n = s.len();
    if n % 2 == 1 {
        s[n / 2]
    } else {
        (s[n / 2 - 1] + s[n / 2]) / 2.0
    }
}

/// Linear-interpolated quantile, `q` in 0..=1.
pub fn quantile(v: &[f64], q: f64) -> f64 {
    let mut s: Vec<f64> = v.iter().copied().filter(|x| x.is_finite()).collect();
    if s.is_empty() {
        return f64::NAN;
    }
    s.sort_by(|a, b| a.partial_cmp(b).unwrap());
    let pos = q.clamp(0.0, 1.0) * (s.len() - 1) as f64;
    let (lo, hi) = (pos.floor() as usize, pos.ceil() as usize);
    s[lo] + (s[hi] - s[lo]) * (pos - lo as f64)
}

/// Median absolute deviation scaled to a standard deviation (1.4826 * MAD).
pub fn mad_sd(v: &[f64]) -> f64 {
    let m = median(v);
    let dev: Vec<f64> = v.iter().map(|x| (x - m).abs()).collect();
    1.4826 * median(&dev)
}

/// Robust z of `x` against `v`: (x - median) / max(1.4826 MAD, `floor`). `floor` keeps a very tight cohort from
/// flagging everything.
pub fn robust_z(x: f64, v: &[f64], floor: f64) -> f64 {
    (x - median(v)) / mad_sd(v).max(floor)
}

/// Adjusted Rand index of two labelings of the same items.
pub fn adjusted_rand(a: &[usize], b: &[usize]) -> f64 {
    assert_eq!(a.len(), b.len());
    let n = a.len();
    if n < 2 {
        return f64::NAN;
    }
    let (ka, kb) = (
        a.iter().copied().max().unwrap() + 1,
        b.iter().copied().max().unwrap() + 1,
    );
    let mut t = vec![0u64; ka * kb];
    let (mut ra, mut rb) = (vec![0u64; ka], vec![0u64; kb]);
    for (&x, &y) in a.iter().zip(b) {
        t[x * kb + y] += 1;
        ra[x] += 1;
        rb[y] += 1;
    }
    let c2 = |x: u64| (x * x.saturating_sub(1) / 2) as f64;
    let sum_t: f64 = t.iter().map(|&x| c2(x)).sum();
    let sum_a: f64 = ra.iter().map(|&x| c2(x)).sum();
    let sum_b: f64 = rb.iter().map(|&x| c2(x)).sum();
    let total = c2(n as u64);
    let expected = sum_a * sum_b / total;
    let max_index = (sum_a + sum_b) / 2.0;
    if (max_index - expected).abs() < 1e-12 {
        return 1.0;
    }
    (sum_t - expected) / (max_index - expected)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn median_and_quantile() {
        assert_eq!(median(&[3.0, 1.0, 2.0]), 2.0);
        assert_eq!(median(&[4.0, 1.0, 2.0, 3.0]), 2.5);
        assert!((quantile(&[0.0, 10.0], 0.9) - 9.0).abs() < 1e-12);
        assert!(median(&[]).is_nan());
    }

    #[test]
    fn mad_ignores_an_outlier() {
        let v = [10.0, 10.1, 9.9, 10.0, 10.2, 9.8, 100.0];
        assert!(mad_sd(&v) < 0.5);
        assert!(robust_z(100.0, &v, 0.01) > 50.0);
        assert!(robust_z(10.0, &v, 0.01).abs() < 0.1);
    }

    #[test]
    fn floor_keeps_a_tight_cohort_from_flagging_everything() {
        let v = [5.0, 5.0, 5.0, 5.0, 5.0];
        assert!(robust_z(5.01, &v, 0.5).abs() < 0.1);
    }

    #[test]
    fn ari_known_values() {
        assert!((adjusted_rand(&[0, 0, 1, 1], &[1, 1, 0, 0]) - 1.0).abs() < 1e-12);
        let r = adjusted_rand(&[0, 0, 0, 1, 1, 1], &[0, 1, 0, 1, 0, 1]);
        assert!(r < 0.2, "{r}");
        // sklearn: adjusted_rand_score([0,0,1,2],[0,0,1,1]) = 0.571428...
        assert!((adjusted_rand(&[0, 0, 1, 2], &[0, 0, 1, 1]) - 4.0 / 7.0).abs() < 1e-9);
    }
}
