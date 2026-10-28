//! Small helpers shared by the chart and page code.

/// HTML/XML escape (also for attribute values).
pub fn esc(s: &str) -> String {
    let mut o = String::with_capacity(s.len() + 8);
    for c in s.chars() {
        match c {
            '&' => o.push_str("&amp;"),
            '<' => o.push_str("&lt;"),
            '>' => o.push_str("&gt;"),
            '"' => o.push_str("&quot;"),
            '\'' => o.push_str("&#x27;"),
            _ => o.push(c),
        }
    }
    o
}

/// Truncate to `n` characters with an ellipsis (character-based, safe for non-ASCII).
pub fn short(s: &str, n: usize) -> String {
    if s.chars().count() <= n {
        s.to_string()
    } else {
        let mut t: String = s.chars().take(n.saturating_sub(1)).collect();
        t.push('…');
        t
    }
}

/// `<text>` element. `extra` is raw attribute text (e.g. `style="font-weight:600"`).
pub fn text(x: f64, y: f64, s: &str, cls: &str, anchor: &str, extra: &str) -> String {
    format!(
        r#"<text x="{x:.1}" y="{y:.1}" class="{cls}" text-anchor="{anchor}" {extra}>{}</text>"#,
        esc(s)
    )
}

/// 1,284 / 12.9K / 3.4M
pub fn compact(n: f64) -> String {
    if n >= 1e6 {
        format!("{:.1}M", n / 1e6).replace(".0M", "M")
    } else if n >= 1e3 {
        format!("{:.1}K", n / 1e3).replace(".0K", "K")
    } else if n.fract() == 0.0 {
        format!("{}", n as i64)
    } else {
        let s = format!("{n:.2}");
        s.trim_end_matches('0').trim_end_matches('.').to_string()
    }
}

/// Round axis ticks (0 .. max) to 1/2/2.5/5 x 10^k steps.
pub fn nice_ticks(vmax: f64, n: usize) -> Vec<f64> {
    if vmax <= 0.0 {
        return vec![0.0, 1.0];
    }
    let raw = vmax / n as f64;
    let mag = 10f64.powf(raw.log10().floor());
    let step = [1.0, 2.0, 2.5, 5.0, 10.0]
        .iter()
        .map(|s| s * mag)
        .find(|s| *s >= raw)
        .unwrap_or(mag);
    let k = (vmax / step + 1.0001) as usize;
    (0..=k).map(|i| i as f64 * step).collect()
}

/// Rounded-end horizontal bar path (4px radius at the data end, square at the baseline).
pub fn hbar_path(x: f64, y: f64, len: f64, thick: f64) -> String {
    let len = len.max(8.0);
    let h = thick - 8.0;
    format!(
        "M{x:.1},{y:.1} h{:.1} a4,4 0 0 1 4,4 v{h:.1} a4,4 0 0 1 -4,4 h-{:.1} z",
        len - 4.0,
        len - 4.0
    )
}

/// Rounded-top column path growing up from the baseline `base`.
pub fn column_path(x: f64, base: f64, width: f64, height: f64) -> String {
    let h = (height - 4.0).max(0.0);
    format!(
        "M{x:.1},{base:.1} v-{h:.1} a4,4 0 0 1 4,-4 h{:.1} a4,4 0 0 1 4,4 v{h:.1} z",
        (width - 8.0).max(0.0)
    )
}

pub fn f(x: f64, d: usize) -> String {
    format!("{x:.d$}")
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn escapes_and_shortens() {
        assert_eq!(
            esc(r#"<a href="x">&'"#),
            "&lt;a href=&quot;x&quot;&gt;&amp;&#x27;"
        );
        assert_eq!(short("abcdef", 4), "abc…");
        assert_eq!(short("αβγδεζ", 4), "αβγ…");
        assert_eq!(short("abc", 4), "abc");
    }

    #[test]
    fn compact_and_ticks() {
        assert_eq!(compact(1500.0), "1.5K");
        assert_eq!(compact(2_000_000.0), "2M");
        assert_eq!(compact(7.0), "7");
        assert_eq!(compact(0.5), "0.5");
        let t = nice_ticks(134.0, 4);
        assert_eq!(t.first(), Some(&0.0));
        assert!(*t.last().unwrap() >= 134.0 && t.len() >= 3);
        assert_eq!(nice_ticks(0.0, 4), vec![0.0, 1.0]);
    }

    #[test]
    fn bar_paths_are_valid_even_for_tiny_values() {
        let p = hbar_path(150.0, 10.0, 0.5, 18.0);
        assert!(!p.contains("--") && !p.contains("NaN"), "{p}");
        let c = column_path(10.0, 100.0, 4.0, 1.0);
        assert!(!c.contains("--") && !c.contains("NaN"), "{c}");
    }
}
