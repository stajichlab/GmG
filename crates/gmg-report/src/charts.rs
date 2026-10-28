//! Inline-SVG charts. Every chart returns the drawing plus a table twin with the same numbers.
//! Rules (dataviz): thin marks with a rounded data end, 2px surface gap between segments, hairline solid grid,
//! text in text tokens never the data colour, values labelled selectively, a tooltip on every mark.

use crate::model::{Hist, Rec};
use crate::svg::{column_path, compact, esc, hbar_path, nice_ticks, short, text};
use std::collections::BTreeMap;

pub struct Chart {
    pub svg: String,
    pub table: String,
}

// ---------------------------------------------------------------- phylum colours
/// Phylum colours follow the entity (the phylum), never its rank in a run: the 8 phyla with the most reference genomes
/// take the 8 validated categorical slots in order; every other phylum shares one neutral grey.
pub const PHYLUM_SLOTS: [&str; 8] = [
    "Ascomycota",
    "Basidiomycota",
    "Mucoromycota",
    "Kickxellomycota",
    "Microsporidia",
    "Mortierellomycota",
    "Chytridiomycota",
    "Glomeromycota",
];

pub fn ph_color(name: &str) -> String {
    match PHYLUM_SLOTS.iter().position(|p| *p == name) {
        Some(i) => format!("var(--ph-{})", i + 1),
        None if name.is_empty() || name == "(not called)" => "var(--nomatch)".into(),
        None => "var(--ph-other)".into(),
    }
}

pub fn ph_order(name: &str) -> usize {
    match PHYLUM_SLOTS.iter().position(|p| *p == name) {
        Some(i) => i,
        None if name.is_empty() || name == "(not called)" => PHYLUM_SLOTS.len() + 1,
        None => PHYLUM_SLOTS.len(),
    }
}

pub fn ph_dot(name: &str) -> String {
    format!(
        r#"<span class="dot" style="background:{}" title="{}"></span>"#,
        ph_color(name),
        esc(if name.is_empty() {
            "phylum not called"
        } else {
            name
        })
    )
}

// ---------------------------------------------------------------- depth bar
#[derive(Default, Clone, Debug)]
pub struct Depth {
    pub species: usize,
    pub genus: usize,
    pub family: usize,
    pub higher: usize,
    pub none: usize,
}

impl Depth {
    pub fn total(&self) -> usize {
        self.species + self.genus + self.family + self.higher + self.none
    }
    pub fn from_recs(recs: &[Rec]) -> Depth {
        let mut d = Depth::default();
        for r in recs.iter().filter(|r| r.status != "error") {
            match r.best_rank {
                Some(5) => d.species += 1,
                Some(4) => d.genus += 1,
                Some(3) => d.family += 1,
                Some(_) => d.higher += 1,
                None => d.none += 1,
            }
        }
        d
    }
}

pub fn depth_chart(d: &Depth) -> Chart {
    let segs_def: [(&str, &str, &str, &str, usize); 5] = [
        ("species", "Species", "var(--ord-1)", "#ffffff", d.species),
        ("genus", "Genus", "var(--ord-2)", "#ffffff", d.genus),
        ("family", "Family", "var(--ord-3)", "#ffffff", d.family),
        (
            "higher",
            "Order or higher",
            "var(--ord-4)",
            "#0b0b0b",
            d.higher,
        ),
        (
            "none",
            "No match",
            "var(--nomatch)",
            "var(--text-primary)",
            d.none,
        ),
    ];
    let total = d.total().max(1) as f64;
    let (w_all, y, h, gap) = (880.0, 14.0, 30.0, 2.0);
    let mut out = String::new();
    let mut x = 0.0;
    let mut segs = Vec::new();
    for (_, label, fill, ink, n) in segs_def {
        if n == 0 {
            continue;
        }
        let w = w_all * n as f64 / total;
        let tip = format!(
            "{label}: {n} of {} inputs ({:.0}%)",
            d.total(),
            100.0 * n as f64 / total
        );
        out.push_str(&format!(
            r#"<rect class="tip" tabindex="0" data-tip="{}" x="{x:.1}" y="{y}" width="{:.1}" height="{h}" rx="4" fill="{fill}"/>"#,
            esc(&tip),
            (w - gap).max(1.0)
        ));
        if w > 70.0 {
            // label only when it fits with padding
            out.push_str(&text(
                x + w / 2.0 - gap / 2.0,
                y + h / 2.0 + 4.5,
                &n.to_string(),
                "t-fill",
                "middle",
                &format!(r#"style="fill:{ink}""#),
            ));
        }
        segs.push((label, n, fill));
        x += w;
    }
    let mut lx = 0.0;
    for (label, n, fill) in &segs {
        out.push_str(&format!(
            r#"<rect x="{lx:.1}" y="{}" width="10" height="10" rx="2" fill="{fill}"/>"#,
            y + h + 16.0
        ));
        let txt = format!("{label} {n}");
        out.push_str(&text(lx + 15.0, y + h + 25.5, &txt, "t-sec", "start", ""));
        lx += 15.0 + txt.chars().count() as f64 * 6.6 + 18.0;
    }
    let mut table = String::from("<table><thead><tr><th>Deepest rank called</th><th>Inputs</th><th>Share</th></tr></thead><tbody>");
    for (label, n, _) in &segs {
        table.push_str(&format!(
            "<tr><td>{}</td><td>{n}</td><td>{:.1}%</td></tr>",
            esc(label),
            100.0 * *n as f64 / total
        ));
    }
    table.push_str("</tbody></table>");
    Chart {
        svg: format!(
            r#"<svg viewBox="0 0 880 84" class="chart" role="img" aria-label="Deepest rank called per input">{out}</svg>"#
        ),
        table,
    }
}

// ---------------------------------------------------------------- single-series horizontal bars
pub fn hbars(items: &[(String, usize)], title_col: &str, colors: Option<&[String]>) -> Chart {
    if items.is_empty() {
        return Chart {
            svg: String::new(),
            table: String::new(),
        };
    }
    let (width, row_h, label_w) = (420.0, 26.0, 150.0);
    let vmax = items.iter().map(|(_, v)| *v).max().unwrap_or(1).max(1) as f64;
    let plot_w = width - label_w - 44.0;
    let hgt = row_h * items.len() as f64 + 6.0;
    let mut out =
        format!(r#"<line x1="{label_w}" x2="{label_w}" y1="0" y2="{hgt}" class="axis"/>"#);
    for (i, (lab, v)) in items.iter().enumerate() {
        let y = i as f64 * row_h + 4.0;
        let bw = (plot_w * *v as f64 / vmax).max(8.0);
        let tip = format!("{lab}: {v} inputs");
        let fill = colors
            .map(|c| c[i].clone())
            .unwrap_or_else(|| "var(--series-1)".into());
        out.push_str(&text(
            label_w - 8.0,
            y + 14.0,
            &short(lab, 22),
            "t-sec",
            "end",
            "",
        ));
        out.push_str(&format!(
            r#"<path class="tip bar" tabindex="0" data-tip="{}" d="{}" fill="{fill}"/>"#,
            esc(&tip),
            hbar_path(label_w, y + 3.0, bw, 18.0)
        ));
        out.push_str(&format!(
            r#"<rect class="tip hit" tabindex="-1" data-tip="{}" x="{label_w}" y="{}" width="{}" height="{}" fill="transparent"/>"#,
            esc(&tip),
            y - 1.0,
            plot_w + 40.0,
            row_h - 2.0
        ));
        out.push_str(&text(
            label_w + bw + 6.0,
            y + 14.5,
            &v.to_string(),
            "t-pri",
            "start",
            "",
        ));
    }
    let mut table = format!(
        "<table><thead><tr><th>{}</th><th>Inputs</th></tr></thead><tbody>",
        esc(title_col)
    );
    for (l, v) in items {
        table.push_str(&format!("<tr><td>{}</td><td>{v}</td></tr>", esc(l)));
    }
    table.push_str("</tbody></table>");
    Chart {
        svg: format!(
            r#"<svg viewBox="0 0 {width} {hgt}" class="chart" role="img" aria-label="{} counts">{out}</svg>"#,
            esc(title_col)
        ),
        table,
    }
}

/// Genus bars grouped under phylum headings and coloured by phylum (heading, legend, tooltip and table all name the phylum,
/// so identity never rests on colour alone). items = (genus, phylum, inputs). Returns (chart, legend html).
pub fn genus_by_phylum(items: &[(String, String, usize)]) -> (Chart, String) {
    if items.is_empty() {
        return (
            Chart {
                svg: String::new(),
                table: String::new(),
            },
            String::new(),
        );
    }
    let (width, row_h, head_h, label_w) = (420.0, 24.0, 26.0, 130.0);
    let mut sorted: Vec<&(String, String, usize)> = items.iter().collect();
    sorted.sort_by(|a, b| {
        ph_order(&a.1)
            .cmp(&ph_order(&b.1))
            .then(a.1.cmp(&b.1))
            .then(b.2.cmp(&a.2))
            .then(a.0.cmp(&b.0))
    });
    let mut groups: Vec<(String, Vec<(&str, usize)>)> = Vec::new();
    for (g, ph, n) in sorted {
        let ph = if ph.is_empty() {
            "(not called)".to_string()
        } else {
            ph.clone()
        };
        match groups.last_mut() {
            Some((p, v)) if *p == ph => v.push((g.as_str(), *n)),
            _ => groups.push((ph, vec![(g.as_str(), *n)])),
        }
    }
    let vmax = items.iter().map(|t| t.2).max().unwrap_or(1).max(1) as f64;
    let plot_w = width - label_w - 40.0;
    let mut y = 2.0;
    let (mut out, mut trs) = (String::new(), String::new());
    for (ph, gl) in &groups {
        out.push_str(&format!(
            r#"<circle cx="{}" cy="{}" r="4.5" fill="{}"/>"#,
            label_w - 117.0,
            y + 11.0,
            ph_color(ph)
        ));
        let name = if PHYLUM_SLOTS.contains(&ph.as_str()) || ph == "(not called)" {
            ph.clone()
        } else {
            format!("{ph} (grey: other phyla)")
        };
        out.push_str(&text(
            label_w - 110.0,
            y + 15.0,
            &name,
            "t-sec",
            "start",
            r#"style="font-weight:600""#,
        ));
        y += head_h;
        for (g, n) in gl {
            let bw = (plot_w * *n as f64 / vmax).max(8.0);
            let tip = format!("{g} ({ph}): {n} inputs");
            out.push_str(&text(
                label_w - 8.0,
                y + 14.0,
                &short(g, 20),
                "t-sec",
                "end",
                "",
            ));
            out.push_str(&format!(
                r#"<path class="tip bar" tabindex="0" data-tip="{}" d="{}" fill="{}"/>"#,
                esc(&tip),
                hbar_path(label_w, y + 3.0, bw, 18.0),
                ph_color(ph)
            ));
            out.push_str(&format!(
                r#"<rect class="tip hit" tabindex="-1" data-tip="{}" x="{label_w}" y="{}" width="{}" height="{}" fill="transparent"/>"#,
                esc(&tip),
                y - 1.0,
                plot_w + 36.0,
                row_h - 2.0
            ));
            out.push_str(&text(
                label_w + bw + 6.0,
                y + 14.5,
                &n.to_string(),
                "t-pri",
                "start",
                "",
            ));
            trs.push_str(&format!(
                "<tr><td>{} {}</td><td>{}</td><td>{n}</td></tr>",
                ph_dot(ph),
                esc(ph),
                esc(g)
            ));
            y += row_h;
        }
        y += 6.0;
    }
    let svg = format!(
        r#"<svg viewBox="0 0 {width} {y}" class="chart" role="img" aria-label="Genus counts grouped by phylum"><line x1="{label_w}" x2="{label_w}" y1="0" y2="{y}" class="axis"/>{out}</svg>"#
    );
    let table = format!("<table><thead><tr><th>Phylum</th><th>Genus</th><th>Inputs</th></tr></thead><tbody>{trs}</tbody></table>");
    let mut legend = String::from(r#"<div class="legend">"#);
    for (ph, _) in groups.iter().filter(|(p, _)| p != "(not called)") {
        legend.push_str(&format!("<span>{}{}</span>", ph_dot(ph), esc(ph)));
    }
    legend.push_str("</div>");
    (Chart { svg, table }, legend)
}

// ---------------------------------------------------------------- ANI histogram
pub fn ani_hist(anis: &[f64], genus_thr: f64, species_thr: f64) -> Chart {
    let (lo, hi, bw): (f64, f64, f64) = (0.70, 1.00, 0.01);
    let nb = ((hi - lo) / bw).round() as usize;
    let mut counts = vec![0usize; nb];
    for &v in anis {
        let i = if v >= lo {
            (((v - lo) / bw) as usize).min(nb - 1)
        } else {
            0
        };
        counts[i] += 1;
    }
    let (w, h, l, r, t, b) = (500.0, 280.0, 38.0, 8.0, 34.0, 34.0);
    let (pw, ph) = (w - l - r, h - t - b);
    let cmax = *counts.iter().max().unwrap_or(&0);
    let ticks = nice_ticks(cmax as f64, 4);
    let ymax = *ticks.last().unwrap_or(&1.0);
    let mut out = String::new();
    for tk in &ticks {
        let y = t + ph - ph * tk / ymax;
        out.push_str(&format!(
            r#"<line x1="{l}" x2="{}" y1="{y:.1}" y2="{y:.1}" class="grid"/>"#,
            w - r
        ));
        out.push_str(&text(l - 8.0, y + 4.0, &compact(*tk), "t-mut", "end", ""));
    }
    let slot = pw / nb as f64;
    let bwid = (slot - 3.0).clamp(6.0, 24.0);
    for (i, &c) in counts.iter().enumerate() {
        let x = l + i as f64 * slot + (slot - bwid) / 2.0;
        let (a, bb) = (lo + i as f64 * bw, lo + (i as f64 + 1.0) * bw);
        let tip = format!("ANI {a:.2}–{bb:.2}: {c} inputs");
        let hgt = ph * c as f64 / ymax;
        if c > 0 {
            out.push_str(&format!(r#"<path class="tip bar" tabindex="0" data-tip="{}" d="{}" fill="var(--series-1)"/>"#, esc(&tip), column_path(x, t + ph, bwid, hgt)));
        }
        out.push_str(&format!(
            r#"<rect class="tip hit" data-tip="{}" x="{:.1}" y="{t}" width="{slot:.1}" height="{ph}" fill="transparent"/>"#,
            esc(&tip),
            l + i as f64 * slot
        ));
        if c == cmax && c > 0 {
            out.push_str(&text(
                x + bwid / 2.0,
                t + ph - hgt - 6.0,
                &c.to_string(),
                "t-pri",
                "middle",
                "",
            )); // label the extreme only
        }
    }
    out.push_str(&format!(
        r#"<line x1="{l}" x2="{}" y1="{}" y2="{}" class="axis"/>"#,
        w - r,
        t + ph,
        t + ph
    ));
    for v in [0.70, 0.80, 0.90, 1.00] {
        let x = l + pw * (v - lo) / (hi - lo);
        out.push_str(&text(
            x,
            t + ph + 18.0,
            &format!("{v:.2}"),
            "t-mut",
            if v >= 1.0 { "end" } else { "middle" },
            "",
        ));
    }
    // threshold hairlines, directly labelled (solid, never dashed)
    for (v, lab) in [(genus_thr, "genus"), (species_thr, "species")] {
        if !(lo..=hi).contains(&v) {
            continue;
        }
        let x = l + pw * (v - lo) / (hi - lo);
        out.push_str(&format!(
            r#"<line x1="{x:.1}" x2="{x:.1}" y1="{}" y2="{}" class="thr"/>"#,
            t - 12.0,
            t + ph
        ));
        let (dx, anchor) = if v > 0.9 {
            (-4.0, "end")
        } else {
            (4.0, "start")
        };
        out.push_str(&text(
            x + dx,
            t - 16.0,
            &format!("{lab} ≥ {v:.3}"),
            "t-sec",
            anchor,
            "",
        ));
    }
    out.push_str(&text(
        l + pw / 2.0,
        h - 2.0,
        "ANI of the nearest reference (containment, k = 21)",
        "t-mut",
        "middle",
        "",
    ));
    let mut table =
        String::from("<table><thead><tr><th>ANI bin</th><th>Inputs</th></tr></thead><tbody>");
    for (i, &c) in counts.iter().enumerate().filter(|(_, c)| **c > 0) {
        table.push_str(&format!(
            "<tr><td>{:.2}–{:.2}</td><td>{c}</td></tr>",
            lo + i as f64 * bw,
            lo + (i as f64 + 1.0) * bw
        ));
    }
    table.push_str("</tbody></table>");
    Chart {
        svg: format!(
            r#"<svg viewBox="0 0 {w} {h}" class="chart" role="img" aria-label="Histogram of nearest-reference ANI">{out}</svg>"#
        ),
        table,
    }
}

// ---------------------------------------------------------------- ANI heatmap among inputs
const HEAT: [&str; 6] = [
    "#cde2fb", "#9ec5f4", "#5598e7", "#2a78d6", "#1c5cab", "#104281",
];
const HEAT_BINS: [f64; 5] = [0.90, 0.95, 0.98, 0.99, 0.995];
const HEAT_LAB: [&str; 6] = [
    "< 0.90",
    "0.90–0.95",
    "0.95–0.98",
    "0.98–0.99",
    "0.99–0.995",
    "≥ 0.995",
];

/// ids, labels (species or genus), ani[i][j]
pub fn heatmap(ids: &[String], labels: &[String], ani: &[Vec<f64>]) -> Chart {
    let n = ids.len();
    let (cell, l, t) = (34.0, 188.0, 26.0);
    let (w, h) = (l + cell * n as f64 + 8.0, t + cell * n as f64 + 62.0);
    let mut out = String::new();
    for i in 0..n {
        out.push_str(&text(
            l - 8.0,
            t + i as f64 * cell + cell / 2.0 + 4.0,
            &format!("{}  {}", i + 1, short(&labels[i], 20)),
            "t-sec",
            "end",
            "",
        ));
        out.push_str(&text(
            l + i as f64 * cell + cell / 2.0,
            t - 8.0,
            &(i + 1).to_string(),
            "t-mut",
            "middle",
            "",
        ));
        for (j, &v) in ani[i].iter().enumerate().take(n) {
            let b = HEAT_BINS.iter().filter(|th| v >= **th).count();
            let tip = format!("{} × {}: ANI {v:.4}", ids[i], ids[j]);
            out.push_str(&format!(
                r#"<rect class="tip" tabindex="0" data-tip="{}" x="{:.1}" y="{:.1}" width="{}" height="{}" rx="3" fill="{}"/>"#,
                esc(&tip),
                l + j as f64 * cell + 1.0,
                t + i as f64 * cell + 1.0,
                cell - 2.0,
                cell - 2.0,
                HEAT[b]
            ));
        }
    }
    let ly = t + cell * n as f64 + 30.0;
    out.push_str(&text(l, ly - 4.0, "Pairwise ANI", "t-mut", "start", ""));
    let mut x = l + 82.0;
    for b in 0..6 {
        out.push_str(&format!(
            r#"<rect x="{x:.1}" y="{}" width="12" height="12" rx="2" fill="{}"/>"#,
            ly - 14.0,
            HEAT[b]
        ));
        out.push_str(&text(x + 17.0, ly - 4.0, HEAT_LAB[b], "t-sec", "start", ""));
        x += 17.0 + HEAT_LAB[b].chars().count() as f64 * 6.4 + 16.0;
    }
    let mut table = String::from("<table class='wide'><thead><tr><th></th>");
    for i in 0..n {
        table.push_str(&format!("<th>{}</th>", i + 1));
    }
    table.push_str("</tr></thead><tbody>");
    for (i, row) in ani.iter().enumerate().take(n) {
        table.push_str(&format!("<tr><th>{} {}</th>", i + 1, esc(&ids[i])));
        for v in row.iter().take(n) {
            table.push_str(&format!("<td>{v:.3}</td>"));
        }
        table.push_str("</tr>");
    }
    table.push_str("</tbody></table>");
    Chart {
        svg: format!(
            r#"<svg viewBox="0 0 {w} {h}" class="chart" role="img" aria-label="Pairwise ANI between inputs">{out}</svg>"#
        ),
        table,
    }
}

// ---------------------------------------------------------------- best matching references (one input)
pub fn tophit_bars(r: &Rec) -> Chart {
    if r.hits.is_empty() {
        return Chart {
            svg: String::new(),
            table: String::new(),
        };
    }
    let (w, rh, lw) = (560.0, 40.0, 218.0);
    let pw = w - lw - 150.0;
    let qh = r.n_hashes.max(1) as f64;
    let mut out = format!(
        r#"<line x1="{lw}" x2="{lw}" y1="0" y2="{}" class="axis"/>"#,
        rh * r.hits.len() as f64
    );
    let mut trs = String::new();
    for (i, hit) in r.hits.iter().enumerate() {
        let c = hit.shared as f64 / qh;
        let ani = c.powf(1.0 / 21.0);
        let y = i as f64 * rh + 4.0;
        let bw = (pw * c).max(8.0);
        let sp = [&hit.lineage[5], &hit.lineage[4]]
            .iter()
            .find(|s| !s.is_empty())
            .map(|s| s.to_string())
            .unwrap_or_else(|| "(unnamed)".into());
        let tip = format!(
            "{} · {sp} · {:.1}% of the query's hashes found (ANI {ani:.4}) · {} shared hashes",
            hit.id,
            100.0 * c,
            hit.shared
        );
        out.push_str(&text(
            lw - 10.0,
            y + 12.0,
            &short(&sp, 25),
            "t-pri",
            "end",
            "",
        ));
        out.push_str(&text(
            lw - 10.0,
            y + 27.0,
            &short(&hit.id, 25),
            "t-mut",
            "end",
            "",
        ));
        let fill = if i == 0 {
            "var(--series-1)"
        } else {
            "var(--deemph)"
        };
        out.push_str(&format!(
            r#"<path class="tip bar" tabindex="0" data-tip="{}" d="{}" fill="{fill}"/>"#,
            esc(&tip),
            hbar_path(lw, y + 6.0, bw, 18.0)
        ));
        out.push_str(&format!(r#"<rect class="tip hit" data-tip="{}" x="{lw}" y="{y}" width="{}" height="{}" fill="transparent"/>"#, esc(&tip), pw + 110.0, rh - 4.0));
        let lab = if i == 0 {
            format!("{:.1}%  ·  ANI {ani:.3}", 100.0 * c)
        } else {
            format!("{:.1}%", 100.0 * c)
        };
        out.push_str(&text(lw + bw + 7.0, y + 21.0, &lab, "t-pri", "start", ""));
        trs.push_str(&format!(
            "<tr><td>{}</td><td>{}</td><td>{}</td><td>{}</td><td>{:.1}%</td><td>{ani:.4}</td></tr>",
            i + 1,
            esc(&hit.id),
            esc(&sp),
            hit.shared,
            100.0 * c
        ));
    }
    let table = format!("<table><thead><tr><th>#</th><th>Reference</th><th>Species label</th><th>Shared hashes</th><th>Query found</th><th>ANI</th></tr></thead><tbody>{trs}</tbody></table>");
    Chart {
        svg: format!(
            r#"<svg viewBox="0 0 {w} {}" class="chart" role="img" aria-label="Best matching reference genomes">{out}</svg>"#,
            rh * r.hits.len() as f64 + 6.0
        ),
        table,
    }
}

// ---------------------------------------------------------------- reads: abundance histogram
pub fn abundance_chart(h: &Hist) -> Chart {
    let (w, hh, l, r, t, b) = (860.0, 210.0, 44.0, 10.0, 14.0, 34.0);
    let (pw, ph) = (w - l - r, hh - t - b);
    let mut vals: Vec<u64> = h.counts.iter().take(12).cloned().collect();
    vals.resize(12, 0);
    vals.push(h.over12);
    let labels: Vec<String> = (1..=12)
        .map(|i| i.to_string())
        .chain(std::iter::once("13+".to_string()))
        .collect();
    let vmax = *vals.iter().max().unwrap_or(&1);
    let ticks = nice_ticks(vmax as f64, 4);
    let ymax = *ticks.last().unwrap_or(&1.0);
    let mut out = String::new();
    for tk in &ticks {
        let y = t + ph - ph * tk / ymax;
        out.push_str(&format!(
            r#"<line x1="{l}" x2="{}" y1="{y:.1}" y2="{y:.1}" class="grid"/>"#,
            w - r
        ));
        out.push_str(&text(l - 8.0, y + 4.0, &compact(*tk), "t-mut", "end", ""));
    }
    let slot = pw / vals.len() as f64;
    let bwid = (slot - 6.0).clamp(6.0, 24.0);
    for (i, &v) in vals.iter().enumerate() {
        let x = l + i as f64 * slot + (slot - bwid) / 2.0;
        let kept = (i as u32 + 1) >= h.min_abund;
        let tip = format!(
            "seen {}×: {v} hashes · {}",
            labels[i],
            if kept {
                "kept".to_string()
            } else {
                format!("removed (below minimum abundance {})", h.min_abund)
            }
        );
        let hgt = ph * v as f64 / ymax;
        let fill = if kept {
            "var(--series-1)"
        } else {
            "var(--deemph)"
        };
        out.push_str(&format!(
            r#"<path class="tip bar" tabindex="0" data-tip="{}" d="{}" fill="{fill}"/>"#,
            esc(&tip),
            column_path(x, t + ph, bwid, hgt)
        ));
        out.push_str(&format!(r#"<rect class="tip hit" data-tip="{}" x="{:.1}" y="{t}" width="{slot:.1}" height="{ph}" fill="transparent"/>"#, esc(&tip), l + i as f64 * slot));
        out.push_str(&text(
            x + bwid / 2.0,
            t + ph + 16.0,
            &labels[i],
            "t-mut",
            "middle",
            "",
        ));
    }
    out.push_str(&format!(
        r#"<line x1="{l}" x2="{}" y1="{}" y2="{}" class="axis"/>"#,
        w - r,
        t + ph,
        t + ph
    ));
    out.push_str(&text(
        l + pw / 2.0,
        hh - 2.0,
        "times a sketch hash was seen in the reads",
        "t-mut",
        "middle",
        "",
    ));
    if h.min_abund > 1 && vals[0] > 0 {
        out.push_str(&text(
            l + slot * 1.15,
            t + 14.0,
            &format!(
                "← {} hashes seen once: mostly sequencing errors, removed",
                compact(vals[0] as f64)
            ),
            "t-sec",
            "start",
            "",
        ));
    }
    let mut table = String::from(
        "<table><thead><tr><th>Seen</th><th>Hashes</th><th>Used</th></tr></thead><tbody>",
    );
    for (i, v) in vals.iter().enumerate() {
        table.push_str(&format!(
            "<tr><td>{}×</td><td>{v}</td><td>{}</td></tr>",
            labels[i],
            if (i as u32 + 1) >= h.min_abund {
                "yes"
            } else {
                "no"
            }
        ));
    }
    table.push_str("</tbody></table>");
    Chart {
        svg: format!(
            r#"<svg viewBox="0 0 {w} {hh}" class="chart" role="img" aria-label="Abundance histogram of read hashes">{out}</svg>"#
        ),
        table,
    }
}

/// phylum -> count, sorted by count (desc) then name; used by the page.
pub fn count_by<'a, I: Iterator<Item = &'a str>>(it: I) -> Vec<(String, usize)> {
    let mut m: BTreeMap<String, usize> = BTreeMap::new();
    for s in it {
        *m.entry(s.to_string()).or_default() += 1;
    }
    let mut v: Vec<(String, usize)> = m.into_iter().collect();
    v.sort_by(|a, b| b.1.cmp(&a.1).then(a.0.cmp(&b.0)));
    v
}
