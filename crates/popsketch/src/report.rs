//! `report.html`: one standalone page with inline SVG. The audit comes first; the population structure follows, always with the
//! resolution limits of the data.

use crate::cluster::Verdict;
use crate::pipeline::{Options, Results, Status};
use crate::tree::{root_tree, Rooted};
use anyhow::Result;
use std::collections::BTreeMap;
use std::fmt::Write as _;
use std::fs;

const PALETTE: [&str; 10] = [
    "#1b6ca8", "#d95f02", "#1b9e77", "#7570b3", "#e7298a", "#66a61e", "#a6761d", "#e6ab02",
    "#666666", "#17becf",
];

fn esc(s: &str) -> String {
    s.replace('&', "&amp;")
        .replace('<', "&lt;")
        .replace('>', "&gt;")
        .replace('"', "&quot;")
}

struct Pt {
    x: f64,
    y: f64,
    color: usize,
    /// 0 circle (included), 1 triangle (marked), 2 diamond (decontaminated), 3 cross (excluded)
    shape: u8,
    title: String,
}

fn marker(p: &Pt, sx: f64, sy: f64) -> String {
    let c = PALETTE[p.color % PALETTE.len()];
    let t = esc(&p.title);
    match p.shape {
        1 => format!("<path d=\"M{:.1},{:.1} l4,7 h-8 z\" fill=\"{c}\" fill-opacity=\".8\"><title>{t}</title></path>", sx, sy - 4.0),
        2 => format!("<path d=\"M{:.1},{:.1} l4,4 l-4,4 l-4,-4 z\" fill=\"{c}\" fill-opacity=\".8\"><title>{t}</title></path>", sx, sy - 4.0),
        3 => format!("<path d=\"M{:.1},{:.1} l6,6 m0,-6 l-6,6\" stroke=\"{c}\" stroke-width=\"1.6\" fill=\"none\"><title>{t}</title></path>", sx - 3.0, sy - 3.0),
        _ => format!("<circle cx=\"{sx:.1}\" cy=\"{sy:.1}\" r=\"3.4\" fill=\"{c}\" fill-opacity=\".75\"><title>{t}</title></circle>"),
    }
}

fn nice(lo: f64, hi: f64) -> (f64, f64) {
    if (hi - lo).abs() < 1e-12 {
        (lo - 1.0, hi + 1.0)
    } else {
        let pad = (hi - lo) * 0.06;
        (lo - pad, hi + pad)
    }
}

/// Scatter plot; `logx` plots log10 of x. Axis labels are the ends of the range.
fn scatter(pts: &[Pt], xl: &str, yl: &str, logx: bool) -> String {
    let (w, h, l, b, t, r) = (520.0, 360.0, 56.0, 40.0, 10.0, 12.0);
    let xs: Vec<f64> = pts
        .iter()
        .map(|p| if logx { p.x.max(1e-3).log10() } else { p.x })
        .collect();
    let ys: Vec<f64> = pts.iter().map(|p| p.y).collect();
    if pts.is_empty() {
        return String::from("<p class=\"note\">nothing to plot</p>");
    }
    let (x0, x1) = nice(
        xs.iter().cloned().fold(f64::INFINITY, f64::min),
        xs.iter().cloned().fold(f64::NEG_INFINITY, f64::max),
    );
    let (y0, y1) = nice(
        ys.iter().cloned().fold(f64::INFINITY, f64::min),
        ys.iter().cloned().fold(f64::NEG_INFINITY, f64::max),
    );
    let px = |x: f64| l + (x - x0) / (x1 - x0) * (w - l - r);
    let py = |y: f64| h - b - (y - y0) / (y1 - y0) * (h - b - t);
    let mut s = format!("<svg viewBox=\"0 0 {w} {h}\" class=\"plot\" role=\"img\"><rect x=\"{l}\" y=\"{t}\" width=\"{}\" height=\"{}\" class=\"frame\"/>", w - l - r, h - b - t);
    for i in 0..=4 {
        let fx = x0 + (x1 - x0) * i as f64 / 4.0;
        let fy = y0 + (y1 - y0) * i as f64 / 4.0;
        let lab = if logx {
            format!("{:.0}", 10f64.powf(fx))
        } else {
            format!("{fx:.3}")
        };
        let _ = write!(
            s,
            "<text x=\"{:.1}\" y=\"{:.1}\" class=\"tick\" text-anchor=\"middle\">{lab}</text>",
            px(fx),
            h - b + 14.0
        );
        let _ = write!(
            s,
            "<text x=\"{:.1}\" y=\"{:.1}\" class=\"tick\" text-anchor=\"end\">{fy:.3}</text>",
            l - 5.0,
            py(fy) + 3.0
        );
    }
    for (p, x) in pts.iter().zip(&xs) {
        s.push_str(&marker(p, px(*x), py(p.y)));
    }
    let _ = write!(
        s,
        "<text x=\"{:.1}\" y=\"{:.1}\" class=\"axis\" text-anchor=\"middle\">{}</text>",
        (l + w - r) / 2.0,
        h - 6.0,
        esc(xl)
    );
    let _ = write!(s, "<text transform=\"translate(14,{:.1}) rotate(-90)\" class=\"axis\" text-anchor=\"middle\">{}</text></svg>", (t + h - b) / 2.0, esc(yl));
    s
}

fn legend(items: &[(String, usize)]) -> String {
    let mut s = String::from("<p class=\"legend\">");
    for (name, c) in items {
        let _ = write!(
            s,
            "<span><i style=\"background:{}\"></i>{}</span> ",
            PALETTE[c % PALETTE.len()],
            esc(name)
        );
    }
    s.push_str("</p>");
    s
}

fn histogram(vals: &[f64], floor: f64) -> String {
    if vals.is_empty() {
        return String::new();
    }
    let (w, h, l, b, t) = (520.0, 220.0, 56.0, 34.0, 10.0);
    let lo = vals
        .iter()
        .cloned()
        .fold(f64::INFINITY, f64::min)
        .min(floor);
    let hi = vals
        .iter()
        .cloned()
        .fold(f64::NEG_INFINITY, f64::max)
        .max(floor)
        + 1e-9;
    let bins = 40usize;
    let mut cnt = vec![0usize; bins];
    for &v in vals {
        cnt[(((v - lo) / (hi - lo)) * bins as f64).min(bins as f64 - 1.0) as usize] += 1;
    }
    let mx = *cnt.iter().max().unwrap() as f64;
    let px = |x: f64| l + (x - lo) / (hi - lo) * (w - l - 12.0);
    let mut s = format!("<svg viewBox=\"0 0 {w} {h}\" class=\"plot\" role=\"img\">");
    for (i, c) in cnt.iter().enumerate() {
        let x = lo + (hi - lo) * i as f64 / bins as f64;
        let bh = (*c as f64 / mx) * (h - b - t);
        let _ = write!(s, "<rect x=\"{:.1}\" y=\"{:.1}\" width=\"{:.1}\" height=\"{bh:.1}\" fill=\"{}\"><title>{c} pairs</title></rect>", px(x), h - b - bh, (w - l - 12.0) / bins as f64 - 0.5, PALETTE[0]);
    }
    let fx = px(floor);
    let _ = write!(s, "<line x1=\"{fx:.1}\" x2=\"{fx:.1}\" y1=\"{t}\" y2=\"{}\" class=\"floor\"/><text x=\"{:.1}\" y=\"{}\" class=\"tick\">noise floor {floor:.4}</text>", h - b, fx + 4.0, t + 10.0);
    for i in 0..=4 {
        let v = lo + (hi - lo) * i as f64 / 4.0;
        let _ = write!(
            s,
            "<text x=\"{:.1}\" y=\"{:.1}\" class=\"tick\" text-anchor=\"middle\">{v:.4}</text>",
            px(v),
            h - b + 14.0
        );
    }
    let _ = write!(s, "<text x=\"{:.1}\" y=\"{}\" class=\"axis\" text-anchor=\"middle\">pairwise distance (1 - ANI), samples in the analysis</text></svg>", (l + w) / 2.0, h - 4.0);
    s
}

fn draw_tree(
    r: &Rooted,
    names: &[String],
    color_of: &dyn Fn(usize) -> usize,
    max_labels: usize,
) -> String {
    let n_tips = r
        .nodes
        .iter()
        .filter(|x| x.tip.is_some() && x.children.is_empty())
        .count();
    // x = distance from the root, y = leaf order
    let mut depth = vec![0.0f64; r.nodes.len()];
    let mut ypos = vec![0.0f64; r.nodes.len()];
    let mut next_leaf = 0.0;
    fn walk(r: &Rooted, i: usize, depth: &mut [f64], ypos: &mut [f64], next: &mut f64) {
        for &c in &r.nodes[i].children {
            depth[c] = depth[i] + r.nodes[c].len;
            walk(r, c, depth, ypos, next);
        }
        ypos[i] = if r.nodes[i].children.is_empty() {
            *next += 1.0;
            *next - 1.0
        } else {
            let ys: Vec<f64> = r.nodes[i].children.iter().map(|&c| ypos[c]).collect();
            (ys.iter().cloned().fold(f64::INFINITY, f64::min)
                + ys.iter().cloned().fold(f64::NEG_INFINITY, f64::max))
                / 2.0
        };
    }
    walk(r, r.root, &mut depth, &mut ypos, &mut next_leaf);
    let maxd = depth.iter().cloned().fold(1e-12, f64::max);
    let (w, rowh) = (720.0, 9.0f64);
    let h = (n_tips as f64 * rowh + 24.0).max(60.0);
    let labels = n_tips <= max_labels;
    let right = if labels { 150.0 } else { 20.0 };
    let px = |d: f64| 10.0 + d / maxd * (w - 10.0 - right);
    let py = |y: f64| 12.0 + y * rowh;
    let mut s = format!("<svg viewBox=\"0 0 {w} {h}\" class=\"plot tree\" role=\"img\">");
    for (i, node) in r.nodes.iter().enumerate() {
        if !node.children.is_empty() {
            let ys: Vec<f64> = node.children.iter().map(|&c| ypos[c]).collect();
            let (a, b) = (
                ys.iter().cloned().fold(f64::INFINITY, f64::min),
                ys.iter().cloned().fold(f64::NEG_INFINITY, f64::max),
            );
            let _ = write!(
                s,
                "<path d=\"M{x:.1},{:.1} V{:.1}\" class=\"br\"/>",
                py(a),
                py(b),
                x = px(depth[i])
            );
        }
        if i != r.root {
            let parent_depth = depth[i] - node.len;
            let _ = write!(
                s,
                "<path d=\"M{:.1},{y:.1} H{:.1}\" class=\"br\"/>",
                px(parent_depth),
                px(depth[i]),
                y = py(ypos[i])
            );
        }
        if let Some(t) = node.tip {
            if node.children.is_empty() {
                let c = PALETTE[color_of(t) % PALETTE.len()];
                let _ = write!(s, "<circle cx=\"{:.1}\" cy=\"{:.1}\" r=\"2.6\" fill=\"{c}\"><title>{}</title></circle>", px(depth[i]), py(ypos[i]), esc(&names[t]));
                if labels {
                    let _ = write!(
                        s,
                        "<text x=\"{:.1}\" y=\"{:.1}\" class=\"tl\">{}</text>",
                        px(depth[i]) + 5.0,
                        py(ypos[i]) + 3.0,
                        esc(&names[t])
                    );
                }
            }
        }
    }
    s.push_str("</svg>");
    s
}

const CSS: &str = r#"
:root{--bg:#fff;--fg:#1d2330;--mute:#5b6577;--line:#d7dce5;--card:#f6f8fb;--warn:#a4560e;--bad:#b3261e}
@media (prefers-color-scheme:dark){:root{--bg:#14171d;--fg:#e6e9ef;--mute:#9aa4b5;--line:#2c3340;--card:#1b2029;--warn:#e0a458;--bad:#f2766d}}
body{background:var(--bg);color:var(--fg);font:15px/1.5 system-ui,sans-serif;margin:0 auto;max-width:980px;padding:16px}
h1{font-size:22px;margin:.2em 0}h2{font-size:18px;margin:1.6em 0 .4em;border-bottom:1px solid var(--line);padding-bottom:4px}
h3{font-size:15px;margin:1.2em 0 .3em}p.note,.note{color:var(--mute);font-size:13px}
table{border-collapse:collapse;width:100%;font-size:13px}th,td{border-bottom:1px solid var(--line);padding:3px 6px;text-align:left;vertical-align:top}
th{color:var(--mute);font-weight:600}code{background:var(--card);padding:0 3px;border-radius:3px}
.card{background:var(--card);border:1px solid var(--line);border-radius:6px;padding:8px 12px;margin:8px 0}
.plot{width:100%;max-width:560px;height:auto}.tree{max-width:760px}.frame{fill:none;stroke:var(--line)}
.tick,.axis,.tl{fill:var(--mute);font-size:10px}.axis{font-size:11px}.floor{stroke:var(--bad);stroke-dasharray:4 3}
.br{stroke:var(--mute);fill:none;stroke-width:1}.legend span{margin-right:12px;font-size:13px}.legend i{display:inline-block;width:10px;height:10px;border-radius:5px;margin-right:4px}
.warn{color:var(--warn)}.bad{color:var(--bad)}.grid{display:flex;flex-wrap:wrap;gap:16px}
"#;

pub fn write(o: &Options, r: &Results) -> Result<()> {
    let mut h = String::new();
    let _ = write!(h, "<!doctype html><html lang=\"en\"><head><meta charset=\"utf-8\"><meta name=\"viewport\" content=\"width=device-width,initial-scale=1\"><title>popsketch report</title><style>{CSS}</style></head><body>");
    let n_total = r.rows.len();
    let count = |s: Status| r.rows.iter().filter(|x| x.status == s).count();
    let _ = write!(
        h,
        "<h1>popsketch: sample audit and coarse structure</h1><p class=\"note\">k={}, scaled={}, {} per sample, flagged samples: {:?}. This tool compares sketches, not alignments: it audits a set of read files and shows coarse structure; it does not build a phylogeny.</p>",
        o.ksize,
        o.scaled,
        if o.max_bases == 0 { "no base limit".to_string() } else { format!("up to {} bases", o.max_bases) },
        o.flagged
    );
    // ---- audit
    let _ = write!(
        h,
        "<h2>1. Audit</h2><div class=\"card\"><b>{n_total}</b> samples: <b>{}</b> in the analysis ({} included, {} marked, {} decontaminated), <b>{}</b> excluded, <b>{}</b> unreadable.</div>",
        r.included.len(), count(Status::Included), count(Status::Marked), count(Status::Decontaminated), count(Status::Excluded), count(Status::Error)
    );
    let mut by_flag: BTreeMap<&str, Vec<&str>> = BTreeMap::new();
    for row in &r.rows {
        for f in row.flags() {
            by_flag.entry(f.code).or_default().push(&row.id);
        }
    }
    if by_flag.is_empty() {
        h.push_str("<p>No sample was flagged.</p>");
    } else {
        h.push_str("<table><tr><th>flag</th><th>samples</th><th>which</th></tr>");
        for (code, ids) in &by_flag {
            let shown: Vec<&str> = ids.iter().take(12).copied().collect();
            let more = if ids.len() > 12 {
                format!(" and {} more", ids.len() - 12)
            } else {
                String::new()
            };
            let _ = write!(
                h,
                "<tr><td><code>{code}</code></td><td>{}</td><td>{}{more}</td></tr>",
                ids.len(),
                esc(&shown.join(", "))
            );
        }
        h.push_str("</table><h3>Flagged samples, with the reason</h3><table><tr><th>sample</th><th>status</th><th>reasons</th></tr>");
        for row in r
            .rows
            .iter()
            .filter(|x| x.qc.as_ref().is_some_and(|q| q.flagged()))
        {
            let why: Vec<String> = row
                .flags()
                .iter()
                .map(|f| format!("<code>{}</code> {}", f.code, esc(&f.detail)))
                .collect();
            let _ = write!(
                h,
                "<tr><td>{}</td><td>{}</td><td>{}</td></tr>",
                esc(&row.id),
                row.status.as_str(),
                why.join("<br>")
            );
        }
        h.push_str("</table>");
    }
    for row in r.rows.iter().filter(|x| x.status == Status::Error) {
        let _ = write!(
            h,
            "<p class=\"bad\">{}: {}</p>",
            esc(&row.id),
            esc(row.message.as_deref().unwrap_or(""))
        );
    }
    h.push_str("<p class=\"note\"><b>What the audit cannot see.</b> Light contamination (about 1 to 7% of the reads) leaves no mark in the abundance spectrum or the size: a contaminant at 1x next to a 40x genome looks like sequencing error. Heavy contamination is caught by size (low coverage) or by an abundance shoulder (a contaminant at several-fold coverage). Giving an index (<code>--index</code>) adds the index's own mixture and organism calls.</p>");
    // size vs coverage
    let pts: Vec<Pt> = r
        .rows
        .iter()
        .filter_map(|row| {
            let q = row.qc.as_ref()?;
            let c = q.kmer_cov?;
            Some(Pt {
                x: c.max(0.1),
                y: q.est_mb,
                color: match row.status {
                    Status::Included => 0,
                    Status::Excluded => 4,
                    _ => 1,
                },
                shape: match row.status {
                    Status::Marked => 1,
                    Status::Decontaminated => 2,
                    Status::Excluded => 3,
                    _ => 0,
                },
                title: format!(
                    "{}: {:.1} Mb at {:.1}x{}",
                    row.id,
                    q.est_mb,
                    c,
                    if q.flagged() {
                        format!(
                            " ({})",
                            row.flags()
                                .iter()
                                .map(|f| f.code)
                                .collect::<Vec<_>>()
                                .join(", ")
                        )
                    } else {
                        String::new()
                    }
                ),
            })
        })
        .collect();
    if !pts.is_empty() {
        h.push_str("<h3>Estimated size against k-mer coverage</h3>");
        h.push_str(&scatter(
            &pts,
            "k-mer coverage (log scale)",
            "estimated size (Mb) = hashes above the cutoff x scaled",
            true,
        ));
        h.push_str(&legend(&[
            ("in the analysis".into(), 0),
            ("excluded (cross)".into(), 4),
            ("marked or decontaminated".into(), 1),
        ]));
    }
    // labels
    for rep in &r.label_reports {
        let _ = write!(h, "<h2>Labels: {}</h2><p class=\"note\">{} labelled samples in the analysis; groups: {}. A blank label means unlabelled.</p>", esc(&rep.name), rep.n_labelled, esc(&rep.groups.iter().map(|(g, n)| format!("{g} ({n})")).collect::<Vec<_>>().join(", ")));
        if rep.conflicts.is_empty() {
            h.push_str("<p>No sample matches another label better than its own.</p>");
        } else {
            h.push_str("<table><tr><th>sample</th><th>label</th><th>median ANI to own label</th><th>better match</th><th>median ANI there</th><th>nearest sample</th></tr>");
            for c in &rep.conflicts {
                let _ = write!(h, "<tr><td>{}</td><td>{}</td><td>{:.4}</td><td><b>{}</b></td><td>{:.4}</td><td>{} ({}, {:.4})</td></tr>", esc(&c.id), esc(&c.label), c.own_median_ani, esc(&c.other_label), c.other_median_ani, esc(&c.nearest), esc(c.nearest_label.as_deref().unwrap_or("unlabelled")), c.nearest_ani);
            }
            h.push_str("</table><p class=\"note\">A conflict means the sample's ANI to the members of another label is higher than to its own label by more than the noise floor. It is evidence to check the sample's metadata, not proof of an error.</p>");
        }
        if !rep.ari.is_empty() {
            h.push_str(
                "<table><tr><th>clusters</th><th>adjusted Rand index vs this label</th></tr>",
            );
            for (k, a) in &rep.ari {
                let _ = write!(h, "<tr><td>{k}</td><td>{a:.3}</td></tr>");
            }
            h.push_str("</table>");
        }
    }
    // ---- structure
    h.push_str("<h2>2. Coarse structure (secondary)</h2>");
    if let (Some(d), Some(ord)) = (&r.dist, &r.ordination) {
        let n = d.n;
        let names: Vec<String> = r.included.iter().map(|&i| r.rows[i].id.clone()).collect();
        // colour by the first label column when present, else by cluster
        let (colors, items): (Vec<usize>, Vec<(String, usize)>) =
            if let (Some(t), Some(_)) = (&r.label_table, r.label_reports.first()) {
                let labs = crate::labels::column_labels(t, &names, 0);
                let mut uniq: Vec<String> = labs.iter().flatten().cloned().collect();
                uniq.sort();
                uniq.dedup();
                let cs = labs
                    .iter()
                    .map(|l| {
                        l.as_ref()
                            .map_or(8, |l| uniq.iter().position(|u| u == l).unwrap() % 8)
                    })
                    .collect();
                let mut it: Vec<(String, usize)> = uniq
                    .iter()
                    .enumerate()
                    .map(|(i, u)| (u.clone(), i % 8))
                    .collect();
                it.push(("unlabelled".into(), 8));
                (cs, it)
            } else if let Some(cl) = &r.clusters {
                let k = cl.iter().copied().max().map_or(0, |m| m + 1);
                (
                    cl.iter().map(|&c| c % 8).collect(),
                    (0..k)
                        .map(|c| (format!("cluster {}", c + 1), c % 8))
                        .collect(),
                )
            } else {
                (vec![0; n], vec![("samples".into(), 0)])
            };
        let color_of = |i: usize| colors[i];
        let by = if r.label_table.is_some() {
            format!(
                "label column {}",
                r.label_reports.first().map_or("", |c| c.name.as_str())
            )
        } else {
            "cluster".to_string()
        };
        let shape = |i: usize| match r.rows[r.included[i]].status {
            Status::Marked => 1,
            Status::Decontaminated => 2,
            _ => 0,
        };
        let pc = |a: usize, b: usize| -> String {
            let pts: Vec<Pt> = (0..n)
                .map(|i| Pt {
                    x: ord.coords[i][a],
                    y: ord.coords[i][b],
                    color: color_of(i),
                    shape: shape(i),
                    title: format!(
                        "{} (cluster {})",
                        names[i],
                        r.clusters
                            .as_ref()
                            .map_or("NA".to_string(), |c| (c[i] + 1).to_string())
                    ),
                })
                .collect();
            scatter(
                &pts,
                &format!("PC{} ({:.1}% of variance)", a + 1, ord.var_frac[a] * 100.0),
                &format!("PC{} ({:.1}%)", b + 1, ord.var_frac[b] * 100.0),
                false,
            )
        };
        let _ = write!(h, "<p class=\"note\">{} on {} {}; coloured by {by}. Triangles are flagged samples kept with <code>--flagged mark</code>, diamonds are decontaminated.</p>", if ord.method == "hash_pca" { "Hash-presence PCA" } else { "PCoA of 1 - ANI" }, ord.n_features, if ord.method == "hash_pca" { "hashes" } else { "samples" });
        h.push_str("<div class=\"grid\">");
        if ord.eigenvalues.len() >= 2 {
            h.push_str(&pc(0, 1));
        }
        if ord.eigenvalues.len() >= 3 {
            h.push_str(&pc(1, 2));
        }
        h.push_str("</div>");
        h.push_str(&legend(&items));
        // clusters
        if let Some(g) = &r.gap {
            let (say, cls) = match &g.verdict {
                Verdict::Clear(k) => (format!("The largest jump in the clustering gives <b>{k} clusters</b>: it is {:.1} times the noise floor. A structureless population would show a gap ratio this large with probability {} (ratio {:.1}), so this is evidence of a split larger than the technical noise, not proof of distinct populations.", g.gap_over_floor, if g.null_p < 0.001 { "under 0.001".to_string() } else { format!("{:.3}", g.null_p) }, g.ratio), ""),
                Verdict::Single => (format!("Everything merges below the noise floor (the last merge is at {:.5}): one population as far as sketches can tell.", g.top_height), ""),
                Verdict::NoClearCut => (format!("<b>No single clear cut.</b> The largest jump is {:.1} times the noise floor (a cut needs to be at least {:.1}). Candidate cuts: {}. Use <code>--k</code> to choose one.", g.gap_over_floor, o.min_floors, g.candidates.iter().map(|c| format!("{} clusters (jump {:.5})", c.k, c.gap)).collect::<Vec<_>>().join("; ")), "warn"),
                Verdict::TooFew => ("Too few samples to cluster.".to_string(), "warn"),
            };
            let forced = if o.k.is_some() || o.cut.is_some() {
                format!(
                    " Clusters shown use {} clusters, as requested.",
                    r.k_used.unwrap_or(0)
                )
            } else {
                String::new()
            };
            let _ = write!(h, "<h3>Clusters</h3><p class=\"{cls}\">{say}{forced}</p>");
        }
        // distance histogram
        let mut vals = Vec::new();
        for i in 0..n {
            for j in (i + 1)..n {
                vals.push(1.0 - d.ani(i, j));
            }
        }
        h.push_str("<h3>Pairwise distances</h3>");
        h.push_str(&histogram(&vals, r.floor.value));
        // tree
        if let Some(nw) = &r.newick {
            let _ = nw;
            let t = {
                let mut tr = crate::tree::neighbor_joining(&d.mash, n, names.clone());
                crate::tree::collapse(&mut tr, o.tree_floor.unwrap_or(r.floor.value));
                root_tree(&tr, o.root.as_deref())
            };
            let _ = write!(h, "<h3>Tree</h3><p class=\"note\">Neighbour-joining on mash distance; {} branches shorter than the noise floor were collapsed, so unresolved parts are drawn as polytomies. {} Tips are coloured as above. This is a first-pass clustering tree, not a phylogeny: recombination, accessory content and uneven coverage distort it.</p>", r.collapsed, if o.root.is_some() { "Rooted on the named sample." } else { "Unrooted (drawn from an arbitrary node)." });
            h.push_str(&draw_tree(&t, &names, &color_of, 120));
            h.push_str("<p class=\"note\">The full tree is in <code>tree.nwk</code>.</p>");
        }
        // advisor
        if let Some(a) = &r.advice {
            h.push_str("<h2>3. Resolution: what these data can and cannot say</h2><ul>");
            for s in &a.statements {
                let _ = write!(h, "<li>{}</li>", esc(s));
            }
            h.push_str("</ul>");
        }
    } else {
        h.push_str("<p class=\"warn\">Too few samples passed the audit for distances, ordination or a tree.</p>");
    }
    if !r.warnings.is_empty() {
        h.push_str("<h2>Warnings</h2><ul>");
        for w in &r.warnings {
            let _ = write!(h, "<li>{}</li>", esc(w));
        }
        h.push_str("</ul>");
    }
    h.push_str("<p class=\"note\">Files: <code>qc.tsv</code> (every sample, flags with reasons), <code>ani.tsv</code>, <code>mash.tsv</code>, <code>pca.tsv</code>, <code>clusters.tsv</code>, <code>tree.nwk</code>, <code>label_check.tsv</code>, <code>popsketch.json</code>.</p></body></html>");
    fs::write(o.out.join("report.html"), h)?;
    Ok(())
}
