//! `report.html`: one standalone page with inline SVG. The audit comes first; the population structure follows, always with the
//! resolution limits of the data.
#![cfg_attr(not(feature = "native"), allow(dead_code))]

use crate::cluster::Verdict;
#[cfg(feature = "native")]
use crate::pipeline::ReportFormat;
use crate::pipeline::{Options, Results, Status};
use crate::tree::{root_tree, Rooted};
#[cfg(feature = "native")]
use anyhow::Result;
use std::collections::BTreeMap;
use std::fmt::Write as _;
#[cfg(feature = "native")]
use std::fs;

/// the deepest merge (1 - ANI) at which the set is said to span more than one species
const DEEP_SPLIT: f64 = 0.10;

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

// ---- a small document model rendered to HTML or Markdown ------------------------------------------------------------

enum In {
    T(String),
    C(String),
    B(String),
    Br,
}

fn t(s: impl Into<String>) -> In {
    In::T(s.into())
}
fn c(s: impl Into<String>) -> In {
    In::C(s.into())
}
fn b(s: impl Into<String>) -> In {
    In::B(s.into())
}

#[derive(Clone, Copy, PartialEq)]
enum Tone {
    Plain,
    Note,
    Warn,
    Bad,
}

enum Block {
    H(u8, String),
    P(Vec<In>, Tone),
    Card(Vec<In>),
    Table(Vec<&'static str>, Vec<Vec<Vec<In>>>),
    List(Vec<Vec<In>>),
    /// one or more SVG plots under one name
    Fig(&'static str, Vec<String>),
    Legend(Vec<(String, usize)>),
}

fn inl_html(v: &[In]) -> String {
    v.iter()
        .map(|x| match x {
            In::T(s) => esc(s),
            In::C(s) => format!("<code>{}</code>", esc(s)),
            In::B(s) => format!("<b>{}</b>", esc(s)),
            In::Br => "<br>".to_string(),
        })
        .collect()
}

fn inl_md(v: &[In]) -> String {
    v.iter()
        .map(|x| match x {
            In::T(s) => s.replace('|', "\\|").replace('\n', " "),
            In::C(s) => format!("`{}`", s.replace('`', "'").replace('|', "\\|")),
            In::B(s) => format!("**{}**", s.replace('|', "\\|")),
            In::Br => "<br>".to_string(),
        })
        .collect()
}

/// The HTML report as one string (no files written).
pub fn html(o: &Options, r: &Results) -> String {
    render_html(&build(o, r))
}

/// The Markdown report and its figures as (file name under `figures/`, SVG) pairs (no files written).
pub fn markdown(o: &Options, r: &Results) -> (String, Vec<(String, String)>) {
    md_text(&build(o, r), false)
}

fn render_html(blocks: &[Block]) -> String {
    let mut h = format!("<!doctype html><html lang=\"en\"><head><meta charset=\"utf-8\"><meta name=\"viewport\" content=\"width=device-width,initial-scale=1\"><title>popsketch report</title><style>{CSS}</style></head><body>");
    for bl in blocks {
        match bl {
            Block::H(l, s) => {
                let _ = write!(h, "<h{l}>{}</h{l}>", esc(s));
            }
            Block::P(v, tone) => {
                let cls = match tone {
                    Tone::Plain => "",
                    Tone::Note => " class=\"note\"",
                    Tone::Warn => " class=\"warn\"",
                    Tone::Bad => " class=\"bad\"",
                };
                let _ = write!(h, "<p{cls}>{}</p>", inl_html(v));
            }
            Block::Card(v) => {
                let _ = write!(h, "<div class=\"card\">{}</div>", inl_html(v));
            }
            Block::Table(head, rows) => {
                h.push_str("<table><tr>");
                for x in head {
                    let _ = write!(h, "<th>{}</th>", esc(x));
                }
                h.push_str("</tr>");
                for row in rows {
                    h.push_str("<tr>");
                    for cell in row {
                        let _ = write!(h, "<td>{}</td>", inl_html(cell));
                    }
                    h.push_str("</tr>");
                }
                h.push_str("</table>");
            }
            Block::List(items) => {
                h.push_str("<ul>");
                for it in items {
                    let _ = write!(h, "<li>{}</li>", inl_html(it));
                }
                h.push_str("</ul>");
            }
            Block::Fig(_, svgs) => {
                if svgs.len() > 1 {
                    h.push_str("<div class=\"grid\">");
                }
                for s in svgs {
                    h.push_str(s);
                }
                if svgs.len() > 1 {
                    h.push_str("</div>");
                }
            }
            Block::Legend(items) => h.push_str(&legend(items)),
        }
    }
    h.push_str("</body></html>");
    h
}

const SVG_CSS: &str = ".frame{fill:none;stroke:#c9cfda}.tick,.axis,.tl{fill:#5b6577;font-family:sans-serif;font-size:10px}.axis{font-size:11px}.floor{stroke:#b3261e;stroke-dasharray:4 3}.br{stroke:#5b6577;fill:none;stroke-width:1}";

/// An SVG that stands on its own: namespace, a width, and the styles the HTML page supplies through its stylesheet.
fn standalone_svg(svg: &str) -> String {
    let width = svg
        .split("viewBox=\"0 0 ")
        .nth(1)
        .and_then(|x| x.split_whitespace().next())
        .unwrap_or("560");
    let svg = svg.replacen(
        "<svg ",
        &format!("<svg xmlns=\"http://www.w3.org/2000/svg\" width=\"{width}\" "),
        1,
    );
    match svg.find('>') {
        Some(i) => format!("{}<style>{SVG_CSS}</style>{}", &svg[..=i], &svg[i + 1..]),
        None => svg,
    }
}

/// The report as Markdown text and its figures (file name, standalone SVG). `pdf` writes the subset the PDF layout understands:
/// `* ` bullets, no italics and no block quotes.
fn md_text(blocks: &[Block], pdf: bool) -> (String, Vec<(String, String)>) {
    let mut m = String::new();
    let mut figs: Vec<(String, String)> = Vec::new();
    for bl in blocks {
        match bl {
            Block::H(l, s) => {
                let _ = write!(m, "\n{} {}\n\n", "#".repeat(*l as usize), s);
            }
            Block::P(v, tone) => {
                let text = inl_md(v);
                let _ = match tone {
                    Tone::Note if !pdf => writeln!(m, "*{}*\n", text.trim()),
                    Tone::Warn | Tone::Bad => writeln!(m, "**Note:** {text}\n"),
                    _ => writeln!(m, "{text}\n"),
                };
            }
            Block::Card(v) => {
                let _ = if pdf {
                    writeln!(m, "{}\n", inl_md(v))
                } else {
                    writeln!(m, "> {}\n", inl_md(v))
                };
            }
            Block::Table(head, rows) => {
                let _ = writeln!(m, "| {} |", head.join(" | "));
                let _ = writeln!(
                    m,
                    "|{}|",
                    head.iter().map(|_| "---").collect::<Vec<_>>().join("|")
                );
                for row in rows {
                    let _ = writeln!(
                        m,
                        "| {} |",
                        row.iter()
                            .map(|c| inl_md(c))
                            .collect::<Vec<_>>()
                            .join(" | ")
                    );
                }
                m.push('\n');
            }
            Block::List(items) => {
                for it in items {
                    let _ = writeln!(m, "{} {}", if pdf { "*" } else { "-" }, inl_md(it));
                }
                m.push('\n');
            }
            Block::Fig(name, svgs) => {
                for (i, s) in svgs.iter().enumerate() {
                    let file = if svgs.len() > 1 {
                        format!("{name}_{}.svg", i + 1)
                    } else {
                        format!("{name}.svg")
                    };
                    figs.push((file.clone(), standalone_svg(s)));
                    let _ = writeln!(m, "![{name}](figures/{file})\n");
                }
            }
            Block::Legend(items) => {
                let _ = writeln!(
                    m,
                    "Colours: {}\n",
                    items
                        .iter()
                        .map(|(n, ci)| format!("{n} ({})", PALETTE[ci % PALETTE.len()]))
                        .collect::<Vec<_>>()
                        .join("; ")
                );
            }
        }
    }
    (m, figs)
}

#[cfg(feature = "native")]
fn render_md(blocks: &[Block], dir: &std::path::Path) -> Result<String> {
    fs::create_dir_all(dir.join("figures"))?;
    let (m, figs) = md_text(blocks, false);
    for (file, svg) in figs {
        fs::write(dir.join("figures").join(file), svg)?;
    }
    Ok(m)
}

#[cfg(feature = "native")]
fn render_pdf(blocks: &[Block]) -> Result<Vec<u8>> {
    let (md, figures) = md_text(blocks, true);
    gmg_pdf::render_pdf(&gmg_report::Markdown { md, figures }, "popsketch report")
        .map_err(|e| anyhow::anyhow!("PDF: {e}"))
}

fn tone_block(v: Vec<In>, tone: Tone) -> Block {
    Block::P(v, tone)
}

/// The report as blocks, shared by the HTML and Markdown renderers.
fn build(o: &Options, r: &Results) -> Vec<Block> {
    let mut d: Vec<Block> = Vec::new();
    let n_total = r.rows.len();
    let count = |s: Status| r.rows.iter().filter(|x| x.status == s).count();
    d.push(Block::H(
        1,
        "popsketch: sample audit and coarse structure".into(),
    ));
    d.push(tone_block(
        vec![t(format!(
            "k={}, scaled={}, {} per sample, flagged samples: {:?}. This tool compares sketches, not alignments: it audits a set of read files and shows coarse structure; it does not build a phylogeny.",
            o.ksize,
            o.scaled,
            if o.max_bases == 0 { "no base limit".to_string() } else { format!("up to {} bases", o.max_bases) },
            o.flagged
        ))],
        Tone::Note,
    ));
    // ---- audit
    d.push(Block::H(2, "1. Audit".into()));
    d.push(Block::Card(vec![
        b(n_total.to_string()),
        t(" samples: "),
        b(r.included.len().to_string()),
        t(format!(
            " in the analysis ({} included, {} marked, {} decontaminated), ",
            count(Status::Included),
            count(Status::Marked),
            count(Status::Decontaminated)
        )),
        b(count(Status::Excluded).to_string()),
        t(" excluded, "),
        b(count(Status::Error).to_string()),
        t(" unreadable."),
    ]));
    let mut by_flag: BTreeMap<&str, Vec<&str>> = BTreeMap::new();
    for row in &r.rows {
        for f in row.flags() {
            by_flag.entry(f.code).or_default().push(&row.id);
        }
    }
    if by_flag.is_empty() {
        d.push(Block::P(vec![t("No sample was flagged.")], Tone::Plain));
    } else {
        let rows = by_flag
            .iter()
            .map(|(code, ids)| {
                let shown: Vec<&str> = ids.iter().take(12).copied().collect();
                let more = if ids.len() > 12 {
                    format!(" and {} more", ids.len() - 12)
                } else {
                    String::new()
                };
                vec![
                    vec![c(*code)],
                    vec![t(ids.len().to_string())],
                    vec![t(format!("{}{more}", shown.join(", ")))],
                ]
            })
            .collect();
        d.push(Block::Table(vec!["flag", "samples", "which"], rows));
        d.push(Block::H(3, "Flagged samples, with the reason".into()));
        let rows = r
            .rows
            .iter()
            .filter(|x| x.qc.as_ref().is_some_and(|q| q.flagged()))
            .map(|row| {
                let mut why: Vec<In> = Vec::new();
                for (i, f) in row.flags().iter().enumerate() {
                    if i > 0 {
                        why.push(In::Br);
                    }
                    why.push(c(f.code));
                    why.push(t(format!(" {}", f.detail)));
                }
                vec![vec![t(row.id.clone())], vec![t(row.status.as_str())], why]
            })
            .collect();
        d.push(Block::Table(vec!["sample", "status", "reasons"], rows));
    }
    for row in r.rows.iter().filter(|x| x.status == Status::Error) {
        d.push(Block::P(
            vec![t(format!(
                "{}: {}",
                row.id,
                row.message.as_deref().unwrap_or("")
            ))],
            Tone::Bad,
        ));
    }
    let noted: Vec<&crate::pipeline::SampleRow> = r
        .rows
        .iter()
        .filter(|x| x.qc.as_ref().is_some_and(|q| !q.notes.is_empty()))
        .collect();
    if !noted.is_empty() {
        d.push(Block::H(3, "Notes (these samples are not excluded)".into()));
        let rows = noted
            .iter()
            .take(40)
            .map(|row| {
                let mut why: Vec<In> = Vec::new();
                for (i, n) in row.qc.iter().flat_map(|q| &q.notes).enumerate() {
                    if i > 0 {
                        why.push(In::Br);
                    }
                    why.push(c(n.code));
                    why.push(t(format!(" {}", n.detail)));
                }
                vec![vec![t(row.id.clone())], why]
            })
            .collect();
        d.push(Block::Table(vec!["sample", "note"], rows));
        if noted.len() > 40 {
            d.push(Block::P(
                vec![t(format!("and {} more, in qc.tsv", noted.len() - 40))],
                Tone::Note,
            ));
        }
    }
    d.push(Block::P(
        vec![
            b("What the audit cannot see."),
            t(" Light contamination (about 1 to 7% of the reads) leaves no mark in the abundance spectrum or the size: a contaminant at 1x next to a 40x genome looks like sequencing error. Heavy contamination is caught by size (low coverage) or by an abundance shoulder (a contaminant at several-fold coverage). Giving an index ("),
            c("--index"),
            t(") adds the index's own mixture and organism calls."),
        ],
        Tone::Note,
    ));
    // size vs coverage
    let pts: Vec<Pt> = r
        .rows
        .iter()
        .filter_map(|row| {
            let q = row.qc.as_ref()?;
            let cov = q.kmer_cov?;
            Some(Pt {
                x: cov.max(0.1),
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
                    cov,
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
        d.push(Block::H(3, "Estimated size against k-mer coverage".into()));
        d.push(Block::Fig(
            "size_vs_coverage",
            vec![scatter(
                &pts,
                "k-mer coverage (log scale)",
                "estimated size (Mb) = hashes above the cutoff x scaled",
                true,
            )],
        ));
        d.push(Block::Legend(vec![
            ("in the analysis".into(), 0),
            ("excluded (cross)".into(), 4),
            ("marked or decontaminated".into(), 1),
        ]));
    }
    // labels
    for rep in &r.label_reports {
        d.push(Block::H(2, format!("Labels: {}", rep.name)));
        d.push(Block::P(
            vec![t(format!(
                "{} labelled samples in the analysis; groups: {}. A blank label means unlabelled.",
                rep.n_labelled,
                rep.groups
                    .iter()
                    .map(|(g, n)| format!("{g} ({n})"))
                    .collect::<Vec<_>>()
                    .join(", ")
            ))],
            Tone::Note,
        ));
        if !rep.unchecked.is_empty() {
            d.push(Block::P(vec![t(format!("Not checked for conflicts, because these labels do not form coherent clusters in the sketches (members are not more similar to each other than to the rest): {}.", rep.unchecked.join(", ")))], Tone::Note));
        }
        if rep.conflicts.is_empty() {
            d.push(Block::P(
                vec![t("No sample matches another label better than its own.")],
                Tone::Plain,
            ));
        } else {
            let rows = rep
                .conflicts
                .iter()
                .map(|x| {
                    vec![
                        vec![t(x.id.clone())],
                        vec![t(x.label.clone())],
                        vec![t(format!("{:.4}", x.own_median_ani))],
                        vec![b(x.other_label.clone())],
                        vec![t(format!("{:.4}", x.other_median_ani))],
                        vec![t(format!(
                            "{} ({}, {:.4})",
                            x.nearest,
                            x.nearest_label.as_deref().unwrap_or("unlabelled"),
                            x.nearest_ani
                        ))],
                    ]
                })
                .collect();
            d.push(Block::Table(
                vec![
                    "sample",
                    "label",
                    "median ANI to own label",
                    "better match",
                    "median ANI there",
                    "nearest sample",
                ],
                rows,
            ));
            d.push(Block::P(vec![t("A conflict means the sample's ANI to the members of another label is higher than to its own label by more than the noise floor. It is evidence to check the sample's metadata, not proof of an error.")], Tone::Note));
        }
        if !rep.ari.is_empty() {
            let rows = rep
                .ari
                .iter()
                .map(|(k, a)| vec![vec![t(k.to_string())], vec![t(format!("{a:.3}"))]])
                .collect();
            d.push(Block::Table(
                vec!["clusters", "adjusted Rand index vs this label"],
                rows,
            ));
        }
    }
    // does each label explain the distances?
    if !r.group_tests.is_empty() {
        d.push(Block::H(2, "Do the labels explain the distances?".into()));
        let mut rows = Vec::new();
        let mut notes: Vec<String> = Vec::new();
        for g in &r.group_tests {
            let fp = |x: Option<f64>| x.map_or("n/a".to_string(), |v| format!("{v:.3}"));
            match &g.result {
                Some(res) => {
                    rows.push(vec![
                        vec![t(g.column.clone())],
                        vec![t(res.n.to_string())],
                        vec![t(res.n_groups.to_string())],
                        vec![t(format!("{:.3}", res.r2))],
                        vec![t(fp(res.p))],
                        vec![t(fp(res.disp_p))],
                        vec![t(g.strata.clone().unwrap_or_else(|| "none".to_string()))],
                    ]);
                    for w in &res.warnings {
                        notes.push(format!("{}: {w}", g.column));
                    }
                }
                None => {
                    rows.push(vec![
                        vec![t(g.column.clone())],
                        vec![t("0")],
                        vec![t("0")],
                        vec![t("n/a")],
                        vec![t("n/a")],
                        vec![t("n/a")],
                        vec![t(g.strata.clone().unwrap_or_else(|| "none".to_string()))],
                    ]);
                    notes.push(format!(
                        "{}: fewer than two groups with two or more samples",
                        g.column
                    ));
                }
            }
        }
        d.push(Block::Table(
            vec![
                "label",
                "samples",
                "groups",
                "R2",
                "PERMANOVA p",
                "PERMDISP p",
                "permuted inside",
            ],
            rows,
        ));
        let perms = r.group_tests.first().map_or(0, |g| g.permutations);
        d.push(Block::P(vec![t(format!("R2 is the share of the squared distances the label explains. PERMANOVA p asks whether that is more than chance when the labels are shuffled ({perms} permutations, among the samples of one stratum when a stratum column is given). PERMDISP p asks whether the groups differ in spread; if it is small, a small PERMANOVA p may reflect spread and not different centres."))], Tone::Note));
        for w in notes {
            d.push(Block::P(vec![t(w)], Tone::Warn));
        }
    }
    // ---- structure
    if let Some(c) = &r.collapse {
        d.push(Block::P(
            vec![t(format!(
                "Near-identical samples were collapsed (rule: {}): {} samples were reduced to {} representatives, and the structure, tree, clusters and label tests below use the representatives (the label tests weight each by the number of samples it stands for). Which sample stands for which is in collapse.tsv.",
                c.rule, c.input, c.kept
            ))],
            Tone::Note,
        ));
    }
    d.push(Block::H(2, "2. Coarse structure (secondary)".into()));
    if let (Some(dm), Some(ord)) = (&r.dist, &r.ordination) {
        let n = dm.n;
        let names: Vec<String> = r.included.iter().map(|&i| r.rows[i].id.clone()).collect();
        // colour by the first label column when present, else by cluster
        let (colors, items): (Vec<usize>, Vec<(String, usize)>) =
            if let (Some(tb), Some(_)) = (&r.label_table, r.label_reports.first()) {
                let labs = crate::labels::column_labels(tb, &names, 0);
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
                    cl.iter().map(|&x| x % 8).collect(),
                    (0..k)
                        .map(|x| (format!("cluster {}", x + 1), x % 8))
                        .collect(),
                )
            } else {
                (vec![0; n], vec![("samples".into(), 0)])
            };
        let color_of = |i: usize| colors[i];
        let by = if r.label_table.is_some() {
            format!(
                "label column {}",
                r.label_reports.first().map_or("", |x| x.name.as_str())
            )
        } else {
            "cluster".to_string()
        };
        let shape = |i: usize| match r.rows[r.included[i]].status {
            Status::Marked => 1,
            Status::Decontaminated => 2,
            _ => 0,
        };
        let pc = |a: usize, bb: usize| -> String {
            let pts: Vec<Pt> = (0..n)
                .map(|i| Pt {
                    x: ord.coords[i][a],
                    y: ord.coords[i][bb],
                    color: color_of(i),
                    shape: shape(i),
                    title: format!(
                        "{} (cluster {})",
                        names[i],
                        r.clusters
                            .as_ref()
                            .map_or("NA".to_string(), |x| (x[i] + 1).to_string())
                    ),
                })
                .collect();
            scatter(
                &pts,
                &format!("PC{} ({:.1}% of variance)", a + 1, ord.var_frac[a] * 100.0),
                &format!("PC{} ({:.1}%)", bb + 1, ord.var_frac[bb] * 100.0),
                false,
            )
        };
        d.push(Block::P(
            vec![
                t(format!(
                    "{} on {} {}; coloured by {by}. Triangles are flagged samples kept with ",
                    if ord.method == "hash_pca" {
                        "Hash-presence PCA"
                    } else {
                        "PCoA of 1 - ANI"
                    },
                    ord.n_features,
                    if ord.method == "hash_pca" {
                        "hashes"
                    } else {
                        "samples"
                    }
                )),
                c("--flagged mark"),
                t(", diamonds are decontaminated."),
            ],
            Tone::Note,
        ));
        let mut figs = Vec::new();
        if ord.eigenvalues.len() >= 2 {
            figs.push(pc(0, 1));
        }
        if ord.eigenvalues.len() >= 3 {
            figs.push(pc(1, 2));
        }
        d.push(Block::Fig("pca", figs));
        d.push(Block::Legend(items));
        // DAPC on the primary label column
        if let Some(dp) = &r.dapc {
            let axes = dp.dapc.eigenvalues.len();
            let pts: Vec<Pt> = (0..n)
                .map(|i| Pt {
                    x: dp.dapc.scores[i][0],
                    y: if axes >= 2 {
                        dp.dapc.scores[i][1]
                    } else {
                        ((i * 37) % 11) as f64 / 11.0
                    },
                    color: color_of(i),
                    shape: shape(i),
                    title: format!(
                        "{} (label {}, assigned {})",
                        names[i],
                        dp.labels[i].as_deref().unwrap_or("none"),
                        dp.dapc.assigned[i]
                            .map_or("n/a".to_string(), |g| dp.dapc.groups[g].clone())
                    ),
                })
                .collect();
            d.push(Block::H(3, format!("DAPC on the labels of {}", dp.column)));
            let elsewhere: Vec<usize> = (0..n)
                .filter(|&i| {
                    dp.dapc.assigned[i].is_some()
                        && dp.labels[i].as_deref()
                            != dp.dapc.assigned[i].map(|g| dp.dapc.groups[g].as_str())
                })
                .collect();
            d.push(Block::P(
                vec![t(format!(
                    "Discriminant analysis of the first {} PCs with the groups of this label column that form coherent clusters in the sketches: {} axes{}. {:.1}% of the samples in those groups are nearest to their own group's centroid; {} are nearest to another group{}",
                    dp.dapc.n_pcs,
                    axes,
                    if axes == 1 { " (the second axis of the plot is jitter)" } else { "" },
                    dp.dapc.accuracy * 100.0,
                    elsewhere.len(),
                    if elsewhere.is_empty() { "." } else { " (a second hint, next to the ANI-based conflict check, that the label may be wrong):" }
                ))],
                Tone::Note,
            ));
            d.push(Block::Fig(
                "dapc",
                vec![scatter(
                    &pts,
                    &format!("LD1 (eigenvalue {:.1})", dp.dapc.eigenvalues[0]),
                    if axes >= 2 { "LD2" } else { "jitter" },
                    false,
                )],
            ));
            if !elsewhere.is_empty() {
                let rows = elsewhere
                    .iter()
                    .take(40)
                    .map(|&i| {
                        vec![
                            vec![t(names[i].clone())],
                            vec![t(dp.labels[i].clone().unwrap_or_default())],
                            vec![b(dp.dapc.assigned[i]
                                .map_or(String::new(), |g| dp.dapc.groups[g].clone()))],
                        ]
                    })
                    .collect();
                d.push(Block::Table(
                    vec!["sample", "label", "nearest group in discriminant space"],
                    rows,
                ));
            }
        }
        // clusters
        if let Some(g) = &r.gap {
            let (say, tone): (Vec<In>, Tone) = match &g.verdict {
                Verdict::Clear(k) => (
                    vec![
                        t("The largest jump in the clustering gives "),
                        b(format!("{k} clusters")),
                        t(format!(
                            ": it is {:.1} times the noise floor. A structureless population would show a gap ratio this large with probability {} (ratio {:.1}), so this is evidence of a split larger than the technical noise, not proof of distinct populations.",
                            g.gap_over_floor,
                            if g.null_p < 0.001 { "under 0.001".to_string() } else { format!("{:.3}", g.null_p) },
                            g.ratio
                        )),
                    ],
                    Tone::Plain,
                ),
                Verdict::Single => (vec![t(format!("Everything merges below the noise floor (the last merge is at {:.5}): one population as far as sketches can tell.", g.top_height))], Tone::Plain),
                Verdict::NoClearCut => (
                    vec![
                        b("No single clear cut."),
                        t(format!(
                            " The largest jump is {:.1} times the noise floor (a cut needs to be at least {:.1}). Candidate cuts: {}. Use ",
                            g.gap_over_floor,
                            o.min_floors,
                            g.candidates.iter().map(|x| format!("{} clusters (jump {:.5})", x.k, x.gap)).collect::<Vec<_>>().join("; ")
                        )),
                        c("--k"),
                        t(" to choose one."),
                    ],
                    Tone::Warn,
                ),
                Verdict::TooFew => (vec![t("Too few samples to cluster.")], Tone::Warn),
            };
            let mut say = say;
            if o.k.is_some() || o.cut.is_some() {
                say.push(t(format!(
                    " Clusters shown use {} clusters, as requested.",
                    r.k_used.unwrap_or(0)
                )));
            }
            d.push(Block::H(3, "Clusters".into()));
            d.push(Block::P(say, tone));
            if g.top_height >= DEEP_SPLIT {
                d.push(Block::P(
                    vec![
                        b("The set spans distant genomes."),
                        t(format!(
                            " The deepest merge is at 1 - ANI = {:.3} (ANI {:.1}%), beyond what one species usually shows. The cohort-relative checks (size, unrelated to the cohort, abundance shoulder) assume most samples are one organism, so a small distant group is not flagged by them, and a large mixed set can have its minority species excluded. Check the clusters against the labels, or run one species or genus at a time.",
                            g.top_height,
                            100.0 * (1.0 - g.top_height)
                        )),
                    ],
                    Tone::Warn,
                ));
            }
        }
        // distance histogram
        let mut vals = Vec::new();
        for i in 0..n {
            for j in (i + 1)..n {
                vals.push(1.0 - dm.ani(i, j));
            }
        }
        d.push(Block::H(3, "Pairwise distances".into()));
        d.push(Block::Fig(
            "distances",
            vec![histogram(&vals, r.floor.value)],
        ));
        // tree
        if r.newick.is_some() {
            let tr = {
                let mut tr = crate::tree::neighbor_joining(&dm.mash, n, names.clone());
                crate::tree::collapse(&mut tr, o.tree_floor.unwrap_or(r.floor.value));
                root_tree(&tr, o.root.as_deref())
            };
            d.push(Block::H(3, "Tree".into()));
            d.push(Block::P(
                vec![t(format!(
                    "Neighbour-joining on mash distance; {} branches shorter than the noise floor were collapsed, so unresolved parts are drawn as polytomies. {} Tips are coloured as above. This is a first-pass clustering tree, not a phylogeny: recombination, accessory content and uneven coverage distort it.",
                    r.collapsed,
                    if o.root.is_some() { "Rooted on the named sample." } else { "Unrooted (drawn from an arbitrary node)." }
                ))],
                Tone::Note,
            ));
            d.push(Block::Fig(
                "tree",
                vec![draw_tree(&tr, &names, &color_of, 120)],
            ));
            d.push(Block::P(
                vec![t("The full tree is in "), c("tree.nwk"), t(".")],
                Tone::Note,
            ));
        }
        // advisor
        if let Some(a) = &r.advice {
            d.push(Block::H(
                2,
                "3. Resolution: what these data can and cannot say".into(),
            ));
            d.push(Block::List(
                a.statements.iter().map(|s| vec![t(s.clone())]).collect(),
            ));
        }
        // comparison with an external tree and PCs
        if r.tree_comparison.is_some() || r.pcs_comparison.is_some() {
            d.push(Block::H(
                2,
                "4. Comparison with the supplied SNP results".into(),
            ));
        }
        if let Some(cmp) = &r.tree_comparison {
            d.push(Block::H(3, "Tree".into()));
            d.push(Block::P(
                vec![t(format!(
                    "{} samples are in both the run and the tree ({} tips). Spearman correlation of the tree's patristic distance with the sketch distance: {:.2} for all pairs, {:.2} for the closest half of the pairs (where sketches have the least to say).",
                    cmp.shared, cmp.tree_tips, cmp.spearman_all, cmp.spearman_closest_half
                ))],
                Tone::Note,
            ));
            if !cmp.clusters.is_empty() {
                let rows = cmp
                    .clusters
                    .iter()
                    .map(|x| {
                        vec![
                            vec![t(x.cluster.to_string())],
                            vec![t(x.n.to_string())],
                            vec![t(if x.exact_clade { "yes" } else { "no" })],
                            vec![t(format!("{:.3}", x.best_f1))],
                        ]
                    })
                    .collect();
                d.push(Block::Table(
                    vec![
                        "cluster",
                        "samples in the tree",
                        "is a clade of the tree",
                        "best match (F1) with any clade or its complement",
                    ],
                    rows,
                ));
            }
        }
        if let Some(cmp) = &r.pcs_comparison {
            let mut s = format!("{} samples are in both. Cosines of the principal angles between the top-m subspaces (1 = the same plane): ", cmp.shared);
            for (m, v) in &cmp.cosines {
                let _ = write!(
                    s,
                    "top {m}: {}; ",
                    v.iter()
                        .map(|x| format!("{x:.2}"))
                        .collect::<Vec<_>>()
                        .join(", ")
                );
            }
            let _ = write!(s, "Procrustes disparity of the top 3 (0 = identical up to rotation): {}. Spearman correlation of pairwise distances in the top-5 PC space: {:.2}.", cmp.procrustes_top3.map_or("n/a".to_string(), |x| format!("{x:.2}")), cmp.distance_spearman_top5);
            d.push(Block::H(3, "Principal components".into()));
            d.push(Block::P(vec![t(s)], Tone::Note));
        }
    } else if o.no_pca && r.dist.is_some() {
        d.push(Block::P(vec![t("The ordination was switched off (--no-pca); distances and clusters are in the output files.")], Tone::Note));
    } else {
        d.push(Block::P(
            vec![t(
                "Too few samples passed the audit for distances, ordination or a tree.",
            )],
            Tone::Warn,
        ));
    }
    if !r.warnings.is_empty() {
        d.push(Block::H(2, "Warnings".into()));
        d.push(Block::List(
            r.warnings.iter().map(|w| vec![t(w.clone())]).collect(),
        ));
    }
    d.push(Block::P(
        vec![
            t("Files: "),
            c("qc.tsv"),
            t(" (every sample, flags with reasons), "),
            c("ani.tsv"),
            t(", "),
            c("mash.tsv"),
            t(", "),
            c("pca.tsv"),
            t(", "),
            c("clusters.tsv"),
            t(", "),
            c("tree.nwk"),
            t(", "),
            c("label_check.tsv"),
            t(", "),
            c("popsketch.json"),
            t("."),
        ],
        Tone::Note,
    ));
    d
}

#[cfg(feature = "native")]
/// Write `report.html`, `report.md` (with `figures/*.svg`) and/or `report.pdf` as `Options::report_format` asks (`All` writes all three).
pub fn write(o: &Options, r: &Results) -> Result<()> {
    let blocks = build(o, r);
    if matches!(o.report_format, ReportFormat::Html | ReportFormat::All) {
        fs::write(o.out.join("report.html"), render_html(&blocks))?;
    }
    if matches!(o.report_format, ReportFormat::Md | ReportFormat::All) {
        fs::write(o.out.join("report.md"), render_md(&blocks, &o.out)?)?;
    }
    if matches!(o.report_format, ReportFormat::Pdf | ReportFormat::All) {
        fs::write(o.out.join("report.pdf"), render_pdf(&blocks)?)?;
    }
    Ok(())
}
