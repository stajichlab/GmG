//! PDF rendering of a GmG report, in pure Rust (krilla + usvg). It lays out the Markdown report
//! (`gmg_report::render_markdown`): headings, paragraphs, bullet lists, tables and the SVG figures, on A4 pages with
//! Inter embedded (SIL OFL, see `assets/fonts/`), so text is selectable and no font is missing.
//!
//! Layout is two-phase: `layout` turns the Markdown into positioned drawing operations per page, `render_pdf` draws them.
//! Only the Markdown subset the report emits is understood: `#`..`###` headings, paragraphs, `* ` bullets, pipe tables,
//! `![alt](figures/x.svg)` images, `**bold**` and `` `code` `` (shown as plain text).

use gmg_report::Markdown;
use krilla::color::rgb;
use krilla::geom::{PathBuilder, Point, Rect, Size, Transform};
use krilla::metadata::Metadata;
use krilla::page::PageSettings;
use krilla::paint::Fill;
use krilla::text::{Font, TextDirection};
use krilla::Document;
use krilla_svg::{SurfaceExt, SvgSettings};
use std::collections::HashMap;
use std::sync::Arc;
use ttf_parser::Face;

const REGULAR: &[u8] = include_bytes!("../assets/fonts/Inter-Regular.ttf");
const BOLD: &[u8] = include_bytes!("../assets/fonts/Inter-Bold.ttf");

const PAGE_W: f32 = 595.0;
const PAGE_H: f32 = 842.0;
const MARGIN: f32 = 42.0;
const CONTENT_W: f32 = PAGE_W - 2.0 * MARGIN;
const BODY: f32 = 9.5;
const LEADING: f32 = 1.4;

type Color = (u8, u8, u8);
const INK: Color = (0x0b, 0x0b, 0x0b);
const MUTED: Color = (0x6b, 0x6a, 0x66);
const RULE: Color = (0xd9, 0xd8, 0xd1);
const HEAD_BG: Color = (0xf0, 0xef, 0xec);

enum Op {
    Text {
        x: f32,
        y: f32,
        size: f32,
        bold: bool,
        color: Color,
        text: String,
    },
    Rect {
        x: f32,
        y: f32,
        w: f32,
        h: f32,
        color: Color,
    },
    Svg {
        name: String,
        x: f32,
        y: f32,
        w: f32,
        h: f32,
    },
}

struct Metrics {
    regular: Face<'static>,
    bold: Face<'static>,
}

impl Metrics {
    fn new() -> Metrics {
        Metrics {
            regular: Face::parse(REGULAR, 0).expect("embedded Inter Regular"),
            bold: Face::parse(BOLD, 0).expect("embedded Inter Bold"),
        }
    }

    fn width(&self, text: &str, size: f32, bold: bool) -> f32 {
        let face = if bold { &self.bold } else { &self.regular };
        let upm = face.units_per_em() as f32;
        text.chars()
            .map(|c| {
                face.glyph_index(c)
                    .and_then(|g| face.glyph_hor_advance(g))
                    .unwrap_or((upm * 0.55) as u16) as f32
            })
            .sum::<f32>()
            * size
            / upm
    }

    /// Greedy word wrap; words longer than the line are broken.
    fn wrap(&self, text: &str, size: f32, bold: bool, max: f32) -> Vec<String> {
        let mut lines = Vec::new();
        let mut cur = String::new();
        for word in text.split_whitespace() {
            let mut word = word.to_string();
            while self.width(&word, size, bold) > max && word.chars().count() > 1 {
                // break an over-long word (a path, an accession) at the line width
                let mut cut = String::new();
                for c in word.chars() {
                    let mut t = cut.clone();
                    t.push(c);
                    if self.width(&t, size, bold) > max {
                        break;
                    }
                    cut = t;
                }
                if !cur.is_empty() {
                    lines.push(std::mem::take(&mut cur));
                }
                let n = cut.chars().count().max(1);
                lines.push(cut);
                word = word.chars().skip(n).collect();
            }
            let candidate = if cur.is_empty() {
                word.clone()
            } else {
                format!("{cur} {word}")
            };
            if self.width(&candidate, size, bold) <= max {
                cur = candidate;
            } else {
                if !cur.is_empty() {
                    lines.push(std::mem::take(&mut cur));
                }
                cur = word;
            }
        }
        if !cur.is_empty() {
            lines.push(cur);
        }
        lines
    }

    /// Shorten to fit one line, ending with an ellipsis.
    fn fit(&self, text: &str, size: f32, bold: bool, max: f32) -> String {
        if self.width(text, size, bold) <= max {
            return text.to_string();
        }
        let mut s: String = text.to_string();
        while !s.is_empty() && self.width(&format!("{s}…"), size, bold) > max {
            s.pop();
        }
        format!("{}…", s.trim_end())
    }
}

struct Layout {
    pages: Vec<Vec<Op>>,
    y: f32,
}

impl Layout {
    fn new() -> Layout {
        Layout {
            pages: vec![Vec::new()],
            y: MARGIN,
        }
    }
    fn page(&mut self) -> &mut Vec<Op> {
        self.pages.last_mut().expect("at least one page")
    }
    fn room(&self, h: f32) -> bool {
        self.y + h <= PAGE_H - MARGIN - 14.0
    }
    fn new_page(&mut self) {
        self.pages.push(Vec::new());
        self.y = MARGIN;
    }
    fn ensure(&mut self, h: f32) {
        if !self.room(h) {
            self.new_page();
        }
    }
    fn text(&mut self, x: f32, baseline: f32, size: f32, bold: bool, color: Color, text: String) {
        self.page().push(Op::Text {
            x,
            y: baseline,
            size,
            bold,
            color,
            text,
        });
    }
}

/// Strip the inline Markdown the report uses (`**bold**`, `` `code` ``).
fn plain(s: &str) -> String {
    s.replace("**", "").replace('`', "")
}

fn table_cells(line: &str) -> Vec<String> {
    let t = line.trim();
    let t = t.strip_prefix('|').unwrap_or(t);
    let t = t.strip_suffix('|').unwrap_or(t);
    let mut cells = Vec::new();
    let mut cur = String::new();
    let mut chars = t.chars().peekable();
    while let Some(c) = chars.next() {
        if c == '\\' && chars.peek() == Some(&'|') {
            cur.push('|');
            chars.next();
        } else if c == '|' {
            cells.push(plain(cur.trim()));
            cur.clear();
        } else {
            cur.push(c);
        }
    }
    cells.push(plain(cur.trim()));
    cells
}

/// A separator row such as `|---|---|`.
fn is_separator(line: &str) -> bool {
    line.trim()
        .chars()
        .all(|c| matches!(c, '|' | '-' | ':' | ' '))
}

fn svg_size(svg: &str) -> Option<(f32, f32)> {
    let vb = svg.split("viewBox=\"").nth(1)?.split('"').next()?;
    let v: Vec<f32> = vb
        .split_whitespace()
        .filter_map(|x| x.parse().ok())
        .collect();
    (v.len() == 4 && v[2] > 0.0 && v[3] > 0.0).then_some((v[2], v[3]))
}

/// Height of the figure that comes next after a heading (0 if the next block is not a figure).
fn follow(lines: &[&str], from: usize, figs: &HashMap<&str, &str>) -> f32 {
    let mut blocks = lines[from..]
        .iter()
        .map(|l| l.trim())
        .filter(|l| !l.is_empty());
    let mut extra = 0.0;
    let mut next = blocks.next();
    // a one-line paragraph between the heading and its figure stays with them
    if let Some(t) = next {
        if !t.starts_with(['#', '|', '!', '*']) {
            extra = BODY * LEADING * 2.0 + 6.0;
            next = blocks.next();
        }
    }
    let name = next
        .and_then(|t| t.strip_prefix("!["))
        .and_then(|r| r.split("](figures/").nth(1))
        .and_then(|r| r.split(')').next());
    name.and_then(|n| figs.get(n))
        .and_then(|svg| svg_size(svg))
        .map(|(w, h)| image_size(w, h).1 + 6.0 + extra)
        .unwrap_or(0.0)
}

fn layout(md: &Markdown, m: &Metrics) -> Layout {
    let figs: HashMap<&str, &str> = md
        .figures
        .iter()
        .map(|(n, s)| (n.as_str(), s.as_str()))
        .collect();
    let mut l = Layout::new();
    let lines: Vec<&str> = md.md.lines().collect();
    let mut i = 0;
    while i < lines.len() {
        let t = lines[i].trim();
        if t.is_empty() {
            i += 1;
        } else if let Some(h) = t.strip_prefix("### ") {
            heading(&mut l, m, h, 11.5, 12.0, follow(&lines, i + 1, &figs));
            i += 1;
        } else if let Some(h) = t.strip_prefix("## ") {
            heading(&mut l, m, h, 14.0, 16.0, follow(&lines, i + 1, &figs));
            i += 1;
        } else if let Some(h) = t.strip_prefix("# ") {
            heading(&mut l, m, h, 21.0, 0.0, 0.0);
            i += 1;
        } else if t.starts_with("![") {
            let name = t
                .split("](figures/")
                .nth(1)
                .and_then(|r| r.split(')').next())
                .unwrap_or("");
            if let Some(svg) = figs.get(name) {
                image(&mut l, name, svg);
            }
            i += 1;
        } else if t.starts_with('|') {
            let start = i;
            while i < lines.len() && lines[i].trim().starts_with('|') {
                i += 1;
            }
            table(&mut l, m, &lines[start..i]);
        } else if t.starts_with("* ") {
            while i < lines.len() && lines[i].trim().starts_with("* ") {
                let item = plain(lines[i].trim().trim_start_matches("* "));
                let wrapped = m.wrap(&item, BODY, false, CONTENT_W - 12.0);
                let h = wrapped.len() as f32 * BODY * LEADING;
                l.ensure(h);
                let y0 = l.y + BODY;
                l.text(MARGIN, y0, BODY, false, MUTED, "•".into());
                for (k, w) in wrapped.into_iter().enumerate() {
                    let y = l.y + BODY + k as f32 * BODY * LEADING;
                    l.text(MARGIN + 12.0, y, BODY, false, INK, w);
                }
                l.y += h + 3.0;
                i += 1;
            }
            l.y += 4.0;
        } else {
            let mut para = String::new();
            while i < lines.len() {
                let u = lines[i].trim();
                if u.is_empty() || u.starts_with(['#', '|', '!']) || u.starts_with("* ") {
                    break;
                }
                para.push_str(u);
                para.push(' ');
                i += 1;
            }
            for w in m.wrap(&plain(&para), BODY, false, CONTENT_W) {
                l.ensure(BODY * LEADING);
                let y = l.y + BODY;
                l.text(MARGIN, y, BODY, false, INK, w);
                l.y += BODY * LEADING;
            }
            l.y += 6.0;
        }
    }
    l
}

fn heading(l: &mut Layout, m: &Metrics, text: &str, size: f32, before: f32, follow: f32) {
    let lines = m.wrap(&plain(text), size, true, CONTENT_W);
    let h = before + lines.len() as f32 * size * 1.25 + 6.0;
    // keep a heading with what follows it (the next figure, or at least a few lines)
    l.ensure(h + follow.max(40.0));
    if l.y > MARGIN {
        l.y += before;
    }
    for w in lines {
        let y = l.y + size;
        l.text(MARGIN, y, size, true, INK, w);
        l.y += size * 1.25;
    }
    l.y += 6.0;
}

/// Charts keep a constant type size (12 px of the drawing = 7.5 pt) and never exceed the text width.
fn image_size(w0: f32, h0: f32) -> (f32, f32) {
    let scale = (CONTENT_W / w0).min(0.62);
    (w0 * scale, h0 * scale)
}

fn image(l: &mut Layout, name: &str, svg: &str) {
    let Some((w0, h0)) = svg_size(svg) else {
        return;
    };
    let (w, h) = image_size(w0, h0);
    l.ensure(h + 6.0);
    let (x, y) = (MARGIN, l.y);
    l.page().push(Op::Svg {
        name: name.to_string(),
        x,
        y,
        w,
        h,
    });
    l.y += h + 10.0;
}

fn table(l: &mut Layout, m: &Metrics, raw: &[&str]) {
    let rows: Vec<Vec<String>> = raw
        .iter()
        .filter(|r| !is_separator(r))
        .map(|r| table_cells(r))
        .collect();
    if rows.is_empty() {
        return;
    }
    let ncol = rows[0].len();
    let size = if ncol > 5 { 8.0 } else { BODY };
    let pad = 4.0;
    // column widths: content width capped, then scaled to the page
    let want: Vec<f32> = (0..ncol)
        .map(|c| {
            rows.iter()
                .enumerate()
                .map(|(ri, r)| m.width(r.get(c).map(String::as_str).unwrap_or(""), size, ri == 0))
                .fold(0.0, f32::max)
                + 2.0 * pad
        })
        .collect();
    let total: f32 = want.iter().sum();
    let widths = if total <= CONTENT_W {
        want
    } else {
        // water-filling: columns that need less than an equal share keep their natural width, the rest share what is left
        let mut widths = vec![0.0f32; ncol];
        let mut open: Vec<usize> = (0..ncol).collect();
        let mut remaining = CONTENT_W;
        while !open.is_empty() {
            let share = remaining / open.len() as f32;
            let (fit, over): (Vec<usize>, Vec<usize>) =
                open.iter().partition(|&&c| want[c] <= share);
            if fit.is_empty() {
                for c in over {
                    widths[c] = share;
                }
                break;
            }
            for c in fit {
                widths[c] = want[c];
                remaining -= want[c];
            }
            open = over;
        }
        widths
    };
    let row_h = size * LEADING + 5.0;
    // a short table stays on one page; a long one needs room for a header and a few rows to start
    let need = if rows.len() <= 12 { rows.len() } else { 4 };
    l.ensure(row_h * need as f32);
    for (ri, r) in rows.iter().enumerate() {
        let head = ri == 0;
        if !l.room(row_h * if head { 2.0 } else { 1.0 }) {
            l.new_page();
            // repeat the header on the next page
            if !head {
                draw_row(l, m, &rows[0], &widths, size, row_h, true, pad);
            }
        }
        draw_row(l, m, r, &widths, size, row_h, head, pad);
    }
    l.y += 10.0;
}

#[allow(clippy::too_many_arguments)]
fn draw_row(
    l: &mut Layout,
    m: &Metrics,
    r: &[String],
    widths: &[f32],
    size: f32,
    row_h: f32,
    head: bool,
    pad: f32,
) {
    let y = l.y;
    let total: f32 = widths.iter().sum();
    if head {
        l.page().push(Op::Rect {
            x: MARGIN,
            y,
            w: total,
            h: row_h,
            color: HEAD_BG,
        });
    }
    let mut x = MARGIN;
    for (c, w) in widths.iter().enumerate() {
        let cell = r.get(c).map(String::as_str).unwrap_or("");
        let shown = m.fit(cell, size, head, w - 2.0 * pad);
        let numeric = !head && cell.parse::<f64>().is_ok();
        let tx = if numeric {
            x + w - pad - m.width(&shown, size, false)
        } else {
            x + pad
        };
        l.text(tx, y + 2.5 + size, size, head, INK, shown);
        x += w;
    }
    l.page().push(Op::Rect {
        x: MARGIN,
        y: y + row_h - 0.5,
        w: total,
        h: 0.5,
        color: RULE,
    });
    l.y += row_h;
}

fn fill((r, g, b): Color) -> Fill {
    Fill {
        paint: rgb::Color::new(r, g, b).into(),
        ..Default::default()
    }
}

/// Render the report as a PDF. Returns the file bytes.
pub fn render_pdf(md: &Markdown, title: &str) -> Result<Vec<u8>, String> {
    let m = Metrics::new();
    let l = layout(md, &m);

    let regular =
        Font::new(Arc::new(REGULAR.to_vec()).into(), 0).ok_or("cannot load Inter Regular")?;
    let bold = Font::new(Arc::new(BOLD.to_vec()).into(), 0).ok_or("cannot load Inter Bold")?;
    let mut db = usvg::fontdb::Database::new();
    db.load_font_data(REGULAR.to_vec());
    db.load_font_data(BOLD.to_vec());
    db.set_sans_serif_family("Inter");
    db.set_serif_family("Inter");
    db.set_monospace_family("Inter");
    let opts = usvg::Options {
        fontdb: Arc::new(db),
        font_family: "Inter".to_string(),
        ..Default::default()
    };
    let mut trees: HashMap<&str, usvg::Tree> = HashMap::new();
    for (name, svg) in &md.figures {
        if let Ok(t) = usvg::Tree::from_str(svg, &opts) {
            trees.insert(name.as_str(), t);
        }
    }

    let mut doc = Document::new();
    doc.set_metadata(
        Metadata::new()
            .title(title.to_string())
            .creator("gmg".to_string()),
    );
    let n = l.pages.len();
    for (pi, ops) in l.pages.iter().enumerate() {
        let mut page =
            doc.start_page_with(PageSettings::from_wh(PAGE_W, PAGE_H).ok_or("bad page size")?);
        let mut s = page.surface();
        for op in ops {
            match op {
                Op::Text {
                    x,
                    y,
                    size,
                    bold: b,
                    color,
                    text,
                } => {
                    s.set_fill(Some(fill(*color)));
                    s.draw_text(
                        Point::from_xy(*x, *y),
                        if *b { bold.clone() } else { regular.clone() },
                        *size,
                        text,
                        false,
                        TextDirection::Auto,
                    );
                }
                Op::Rect { x, y, w, h, color } => {
                    if let Some(r) = Rect::from_xywh(*x, *y, *w, *h) {
                        let mut pb = PathBuilder::new();
                        pb.push_rect(r);
                        if let Some(p) = pb.finish() {
                            s.set_fill(Some(fill(*color)));
                            s.draw_path(&p);
                        }
                    }
                }
                Op::Svg { name, x, y, w, h } => {
                    if let (Some(tree), Some(size)) =
                        (trees.get(name.as_str()), Size::from_wh(*w, *h))
                    {
                        s.push_transform(&Transform::from_translate(*x, *y));
                        s.draw_svg(tree, size, SvgSettings::default());
                        s.pop();
                    }
                }
            }
        }
        s.set_fill(Some(fill(MUTED)));
        s.draw_text(
            Point::from_xy(MARGIN, PAGE_H - MARGIN + 14.0),
            regular.clone(),
            7.5,
            &format!("GmG · Guess My Genome · page {} of {n}", pi + 1),
            false,
            TextDirection::Auto,
        );
        s.finish();
        page.finish();
    }
    doc.finish().map_err(|e| format!("{e:?}"))
}

/// Number of pages the report needs (for tests and size checks).
pub fn page_count(md: &Markdown) -> usize {
    layout(md, &Metrics::new()).pages.len()
}
