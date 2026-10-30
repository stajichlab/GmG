//! Validation report (#25): scores the rows of a `gmg run` results.tsv against a truth table and renders the standard tables and
//! figures (Markdown, long-format TSV, one HTML page). Pure functions, no I/O. The definitions are those of
//! `scripts/validate_run.py` (which it replaces for new work):
//! * `n` = inputs with a truth value at the rank (blank truth is not scored); `called` = share of those with a call;
//! * `precision` = correct among the calls made; `recall` = correct calls / `n` (an input without a call counts as missed);
//! * inputs with status `error` are not scored; results without a truth row are counted and dropped.

use crate::svg::{compact, esc, f, text};
use std::collections::{BTreeMap, HashMap};

pub type Row = HashMap<String, String>;

const RANKS: [&str; 6] = ["phylum", "class", "order", "family", "genus", "species"];
/// (rank, [(ANI bin label, n, precision)]) for the calibration chart
type Calib = Vec<(String, Vec<(String, usize, f64)>)>;
type Keyed = Vec<(String, Box<dyn Fn(&J) -> String>)>;
const CSS: &str = include_str!("../assets/report.css");

pub struct Options {
    /// truth/result columns to break the metrics down by (those that exist; blank = a group of its own)
    pub by: Vec<String>,
    pub title: String,
}

impl Default for Options {
    fn default() -> Self {
        Options {
            by: ["species_in_bfd", "genus_in_bfd", "level", "phylum"]
                .iter()
                .map(|s| s.to_string())
                .collect(),
            title: "Validation".into(),
        }
    }
}

pub struct Output {
    pub md: String,
    /// long format: group, rank, n, called, precision, recall
    pub tsv: String,
    pub html: String,
}

/// One result joined with its truth row.
struct J<'a> {
    res: &'a Row,
    tru: &'a Row,
}

impl J<'_> {
    /// results columns win over truth columns of the same name (as in the pandas merge of validate_run.py)
    fn get(&self, c: &str) -> &str {
        self.res
            .get(c)
            .or_else(|| self.tru.get(c))
            .map_or("", |s| s.as_str())
    }
    fn num(&self, c: &str) -> Option<f64> {
        self.get(c)
            .trim()
            .parse::<f64>()
            .ok()
            .filter(|v| v.is_finite())
    }
    fn truth(&self, rank: &str) -> &str {
        self.tru.get(rank).map_or("", |s| s.as_str())
    }
    fn call(&self, rank: &str) -> &str {
        self.res
            .get(&format!("call_{rank}"))
            .map_or("", |s| s.as_str())
    }
    fn has(&self, rank: &str) -> bool {
        !self.truth(rank).is_empty()
    }
    fn called(&self, rank: &str) -> bool {
        !self.call(rank).is_empty()
    }
    fn ok(&self, rank: &str) -> bool {
        self.call(rank) == self.truth(rank)
    }
}

#[derive(Clone, Debug)]
pub struct Metric {
    pub group: String,
    pub rank: &'static str,
    pub n: usize,
    pub called: f64,
    pub n_called: usize,
    pub precision: Option<f64>,
    pub recall: Option<f64>,
}

fn metrics(g: &[&J], rank: &'static str, label: &str) -> Option<Metric> {
    let t: Vec<&&J> = g.iter().filter(|j| j.has(rank)).collect();
    if t.is_empty() {
        return None;
    }
    let called: Vec<&&&J> = t.iter().filter(|j| j.called(rank)).collect();
    let right = called.iter().filter(|j| j.ok(rank)).count();
    Some(Metric {
        group: label.to_string(),
        rank,
        n: t.len(),
        called: called.len() as f64 / t.len() as f64,
        n_called: called.len(),
        precision: (!called.is_empty()).then(|| right as f64 / called.len() as f64),
        recall: Some(t.iter().filter(|j| j.ok(rank)).count() as f64 / t.len() as f64),
    })
}

fn table_of(g: &[&J], label: &str) -> Vec<Metric> {
    RANKS.iter().filter_map(|r| metrics(g, r, label)).collect()
}

/// Wilson 95% interval of a proportion `p` out of `n`.
pub fn wilson(p: f64, n: usize) -> (f64, f64) {
    if n == 0 {
        return (f64::NAN, f64::NAN);
    }
    let (n, z) = (n as f64, 1.96f64);
    let c = (p + z * z / (2.0 * n)) / (1.0 + z * z / n);
    let h = z * (p * (1.0 - p) / n + z * z / (4.0 * n * n)).sqrt() / (1.0 + z * z / n);
    (c - h, c + h)
}

/// A titled block of the report: some text and optionally a table.
struct Section {
    title: String,
    text: String,
    head: Vec<String>,
    rows: Vec<Vec<String>>,
}

fn rank_rows(ms: &[Metric]) -> Vec<Vec<String>> {
    ms.iter()
        .map(|m| {
            let p = match m.precision {
                Some(p) => {
                    let (lo, hi) = wilson(p, m.n_called);
                    format!("{p:.3} [{lo:.3}, {hi:.3}]")
                }
                None => "-".into(),
            };
            vec![
                m.rank.to_string(),
                m.n.to_string(),
                format!("{:.3}", m.called),
                p,
                format!("{:.3}", m.recall.unwrap_or(f64::NAN)),
            ]
        })
        .collect()
}

const RANK_HEAD: [&str; 5] = ["rank", "n", "called", "precision [95% CI]", "recall"];

fn head(cols: &[&str]) -> Vec<String> {
    cols.iter().map(|s| s.to_string()).collect()
}

/// Bin of a value: `edges` are upper bounds (a value `v` is in the first bin with `v <= edge`, like pandas `cut`).
fn cut_right(v: Option<f64>, edges: &[f64], labels: &[&str]) -> String {
    match v {
        None => "nan".into(),
        Some(v) => edges
            .iter()
            .position(|e| v <= *e)
            .map_or("nan", |i| labels[i])
            .to_string(),
    }
}

/// Group key list of a breakdown column (sorted, like a pandas groupby on strings).
fn groups<'a, 'b>(
    rows: &'b [&'b J<'a>],
    key: &dyn Fn(&J) -> String,
) -> BTreeMap<String, Vec<&'b J<'a>>> {
    let mut m: BTreeMap<String, Vec<&J>> = BTreeMap::new();
    for j in rows {
        m.entry(key(j)).or_default().push(*j);
    }
    m
}

pub fn score(results: &[Row], truth: &[Row], opts: &Options) -> Output {
    let tmap: HashMap<&str, &Row> = {
        let mut m = HashMap::new();
        for t in truth {
            if let Some(id) = t.get("id") {
                m.entry(id.as_str()).or_insert(t);
            }
        }
        m
    };
    let joined: Vec<J> = results
        .iter()
        .filter_map(|r| {
            r.get("query_id")
                .and_then(|q| tmap.get(q.as_str()))
                .map(|t| J { res: r, tru: t })
        })
        .collect();
    let unmatched = results.len() - joined.len();
    let scored: Vec<&J> = joined
        .iter()
        .filter(|j| j.get("status") != "error")
        .collect();
    let has_col = |c: &str| {
        joined
            .iter()
            .any(|j| j.res.contains_key(c) || j.tru.contains_key(c))
    };

    let mut status: BTreeMap<&str, usize> = BTreeMap::new();
    for r in results {
        *status
            .entry(r.get("status").map_or("", |s| s.as_str()))
            .or_default() += 1;
    }
    let mut by_count: Vec<(&&str, &usize)> = status.iter().collect();
    by_count.sort_by(|a, b| b.1.cmp(a.1).then(a.0.cmp(b.0)));
    let status_txt = by_count
        .iter()
        .map(|(k, v)| format!("{k} {v}"))
        .collect::<Vec<_>>()
        .join(", ");

    let mut all: Vec<Metric> = table_of(&scored, "all");
    let overall = all.clone();
    let mut sections: Vec<Section> = vec![Section {
        title: "All inputs".into(),
        text: String::new(),
        head: head(&RANK_HEAD),
        rows: rank_rows(&all),
    }];

    // breakdowns
    let mut cols: Keyed = Vec::new();
    for c in &opts.by {
        if has_col(c) {
            let c2 = c.clone();
            cols.push((c.clone(), Box::new(move |j: &J| j.get(&c2).to_string())));
        }
    }
    if has_col("n50") {
        cols.push((
            "n50_bin".into(),
            Box::new(|j: &J| {
                cut_right(
                    j.num("n50"),
                    &[1e4, 1e5, 1e6, 1e12],
                    &["<10 kb", "10-100 kb", "0.1-1 Mb", ">= 1 Mb"],
                )
            }),
        ));
    }
    if has_col("total_len") {
        cols.push((
            "size_bin".into(),
            Box::new(|j: &J| {
                cut_right(
                    j.num("total_len"),
                    &[1e7, 2.5e7, 5e7, 1e12],
                    &["<10 Mb", "10-25 Mb", "25-50 Mb", ">= 50 Mb"],
                )
            }),
        ));
    }
    for (name, key) in &cols {
        for (val, g) in groups(&scored, key.as_ref()) {
            let label = format!("{name}={val}");
            let ms = table_of(&g, &label);
            if !ms.is_empty() && g.len() >= 5 {
                sections.push(Section {
                    title: format!("{name} = {val} (n = {})", g.len()),
                    text: String::new(),
                    head: head(&RANK_HEAD),
                    rows: rank_rows(&ms),
                });
            }
            all.extend(ms);
        }
    }

    // taxa that are not in the reference: any call at the rank is wrong by construction
    if has_col("species_in_bfd") {
        let nov: Vec<&&J> = scored
            .iter()
            .filter(|j| j.get("species_in_bfd") == "False" && j.has("species"))
            .collect();
        if !nov.is_empty() {
            let cg = nov.iter().filter(|j| j.called("genus")).count();
            let og = nov
                .iter()
                .filter(|j| j.called("genus") && j.ok("genus"))
                .count();
            sections.push(Section {
                title: format!("Named species NOT in the reference (n = {})", nov.len()),
                text: format!(
                    "Species called for {} of {} (every one of those is wrong by construction: the species is not in the reference); genus called for {cg}, correct {og}.",
                    nov.iter().filter(|j| j.called("species")).count(),
                    nov.len()
                ),
                head: vec![],
                rows: vec![],
            });
        }
    }
    if has_col("genus_in_bfd") {
        let ng: Vec<&&J> = scored
            .iter()
            .filter(|j| j.get("genus_in_bfd") == "False" && j.has("genus"))
            .collect();
        if !ng.is_empty() {
            let part = |r: &str| {
                (
                    ng.iter().filter(|j| j.called(r)).count(),
                    ng.iter().filter(|j| j.called(r) && j.ok(r)).count(),
                )
            };
            let (f_, o_, p_) = (part("family"), part("order"), part("phylum"));
            sections.push(Section {
                title: format!("Genus NOT in the reference (n = {})", ng.len()),
                text: format!(
                    "Genus called for {} (wrong by construction); family called {}, correct {}; order called {}, correct {}; phylum called {}, correct {}.",
                    ng.iter().filter(|j| j.called("genus")).count(), f_.0, f_.1, o_.0, o_.1, p_.0, p_.1
                ),
                head: vec![],
                rows: vec![],
            });
        }
    }

    // confusion: the most frequent wrong calls (truth -> call)
    for r in ["genus", "species"] {
        let w: Vec<&&J> = scored
            .iter()
            .filter(|j| j.has(r) && j.called(r) && !j.ok(r))
            .collect();
        if w.is_empty() {
            continue;
        }
        let mut c: HashMap<(&str, &str), usize> = HashMap::new();
        for j in &w {
            *c.entry((j.truth(r), j.call(r))).or_default() += 1;
        }
        let mut c: Vec<_> = c.into_iter().collect();
        c.sort_by(|a, b| b.1.cmp(&a.1).then(a.0.cmp(&b.0)));
        let called_n = scored.iter().filter(|j| j.has(r) && j.called(r)).count();
        sections.push(Section {
            title: format!(
                "Most frequent wrong {r} calls ({} wrong of {called_n} called)",
                w.len()
            ),
            text: String::new(),
            head: head(&["truth", "called", "n"]),
            rows: c
                .iter()
                .take(15)
                .map(|((t, k), n)| vec![t.to_string(), k.to_string(), n.to_string()])
                .collect(),
        });
    }

    // calibration: precision by reported ANI and by vote consistency
    let mut calib: Calib = Vec::new();
    for r in ["genus", "species"] {
        let t: Vec<&&J> = scored.iter().filter(|j| j.has(r) && j.called(r)).collect();
        let cfg: [(String, Vec<f64>, String); 2] = [
            (
                "ani".into(),
                vec![0.0, 0.90, 0.95, 0.98, 0.995, 1.01],
                "top-hit ANI".into(),
            ),
            (
                format!("conf_{r}"),
                vec![-1.0, 0.6, 0.8, 0.9, 0.95, 1.01],
                format!("{r} vote consistency"),
            ),
        ];
        for (col, edges, lab) in cfg {
            if t.is_empty() || !has_col(&col) {
                continue;
            }
            let mut bins: Vec<(usize, usize)> = vec![(0, 0); edges.len() - 1];
            for j in &t {
                if let Some(v) = j.num(&col) {
                    if let Some(i) =
                        (0..edges.len() - 1).find(|&i| v >= edges[i] && v < edges[i + 1])
                    {
                        bins[i].0 += 1;
                        bins[i].1 += usize::from(j.ok(r));
                    }
                }
            }
            let mut rows = Vec::new();
            let mut pts = Vec::new();
            for (i, (n, k)) in bins.iter().enumerate() {
                if *n == 0 {
                    continue;
                }
                let p = *k as f64 / *n as f64;
                let b = format!("[{:?}, {:?})", edges[i], edges[i + 1]);
                rows.push(vec![b.clone(), n.to_string(), format!("{p:.3}")]);
                all.push(Metric {
                    group: format!("calib_{col}={b}"),
                    rank: if r == "genus" { "genus" } else { "species" },
                    n: *n,
                    called: 1.0,
                    n_called: *n,
                    precision: Some(p),
                    recall: None,
                });
                pts.push((b, *n, p));
            }
            if !rows.is_empty() {
                sections.push(Section {
                    title: format!("Calibration: {r} precision by {lab}"),
                    text: String::new(),
                    head: head(&[&lab, "n", "precision"]),
                    rows,
                });
                if col == "ani" {
                    calib.push((r.to_string(), pts));
                }
            }
        }
    }

    // reads: accuracy against the amount of sequence
    let reads: Vec<&&J> = scored.iter().filter(|j| j.get("kind") == "reads").collect();
    if !reads.is_empty() && has_col("total_bases") {
        let rd: Vec<&J> = reads.iter().map(|j| **j).collect();
        let key = |j: &J| {
            cut_right(
                j.num("total_bases"),
                &[1e7, 1e8, 1e9, 1e12],
                &["<10 Mb", "10-100 Mb", "0.1-1 Gb", ">= 1 Gb"],
            )
        };
        for (val, g) in groups(&rd, &key) {
            let ms = table_of(&g, &format!("reads_bases={val}"));
            if !ms.is_empty() && g.len() >= 5 {
                sections.push(Section {
                    title: format!("reads, {val} of sequence (n = {})", g.len()),
                    text: String::new(),
                    head: head(&RANK_HEAD),
                    rows: rank_rows(&ms),
                });
            }
            all.extend(ms);
        }
    }

    let intro = format!(
        "{} inputs scored ({unmatched} without a truth row); status: {status_txt}. Precision = correct among the calls made; called = share of inputs with a call at that rank; recall = correct calls / inputs with truth. Errors are excluded. Truth = the lineage in the truth table (NCBI taxonomy or curated); where it is blank the rank is not scored.",
        results.len()
    );

    // outputs
    let mut md = format!("# {}\n\n{intro}\n", opts.title);
    for s in &sections {
        md.push_str(&format!("\n## {}\n\n", s.title));
        if !s.text.is_empty() {
            md.push_str(&format!("{}\n", s.text));
        }
        if !s.head.is_empty() {
            md.push_str(&format!(
                "\n| {} |\n|{}\n",
                s.head.join(" | "),
                "---|".repeat(s.head.len())
            ));
            for r in &s.rows {
                md.push_str(&format!("| {} |\n", r.join(" | ")));
            }
        }
    }
    let fmt = |v: Option<f64>| v.map_or(String::new(), |x| format!("{x:.4}"));
    let mut tsv = String::from("group\trank\tn\tcalled\tprecision\trecall\n");
    for m in &all {
        tsv.push_str(&format!(
            "{}\t{}\t{}\t{:.4}\t{}\t{}\n",
            m.group,
            m.rank,
            m.n,
            m.called,
            fmt(m.precision),
            fmt(m.recall)
        ));
    }
    let html = html_page(&opts.title, &intro, &overall, &calib, &sections);
    Output { md, tsv, html }
}

// ---------------------------------------------------------------- HTML and figures

/// Per rank: share of inputs called (light) and precision of the calls made (dark) with its 95% interval.
fn rank_chart(ms: &[Metric]) -> String {
    let (w, rowh, l, r, t) = (520.0, 44.0, 70.0, 20.0, 26.0);
    let h = t + rowh * ms.len() as f64 + 28.0;
    let pw = w - l - r;
    let x = |v: f64| l + pw * v;
    let mut o = format!(
        r#"<svg viewBox="0 0 {w} {h}" role="img" aria-label="Precision and called share per rank">"#
    );
    for tk in [0.0, 0.25, 0.5, 0.75, 1.0] {
        o.push_str(&format!(
            r#"<line x1="{0:.1}" x2="{0:.1}" y1="{t}" y2="{1:.1}" class="grid"/>"#,
            x(tk),
            h - 22.0
        ));
        o.push_str(&text(
            x(tk),
            h - 8.0,
            &format!("{tk:.2}"),
            "t-mut",
            "middle",
            "",
        ));
    }
    for (i, m) in ms.iter().enumerate() {
        let y = t + i as f64 * rowh;
        o.push_str(&text(
            l - 8.0,
            y + rowh / 2.0 + 4.0,
            m.rank,
            "t-sec",
            "end",
            "",
        ));
        let tip = format!(
            "{}: n {}, called {:.3}, precision {}",
            m.rank,
            m.n,
            m.called,
            m.precision.map_or("-".into(), |p| format!("{p:.3}"))
        );
        o.push_str(&format!(
            r#"<rect class="tip bar" tabindex="0" data-tip="{}" x="{l}" y="{:.1}" width="{:.1}" height="12" rx="3" fill="var(--ord-4)"/>"#,
            esc(&tip),
            y + 6.0,
            (pw * m.called).max(1.0)
        ));
        if let Some(p) = m.precision {
            o.push_str(&format!(
                r#"<rect class="tip bar" tabindex="0" data-tip="{}" x="{l}" y="{:.1}" width="{:.1}" height="12" rx="3" fill="var(--ord-2)"/>"#,
                esc(&tip),
                y + 22.0,
                (pw * p).max(1.0)
            ));
            let (lo, hi) = wilson(p, m.n_called);
            let yy = y + 28.0;
            o.push_str(&format!(
                r#"<line x1="{:.1}" x2="{:.1}" y1="{yy:.1}" y2="{yy:.1}" stroke="var(--text-primary)" stroke-width="1.5"/>"#,
                x(lo),
                x(hi)
            ));
        }
    }
    o.push_str(&text(
        l,
        14.0,
        "light = share of inputs called · dark = precision of the calls (95% interval)",
        "t-mut",
        "start",
        "",
    ));
    o.push_str("</svg>");
    o
}

/// Precision of genus/species calls by ANI of the top hit (a calibration chart: should rise with ANI).
fn calib_chart(series: &Calib) -> String {
    let (w, h, l, r, t, b) = (520.0, 260.0, 40.0, 12.0, 28.0, 46.0);
    let (pw, ph) = (w - l - r, h - t - b);
    let nb = series.iter().map(|s| s.1.len()).max().unwrap_or(1).max(1);
    let slot = pw / nb as f64;
    let mut o =
        format!(r#"<svg viewBox="0 0 {w} {h}" role="img" aria-label="Precision by top-hit ANI">"#);
    for tk in [0.0, 0.25, 0.5, 0.75, 1.0] {
        let y = t + ph - ph * tk;
        o.push_str(&format!(
            r#"<line x1="{l}" x2="{}" y1="{y:.1}" y2="{y:.1}" class="grid"/>"#,
            w - r
        ));
        o.push_str(&text(l - 6.0, y + 4.0, &f(tk, 2), "t-mut", "end", ""));
    }
    let colors = ["var(--ord-2)", "var(--ph-2)"];
    for (si, (name, pts)) in series.iter().enumerate() {
        let n_s = series.len() as f64;
        let bw = ((slot - 6.0) / n_s).clamp(6.0, 28.0);
        for (i, (lab, n, p)) in pts.iter().enumerate() {
            let x0 = l + i as f64 * slot + (slot - bw * n_s) / 2.0 + si as f64 * bw;
            let hh = ph * p;
            let tip = format!(
                "{name}, ANI {lab}: precision {p:.3} (n = {})",
                compact(*n as f64)
            );
            o.push_str(&format!(
                r#"<rect class="tip bar" tabindex="0" data-tip="{}" x="{x0:.1}" y="{:.1}" width="{:.1}" height="{hh:.1}" rx="3" fill="{}"/>"#,
                esc(&tip),
                t + ph - hh,
                bw - 2.0,
                colors[si % 2]
            ));
        }
    }
    if let Some((_, pts)) = series.iter().max_by_key(|s| s.1.len()) {
        for (i, (lab, _, _)) in pts.iter().enumerate() {
            o.push_str(&text(
                l + (i as f64 + 0.5) * slot,
                h - 28.0,
                lab,
                "t-mut",
                "middle",
                "",
            ));
        }
    }
    for (si, (name, _)) in series.iter().enumerate() {
        o.push_str(&format!(
            r#"<rect x="{}" y="{}" width="10" height="10" rx="2" fill="{}"/>"#,
            l + si as f64 * 90.0,
            h - 14.0,
            colors[si % 2]
        ));
        o.push_str(&text(
            l + si as f64 * 90.0 + 14.0,
            h - 5.0,
            name,
            "t-sec",
            "start",
            "",
        ));
    }
    o.push_str("</svg>");
    o
}

fn html_table(s: &Section) -> String {
    if s.head.is_empty() {
        return String::new();
    }
    let th: String = s
        .head
        .iter()
        .map(|h| format!("<th>{}</th>", esc(h)))
        .collect();
    let body: String = s
        .rows
        .iter()
        .map(|r| {
            format!(
                "<tr>{}</tr>",
                r.iter()
                    .map(|c| format!("<td>{}</td>", esc(c)))
                    .collect::<String>()
            )
        })
        .collect();
    format!(
        r#"<div class="scroll"><table><thead><tr>{th}</tr></thead><tbody>{body}</tbody></table></div>"#
    )
}

fn html_page(
    title: &str,
    intro: &str,
    overall: &[Metric],
    calib: &Calib,
    sections: &[Section],
) -> String {
    let tile = |label: &str, rank: &str, pick: fn(&Metric) -> Option<f64>| {
        let v = overall
            .iter()
            .find(|m| m.rank == rank)
            .and_then(pick)
            .map_or("-".into(), |v| format!("{v:.3}"));
        format!(
            r#"<div class="tile"><div class="tile-l">{}</div><div class="tile-v">{v}</div></div>"#,
            esc(label)
        )
    };
    let tiles = format!(
        "{}{}{}{}",
        tile("Phylum precision", "phylum", |m| m.precision),
        tile("Genus precision", "genus", |m| m.precision),
        tile("Genus called", "genus", |m| Some(m.called)),
        tile("Species precision", "species", |m| m.precision),
    );
    let figs = format!(
        r#"<section class="card2"><h2>Precision and calls per rank</h2><figure class="fig">{}</figure></section>{}"#,
        rank_chart(overall),
        if calib.is_empty() {
            String::new()
        } else {
            format!(
                r#"<section class="card2"><h2>Calibration</h2><p class="lede">Precision of genus and species calls by ANI of the top hit; it should rise with ANI.</p><figure class="fig">{}</figure></section>"#,
                calib_chart(calib)
            )
        }
    );
    let blocks: String = sections
        .iter()
        .map(|s| {
            format!(
                r#"<section class="card2"><h2>{}</h2>{}{}</section>"#,
                esc(&s.title),
                if s.text.is_empty() {
                    String::new()
                } else {
                    format!(r#"<p class="lede">{}</p>"#, esc(&s.text))
                },
                html_table(s)
            )
        })
        .collect();
    format!(
        r##"<!doctype html>
<html lang="en"><head><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1">
<title>{t}</title><style>{CSS}.scroll{{overflow-x:auto}}table{{border-collapse:collapse;width:100%;font-size:14px}}th,td{{text-align:left;padding:6px 10px;border-bottom:1px solid var(--border)}}</style></head>
<body><main class="page viz-root">
<header class="top"><div><div class="eyebrow">GmG · Guess My Genome</div><h1>{t}</h1><p class="sub">{i}</p></div></header>
<div class="tiles" style="margin-bottom:16px">{tiles}</div>
{figs}{blocks}
</main></body></html>
"##,
        t = esc(title),
        i = esc(intro)
    )
}

#[cfg(test)]
mod tests {
    use super::*;

    fn row(kv: &[(&str, &str)]) -> Row {
        kv.iter()
            .map(|(k, v)| (k.to_string(), v.to_string()))
            .collect()
    }

    #[test]
    fn definitions_match_validate_run_py() {
        // 4 inputs with a genus truth; 3 get a genus call, 2 of those right: precision 2/3, called 3/4, recall 2/4
        let res: Vec<Row> = [
            ("a", "ok", "Aspergillus", "0.99"),
            ("b", "ok", "Fusarium", "0.97"),
            ("c", "ok", "Fusarium", "0.91"),
            ("d", "low_confidence", "", "0.80"),
            ("e", "error", "", ""),
        ]
        .iter()
        .map(|(id, st, g, ani)| {
            row(&[
                ("query_id", id),
                ("status", st),
                ("call_genus", g),
                ("ani", ani),
            ])
        })
        .collect();
        let tru: Vec<Row> = [
            ("a", "Aspergillus"),
            ("b", "Fusarium"),
            ("c", "Aspergillus"),
            ("d", "Candida"),
            ("e", "Candida"),
        ]
        .iter()
        .map(|(id, g)| row(&[("id", id), ("genus", g)]))
        .collect();
        let o = score(&res, &tru, &Options::default());
        let line = o.tsv.lines().find(|l| l.starts_with("all\tgenus")).unwrap();
        assert_eq!(line, "all\tgenus\t4\t0.7500\t0.6667\t0.5000");
        assert!(o.md.contains("| genus | 4 | 0.750 | 0.667 ["), "{}", o.md);
        // the wrong call is listed: truth Aspergillus, called Fusarium
        assert!(o.md.contains("| Aspergillus | Fusarium | 1 |"));
        assert!(o.html.contains("<svg") && o.html.contains("Genus precision"));
        // calibration: the ANI >= 0.98 bin has the one correct call out of one
        assert!(
            o.tsv
                .contains("calib_ani=[0.98, 0.995)\tgenus\t1\t1.0000\t1.0000"),
            "{}",
            o.tsv
        );
    }

    #[test]
    fn unmatched_and_blank_truth_are_not_scored() {
        let res = vec![row(&[
            ("query_id", "x"),
            ("status", "ok"),
            ("call_genus", "G"),
        ])];
        let tru = vec![row(&[("id", "y"), ("genus", "G")])];
        let o = score(&res, &tru, &Options::default());
        assert!(o.md.contains("1 inputs scored (1 without a truth row)"));
        assert!(!o.tsv.contains("all\tgenus"));
    }
}
