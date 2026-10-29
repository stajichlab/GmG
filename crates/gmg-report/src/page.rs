//! Page components and assembly. The look is fixed by `assets/report.css` (tokens from the dataviz palette).

use crate::charts::{self, ph_dot, Chart, Depth};
use crate::model::{Rec, RANKS};
use crate::svg::{compact, esc, f, short};
use crate::Options;
use serde_json::Value;

const CSS: &str = include_str!("../assets/report.css");
const JS: &str = include_str!("../assets/report.js");

fn icon(status: &str) -> &'static str {
    match status {
        "ok" => {
            r##"<svg class="ic" viewBox="0 0 16 16" aria-hidden="true"><circle cx="8" cy="8" r="7" fill="var(--st-good)"/><path d="M4.6 8.3l2.2 2.2 4.5-4.8" fill="none" stroke="#fff" stroke-width="1.8" stroke-linecap="round" stroke-linejoin="round"/></svg>"##
        }
        "low_confidence" => {
            r##"<svg class="ic" viewBox="0 0 16 16" aria-hidden="true"><path d="M8 1.5l7 12.5H1z" fill="var(--st-warn)"/><path d="M8 6v3.6M8 11.6v.4" stroke="#0b0b0b" stroke-width="1.7" stroke-linecap="round"/></svg>"##
        }
        "no_match" => {
            r##"<svg class="ic" viewBox="0 0 16 16" aria-hidden="true"><circle cx="8" cy="8" r="7" fill="var(--st-serious)"/><path d="M5 11L11 5" stroke="#0b0b0b" stroke-width="1.8" stroke-linecap="round"/></svg>"##
        }
        _ => {
            r##"<svg class="ic" viewBox="0 0 16 16" aria-hidden="true"><circle cx="8" cy="8" r="7" fill="var(--st-crit)"/><path d="M5.3 5.3l5.4 5.4M10.7 5.3l-5.4 5.4" stroke="#fff" stroke-width="1.8" stroke-linecap="round"/></svg>"##
        }
    }
}

pub(crate) fn status_label(s: &str) -> &'static str {
    match s {
        "ok" => "classified",
        "low_confidence" => "low confidence",
        "no_match" => "no match",
        _ => "error",
    }
}

/// Icon + label always travel together: colour is never the only channel for status.
fn chip(status: &str) -> String {
    format!(
        r#"<span class="chip">{}{}</span>"#,
        icon(status),
        status_label(status)
    )
}

/// Thresholds used for the run (from run.json), for the "why was this rank not called" text.
#[derive(Clone)]
pub struct Th {
    pub min_ani: [f64; 6],
    pub min_cons: [f64; 6],
    pub tentative_genus_min_ani: f64,
}

impl Th {
    pub fn from_run(run: &Value) -> Th {
        let v = &run["thresholds"]["values"];
        let mut t = Th {
            min_ani: [0.7109, 0.7109, 0.7109, 0.7715, 0.90, 0.9938],
            min_cons: [0.5, 0.5, 0.6, 0.9, 0.5, 0.9],
            tentative_genus_min_ani: 0.7765,
        };
        for i in 0..6 {
            if let Some(x) = v["min_ani"][i].as_f64() {
                t.min_ani[i] = x;
            }
            if let Some(x) = v["min_consistency"][i].as_f64() {
                t.min_cons[i] = x;
            }
        }
        if let Some(x) = v["tentative_genus_min_ani"].as_f64() {
            t.tentative_genus_min_ani = x;
        }
        t
    }
}

/// Why was the first uncalled rank not called? Mirrors the calling logic. Returns (rank index, sentence).
pub(crate) fn reasons(r: &Rec, th: &Th) -> Vec<(usize, String)> {
    let best = r.best_rank.map(|b| b as i32).unwrap_or(-1);
    let ani = r.ani().unwrap_or(0.0);
    for (i, k) in RANKS.iter().enumerate() {
        if i as i32 <= best {
            continue;
        }
        let cons = r.calls[i].consistency;
        let text = if r.calls[i].vote.is_none() {
            format!("The nearest reference genomes have no {k} recorded in the reference metadata, so there is nothing to call. This reflects the reference annotation, not weak evidence")
        } else if cons < th.min_cons[i] {
            format!("Vote consistency {cons:.3} is below the {:.3} required: the near-best references disagree", th.min_cons[i])
        } else if ani < th.min_ani[i] {
            match (&r.tentative_genus, i) {
                (Some((t, _, _)), 4) => format!(
                    "ANI {ani:.4} is below the {:.4} required for a genus call. The evidence points to {t} as the nearest genus, shown as tentative: in this range a named genus is right only about 80% of the time even when the genus is in the reference, and it may be a genus that is not in the reference",
                    th.min_ani[i]
                ),
                _ => format!("ANI {ani:.4} is below the {:.4} required for a {k} call", th.min_ani[i]),
            }
        } else {
            "A higher rank was not called, so deeper ranks are not called".to_string()
        };
        return vec![(i, text)];
    }
    vec![]
}

fn ladder(r: &Rec) -> String {
    let best = r.best_rank.map(|b| b as i32).unwrap_or(-1);
    let mut cells = String::new();
    for (i, k) in RANKS.iter().enumerate() {
        let called = i as i32 <= best && !r.call(i).is_empty();
        let cons = r.calls[i].consistency;
        let gap = !called && (i as i32) < best; // a deeper rank was called; the reference metadata has no value here
        let tent = if i == 4 && !called {
            r.tentative_genus.as_ref()
        } else {
            None
        };
        let name = if let Some((t, _, _)) = tent {
            format!("{t} (tentative)")
        } else if called {
            r.call(i).to_string()
        } else if gap {
            "not in reference".into()
        } else {
            "withheld".into()
        };
        let tip = if let Some((t, a, c)) = tent {
            format!("genus: {t} is the nearest genus but is only tentative (ANI {a:.3} is below the genus cutoff; vote consistency {c:.2})")
        } else if called {
            format!("{k}: {} · vote consistency {cons:.2}", r.call(i))
        } else if gap {
            format!("{k}: the reference metadata has no {k} for these genomes")
        } else {
            format!("{k}: not called")
        };
        let class = if tent.is_some() {
            "tent"
        } else if called {
            "on"
        } else if gap {
            "na"
        } else {
            "off"
        };
        let style = if *k == "phylum" && called {
            format!(
                r#" style="border-left-color:{}""#,
                charts::ph_color(r.call(0))
            )
        } else {
            String::new()
        };
        cells.push_str(&format!(
            r#"<div class="rk {class} tip" tabindex="0" data-tip="{}"{style}><div class="rk-n">{k}</div><div class="rk-v">{}</div><div class="rk-m"><span style="width:{:.0}%"></span></div></div>"#,
            esc(&tip),
            esc(&name),
            if called {
                100.0 * cons
            } else if let Some((_, _, c)) = tent {
                100.0 * c
            } else {
                0.0
            }
        ));
    }
    format!(r#"<div class="ladder">{cells}</div>"#)
}

fn figtools(c: &Chart) -> String {
    format!(
        r#"<div class="figbtn"><button class="tbtn" data-table-toggle aria-pressed="false">Table view</button></div><div class="chartbox">{}</div><div class="tableview" hidden>{}</div>"#,
        c.svg, c.table
    )
}

fn qc_panel(r: &Rec) -> String {
    let bases = if r.total_bases >= 1e6 {
        format!("{:.1} Mb", r.total_bases / 1e6)
    } else {
        compact(r.total_bases)
    };
    let n50 = if r.kind == "reads" || r.n50 == 0.0 {
        "n/a (reads)".to_string()
    } else {
        compact(r.n50)
    };
    let flags = if r.flags.is_empty() {
        r#"<span class="muted">none</span>"#.to_string()
    } else {
        r.flags
            .iter()
            .map(|x| format!(r#"<span class="flag">{}</span>"#, esc(x)))
            .collect::<String>()
    };
    format!(
        r#"<div class="qc"><div><span class="muted">Bases</span><b>{}</b></div><div><span class="muted">N50</span><b>{}</b></div><div><span class="muted">Sketch hashes</span><b>{}</b></div><div><span class="muted">QC flags</span>{flags}</div></div>"#,
        esc(&bases),
        esc(&n50),
        r.n_hashes
    )
}

fn card(r: &Rec, th: &Th, open: bool) -> String {
    let id = esc(&r.id);
    let dot = if r.phylum().is_empty() {
        String::new()
    } else {
        ph_dot(r.phylum())
    };
    let call = r.call(5).to_string();
    let call = if call.is_empty() {
        r.call(4).to_string()
    } else {
        call
    };
    let call = match (&r.tentative_genus, call.is_empty()) {
        (Some((t, _, _)), true) => format!("{t} (tentative)"),
        _ => call,
    };
    let head = format!(
        r#"<summary><span class="card-id">{id}</span>{}{dot}<span class="card-call">{}</span></summary>"#,
        chip(&r.status),
        esc(&call)
    );
    if r.status == "error" {
        return format!(
            r#"<details class="card" id="q-{id}" open>{head}<div class="card-b"><p class="why">This input could not be processed{}: <b>{}</b>{} The run continued; this row stays in results.tsv with status <code>error</code>.</p></div></details>"#,
            if r.error_code.is_empty() {
                String::new()
            } else {
                format!(" (<code>{}</code>)", esc(&r.error_code))
            },
            esc(&r.message),
            if r.error_action.is_empty() {
                ".".to_string()
            } else {
                format!(" {}", esc(&r.error_action))
            }
        );
    }
    let bars = charts::tophit_bars(r);
    let reasons = reasons(r, th);
    let ev = if let Some(note) = little_data_note(r) {
        format!("<li><b>Too little data.</b> {note}</li>")
    } else if r.status == "no_match" {
        "<li><b>No rank called.</b> The nearest reference shares too little with this input to support even a phylum call. Not a fungus, or a lineage the reference does not cover.</li>".to_string()
    } else {
        reasons
            .iter()
            .map(|(i, t)| {
                format!(
                    "<li><b>{} not called.</b> {}.</li>",
                    capitalize(RANKS[*i]),
                    esc(t)
                )
            })
            .collect::<String>()
    };
    let ev_html = if ev.is_empty() {
        r#"<p class="why ok">Every rank down to species passed its threshold.</p>"#.to_string()
    } else {
        format!(r#"<ul class="why">{ev}</ul>"#)
    };
    let top = r.top.clone().unwrap_or_default();
    let clus = if top.cluster.is_empty() {
        String::new()
    } else {
        format!(
            r#"<div><span class="muted">ANI cluster of the top hit</span><b class="mono">{}</b></div>"#,
            esc(&short(&top.cluster, 40))
        )
    };
    let reads = match (&r.hist, r.kind.as_str()) {
        (Some(h), "reads") => {
            let c = charts::abundance_chart(h);
            format!(
                r#"<h4>Reads: what was kept</h4><figure class="fig">{}<figcaption>Hashes seen once are mostly sequencing errors and are dropped (minimum abundance {}). The peak of the kept hashes (median {}×) estimates k-mer coverage; {:.1}% of the kept hashes are found in the best reference.</figcaption></figure>"#,
                figtools(&c),
                h.min_abund,
                h.median_kept,
                100.0 * top.containment_query
            )
        }
        _ => String::new(),
    };
    format!(
        r#"<details class="card" id="q-{id}"{}>{head}<div class="card-b"><h4>Call by rank</h4>{}{ev_html}<div class="two"><div><h4>Best matching references</h4><figure class="fig">{}<figcaption>Bar = share of this input's sketch hashes found in the reference; the first bar is the top hit.</figcaption></figure></div><div><h4>Input</h4>{}<div class="qc"><div><span class="muted">ANI to top hit</span><b>{}</b></div><div><span class="muted">Shared hashes</span><b>{}</b></div>{clus}</div></div></div>{reads}</div></details>"#,
        if open { " open" } else { "" },
        ladder(r),
        figtools(&bars),
        qc_panel(r),
        f(top.ani, 4),
        top.shared
    )
}

/// Explanation for an unresolved input that is simply too small or fragmented to classify, if the QC says so.
pub(crate) fn little_data_note(r: &Rec) -> Option<String> {
    if r.status != "no_match" && r.status != "low_confidence" {
        return None;
    }
    let small = r
        .flags
        .iter()
        .any(|f| matches!(f.as_str(), "short_genome" | "few_hashes" | "fragmented"));
    if !small {
        return None;
    }
    let size = if r.total_bases >= 1e6 {
        format!("{:.1} Mb", r.total_bases / 1e6)
    } else {
        format!("{:.0} kb", r.total_bases / 1e3)
    };
    Some(format!(
        "{size} and {} sketch hashes ({}): too little of a genome to classify. This says nothing about the organism; a more complete assembly or more reads may resolve it.",
        r.n_hashes,
        r.flags.join(", ")
    ))
}

fn capitalize(s: &str) -> String {
    let mut c = s.chars();
    match c.next() {
        Some(h) => h.to_uppercase().collect::<String>() + c.as_str(),
        None => String::new(),
    }
}

fn stat(label: &str, value: usize, status: Option<&str>) -> String {
    let ic = status.map(icon).unwrap_or("");
    format!(
        r#"<div class="tile"><div class="tile-l">{ic}{}</div><div class="tile-v">{value}</div></div>"#,
        esc(label)
    )
}

fn trow(r: &Rec, has_card: bool) -> String {
    let call = r.call_name();
    let link = if has_card {
        format!(r##"<a href="#q-{}">{}</a>"##, esc(&r.id), esc(&r.id))
    } else {
        esc(&r.id)
    };
    let best = r.best_rank.map(|b| RANKS[b]).unwrap_or("—");
    let flags: String = r
        .flags
        .iter()
        .map(|x| format!(r#"<span class="flag">{}</span>"#, esc(x)))
        .collect();
    let ph = r.phylum();
    let ani = r.ani();
    format!(
        r#"<tr data-status="{st}" data-q="{q}"><td data-v="{idv}">{link}</td><td data-v="{st}">{}</td><td data-v="{phv}">{}{}</td><td data-v="{rk}">{best}</td><td data-v="{callv}">{}</td><td data-v="{anv}" class="num">{}</td><td>{flags}</td></tr>"#,
        chip(&r.status),
        if ph.is_empty() {
            String::new()
        } else {
            ph_dot(ph)
        },
        if ph.is_empty() {
            "—".to_string()
        } else {
            esc(ph)
        },
        {
            let tent = r
                .tentative_genus
                .as_ref()
                .map(|(t, _, _)| {
                    format!(
                        r#" <span class="muted">· tentative genus: {}</span>"#,
                        esc(t)
                    )
                })
                .unwrap_or_default();
            if call.is_empty() && tent.is_empty() {
                "—".to_string()
            } else {
                format!("{}{tent}", esc(call))
            }
        },
        ani.map(|a| f(a, 3)).unwrap_or_else(|| "—".into()),
        st = r.status,
        q = esc(&format!("{} {} {}", r.id, call, ph).to_lowercase()),
        idv = esc(&r.id),
        phv = esc(ph),
        rk = r.best_rank.map(|b| b as i32).unwrap_or(-1),
        callv = esc(call),
        anv = ani.unwrap_or(0.0)
    )
}

/// Inputs that deserve a look first: errors, no match, low confidence, then the lowest ANI.
fn priority_order(recs: &[Rec]) -> Vec<usize> {
    let mut idx: Vec<usize> = (0..recs.len()).collect();
    let rank = |r: &Rec| match r.status.as_str() {
        "error" => 0,
        "no_match" => 1,
        "low_confidence" => 2,
        _ => 3,
    };
    idx.sort_by(|&a, &b| {
        rank(&recs[a])
            .cmp(&rank(&recs[b]))
            .then(
                recs[a]
                    .ani()
                    .unwrap_or(1.0)
                    .partial_cmp(&recs[b].ani().unwrap_or(1.0))
                    .unwrap_or(std::cmp::Ordering::Equal),
            )
            .then(a.cmp(&b))
    });
    idx
}

pub fn render(
    results: &[Value],
    run: &Value,
    relationships: Option<&Value>,
    opts: &Options,
) -> String {
    let recs: Vec<Rec> = results.iter().map(Rec::from_value).collect();
    let th = Th::from_run(run);
    let n_in = recs.len();
    let errors = recs.iter().filter(|r| r.status == "error").count();
    let ok_total = n_in - errors;
    let depth = Depth::from_recs(&recs);
    let classified = depth.species + depth.genus;
    let pct = if ok_total > 0 {
        (100.0 * classified as f64 / ok_total as f64).round() as usize
    } else {
        0
    };

    let crate::figs::Figures {
        depth_chart,
        ph_chart,
        genus_chart,
        genus_legend,
        n_genus,
        ani_chart,
        low_tbl,
        heat,
    } = crate::figs::figures(&recs, &depth, &th, relationships);

    // cards: problems first, capped
    let order = priority_order(&recs);
    let card_idx: Vec<usize> = order.into_iter().take(opts.max_cards).collect();
    let carded: std::collections::HashSet<usize> = card_idx.iter().cloned().collect();
    let open_all = n_in <= 3;
    let cards: String = card_idx
        .iter()
        .map(|&i| card(&recs[i], &th, open_all))
        .collect();
    let cards_note = if card_idx.len() < n_in {
        format!(
            r#"<p class="muted small">Detail cards are shown for the {} inputs that need the most attention; all {} inputs are in the table above and in results.tsv.</p>"#,
            card_idx.len(),
            n_in
        )
    } else {
        String::new()
    };

    let shown = recs.len().min(opts.max_rows);
    let rows: String = recs
        .iter()
        .enumerate()
        .take(shown)
        .map(|(i, r)| trow(r, carded.contains(&i)))
        .collect();
    let rows_note = if shown < n_in {
        format!(
            r#"<p class="muted small">The table shows the first {shown} of {n_in} inputs; all are in results.tsv.</p>"#
        )
    } else {
        String::new()
    };

    let idx = &run["index"];
    let th_name = run["thresholds"]["name"].as_str().unwrap_or("known");
    let thr_rows: String = RANKS
        .iter()
        .enumerate()
        .map(|(i, k)| {
            let a = if th.min_ani[i] > 1.0 {
                "never called".to_string()
            } else {
                format!("{:.4}", th.min_ani[i])
            };
            format!(
                "<tr><td>{k}</td><td>{a}</td><td>{:.2}</td></tr>",
                th.min_cons[i]
            )
        })
        .collect();
    let started = run["started_utc"].as_str().unwrap_or("");
    let commit = run["commit"].as_str().unwrap_or("");
    let tsv = opts
        .tsv_link
        .as_deref()
        .map(|l| {
            format!(
                r#"<a class="dl" href="{}" download>Download results.tsv</a>"#,
                esc(l)
            )
        })
        .unwrap_or_default();
    let banner = opts
        .banner
        .as_deref()
        .map(|b| format!(r#"<div class="mockbar">{}</div>"#, esc(b)))
        .unwrap_or_default();
    let heat_html = heat
        .map(|(c, n)| {
            format!(
                r#"<section class="card2"><h2>Relationships among inputs <span class="opt">optional · up to 60 inputs</span></h2><p class="lede">Pairwise ANI between the {n} assemblies, computed from the same sketches. Blocks along the diagonal are species and species complexes.</p><figure class="fig">{}</figure></section>"#,
                figtools(&c)
            )
        })
        .unwrap_or_default();
    let depth_html = if depth.total() == 0 {
        String::from(r#"<p class="lede">No input could be classified.</p>"#)
    } else {
        format!(r#"<figure class="fig">{}</figure>"#, figtools(&depth_chart))
    };

    format!(
        r##"<!doctype html>
<html lang="en"><head><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1">
<title>{title}</title><style>{CSS}</style></head>
<body>
{banner}
<main class="page viz-root">
<header class="top">
  <div><div class="eyebrow">GmG · Guess My Genome</div><h1>Classification report</h1>
  <p class="sub">Run {started_day} · {n_in} inputs · gmg {version} ({commit9}) · index of {ngen} genomes (k = {k}, scaled = {scaled}) · thresholds: {th_name}</p></div>
  <button id="theme" class="tbtn" aria-label="Toggle dark mode">Dark / light</button>
</header>

<section class="summary" aria-label="Summary">
  <div class="hero"><div class="hero-v">{pct}%</div><div class="hero-l">of the {ok_total} processed inputs were classified to genus or species</div></div>
  <div class="tiles">{t_species}{t_genus}{t_low}{t_none}{t_err}</div>
</section>

<section class="card2"><h2>How deep each call went</h2>
<p class="lede">A rank is called only when its evidence passes the calibrated threshold; otherwise the report stops higher up and says why.</p>
{depth_html}</section>

<section class="grid2">
  <div class="card2"><h2>What the run contains</h2><h3>Phylum</h3><figure class="fig">{ph_fig}<figcaption>Colours are fixed per phylum (the eight largest phyla in the reference). Other phyla share grey; the table view names every phylum.</figcaption></figure>
  <h3>Genus (top {n_genus}), coloured by phylum</h3><figure class="fig">{genus_fig}{genus_legend}</figure></div>
  <div class="card2"><h2>How close the nearest reference is</h2>
  <p class="lede">Calls to species need ANI ≥ {sp_thr}; lower values mean less well covered lineages.</p>
  <figure class="fig">{ani_fig}</figure>
  {low_tbl}</div>
</section>

{heat_html}

<section class="card2" id="results"><h2>All inputs</h2>
<div class="filters"><label>Status <select id="f-status"><option value="">all</option><option value="ok">classified</option><option value="low_confidence">low confidence</option><option value="no_match">no match</option><option value="error">error</option></select></label>
<label>Search <input id="f-q" type="search" placeholder="input or taxon"></label><span id="f-n" class="muted"></span>
{tsv}</div>
<div class="tscroll"><table class="results"><thead><tr><th data-sort>Input</th><th data-sort>Status</th><th data-sort>Phylum</th><th data-sort>Deepest rank</th><th data-sort>Call</th><th data-sort class="num">ANI</th><th>QC flags</th></tr></thead>
<tbody>{rows}</tbody></table></div>{rows_note}</section>

<section><h2 class="cards-h">Per-input detail</h2>{cards}{cards_note}</section>

<section class="card2 meth"><h2>How to read this report</h2>
<ul>
<li><b>Call by rank.</b> Phylum to species. A rank is called when the nearest reference genomes agree on it (vote consistency) and the input is similar enough to them (ANI). Higher ranks are easier than lower ones, so the ladder can stop part-way.</li>
<li><b>Tentative genus.</b> A genus is called only at ANI of at least the cutoff shown in the thresholds table (0.90 by default). Between that and ANI 0.7765 the nearest genus is shown as tentative: it is right only about 80% of the time in that range even when the genus is in the reference, and wrong when the genus is not in the reference. It never counts as a call.</li>
<li><b>Vote consistency</b> is the share of near-best references agreeing on the taxon, an evidence value, not a probability.</li>
<li><b>ANI</b> here is containment ANI: the k = 21 sketch fraction of the input found in a reference, raised to 1/21. From a very small overlap it is not meaningful: look at the percentage of the input found.</li>
<li><b>Species labels follow the reference metadata.</b> Where references disagree about a name (hybrids, species complexes) the genus is usually right and the species is withheld; the ANI cluster is the safer unit.</li>
<li>Thresholds are calibrated to about 95% precision on held-out genomes of the reference (<code>known</code> set). The genus cutoff is raised to ANI 0.90 because with the earlier 0.7765 a genus was named for 58% of queries whose genus is not in the reference. Genomes from a genus the reference lacks may need the stricter <code>novel</code> set.</li>
</ul>
<h3>Thresholds used ({th_name})</h3>
<table class="small"><thead><tr><th>Rank</th><th>Minimum ANI</th><th>Minimum vote consistency</th></tr></thead><tbody>{thr_rows}</tbody></table>
<h3>Reproducibility</h3>
<p class="small mono">gmg {version} ({commit}) · index {index_name} ({index_bytes} bytes) · {elapsed} s · results.tsv · hits.tsv · run.json</p>
</section>
</main>
<div id="tooltip" role="tooltip" hidden></div>
<script>{JS}</script>
</body></html>"##,
        title = esc(&opts.title),
        started_day = esc(started.get(..10).unwrap_or(started)),
        version = esc(run["version"].as_str().unwrap_or("")),
        commit9 = esc(commit.get(..9).unwrap_or(commit)),
        commit = esc(commit),
        ngen = idx["n_genomes"].as_u64().unwrap_or(0),
        k = idx["k"].as_u64().unwrap_or(0),
        scaled = idx["scaled"].as_u64().unwrap_or(0),
        th_name = esc(th_name),
        t_species = stat("species", depth.species, None),
        t_genus = stat("genus", depth.genus, None),
        t_low = stat(
            "low confidence",
            depth.family + depth.higher,
            Some("low_confidence")
        ),
        t_none = stat("no match", depth.none, Some("no_match")),
        t_err = stat("errors", errors, Some("error")),
        ph_fig = figtools(&ph_chart),
        genus_fig = figtools(&genus_chart),
        n_genus = n_genus,
        ani_fig = figtools(&ani_chart),
        sp_thr = if th.min_ani[5] > 1.0 {
            "never (this threshold set does not call species)".to_string()
        } else {
            format!("{:.3}", th.min_ani[5])
        },
        index_name = esc(&crate::model::file_name(idx["path"].as_str().unwrap_or(""))),
        index_bytes = idx["bytes"].as_u64().unwrap_or(0),
        elapsed = run["elapsed_seconds"].as_f64().unwrap_or(0.0),
    )
}
