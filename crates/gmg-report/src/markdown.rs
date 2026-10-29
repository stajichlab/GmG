//! Markdown report: `report.md` plus standalone SVG figures, so the report renders on GitHub or in a notebook.
//! Same numbers as the HTML page (they share `figs`); no scripts, no dark mode, one table row per input.
use crate::charts::Depth;
use crate::model::{Rec, RANKS};
use crate::page::{little_data_note, reasons, status_label, Th, NOT_FUNGAL_NOTE};
use crate::Options;
use serde_json::Value;

pub struct Markdown {
    /// the document; figures are referenced as `figures/<name>`
    pub md: String,
    /// (file name, standalone SVG text)
    pub figures: Vec<(String, String)>,
}

const TOKENS: &str = ":root{--ph-1:#2a78d6;--ph-2:#eb6834;--ph-3:#1baf7a;--ph-4:#eda100;--ph-5:#e87ba4;--ph-6:#008300;--ph-7:#6250d6;--ph-8:#e34948;--ph-other:#a8a7a0;\
--text-primary:#0b0b0b;--text-secondary:#52514e;--text-muted:#898781;--grid:#e1e0d9;--axis:#c3c2b7;--series-1:#2a78d6;--deemph:#c3c2b7;--nomatch:#d9d8d1;\
--ord-1:#184f95;--ord-2:#256abf;--ord-3:#3987e5;--ord-4:#86b6ef;--page:#fff;--surface-1:#fcfcfb;--border:rgba(11,11,11,.1);\
--st-good:#0ca30c;--st-warn:#fab219;--st-serious:#ec835a;--st-crit:#d03b3b;--good-text:#006300}";
/// Chart text and line classes, with the colours written out.
const RULES: &str = "svg text{font:12px Inter,system-ui,-apple-system,'Segoe UI',Helvetica,Arial,sans-serif}.t-pri{fill:#0b0b0b}.t-sec{fill:#52514e}.t-mut{fill:#898781}.t-fill{font-weight:600}\
.grid{stroke:#e1e0d9;stroke-width:1}.axis{stroke:#c3c2b7;stroke-width:1}.thr{stroke:#52514e;stroke-width:1}";

/// Make a report chart usable as a file: namespace, white background, and every `var(--token)` replaced by its colour
/// (renderers such as PDF converters do not implement CSS variables).
pub fn standalone_svg(svg: &str) -> String {
    let Some(end) = svg.find('>') else {
        return svg.to_string();
    };
    let (open, rest) = svg.split_at(end + 1);
    let open = open.replacen("<svg ", "<svg xmlns=\"http://www.w3.org/2000/svg\" ", 1);
    let mut body = rest.to_string();
    for (name, value) in token_values() {
        body = body.replace(&format!("var(--{name})"), value);
    }
    format!(
        "<?xml version=\"1.0\" encoding=\"UTF-8\"?>\n{open}<style>{RULES}</style><rect width=\"100%\" height=\"100%\" fill=\"#fff\"/>{body}"
    )
}

/// `(name, value)` of every colour token in `TOKENS`.
fn token_values() -> Vec<(&'static str, &'static str)> {
    let inner = TOKENS
        .split_once(":root{")
        .and_then(|(_, r)| r.split_once('}'))
        .map(|(t, _)| t)
        .unwrap_or("");
    inner
        .split(';')
        .filter_map(|d| {
            let d = d.trim().strip_prefix("--")?;
            d.split_once(':')
        })
        .collect()
}

fn cell(s: &str) -> String {
    s.replace('|', "\\|").replace(['\n', '\r'], " ")
}

pub fn render(
    results: &[Value],
    run: &Value,
    relationships: Option<&Value>,
    opts: &Options,
) -> Markdown {
    let recs: Vec<Rec> = results.iter().map(Rec::from_value).collect();
    let th = Th::from_run(run);
    let depth = Depth::from_recs(&recs);
    let f = crate::figs::figures(&recs, &depth, &th, relationships);
    let n = recs.len();
    let count = |st: &str| recs.iter().filter(|r| r.status == st).count();
    let mut figures: Vec<(String, String)> = Vec::new();
    let mut fig = |name: &str, alt: &str, svg: &str, out: &mut String| {
        figures.push((name.to_string(), standalone_svg(svg)));
        out.push_str(&format!("![{alt}](figures/{name})\n\n"));
    };

    let mut md = String::new();
    md.push_str(&format!("# {}\n\n", opts.title));
    md.push_str(&format!(
        "Run {} · {} inputs · gmg {} ({}) · index of {} genomes (k = {}, scaled = {}) · thresholds: {}\n\n",
        run["started_utc"].as_str().unwrap_or("").get(..10).unwrap_or(""),
        n,
        run["version"].as_str().unwrap_or(""),
        run["commit"].as_str().unwrap_or(""),
        run["index"]["n_genomes"],
        run["index"]["k"],
        run["index"]["scaled"],
        run["thresholds"]["name"].as_str().unwrap_or("known"),
    ));
    md.push_str("## Summary\n\n| | inputs |\n|---|---|\n");
    md.push_str(&format!(
        "| classified (genus or species) | {} |\n| classified to a higher rank only | {} |\n| low confidence | {} |\n| no match | {} |\n| error | {} |\n\n",
        depth.species + depth.genus,
        depth.family + depth.higher,
        count("low_confidence"),
        count("no_match"),
        count("error"),
    ));
    if depth.total() > 0 {
        md.push_str("### Deepest rank called\n\n");
        fig(
            "depth.svg",
            "Deepest rank called per input",
            &f.depth_chart.svg,
            &mut md,
        );
        md.push_str("### Phyla\n\n");
        fig("phylum.svg", "Inputs per phylum", &f.ph_chart.svg, &mut md);
        md.push_str("### Genera, coloured by phylum\n\n");
        fig(
            "genus.svg",
            "Genus counts grouped by phylum",
            &f.genus_chart.svg,
            &mut md,
        );
        md.push_str("### ANI to the nearest reference\n\n");
        fig(
            "ani.svg",
            "Distribution of ANI to the top reference hit",
            &f.ani_chart.svg,
            &mut md,
        );
    } else {
        md.push_str("No input could be classified.\n\n");
    }
    if let Some((c, k)) = &f.heat {
        md.push_str(&format!(
            "### Relationships among {k} assemblies\n\nPairwise ANI from the same sketches.\n\n"
        ));
        fig("heatmap.svg", "Pairwise ANI heatmap", &c.svg, &mut md);
    }

    // inputs that need attention
    let attention: Vec<&Rec> = recs
        .iter()
        .filter(|r| r.status != "ok" || r.flags.iter().any(|f| f == "not_fungal"))
        .collect();
    if !attention.is_empty() {
        md.push_str("## Needs attention\n\n");
        for r in attention.iter().take(opts.max_cards) {
            let why = if r.flags.iter().any(|f| f == "not_fungal") {
                format!("Probably not a fungus. {NOT_FUNGAL_NOTE}")
            } else if r.status == "error" {
                format!("`{}`: {} {}", r.error_code, r.message, r.error_action)
            } else if let Some(note) = little_data_note(r) {
                note
            } else if r.status == "no_match" {
                "No rank called: the nearest reference shares too little with this input."
                    .to_string()
            } else {
                reasons(r, &th)
                    .first()
                    .map(|(i, t)| format!("{} not called: {t}.", RANKS[*i]))
                    .unwrap_or_default()
            };
            md.push_str(&format!(
                "* **{}** ({}): {}\n",
                cell(&r.id),
                status_label(&r.status),
                cell(&why)
            ));
        }
        md.push('\n');
    }

    md.push_str("## All inputs\n\n| input | status | phylum | deepest rank | call | ANI | flags |\n|---|---|---|---|---|---|---|\n");
    for r in recs.iter().take(opts.max_rows) {
        let call = match (r.call_name(), &r.tentative_genus) {
            ("", Some((t, _, _))) => format!("{t} (tentative genus)"),
            (c, _) => c.to_string(),
        };
        md.push_str(&format!(
            "| {} | {} | {} | {} | {} | {} | {} |\n",
            cell(&r.id),
            status_label(&r.status),
            cell(if r.phylum().is_empty() {
                "—"
            } else {
                r.phylum()
            }),
            r.best_rank.map(|b| RANKS[b]).unwrap_or("—"),
            cell(if call.is_empty() { "—" } else { &call }),
            r.ani()
                .map(|a| format!("{a:.3}"))
                .unwrap_or_else(|| "—".into()),
            cell(&r.flags.join(", ")),
        ));
    }
    if recs.len() > opts.max_rows {
        md.push_str(&format!(
            "\nThe table shows the first {} of {} inputs; all are in results.tsv.\n",
            opts.max_rows,
            recs.len()
        ));
    }
    md.push_str("\n## Methods\n\nGmG sketches each input (FracMinHash, k = 21), finds the nearest reference genomes by shared hashes and calls each taxonomic rank by a vote among the near-best hits. A rank is called only when the containment ANI and the agreement of the hits pass that rank's threshold.\n\n| rank | min ANI | min agreement |\n|---|---|---|\n");
    for (i, k) in RANKS.iter().enumerate() {
        let a = if th.min_ani[i] > 1.0 {
            "never called".to_string()
        } else {
            format!("{:.4}", th.min_ani[i])
        };
        md.push_str(&format!("| {k} | {a} | {:.2} |\n", th.min_cons[i]));
    }
    md.push_str("\nA tentative genus names the nearest genus when the ANI is too low to call it; it is never a call. Every number here is a column of `results.tsv`.\n");
    Markdown { md, figures }
}
