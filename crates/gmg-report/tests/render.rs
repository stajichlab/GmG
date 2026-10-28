//! Report engine tests on a fixture of real `gmg run` records (34 of the 148 inputs of the design-mockup run, paths removed).

use gmg_report::{parse_results_jsonl, render_html, Options};
use serde_json::{json, Value};
use std::collections::HashSet;
use std::fs;

fn fixture() -> (Vec<Value>, Value, Value) {
    let d = concat!(env!("CARGO_MANIFEST_DIR"), "/tests/fixtures/");
    let results = parse_results_jsonl(&fs::read_to_string(format!("{d}results.jsonl")).unwrap());
    let run: Value =
        serde_json::from_str(&fs::read_to_string(format!("{d}run.json")).unwrap()).unwrap();
    let rel: Value =
        serde_json::from_str(&fs::read_to_string(format!("{d}relationships.json")).unwrap())
            .unwrap();
    (results, run, rel)
}

fn svgs(html: &str) -> Vec<&str> {
    let mut out = Vec::new();
    let mut rest = html;
    while let Some(a) = rest.find("<svg") {
        let b = rest[a..].find("</svg>").expect("unterminated svg") + a + 6;
        out.push(&rest[a..b]);
        rest = &rest[b..];
    }
    out
}

#[test]
fn page_has_every_section_and_valid_svg() {
    let (results, run, rel) = fixture();
    let html = render_html(&results, &run, Some(&rel), &Options::default());
    for needle in [
        "<h1>Classification report</h1>",
        "How deep each call went",
        "What the run contains",
        "How close the nearest reference is",
        "Relationships among inputs",
        "All inputs",
        "Per-input detail",
        "How to read this report",
        "Thresholds used",
    ] {
        assert!(html.contains(needle), "missing section: {needle}");
    }
    let charts = svgs(&html);
    assert!(
        charts.len() > 20,
        "expected charts and icons, got {}",
        charts.len()
    );
    for s in charts {
        roxmltree::Document::parse(s)
            .unwrap_or_else(|e| panic!("invalid SVG: {e}\n{}", &s[..s.len().min(300)]));
    }
}

#[test]
fn no_invalid_path_data_no_nan_no_local_paths() {
    let (results, run, rel) = fixture();
    let html = render_html(&results, &run, Some(&rel), &Options::default());
    for chunk in html.split(" d=\"").skip(1) {
        let d = &chunk[..chunk.find('"').unwrap()];
        assert!(
            !d.contains("--") && !d.contains("NaN") && !d.contains("inf"),
            "bad path data: {d}"
        );
    }
    let markup = &html[..html.find("<script>").expect("script block")]; // the embedded JS legitimately says isNaN
    for bad in ["NaN", ">inf<", "=\"inf\"", "-inf"] {
        if let Some(i) = markup.find(bad) {
            panic!(
                "{bad} in output: …{}…",
                &markup[i.saturating_sub(100)..(i + 60).min(markup.len())]
            );
        }
    }
    for p in ["/Users/", "/private/", "/bigdata/", "/home/", "/tmp/"] {
        if let Some(i) = html.find(p) {
            panic!(
                "local path {p} leaked into the report: …{}…",
                &html[i.saturating_sub(80)..(i + 80).min(html.len())]
            );
        }
    }
}

#[test]
fn cards_and_links_are_consistent() {
    let (results, run, rel) = fixture();
    let html = render_html(&results, &run, Some(&rel), &Options::default());
    let ids: Vec<&str> = html
        .split(" id=\"q-")
        .skip(1)
        .map(|c| &c[..c.find('"').unwrap()])
        .collect();
    let uniq: HashSet<&str> = ids.iter().cloned().collect();
    assert_eq!(ids.len(), uniq.len(), "duplicate card ids");
    assert_eq!(
        ids.len(),
        results.len(),
        "one card per input when under the cap"
    );
    for l in html.split("href=\"#q-").skip(1) {
        let id = &l[..l.find('"').unwrap()];
        assert!(uniq.contains(id), "link to missing card {id}");
    }
    // every bar and cell has a tooltip
    assert!(html.matches("class=\"tip bar\"").count() > 50);
    assert_eq!(
        html.matches("class=\"tip bar\" tabindex=\"0\" data-tip=\"")
            .count(),
        html.matches("class=\"tip bar\"").count()
    );
}

#[test]
fn counts_match_the_records() {
    let (results, run, rel) = fixture();
    let html = render_html(&results, &run, Some(&rel), &Options::default());
    let n = |st: &str| results.iter().filter(|r| r["status"] == st).count();
    assert_eq!(n("error"), 1);
    assert!(html.contains(&format!("<div class=\"tile-v\">{}</div>", n("error"))));
    assert!(html.contains(&format!("{} inputs ·", results.len())));
    // the error card shows the file name, never a path, and says the run continued
    assert!(html.contains("empty_file.fa"));
    assert!(html.contains("The run continued"));
}

#[test]
fn deterministic_and_capped() {
    let (results, run, rel) = fixture();
    let a = render_html(&results, &run, Some(&rel), &Options::default());
    let b = render_html(&results, &run, Some(&rel), &Options::default());
    assert_eq!(a, b);
    let small = render_html(
        &results,
        &run,
        Some(&rel),
        &Options {
            max_cards: 3,
            max_rows: 10,
            ..Options::default()
        },
    );
    assert_eq!(small.matches("<details class=\"card\"").count(), 3);
    assert!(small.contains("The table shows the first 10 of"));
    assert!(small.contains("Detail cards are shown for the 3 inputs that need the most attention"));
    // problem inputs come first: the error and the no-match inputs get the cards
    assert!(small.contains("id=\"q-empty_file\""));
}

#[test]
fn reads_card_has_the_abundance_panel() {
    let (results, run, _) = fixture();
    let html = render_html(&results, &run, None, &Options::default());
    assert!(html.contains("Reads: what was kept"));
    assert!(html.contains("hashes seen once: mostly sequencing errors, removed"));
    assert!(
        !html.contains("Relationships among inputs"),
        "no relationships data, no heatmap"
    );
}

#[test]
fn edge_cases_do_not_panic() {
    let run = json!({"version": "0", "commit": "abc", "index": {"n_genomes": 1, "k": 21, "scaled": 1000, "path": "/x/y.gmgi", "bytes": 1},
                     "thresholds": {"name": "novel", "values": {"min_ani": [0.8, 0.9, 0.9, 1.01, 1.01, 1.01], "min_consistency": [0.5, 0.5, 0.5, 0.5, 0.5, 0.5]}}});
    // no inputs
    let html = render_html(&[], &run, None, &Options::default());
    assert!(html.contains("No input could be classified"));
    // only errors
    let err = json!({"id": "bad", "path": "/p/bad.fa", "kind": "assembly", "status": "error", "message": "opening /p/bad.fa: nope", "qc_flags": []});
    let html = render_html(&[err], &run, None, &Options::default());
    assert!(html.contains("bad.fa") && !html.contains("/p/bad.fa"));
    // novel thresholds: species never called
    assert!(
        html.contains("never (this threshold set does not call species)")
            || html.contains("never called")
    );
    // reads only, no ANI
    let reads = json!({"id": "r1", "path": "r1.fastq.gz", "kind": "reads", "status": "no_match", "best_rank": null, "qc": {"n_hashes": 5}, "qc_flags": ["few_hashes"],
                       "guess": {"calls": {}, "top": null, "hits": []}});
    let html = render_html(&[reads], &run, None, &Options::default());
    assert!(html.contains("no match"));
    for s in svgs(&html) {
        roxmltree::Document::parse(s).unwrap();
    }
}

/// Labels must not run past the left/right edge of their chart. Width is estimated at 7px per character (12px system font).
#[test]
fn svg_text_stays_inside_the_viewbox() {
    let (results, run, rel) = fixture();
    let html = render_html(&results, &run, Some(&rel), &Options::default());
    let mut bad = Vec::new();
    for s in svgs(&html) {
        let doc = roxmltree::Document::parse(s).unwrap();
        let vb: Vec<f64> = match doc.root_element().attribute("viewBox") {
            Some(v) => v.split_whitespace().map(|x| x.parse().unwrap()).collect(),
            None => continue,
        };
        for t in doc.descendants().filter(|n| n.has_tag_name("text")) {
            let x: f64 = t.attribute("x").unwrap().parse().unwrap();
            let w = t.text().unwrap_or("").chars().count() as f64 * 7.0;
            let (lo, hi) = match t.attribute("text-anchor") {
                Some("end") => (x - w, x),
                Some("middle") => (x - w / 2.0, x + w / 2.0),
                _ => (x, x + w),
            };
            if lo < -1.0 || hi > vb[2] + 1.0 {
                bad.push(format!(
                    "{} ({lo:.0}..{hi:.0} of {})",
                    t.text().unwrap_or(""),
                    vb[2]
                ));
            }
        }
    }
    assert!(
        bad.is_empty(),
        "text outside its chart: {:?}",
        &bad[..bad.len().min(5)]
    );
}

#[test]
fn tentative_genus_is_shown_as_tentative_everywhere() {
    let (results, run, _) = fixture();
    let rec = json!({
        "id": "tent1", "path": "tent1.fa", "kind": "assembly", "status": "low_confidence", "best_rank": "family",
        "qc": {"total_bases": 3e7, "n50": 1e6, "n_hashes": 3000}, "qc_flags": [],
        "guess": {
            "calls": {
                "phylum": {"taxon": "Ascomycota", "vote": "Ascomycota", "consistency": 1.0},
                "class": {"taxon": "Sordariomycetes", "vote": "Sordariomycetes", "consistency": 1.0},
                "order": {"taxon": "Hypocreales", "vote": "Hypocreales", "consistency": 1.0},
                "family": {"taxon": "Nectriaceae", "vote": "Nectriaceae", "consistency": 1.0},
                "genus": {"taxon": null, "vote": "Fusarium", "consistency": 0.8},
                "species": {"taxon": null, "vote": "Fusarium solani", "consistency": 0.5}
            },
            "tentative_genus": {"taxon": "Fusarium", "ani": 0.85, "consistency": 0.8},
            "top": {"id": "r1", "cluster": "", "shared": 100, "ani": 0.85, "containment_query": 0.033, "containment_ref": 0.03},
            "hits": [{"id": "r1", "taxid": 1, "shared": 100, "lineage": ["Ascomycota", "Sordariomycetes", "Hypocreales", "Nectriaceae", "Fusarium", "Fusarium solani"], "cluster": ""}]
        }
    });
    let mut all = results;
    all.push(rec);
    let html = render_html(&all, &run, None, &Options::default());
    assert!(html.contains("Fusarium (tentative)"), "ladder / header");
    assert!(html.contains("tentative genus: Fusarium"), "table");
    assert!(html.contains("shown as tentative"), "evidence sentence");
    assert!(html.contains("Tentative genus.</b>"), "how-to-read note");
    for s in svgs(&html) {
        roxmltree::Document::parse(s).unwrap();
    }
}
