//! GmG report engine: turns the records of `gmg run` (`results.jsonl`, `run.json`, optional relationships) into one standalone
//! HTML page with inline SVG charts. Pure functions, no I/O: the same code runs in the CLI and in the browser (WASM).
//!
//! Design: `docs/report-mockup/` (issue #22). Every number in the report is a column of results.tsv.

mod charts;
mod model;
mod page;
mod svg;

pub use model::{file_names_only, Rec};

use serde_json::Value;

#[derive(Clone, Debug)]
pub struct Options {
    /// Detail cards are written for at most this many inputs (problems first).
    pub max_cards: usize,
    /// The table lists at most this many inputs (all are in results.tsv).
    pub max_rows: usize,
    /// Relative link for the "Download results.tsv" button; None hides it (e.g. in the browser).
    pub tsv_link: Option<String>,
    /// A line shown above the page (used by the design mockup); None in real reports.
    pub banner: Option<String>,
    pub title: String,
}

impl Default for Options {
    fn default() -> Self {
        Options {
            max_cards: 400,
            max_rows: 5000,
            tsv_link: Some("results.tsv".into()),
            banner: None,
            title: "GmG classification report".into(),
        }
    }
}

/// `results`: one JSON value per input (lines of results.jsonl); `run`: run.json; `relationships`: optional
/// `{ids, labels, ani}` (pairwise ANI between up to 60 assemblies).
pub fn render_html(
    results: &[Value],
    run: &Value,
    relationships: Option<&Value>,
    opts: &Options,
) -> String {
    page::render(results, run, relationships, opts)
}

/// Parse `results.jsonl` text (one JSON object per line; blank/invalid lines are skipped).
pub fn parse_results_jsonl(text: &str) -> Vec<Value> {
    text.lines()
        .filter_map(|l| serde_json::from_str::<Value>(l).ok())
        .filter(|v| v.get("id").is_some())
        .collect()
}
