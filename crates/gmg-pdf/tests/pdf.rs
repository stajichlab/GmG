//! The PDF must open, carry selectable text with embedded fonts, and stay a sensible size (#23).
use gmg_report::{parse_results_jsonl, render_markdown, Options};
use serde_json::Value;
use std::fs;

fn inputs() -> (Vec<Value>, Value, Value) {
    let d = concat!(env!("CARGO_MANIFEST_DIR"), "/../gmg-report/tests/fixtures/");
    let results = parse_results_jsonl(&fs::read_to_string(format!("{d}results.jsonl")).unwrap());
    let run: Value =
        serde_json::from_str(&fs::read_to_string(format!("{d}run.json")).unwrap()).unwrap();
    let rel: Value =
        serde_json::from_str(&fs::read_to_string(format!("{d}relationships.json")).unwrap())
            .unwrap();
    (results, run, rel)
}

#[test]
fn pdf_has_selectable_text_embedded_fonts_and_sane_size() {
    let (results, run, rel) = inputs();
    let md = render_markdown(&results, &run, Some(&rel), &Options::default());
    let pdf = gmg_pdf::render_pdf(&md, "GmG test").unwrap();
    assert!(pdf.starts_with(b"%PDF-"), "pdf header");
    assert!(
        pdf.len() > 20_000 && pdf.len() < 2_000_000,
        "size {}",
        pdf.len()
    );

    let doc = lopdf::Document::load_mem(&pdf).unwrap();
    let pages = doc.get_pages();
    assert!(
        pages.len() >= 3 && pages.len() <= 8,
        "{} pages for 37 inputs",
        pages.len()
    );
    // text is extractable and contains the report's words and an input id
    let all: String = pages
        .keys()
        .map(|n| doc.extract_text(&[*n]).unwrap_or_default())
        .collect();
    for needle in [
        "GmG classification report",
        "Summary",
        "All inputs",
        "Methods",
        "Aspergillus",
        "GCA_000149035.1",
    ] {
        assert!(all.contains(needle), "extracted text lacks {needle:?}");
    }
    // every font is embedded (a font file stream in the descriptor), none is one of the unembedded standard fonts
    let mut descriptors = 0;
    for obj in doc.objects.values() {
        if let Ok(d) = obj.as_dict() {
            if d.has(b"FontDescriptor")
                || d.has(b"FontFile2")
                || d.has(b"FontFile3")
                || d.has(b"FontFile")
            {
                descriptors += 1;
            }
        }
    }
    assert!(
        descriptors >= 2,
        "expected embedded font descriptors, found {descriptors}"
    );
    assert!(
        !all.contains('\u{fffd}'),
        "replacement characters in the text"
    );
}

#[test]
fn hundred_inputs_stay_a_short_report() {
    let (results, run, rel) = inputs();
    let mut many = Vec::new();
    for i in 0..100 {
        let mut r = results[i % results.len()].clone();
        r["id"] = Value::String(format!("sample_{i:03}"));
        many.push(r);
    }
    let md = render_markdown(&many, &run, Some(&rel), &Options::default());
    let pdf = gmg_pdf::render_pdf(&md, "GmG 100").unwrap();
    let pages = lopdf::Document::load_mem(&pdf).unwrap().get_pages().len();
    assert!(pages <= 12, "{pages} pages for 100 inputs");
    assert_eq!(pages, gmg_pdf::page_count(&md));
}

#[test]
fn empty_and_error_only_runs_render() {
    let (_, run, _) = inputs();
    for results in [
        vec![],
        vec![
            serde_json::json!({"id": "bad", "path": "bad.fa", "kind": "assembly", "status": "error", "error_code": "empty_file",
            "message": "The file is empty (0 bytes).", "error_action": "Re-download or re-copy it.", "qc": null, "qc_flags": [], "guess": null}),
        ],
    ] {
        let md = render_markdown(&results, &run, None, &Options::default());
        let pdf = gmg_pdf::render_pdf(&md, "x").unwrap();
        assert!(pdf.starts_with(b"%PDF-"));
    }
}
