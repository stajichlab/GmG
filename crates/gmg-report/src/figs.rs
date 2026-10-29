//! The summary figures shared by the HTML page and the Markdown report.
use crate::charts::{self, Chart, Depth};
use crate::model::Rec;
use crate::page::Th;
use crate::svg::{esc, short};
use serde_json::Value;

pub struct Figures {
    pub depth_chart: Chart,
    pub ph_chart: Chart,
    pub genus_chart: Chart,
    pub genus_legend: String,
    pub n_genus: usize,
    pub ani_chart: Chart,
    pub low_tbl: String,
    /// (chart, number of assemblies) when relationships were supplied
    pub heat: Option<(Chart, usize)>,
}

pub fn figures(recs: &[Rec], depth: &Depth, th: &Th, relationships: Option<&Value>) -> Figures {
    let depth_chart = charts::depth_chart(depth);
    let ph_items = charts::count_by(recs.iter().filter(|r| r.status != "error").map(|r| {
        if r.phylum().is_empty() {
            "(not called)"
        } else {
            r.phylum()
        }
    }));
    let ph_colors: Vec<String> = ph_items.iter().map(|(n, _)| charts::ph_color(n)).collect();
    let ph_chart = charts::hbars(&ph_items, "Phylum", Some(&ph_colors));
    let mut genus_ph: std::collections::BTreeMap<
        String,
        std::collections::BTreeMap<String, usize>,
    > = Default::default();
    for r in recs
        .iter()
        .filter(|r| r.status != "error" && !r.call(4).is_empty())
    {
        *genus_ph
            .entry(r.call(4).to_string())
            .or_default()
            .entry(if r.phylum().is_empty() {
                "(not called)".into()
            } else {
                r.phylum().to_string()
            })
            .or_default() += 1;
    }
    let mut genus_items: Vec<(String, String, usize)> = genus_ph
        .into_iter()
        .map(|(g, m)| {
            let n: usize = m.values().sum();
            let ph = m
                .iter()
                .max_by(|a, b| a.1.cmp(b.1).then(b.0.cmp(a.0)))
                .map(|(p, _)| p.clone())
                .unwrap_or_default();
            (g, ph, n)
        })
        .collect();
    genus_items.sort_by(|a, b| b.2.cmp(&a.2).then(a.0.cmp(&b.0)));
    genus_items.truncate(12);
    let (genus_chart, genus_legend) = charts::genus_by_phylum(&genus_items);

    let mut anis: Vec<f64> = recs
        .iter()
        .filter(|r| r.kind == "assembly")
        .filter_map(|r| r.ani())
        .collect();
    if anis.is_empty() {
        anis = recs.iter().filter_map(|r| r.ani()).collect();
    }
    let ani_chart = charts::ani_hist(&anis, th.min_ani[4], th.min_ani[5]);
    let mut low: Vec<&Rec> = recs
        .iter()
        .filter(|r| r.status != "error" && r.ani().is_some())
        .collect();
    low.sort_by(|a, b| {
        a.ani()
            .partial_cmp(&b.ani())
            .unwrap_or(std::cmp::Ordering::Equal)
    });
    let low_tbl: String = if low.is_empty() {
        String::new()
    } else {
        let rows: String = low
            .iter()
            .take(5)
            .map(|r| {
                format!(
                    r#"<tr><td class="mono">{}</td><td>{}</td><td class="num">{:.3}</td></tr>"#,
                    esc(&short(&r.id, 22)),
                    if r.call_name().is_empty() {
                        "no call".into()
                    } else {
                        esc(r.call_name())
                    },
                    r.ani().unwrap_or(0.0)
                )
            })
            .collect();
        format!(
            r#"<h3>Lowest ANI to a reference: look at these first</h3><table class="small full"><thead><tr><th>Input</th><th>Call</th><th class="num">ANI</th></tr></thead><tbody>{rows}</tbody></table>"#
        )
    };

    let heat = relationships.and_then(|rel| {
        let ids: Vec<String> = rel["ids"]
            .as_array()?
            .iter()
            .map(|v| v.as_str().unwrap_or("").to_string())
            .collect();
        let labels: Vec<String> = rel["labels"]
            .as_array()?
            .iter()
            .map(|v| v.as_str().unwrap_or("").to_string())
            .collect();
        let ani: Vec<Vec<f64>> = rel["ani"]
            .as_array()?
            .iter()
            .map(|row| {
                row.as_array()
                    .map(|r| r.iter().map(|x| x.as_f64().unwrap_or(0.0)).collect())
                    .unwrap_or_default()
            })
            .collect();
        if ids.len() < 2 || ids.len() != ani.len() || ani.iter().any(|r| r.len() != ids.len()) {
            return None;
        }
        Some((charts::heatmap(&ids, &labels, &ani), ids.len()))
    });

    Figures {
        depth_chart,
        ph_chart,
        genus_chart,
        genus_legend,
        n_genus: genus_items.len(),
        ani_chart,
        low_tbl,
        heat,
    }
}
