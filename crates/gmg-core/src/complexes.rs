//! Species complexes: named groups of species (e.g. the *Fusarium oxysporum* species complex, FOSC) that gmg reports
//! next to the species, and for some (`collapse`) instead of it.
//!
//! The table is embedded (`data/species_complexes.tsv`, written by `scripts/make_species_complexes.py` from the FusariumID
//! database) so the CLI, Python and wasm builds behave the same. Why FOSC collapses: FusariumID species inside FOSC are
//! 99.6-99.99% genome ANI apart and 38% of FOSC genomes have a genome of another FOSC species above the species ANI cutoff
//! (`bench/F1_fusariumid/README.md`), so a species name from a sketch would be a coin toss between close names.

use std::collections::HashMap;
use std::sync::OnceLock;

const TABLE: &str = include_str!("../data/species_complexes.tsv");

/// A complex is reported only when the top hit is at least this close (ANI). FOSC genomes sit >= 0.99 from their nearest
/// FOSC reference; other complexes of the genus are ~0.90-0.95 away, so 0.98 keeps an unknown relative from inheriting a
/// complex it may not belong to.
pub const COMPLEX_MIN_ANI: f64 = 0.98;

#[derive(Clone, Debug, PartialEq)]
pub struct ComplexEntry {
    pub complex: String,
    pub collapse: bool,
}

fn table() -> &'static HashMap<String, ComplexEntry> {
    static T: OnceLock<HashMap<String, ComplexEntry>> = OnceLock::new();
    T.get_or_init(|| parse(TABLE))
}

/// Parse `species<TAB>complex<TAB>collapse` lines (header and `#` comments skipped).
pub fn parse(text: &str) -> HashMap<String, ComplexEntry> {
    text.lines()
        .filter(|l| !l.starts_with('#') && !l.starts_with("species\t") && !l.trim().is_empty())
        .filter_map(|l| {
            let mut f = l.split('\t');
            let sp = f.next()?.trim();
            let cx = f.next()?.trim();
            let col = f.next().map(|c| c.trim() == "1").unwrap_or(false);
            (!sp.is_empty() && !cx.is_empty()).then(|| {
                (
                    sp.to_string(),
                    ComplexEntry {
                        complex: cx.to_string(),
                        collapse: col,
                    },
                )
            })
        })
        .collect()
}

/// The complex of a species name, if the embedded table lists it.
pub fn lookup(species: &str) -> Option<&'static ComplexEntry> {
    table().get(species)
}

/// Number of species in the embedded table (for `gmg info` and tests).
pub fn n_species() -> usize {
    table().len()
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn embedded_table_has_fosc_collapsed() {
        let e = lookup("Fusarium oxysporum").expect("F. oxysporum listed");
        assert!(e.collapse);
        assert!(e.complex.contains("FOSC"));
        let g = lookup("Fusarium graminearum").expect("F. graminearum listed");
        assert!(!g.collapse);
        assert!(lookup("Aspergillus fumigatus").is_none());
        assert!(n_species() > 500);
    }

    #[test]
    fn parse_skips_comments_and_header() {
        let t = parse("# c\nspecies\tcomplex\tcollapse\nA b\tX complex\t1\nC d\tY\t0\nbad\n");
        assert_eq!(t.len(), 2);
        assert!(t["A b"].collapse && !t["C d"].collapse);
    }
}
