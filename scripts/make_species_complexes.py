#!/usr/bin/env python3
"""Write gmg-core's species -> species-complex table from FusariumID records.

Input: records.tsv from bench/F1_fusariumid/parse_fusariumid.py (columns species, genus, complex, ...).
Output (stdout): TSV `species complex collapse` with a comment header, embedded by gmg-core
(crates/gmg-core/data/species_complexes.tsv).

* complex: the FusariumID group with `F.` expanded and the abbreviation kept, e.g.
  `Fusarium oxysporum species complex (FOSC)`; genera that FusariumID lists as a group (Neocosmospora, Bisifusarium, ...)
  keep the genus name.
* collapse: 1 for complexes whose species gmg must not name from a genome (--collapse, default FOSC: FusariumID species
  in FOSC are 99.6-99.99% genome ANI apart, bench/F1_fusariumid/README.md); gmg then reports the complex instead.
* Neocosmospora species also get a `Fusarium <epithet>` row when no FusariumID species has that name, because NCBI files
  most of them under Fusarium (e.g. Fusarium solani).
* A species listed under several groups keeps the group of most of its records.

usage: make_species_complexes.py records.tsv [--collapse FOSC] [--source TEXT] > species_complexes.tsv
"""
import argparse
import collections
import re
import sys

ABBR = re.compile(r"\(([A-Z][A-Za-z]*SC)\b")
JUNK = ("not assigned", "Art.", "Arts ")


def clean_complex(raw):
    """`F. oxysporum species complex (FOSC)` -> (`Fusarium oxysporum species complex (FOSC)`, `FOSC`)."""
    raw = raw.strip()
    if not raw or any(j in raw for j in JUNK):
        return None, None
    head = raw.split(" (")[0].strip()
    head = re.sub(r"^F\. ", "Fusarium ", head)
    m = ABBR.search(raw)
    abbr = m.group(1) if m else ""
    if "species complex" in head and abbr:
        return f"{head} ({abbr})", abbr
    prev = re.search(r"previously (?:the |part of the )?F\. ([a-z-]+) species complex", raw)
    if prev and " " not in head:
        return f"{head} (formerly Fusarium {prev.group(1)} species complex)", abbr
    return head, abbr


def build(rows, collapse):
    votes = collections.defaultdict(collections.Counter)
    for r in rows:
        name, abbr = clean_complex(r["complex"])
        sp = r["species"]
        # a group that is just the species' own genus (e.g. Atractium) says nothing beyond the genus call
        # informal taxa (`Fusarium sp. 3 in FOSC (Crous et al. 2021)`) stay in: they are labels in FusariumID-labelled
        # references, and leaving them out sent FOSC genomes to "not in a complex" (F1 complex_effect check)
        if name and name != sp.split()[0] and len(sp.split()) >= 2:
            votes[sp][(name, abbr)] += 1
    table = {}
    for sp, c in votes.items():
        (name, abbr), _ = c.most_common(1)[0]
        table[sp] = (name, int(abbr in collapse))
    for sp in list(table):
        g, ep = sp.split()[0], sp.split()[1]
        alias = f"Fusarium {ep}"
        if g == "Neocosmospora" and ep != "sp." and alias not in table:
            table[alias] = table[sp]
    return table


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("records")
    ap.add_argument("--collapse", default="FOSC")
    ap.add_argument("--source", default="FusariumID (fusarium.org) sequence download 2026-06-23")
    a = ap.parse_args()
    with open(a.records) as f:
        cols = f.readline().rstrip("\n").split("\t")
        rows = [dict(zip(cols, l.rstrip("\n").split("\t"))) for l in f]
    table = build(rows, set(a.collapse.split(",")))
    out = sys.stdout
    out.write(f"# species -> species complex; source: {a.source}; written by scripts/make_species_complexes.py\n")
    out.write("# collapse=1: gmg reports the complex instead of a species name (genome sketches cannot separate its species)\n")
    out.write("species\tcomplex\tcollapse\n")
    for sp in sorted(table):
        name, col = table[sp]
        out.write(f"{sp}\t{name}\t{col}\n")
    print(f"{len(table)} species, {sum(v[1] for v in table.values())} collapsed", file=sys.stderr)


if __name__ == "__main__":
    main()
