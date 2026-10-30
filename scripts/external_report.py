#!/usr/bin/env python3
"""Aggregate report for one external validation set (#40), written to REPORT_DIR/external_<set>.md (aggregates only).

usage: external_report.py SET EXTERNAL_ROOT REPORT_DIR
Reads EXTERNAL_ROOT/<set>/{results.tsv,labels.tsv,score.tsv} (made by scripts/hpcc/20_external_set.sbatch).
Sets: kasson, zygolife (provider names: genus precision from score_by_name.py plus species-name agreement),
antarctica (classification only: what the bins look like), ccfee (aaftf vs shovill agreement of the same strains).
"""
import collections, csv, sys

name, root, outdir = sys.argv[1:4]
RANKS = ["phylum", "class", "order", "family", "genus", "species"]


def rd(p):
    return list(csv.DictReader(open(p), delimiter="\t"))


def counts(rows, col):
    c = collections.Counter(r[col] or "(none)" for r in rows)
    return ", ".join(f"{k} {v}" for k, v in c.most_common())


def flagshare(rows):
    c = collections.Counter(f for r in rows for f in (r["qc_flags"].split(",") if r["qc_flags"] else []))
    return ", ".join(f"{k} {v}" for k, v in c.most_common(12))


def overview(rows, title):
    n = len(rows)
    md = [f"## {title}\n", f"{n} inputs. Status: {counts(rows, 'status')}.\n",
          f"Deepest rank called: {counts(rows, 'best_rank')}.\n",
          "| rank | called | share |", "|---|---|---|"]
    for k in RANKS:
        c = sum(bool(r[f"call_{k}"]) for r in rows)
        md.append(f"| {k} | {c} | {c / n:.3f} |")
    md.append(f"\nPhyla called: {counts(rows, 'call_phylum')}.\n")
    md.append(f"QC flags: {flagshare(rows) or 'none'}.\n")
    ident = sum("identical_to_reference" in r["qc_flags"] for r in rows)
    md.append(f"`identical_to_reference` (the input, or a copy, is in the reference): {ident} of {n}.\n")
    return md


def norm(s):
    return " ".join(s.replace("[", "").replace("]", "").split()).lower()


def species_agreement(res, labels):
    """Exact species-name agreement: the provider name has genus + epithet (not sp./cf.), the call is that species."""
    lab = {r["sample"]: r["name"] for r in labels}
    tot = called = ok = 0
    wrong = []
    for r in res:
        w = norm(lab.get(r["query_id"], "")).split()
        if len(w) < 2 or w[1] in ("sp.", "sp", "cf.", "aff.", "spp.") or "sp." in w[:3]:
            continue
        tot += 1
        if r["call_species"]:
            called += 1
            if norm(r["call_species"]).split()[:2] == w[:2]:
                ok += 1
            else:
                wrong.append(f"{' '.join(w[:2])} -> {r['call_species']}")
    return tot, called, ok, wrong


def scored(set_dir):
    try:
        return open(f"{set_dir}/score.md").read().split("\n", 1)[1]
    except OSError:
        return ""


sd = f"{root}/{name}"
md = [f"# External set: {name}\n"]
if name == "ccfee":
    a = {r["query_id"]: r for r in rd(f"{root}/ccfee_aaftf/results.tsv")}
    s = {r["query_id"]: r for r in rd(f"{root}/ccfee_shovill/results.tsv")}
    both = sorted(set(a) & set(s))
    md.append(f"Two assemblies of the same strains (AAFTF and shovill): {len(a)} and {len(s)} inputs, {len(both)} strains with both.\n")
    md += overview(list(a.values()), "AAFTF assemblies")
    md += overview(list(s.values()), "shovill assemblies")
    md += ["## Agreement of the two assemblies of one strain\n", "| rank | both call | same call | disagree | one only |", "|---|---|---|---|---|"]
    for k in RANKS:
        bc = [x for x in both if a[x][f"call_{k}"] and s[x][f"call_{k}"]]
        same = sum(a[x][f"call_{k}"] == s[x][f"call_{k}"] for x in bc)
        one = sum(bool(a[x][f"call_{k}"]) != bool(s[x][f"call_{k}"]) for x in both)
        md.append(f"| {k} | {len(bc)} | {same} | {len(bc) - same} | {one} |")
    dis = [(x, a[x]["call_genus"], s[x]["call_genus"]) for x in both if a[x]["call_genus"] and s[x]["call_genus"] and a[x]["call_genus"] != s[x]["call_genus"]]
    md.append(f"\nGenus disagreements between the assemblers: {len(dis)}" + (": " + ", ".join(f"{x} ({p} vs {q})" for x, p, q in dis) if dis else "") + ".\n")
    for tag, rows in (("AAFTF", a), ("shovill", s)):
        lab = {r["sample"]: r["note"] for r in rd(f"{root}/ccfee_{tag.lower()}/labels.tsv")}
        ok = tot = 0
        for x, r in rows.items():
            ph = lab.get(x, "").replace("phylum=", "")
            if ph and r["call_phylum"]:
                tot += 1
                ok += ph == r["call_phylum"]
        md.append(f"Phylum call vs the samples.csv Phylum column, {tag}: {ok} of {tot} called agree.\n")
else:
    res, labels = rd(f"{sd}/results.tsv"), rd(f"{sd}/labels.tsv")
    md += overview(res, "Classification")
    if any(r["name"] for r in labels):
        md.append("## Against the provider names\n")
        md.append("Truth = the NCBI lineage of the provider genus (names.dmp scientific name or synonym, pinned taxdump); "
                  "see scripts/score_by_name.py.\n")
        md.append(scored(sd))
        t, c, ok, wrong = species_agreement(res, labels)
        if t:
            md.append(f"\nSpecies, exact name (names with an epithet only): {t} inputs, {c} called, {ok} equal to the provider species"
                      + (f" ({ok / c:.3f} of called)" if c else "") + ".")
            if wrong:
                md.append("Species calls that differ: " + "; ".join(wrong[:20]) + ".")
        notes = collections.Counter(r["note"] for r in labels if "metagenome" in r["note"].lower())
        if notes:
            md.append("")
            meta = {r["sample"] for r in labels if "metagenome" in r["note"].lower()}
            mm = [r for r in res if r["query_id"] in meta]
            md.append(f"Metagenome-derived assemblies (herbarium or mixed DNA): {len(mm)} inputs, status {counts(mm, 'status')}; "
                      f"genus called for {sum(bool(r['call_genus']) for r in mm)}.")
    if name == "antarctica":
        md.append("\n## What the bins look like\n")
        for flag in ("not_fungal", "distant_non_fungal_hit"):
            ids = [r["query_id"] for r in res if flag in r["qc_flags"]]
            md.append(f"* `{flag}`: {len(ids)} bins" + (": " + ", ".join(ids[:40]) if ids else "") + ".")
        nm = [r["query_id"] for r in res if r["status"] == "no_match"]
        md.append(f"* `no_match`: {len(nm)} bins" + (": " + ", ".join(nm[:40]) if nm else "") + ".")
        lc = [r for r in res if r["status"] == "low_confidence"]
        md.append(f"* `low_confidence`: {len(lc)} bins; phyla: {counts(lc, 'call_phylum')}.")
        md.append("\nNo names are known for these bins, so this is a classification report, not a validation: a call is a "
                  "hypothesis. Bins with few hashes (small or fragmentary) are expected to get no call.")
open(f"{outdir}/external_{name}.md", "w").write("\n".join(md) + "\n")
print("\n".join(md))
