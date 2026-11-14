#!/usr/bin/env python3
"""Choose a diverse, reproducible set of runs from a PopSeq candidate table: spread over geography and, for yeast, over how and
where the strain was isolated (wine, beer, sake, soy sauce, bread, clinical, wild ...).

usage:
  popset_select.py CANDIDATES.csv --biosamples ATTRS.csv -n 100 --out SELECTED.csv [--by cultivation|region|species]
                   [--max-per-study N] [--lab exclude|keep] [--quota CATEGORY=N ...] [--seed S] [--summary SUMMARY.tsv]

CANDIDATES.csv is the output of `popseq runs --species ... --format csv` (comments starting with # are skipped). ATTRS.csv has
one row per BioSample (columns biosample, isolation_source, env_broad_scale, env_local_scale, env_medium, geo_loc_name, lat_lon,
host, isolate, strain, culture_collection, sample_name, ...), as harvested by PopSeq into its biosample cache; the join is on the
run's sample_accession. Nothing here queries an archive.

How runs are chosen
* One run per BioSample (the deepest), no run of an unrelated organism: whatever PopSeq's `exclusion` column marks (mutant,
  evolution, cross, not WGS) is left out by `popseq runs` already; here, in addition, samples whose text says laboratory strain,
  evolution experiment, chemostat, FACS, spore, tetrad, zygote, plate, mutant and so on are left out with `--lab exclude`
  (the default), because a lab strain is not a wild or domesticated isolate.
* Each run gets a category (`--by cultivation`: wine, beer, sake_rice, soy_other_ferment, bread_commercial, spirits_ethanol, dairy, clinical,
  wild_nature, other_known, unannotated; `--by region`: continent; `--by species`: the species name) from the keywords in the
  BioSample text, and a region from its country.
* Quotas per category (`--quota wine=60 ...`, default: an equal share per category that has runs, with the shortfall of small
  categories given to the others) and, inside a category, a greedy pick that always takes the run from the least represented
  (country, region, study) so far, so geography and studies are spread; ties are broken by a hash of the run accession and
  `--seed`, so the same input and seed always give the same set. `--max-per-study` caps any one study.
* The classification is by keywords and is a first pass: the summary says how many runs of each category had any text at all
  (only about a quarter of S. cerevisiae BioSamples have an isolation source), and `unannotated` is its own category.
"""
import argparse
import csv
import hashlib
import re
import sys
from collections import Counter, defaultdict

# ---------------------------------------------------------------------------------------------------------------- geography
# ISO 3166 alpha-2 code -> continent (UN geoscheme, simplified), for the countries that occur in fungal isolate collections; any
# other code is "unknown_region".
_CONTINENTS = {
    "Europe": "AD AL AT BA BE BG BY CH CY CZ DE DK EE ES FI FO FR GB GG GI GR HR HU IE IM IS IT JE LI LT LU LV MC MD ME MK MT NL NO PL PT RO RS RU SE SI SK SM UA VA XK",
    "Asia": "AE AF AM AZ BD BH BN BT CN GE HK ID IL IN IQ IR JO JP KG KH KP KR KW KZ LA LB LK MM MN MO MV MY NP OM PH PK PS QA SA SG SY TH TJ TL TM TR TW UZ VN YE",
    "Africa": "AO BF BI BJ BW CD CF CG CI CM CV DJ DZ EG ER ET GA GH GM GN GQ GW KE KM LR LS LY MA MG ML MR MU MW MZ NA NE NG RE RW SC SD SL SN SO SS ST SZ TD TG TN TZ UG ZA ZM ZW",
    "North America": "BM BS BZ CA CR CU DO GL GT HN HT JM MX NI PA PR SV TT US",
    "South America": "AR BO BR CL CO EC FK GF GY PE PY SR UY VE",
    "Oceania": "AU FJ GU NC NZ PF PG SB TO VU WS",
    "Antarctica": "AQ",
}
CONTINENT = {code: cont for cont, codes in _CONTINENTS.items() for code in codes.split()}


def region_of(iso):
    return CONTINENT.get((iso or "").strip().upper(), "unknown_region") if (iso or "").strip() else "unknown_region"


# ---------------------------------------------------------------------------------------------------------------- categories
# Words that mark a laboratory or experimental sample (whole-word match on the lower-cased text).
LAB_WORDS = [
    r"laborator\w*", r"lab", r"s288c", r"by4741", r"by4742", r"by4716", r"w303", r"cen\.pk\w*", r"sk1", r"fy\d+",
    r"evolution", r"evolved", r"chemostat", r"facs\w*", r"spore\w*", r"tetrad\w*", r"zygote\w*", r"haploid cells?",
    r"mitotic", r"meiotic", r"mutant\w*", r"deletion", r"knock-?out", r"transformant\w*", r"plate", r"96-well", r"ypd",
    r"cell culture", r"liquid culture", r"batch culture", r"crossing", r"cross", r"segregant\w*", r"pooled", r"pool",
    r"reference strain", r"\[todo\]", r"yeast extract", r"engineered", r"plasmid\w*", r"pyes\w*", r"sc-(?:his|ura|leu|lys|trp|complete)\w*", r"\w+-derived",
    r"levy et al", r"barcode\w*", r"genomic dna", r"mnase\w*", r"\w+ x \w+",
]
LAB_RE = re.compile(r"\b(?:" + "|".join(LAB_WORDS) + r")\b", re.I)

# category -> keywords (whole words or prefixes, lower case). The first category that matches wins, in this order: the more
# specific uses before the general ones.
CULTIVATION = [
    ("sake_rice", [r"sake", r"shochu", r"makgeolli", r"rice wine", r"koji", r"awamori", r"baijiu", r"huangjiu", r"mijiu", r"rice"]),
    ("soy_other_ferment", [r"soy\w*", r"miso", r"kimchi", r"fermented\w*", r"tempeh", r"cocoa", r"cacao", r"coffee",
                           r"palm wine", r"pulque", r"kefir", r"kombucha", r"cassava", r"vinegar", r"pickle\w*", r"sourdough starter",
                           r"sorghum", r"olive brine", r"alpechin"]),
    ("beer", [r"beer", r"brew\w*", r"\bale\b", r"lager", r"stout", r"porter", r"saison", r"wort", r"hops?", r"craft"]),
    ("wine", [r"wine\w*", r"vineyard\w*", r"grape\w*", r"must", r"winery", r"viticol\w*", r"vin", r"cider", r"fermenting juice",
              r"juice", r"flor", r"sherry"]),
    ("bread_commercial", [r"bread", r"dough", r"bakery", r"baker\w*", r"sourdough", r"active dr\w+ yeast", r"dry yeast", r"liquid yeast",
                          r"baking", r"commer\w+"]),
    ("spirits_ethanol", [r"whisk\w*", r"rum", r"tequila", r"mezcal", r"agave", r"cacha\w+", r"distill\w*", r"ethanol", r"bioethanol",
                         r"sugar ?cane", r"molasses", r"spirits?", r"cognac", r"brandy", r"industrial"]),
    ("dairy", [r"cheese", r"milk", r"dairy", r"yogh?urt"]),
    ("clinical", [r"clinical", r"patient", r"blood", r"sputum", r"stool", r"faec\w*|fec\w*", r"vagin\w*", r"hospital", r"urine",
                  r"human", r"homo sapiens", r"throat", r"respiratory", r"bronch\w*", r"wound", r"infect\w*", r"mucosa\w*"]),
    ("wild_nature", [r"oak", r"quercus", r"bark", r"soil", r"tree", r"forest", r"woodland", r"insect", r"fly", r"drosophila",
                     r"wasp", r"bee", r"beetle", r"fruit\w*", r"flower", r"nectar", r"natural", r"wild", r"water", r"lake", r"river",
                     r"sea", r"leaf", r"leaves", r"environment\w*", r"plant", r"rotting", r"decay\w*", r"cactus", r"sap", r"exudate",
                     r"termite", r"ant", r"bird", r"feather", r"cave", r"glacier", r"ice", r"sewage"]),
]
def _word(p):
    """A keyword matches as a whole word unless it ends in `\\w*` or `\\w+` (a prefix, like wine\\w* for wines)."""
    return p if p.endswith(("\\w*", "\\w+")) else p + r"\b"


_CULT_RE = [(c, re.compile(r"\b(?:" + "|".join(_word(w) for w in ws) + r")", re.I)) for c, ws in CULTIVATION]
TEXT_FIELDS = ["isolation_source", "env_broad_scale", "env_local_scale", "env_medium", "host", "sample_name", "isolate", "strain",
               "culture_collection", "geo_loc_name"]


def text_of(attrs):
    return " ; ".join(attrs.get(f, "") for f in TEXT_FIELDS if attrs.get(f, "")).lower()


def is_lab(text):
    return bool(LAB_RE.search(text))


def cultivation_of(text, has_text):
    """(category, matched keyword). `unannotated` when the BioSample has no isolation text at all, `other_known` when it has
    some that no rule matches."""
    if not has_text:
        return "unannotated", ""
    for cat, rx in _CULT_RE:
        m = rx.search(text)
        if m:
            return cat, m.group(0)
    return "other_known", ""


# ---------------------------------------------------------------------------------------------------------------- selection
def stable_hash(*parts):
    return int(hashlib.sha256("|".join(str(p) for p in parts).encode()).hexdigest()[:12], 16)


def read_candidates(path):
    with open(path, newline="") as fh:
        lines = [ln for ln in fh if not ln.startswith("#")]
    return list(csv.DictReader(lines))


def read_attrs(path):
    with open(path, newline="") as fh:
        return {r["biosample"]: r for r in csv.DictReader(fh)}


def annotate(rows, attrs, by, lab):
    """The candidate rows with `category`, `region`, `text`, `keyword` added; one run (the deepest) per BioSample; lab
    samples removed when lab == "exclude". Returns (rows, counts of what was dropped)."""
    dropped = Counter()
    best = {}
    for r in rows:
        k = r["sample_accession"]
        if k not in best or float(r.get("est_depth_x") or 0) > float(best[k].get("est_depth_x") or 0):
            best[k] = r
    dropped["second run of a sample"] = len(rows) - len(best)
    out = []
    for r in best.values():
        a = attrs.get(r["sample_accession"], {})
        text = text_of(a)
        has_text = bool((a.get("isolation_source") or "").strip() or (a.get("env_broad_scale") or "").strip()
                        or (a.get("env_medium") or "").strip())
        if lab == "exclude" and is_lab(text):
            dropped["laboratory or experimental"] += 1
            continue
        iso = r.get("country_iso", "")
        if not iso and a.get("geo_loc_name"):
            iso = ""  # the country name in geo_loc_name is not mapped to a code here; the region stays unknown
        cat, kw = cultivation_of(text, has_text)
        r = dict(r)
        r["region"] = region_of(iso)
        r["cultivation"] = cat
        r["keyword"] = kw
        r["isolation_source"] = a.get("isolation_source", "")
        r["geo_loc_name"] = a.get("geo_loc_name", "")
        r["category"] = {"cultivation": cat, "region": r["region"], "species": r.get("species_name", "")}[by]
        out.append(r)
    return out, dropped


def allocate(counts, n, quotas):
    """Seats per category: explicit quotas first (capped by availability), then an equal share of what is left among the
    categories that still have runs (quota categories included, so a quota a small category cannot fill goes to the others)."""
    seats = {}
    left = n
    for c, q in quotas.items():
        s = min(q, counts.get(c, 0), left)
        seats[c] = s
        left -= s
    # what the quotas leave (a small category cannot fill its quota) is shared by every category that still has runs
    open_cats = [c for c in sorted(counts) if counts[c] - seats.get(c, 0) > 0]
    while left > 0 and open_cats:
        share = max(1, left // len(open_cats))
        progressed = False
        # a remainder smaller than one seat per category goes to the categories with the most runs left
        for c in sorted(open_cats, key=lambda c: (-(counts[c] - seats.get(c, 0)), c)):
            room = counts[c] - seats.get(c, 0)
            take = min(share, room, left)
            if take > 0:
                seats[c] = seats.get(c, 0) + take
                left -= take
                progressed = True
            if counts[c] - seats.get(c, 0) <= 0:
                open_cats.remove(c)
            if left == 0:
                break
        if not progressed:
            break
    return seats


def place_of(r, use_place):
    """The place the greedy pick spreads over: the country code, or with `--place` the BioSample's geo_loc_name (for example
    `USA: Arizona`, so that a one-country species is spread over its states)."""
    if use_place and r.get("geo_loc_name"):
        return re.sub(r"\s*:\s*", ": ", r["geo_loc_name"].strip().lower())  # "USA:Washington" and "USA: Washington" are one place
    return r["country_iso"] or ""


def pick(rows, seats, seed, max_per_study, use_place=False):
    """Greedy: in each category take, seat by seat, the run from the least represented (country, region, study)."""
    chosen = []
    seen_country, seen_region, seen_study = Counter(), Counter(), Counter()
    by_cat = defaultdict(list)
    for r in rows:
        by_cat[r["category"]].append(r)
    # categories with fewer seats first, so a scarce category gets the first pick of the shared countries and studies
    for cat in sorted(seats, key=lambda c: (seats[c], c)):
        pool = by_cat.get(cat, [])
        for _ in range(seats[cat]):
            best = None
            for r in pool:
                if r.get("_taken"):
                    continue
                if max_per_study and seen_study[r["study_accession"]] >= max_per_study:
                    continue
                place = place_of(r, use_place)
                key = (
                    seen_country[place] if place else 10**6,
                    seen_region[r["region"]],
                    seen_study[r["study_accession"]],
                    float(r.get("est_depth_x") or 0) // 50,  # prefer shallower runs (above --min-depth), in coarse steps: they download less
                    stable_hash(seed, r["run_accession"]),
                )
                if best is None or key < best[0]:
                    best = (key, r)
            if best is None:
                break
            r = best[1]
            r["_taken"] = True
            chosen.append(r)
            seen_country[place_of(r, use_place)] += 1
            seen_region[r["region"]] += 1
            seen_study[r["study_accession"]] += 1
    return chosen


OUT_COLS = ["run_accession", "sample_accession", "study_accession", "species_name", "strain", "category", "cultivation", "keyword",
            "region", "country_iso", "country", "geo_loc_name", "isolation_source", "collection_date", "est_depth_x", "base_count",
            "library_layout", "fastq_ftp"]


def summarize(annotated, chosen, by):
    lines = ["what\tvalue\tcandidates\tselected"]
    cc, cs = Counter(r["category"] for r in annotated), Counter(r["category"] for r in chosen)
    for k in sorted(cc):
        lines.append(f"{by}\t{k}\t{cc[k]}\t{cs.get(k, 0)}")
    rc, rs = Counter(r["region"] for r in annotated), Counter(r["region"] for r in chosen)
    for k in sorted(rc):
        lines.append(f"region\t{k}\t{rc[k]}\t{rs.get(k, 0)}")
    tc, ts = Counter(r["country_iso"] or "none" for r in annotated), Counter(r["country_iso"] or "none" for r in chosen)
    lines.append(f"countries\tdistinct\t{len(tc)}\t{len(ts)}")
    pc, ps_ = Counter(place_of(r, True) or "none" for r in annotated), Counter(place_of(r, True) or "none" for r in chosen)
    lines.append(f"places (geo_loc_name)\tdistinct\t{len(pc)}\t{len(ps_)}")
    sc, ss = Counter(r["study_accession"] for r in annotated), Counter(r["study_accession"] for r in chosen)
    lines.append(f"studies\tdistinct\t{len(sc)}\t{len(ss)}")
    lines.append(f"studies\tlargest share of the selection\t\t{max(ss.values()) / len(chosen):.3f}" if chosen else "studies\tlargest\t\t")
    return "\n".join(lines) + "\n"


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("candidates")
    ap.add_argument("--biosamples", required=True)
    ap.add_argument("-n", type=int, required=True, help="how many runs to choose")
    ap.add_argument("--out", required=True)
    ap.add_argument("--by", choices=["cultivation", "region", "species"], default="region")
    ap.add_argument("--lab", choices=["exclude", "keep"], default="exclude")
    ap.add_argument("--max-per-study", type=int, default=0, help="at most this many runs from one study (0: no cap)")
    ap.add_argument("--quota", action="append", default=[], metavar="CATEGORY=N")
    ap.add_argument("--place", action="store_true", help="spread over the BioSample's geo_loc_name (state or province level) instead of the country")
    ap.add_argument("--min-depth", type=float, default=0.0, help="drop runs with an estimated depth below this")
    ap.add_argument("--seed", type=int, default=1)
    ap.add_argument("--summary")
    a = ap.parse_args(argv)
    quotas = {}
    for q in a.quota:
        k, _, v = q.partition("=")
        quotas[k] = int(v)
    rows = [r for r in read_candidates(a.candidates) if float(r.get("est_depth_x") or 0) >= a.min_depth]
    annotated, dropped = annotate(rows, read_attrs(a.biosamples), a.by, a.lab)
    counts = Counter(r["category"] for r in annotated)
    seats = allocate(counts, a.n, quotas)
    chosen = pick(annotated, seats, a.seed, a.max_per_study, a.place)
    chosen.sort(key=lambda r: (r["category"], r["country_iso"], r["run_accession"]))
    with open(a.out, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=OUT_COLS, extrasaction="ignore", lineterminator="\n")
        w.writeheader()
        w.writerows(chosen)
    text = summarize(annotated, chosen, a.by)
    print(f"{len(rows)} candidates, dropped: {dict(dropped)}; {len(annotated)} usable; selected {len(chosen)} of {a.n}", file=sys.stderr)
    if a.summary:
        with open(a.summary, "w") as fh:
            fh.write(text)
    else:
        print(text, file=sys.stderr)
    return 0 if len(chosen) == a.n else 3


if __name__ == "__main__":
    raise SystemExit(main())
