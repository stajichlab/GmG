#!/usr/bin/env python3
"""Small and non-fungal inputs for the stress tests (#41): usage fetch_small_inputs.py OUT_DIR

Fungal mitochondrial genomes and plasmids (RefSeq, found with NCBI E-utilities) and viruses (RefSeq nucleotide accessions).
Writes OUT_DIR/<kind>__<accession>.fa and OUT_DIR/small.tsv (id, kind, accession, length, title). Anything that cannot be found is reported
and skipped, never faked.
"""
import csv, json, os, sys, time, urllib.parse, urllib.request

out = sys.argv[1]
os.makedirs(out, exist_ok=True)
EU = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/"


def get(url):
    for t in range(4):
        try:
            return urllib.request.urlopen(url, timeout=60).read().decode()
        except Exception:
            time.sleep(2 * (t + 1))
    return ""


def search(term, n=1):
    time.sleep(0.4)
    r = get(EU + "esearch.fcgi?" + urllib.parse.urlencode({"db": "nuccore", "term": term, "retmax": n, "retmode": "json"}))
    try:
        return json.loads(r)["esearchresult"]["idlist"]
    except Exception:
        return []


def fetch(uid):
    time.sleep(0.4)
    return get(EU + "efetch.fcgi?" + urllib.parse.urlencode({"db": "nuccore", "id": uid, "rettype": "fasta", "retmode": "text"}))


rows = []


def save(kind, text):
    if not text.startswith(">"):
        return
    head = text.split("\n", 1)[0]
    acc = head[1:].split()[0]
    seq = "".join(text.split("\n")[1:])
    name = f"{kind}__{acc}"
    open(f"{out}/{name}.fa", "w").write(text)
    rows.append([name, kind, acc, len(seq), head[1:][:120]])


for org in ["Saccharomyces cerevisiae", "Neurospora crassa", "Aspergillus nidulans", "Candida albicans", "Schizosaccharomyces pombe",
            "Fusarium graminearum", "Cryptococcus neoformans", "Rhizopus delemar", "Agaricus bisporus", "Batrachochytrium dendrobatidis",
            "Trichoderma reesei", "Magnaporthe oryzae"]:
    ids = search(f'"{org}"[Organism] AND mitochondrion[Title] AND srcdb_refseq[PROP] AND complete genome[Title]')
    if not ids:
        ids = search(f'"{org}"[Organism] AND mitochondrion[Title] AND complete genome[Title]')
    if ids:
        save("mito", fetch(ids[0]))
    else:
        print("no mitochondrial genome found for", org, file=sys.stderr)
for uid in search('plasmid[Title] AND (Saccharomyces[Organism] OR Kluyveromyces[Organism] OR Zygosaccharomyces[Organism]) AND srcdb_refseq[PROP]', 6):
    save("plasmid", fetch(uid))
for acc in ["NC_001416.1", "NC_000866.4", "NC_014649.1", "NC_000852.5", "NC_006998.1", "NC_001806.2", "NC_001422.1"]:
    t = fetch(acc)
    if t.startswith(">"):
        save("virus", t)
    else:
        print("virus", acc, "not fetched", file=sys.stderr)
with open(f"{out}/small.tsv", "w") as f:
    w = csv.writer(f, delimiter="\t", lineterminator="\n")
    w.writerow(["id", "kind", "accession", "length", "title"])
    w.writerows(rows)
print(len(rows), "sequences", file=sys.stderr)
