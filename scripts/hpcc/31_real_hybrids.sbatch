#!/bin/bash -l
#SBATCH -J gmg_hybrids
#SBATCH -p short
#SBATCH -c 4
#SBATCH --mem=16G
#SBATCH -t 01:00:00
#SBATCH -o /bigdata/stajichlab/shared/projects/GmG_runs/logs/31_real_hybrids.%j.log
# Real interspecies hybrids (#41): assemblies of known hybrids from NCBI, classified against the release index; truth = the genus of the parents.
set -uo pipefail
source $HOME/projects/GmG/scripts/hpcc/env.sh
GMG_REPO=${STRESS2_REPO:-$GMG_REPO}
module load ncbi_datasets/18.30.1
INDEX=${INDEX:-$(ls -t $GR/release/2026.10.4/bfd-2026.10.4-s10000.gmgi $GR/release/2026.10.3/bfd-2026.10.3-s10000.gmgi 2>/dev/null | head -1)}
echo "index: $INDEX"; echo "gmg: $(readlink -f $GMG)"
R=$GR/stress2; W=$R/work_hybrids; rm -rf $W; mkdir -p $W/g; cd $W
: > acc.tsv
# name -> expected genus; the best assembly (reference or chromosome level first) of each
while IFS='|' read -r name genus; do
  datasets summary genome taxon "$name" --as-json-lines --limit 5 2>/dev/null > s.jsonl || continue
  python3 - "$name" "$genus" >> acc.tsv <<'PY'
import json, sys
name, genus = sys.argv[1:3]
rows = [json.loads(l) for l in open("s.jsonl") if l.strip()]
rank = {"Complete Genome": 0, "Chromosome": 1, "Scaffold": 2, "Contig": 3}
rows.sort(key=lambda r: (rank.get(r["assembly_info"].get("assembly_level"), 9), -int(r.get("assembly_stats", {}).get("total_sequence_length", 0) or 0)))
if rows:
    r = rows[0]
    print("\t".join([r["accession"], genus, name, r["assembly_info"].get("assembly_level", ""), r.get("organism", {}).get("organism_name", "")]))
PY
done <<'LIST'
Saccharomyces pastorianus|Saccharomyces
Saccharomyces bayanus|Saccharomyces
Saccharomyces cerevisiae x Saccharomyces kudriavzevii|Saccharomyces
Candida orthopsilosis|Candida
Cryptococcus neoformans x Cryptococcus deneoformans hybrid|Cryptococcus
Zygosaccharomyces parabailii|Zygosaccharomyces
Verticillium longisporum|Verticillium
Saccharomyces eubayanus|Saccharomyces
LIST
cat acc.tsv
cut -f1 acc.tsv | sort -u > acc.txt
[ -s acc.txt ] || { echo "no hybrid assemblies found"; exit 0; }
datasets download genome accession --inputfile acc.txt --include genome --filename h.zip --no-progressbar || exit 1
unzip -oq h.zip -d hx
for f in hx/ncbi_dataset/data/GC*/*.fna; do mv "$f" g/$(basename $(dirname "$f")).fna; done
$GMG --threads ${SLURM_CPUS_PER_TASK} run g --index $INDEX -o $R/out_hybrids --no-report --top 10
cp acc.tsv $R/hybrids_acc.tsv
# most of these assemblies are in the reference: hide each one's own species, then its genus (#85), for a test that is not a lookup
META=$(dirname $INDEX)/meta_full.tsv
python3 - $META <<PY
import csv, sys
ids = [r["id"] for r in csv.DictReader(open(sys.argv[1]), delimiter="\t")]
with open("sheet.tsv", "w") as o:
    o.write("id\tpath\texpected\tref\n")
    for l in open("acc.tsv"):
        acc, genus, name, level, org = l.rstrip("\n").split("\t")
        ref = next((i for i in ids if i.startswith(acc + "_") or i == acc), "")
        o.write(f"{acc}\tg/{acc}.fna\t{genus}\t{ref}\n")
PY
cat sheet.tsv | cut -c1-120
for lv in species genus; do
  $GMG --threads ${SLURM_CPUS_PER_TASK} run --sample-sheet sheet.tsv --index $INDEX -o $R/out_hybrids_$lv --no-report --top 10 --exclude-level $lv -q
done
python3 - <<PY
import csv
for lv in ("species", "genus"):
    print("--- own", lv, "hidden")
    print("accession | expected genus | status | call (deepest) | genus call | top hit species | ANI | flags")
    for r in csv.DictReader(open("$R/out_hybrids_%s/results.tsv" % lv), delimiter="\t"):
        deepest = next((r["call_" + x] for x in ("species", "genus", "family", "order", "class", "phylum") if r["call_" + x]), "")
        print(" | ".join([r["query_id"], r.get("expected", ""), r["status"], deepest, r["call_genus"], r["top_hit_species"], r["ani"][:6], r["qc_flags"]]))
PY
python3 - <<PY
import csv
res = {r["query_id"]: r for r in csv.DictReader(open("$R/out_hybrids/results.tsv"), delimiter="\t")}
print("accession | named as | level | expected genus | status | genus call | species call | ANI | flags")
for l in open("acc.tsv"):
    acc, genus, name, level, org = l.rstrip("\n").split("\t")
    r = res.get(acc)
    if r:
        print(" | ".join([acc, org or name, level, genus, r["status"], r["call_genus"], r["call_species"], r["ani"][:6], r["qc_flags"]]))
PY
rm -rf $W
