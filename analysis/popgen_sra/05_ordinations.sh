#!/bin/bash
# usage: 05_ordinations.sh SLUG   (from the popsra folder)
# Runs popsketch twice on the finished sketches, hash PCA and PCoA, with labels (country first, then study) from runs.tsv,
# and writes the Markdown report with its SVG figures (PCA/PCoA plots and tree) beside the HTML.
set -euo pipefail
SLUG=${1:?slug}; BIN=${BIN:-../popsketch_build/target/release/popsketch}; T=${THREADS:-8}
python3 - "$SLUG" <<'PY'
import csv, sys, os
s = sys.argv[1]
have = {f[:-5] for f in os.listdir(f"{s}/sketches") if f.endswith(".json")}
rows = list(csv.DictReader(open(f"{s}/runs.tsv"), delimiter="\t"))
with open(f"{s}_labels.tsv", "w") as o:
    o.write("id\tcountry\tstudy\n")
    for r in rows:
        if r["run_accession"] in have:
            c = (r.get("country") or "").split(":")[0].strip()
            o.write(f'{r["run_accession"]}\t{c}\t{r["study_accession"]}\n')
PY
for ord in hash-pca pcoa; do
  $BIN run --from-sketches $SLUG/sketches -o ${SLUG}_${ord//-/} --labels ${SLUG}_labels.tsv \
     --ordination $ord --report-format all --threads $T -q
done
