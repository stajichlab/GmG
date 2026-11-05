#!/bin/bash
# sketch the Coccidioides reference genomes at scaled=200 and build the reference hash union
set -euo pipefail
cd /bigdata/stajichlab/jstajich/projects/GmG_popgen/popgen_cocci
B=/bigdata/stajichlab/shared/projects/BFD/Fungi_BFD_runs
GMG=/bigdata/stajichlab/shared/projects/GmG_runs/bin/gmg
python3 - <<PY > cocci_ref_ids.txt
import csv
for r in csv.DictReader(open("$B/samples.csv")):
    if r["GENUS"] == "Coccidioides": print(r["ASMID"])
PY
mkdir -p refsk200
$GMG sketch -q --threads 4 -k 21 -s 200 -o refsk200 $(sed "s#^#$B/input_clean_genomes/#; s#\$#.fa.gz#" cocci_ref_ids.txt)
python3 02a_refhashes.py $B/samples.csv refsk200 refs_cocci_s200.npy
