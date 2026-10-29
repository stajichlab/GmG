#!/bin/bash -l
#SBATCH -J gmg_release
#SBATCH -p short
#SBATCH -c 8
#SBATCH --mem=48G
#SBATCH -t 02:00:00
#SBATCH -o /bigdata/stajichlab/shared/projects/GmG_runs/logs/15_release.%j.log
# Production indexes for a GitHub release (docs/index-distribution.md): full BFD s=10000 and the representatives-only s=20000,
# both with the ANI-cluster column and a manifest, in the compact v1 format.
set -euo pipefail
source $HOME/projects/GmG/scripts/hpcc/env.sh
VER=${VER:-2026.10}
TAXV="BFD samples.csv lineages (Fungi_BFD_runs), NCBI taxon ids; snapshot 2026-10-27"
R=$GR/release/$VER; mkdir -p $R; cd $R
I=$GR/index
# meta + cluster column (greedy centroid, containment-ANI >= 0.95; gmg cluster)
python3 - <<PY
import csv
cl = {r["id"]: r["cluster"] for r in csv.DictReader(open("$I/clusters_ani0.95.tsv"), delimiter="\t")}
for src, dst in (("$I/meta.tsv", "meta_full.tsv"), ("$I/reps.meta_reps1.tsv", "meta_reps.tsv")):
    rows = list(csv.reader(open(src), delimiter="\t"))
    with open(dst, "w") as o:
        o.write("\t".join(rows[0] + ["cluster"]) + "\n")
        for r in rows[1:]:
            o.write("\t".join(r + [cl.get(r[0], "")]) + "\n")
PY
$GMG --threads 8 build --meta meta_full.tsv --sketches $GR/sketch --scaled 10000 -o full_raw.gmgi \
  --taxonomy-version "$TAXV" --note "full reference; ANI clusters at 0.95; includes 43 genomes flagged as possible mislabels (reports/label_audit/exclude_ids.txt)"
$GMG convert full_raw.gmgi -o bfd-$VER-s10000.gmgi
$GMG --threads 8 build --meta meta_reps.tsv --sketches $GR/sketch --scaled 20000 -o reps_raw.gmgi \
  --taxonomy-version "$TAXV" --note "representatives only (best genome per species, BFD rule); ANI clusters at 0.95 where the centroid is in the full index"
$GMG convert reps_raw.gmgi -o bfd-reps-$VER-s20000.gmgi
rm -f full_raw.gmgi reps_raw.gmgi
bash $GMG_REPO/scripts/make_index_release.sh bfd-$VER-s10000.gmgi bfd-reps-$VER-s20000.gmgi
ls -la $R
