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
VER=${VER:-2026.10.1}
R=$GR/release/$VER; mkdir -p $R; cd $R
I=$GR/index
# NCBI taxonomy: lineages and per-rank taxids come from the BFD NCBI taxids (scripts/normalize_taxonomy.py)
TD=$GR/taxdump_$(date +%Y%m%d)
if [ ! -s $TD/nodes.dmp ]; then
  mkdir -p $TD && cd $TD
  curl -sLO https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/new_taxdump/new_taxdump.tar.gz
  curl -sL https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/new_taxdump/new_taxdump.tar.gz.md5 > md5.txt
  tar xzf new_taxdump.tar.gz nodes.dmp names.dmp merged.dmp delnodes.dmp
  cd $R
fi
MD5=$(md5sum $TD/new_taxdump.tar.gz 2>/dev/null | cut -d" " -f1); [ -n "$MD5" ] || MD5=$(cut -d" " -f1 $TD/md5.txt)
TAXV="NCBI taxonomy new_taxdump md5 $MD5 ($(basename $TD | sed s/taxdump_//)); BFD NCBI taxids; original BFD species string kept as label"
python3 $GMG_REPO/scripts/normalize_taxonomy.py $I/meta.tsv $TD meta_full.tsv --clusters $I/clusters_ani0.95.tsv \
  --report taxonomy_normalization.md --taxdump-label "new_taxdump md5 $MD5"
# representatives: the same normalized rows, restricted to the representative ids
python3 - <<PY
import csv
keep = {r["id"] for r in csv.DictReader(open("$I/reps.meta_reps1.tsv"), delimiter="\t")}
rows = list(csv.reader(open("meta_full.tsv"), delimiter="\t"))
with open("meta_reps.tsv", "w") as o:
    o.write("\t".join(rows[0]) + "\n")
    for r in rows[1:]:
        if r[0] in keep: o.write("\t".join(r) + "\n")
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
