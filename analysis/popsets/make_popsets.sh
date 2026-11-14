#!/bin/bash
# Regenerate the four selections in this folder. Run on the HPCC (rule: processing there), from a clone of GmG, with a PopSeq
# checkout (POPSEQ_ROOT, with its venv) whose biosample cache is at /bigdata/stajichlab/shared/projects/Fungi/PopSeq/cache.
#   POPSEQ_ROOT=~/projects/PopSeq_prof_test PYTHON=~/projects/PopSeq/.venv/bin/python bash analysis/popsets/make_popsets.sh WORKDIR
# WORKDIR (purgeable) receives the candidate tables (`popseq runs`, not committed: 76,000 rows for S. cerevisiae) and the
# BioSample attributes joined to them; the selections are written next to this script.
set -euo pipefail
HERE=$(cd "$(dirname "$0")" && pwd)
W=${1:?usage: make_popsets.sh WORKDIR}
: "${POPSEQ_ROOT:?}" "${PYTHON:?}"
CACHE=${BIOSAMPLE_SHARDS:-/bigdata/stajichlab/shared/projects/Fungi/PopSeq/cache/_biosample_shards}
mkdir -p "$W"
cd "$POPSEQ_ROOT"
for sp in "Saccharomyces cerevisiae" "Aspergillus fumigatus" "Coccidioides immitis" "Coccidioides posadasii"; do
  PYTHONPATH=. "$PYTHON" -m popseq.cli runs --species "$sp" --layout PAIRED --min-depth 15 --format csv 2>/dev/null \
    | grep -v '^#' > "$W/cand_${sp// /_}.csv"
done
"$PYTHON" "$HERE/popset_attrs.py" "$W" "$CACHE" "$W/biosample_attrs.csv"
S="$HERE/../../scripts/popset_select.py"; A="$W/biosample_attrs.csv"
"$PYTHON" "$S" "$W/cand_Saccharomyces_cerevisiae.csv" --biosamples "$A" -n 500 --by cultivation --min-depth 25 --max-per-study 40 \
  --quota wine=70 --quota beer=60 --quota sake_rice=40 --quota soy_other_ferment=40 --quota bread_commercial=50 \
  --quota spirits_ethanol=50 --quota dairy=25 --quota clinical=60 --quota wild_nature=70 --quota other_known=20 --quota unannotated=15 \
  --out "$HERE/saccharomyces_500.csv" --summary "$HERE/saccharomyces_500.summary.tsv"
"$PYTHON" "$S" "$W/cand_Aspergillus_fumigatus.csv" --biosamples "$A" -n 150 --by region --min-depth 25 --max-per-study 30 \
  --out "$HERE/afumigatus_pool150.csv" --summary "$HERE/afumigatus_pool150.summary.tsv"
for sp in immitis posadasii; do
  "$PYTHON" "$S" "$W/cand_Coccidioides_$sp.csv" --biosamples "$A" -n 50 --by species --place --min-depth 25 --max-per-study 12 \
    --out "$HERE/coccidioides_${sp}_50.csv" --summary "$HERE/coccidioides_${sp}_50.summary.tsv"
done
