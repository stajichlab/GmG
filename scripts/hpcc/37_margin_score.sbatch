#!/bin/bash -l
#SBATCH -J gmg_margin_score
#SBATCH -p short
#SBATCH -c 2
#SBATCH --mem=16G
#SBATCH -t 01:00:00
#SBATCH -o /bigdata/stajichlab/shared/projects/GmG_runs/logs/37_margin_score.%j.log
# #58: score the before (default thresholds, index 2026.10.3) and after (margin thresholds, TAG m58) runs with the same scripts:
#   truth set loo/los/log/indep (RefSeq lineage, decoy calls blanked as in 24_truth_score), the 1,257 post-freeze assemblies, IndigoAg.
# usage: REPO=<dir with scripts/ and data/validation, data/decoys> sbatch 37_margin_score.sbatch   (needs pandas)
set -uo pipefail
source $HOME/projects/GmG/scripts/hpcc/env.sh
REPO=${REPO:-$GR/margin/repo}
O=$GR/margin/score; mkdir -p $O
T=$GR/truthset; REL=$GR/release/2026.10.3
for v in before after; do
  for m in loo los log indep; do
    d=$T/runs/$m.known.i3; [ $v = after ] && d=$T/runs/$m.known.m58
    out=$O/$m.$v.results.tsv
    first=$(ls $d/[0-9]*/results.tsv | head -1); head -1 $first > $out
    for f in $d/[0-9]*/results.tsv; do tail -n +2 $f >> $out; done
    python3 $REPO/scripts/blank_decoy_calls.py $out $REL/meta_full.tsv $REPO/data/decoys/decoys.tsv $O/$m.$v.aware.tsv
    python3 $REPO/scripts/validate_run.py $O/$m.$v.aware.tsv $T/sheets/truth_$m.tsv $O/truth_$m.$v \
      --by evidence level species_in_bfd genus_in_bfd phylum --title "Truth set $m, $v the margin (aware)" > /dev/null
  done
  # post-freeze: the 51 batches of external_big
  R=$GR/external_big/rel_2026.10.3; [ $v = after ] && R=$GR/external_big/rel_2026.10.3.m58
  out=$O/post.$v.results.tsv; first=$(ls $R/b*/results.tsv | head -1); head -1 $first > $out
  for f in $R/b*/results.tsv; do tail -n +2 $f >> $out; done
  echo "post-freeze $v: $(($(wc -l < $out)-1)) rows"
  python3 $REPO/scripts/validate_run.py $out $REPO/data/validation/post_freeze.meta.tsv $O/post.$v \
    --by species_in_bfd genus_in_bfd level phylum --title "Post-freeze, $v the margin" > /dev/null
  R=$GR/relstab/indigo_2026.10.3; [ $v = after ] && R=$GR/relstab/indigo_2026.10.3.m58
  [ -s $R/results.tsv ] && python3 $REPO/scripts/score_by_name.py $R/results.tsv $GR/external/indigo_labels.tsv $GR/taxdump_20261029 $O/indigo.$v
done
