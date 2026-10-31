#!/bin/bash -l
#SBATCH -J gmg_truth_score_i5
#SBATCH -p short
#SBATCH -c 2
#SBATCH --mem=24G
#SBATCH -t 01:00:00
#SBATCH -o /bigdata/stajichlab/shared/projects/GmG_runs/logs/truth_score_i5.%j.log
# score the 2026.10.5 truth-set runs (TAG i5, default thresholds) as 24_truth_score.sbatch does (RefSeq lineage truth; raw tables)
set -uo pipefail
source $HOME/projects/GmG/scripts/hpcc/env.sh
REPO=$GR/precision/repo; T=$GR/truthset; cd $T
REL=$GR/release/2026.10.5; mkdir -p score_i5; v=known.i5
for m in loo los log indep; do
  d=runs/$m.$v; first=$(ls $d/[0-9]*/results.tsv | head -1)
  head -1 $first > $d.results.tsv; for f in $d/[0-9]*/results.tsv; do tail -n +2 $f >> $d.results.tsv; done
  # GCF_011074885.1 is one of the 72 genomes left out of the reference (it is not in the index), so it is not a leave-out query: dropped, as in the 2026.10.3 tables
  grep -v -P "^GCF_011074885.1\t" $d.results.tsv > $d.tmp && mv $d.tmp $d.results.tsv
  echo "== $m: $(($(wc -l < $d.results.tsv)-1)) rows"
  python3 $REPO/scripts/blank_decoy_calls.py $d.results.tsv $REL/meta_full.tsv $REPO/data/decoys/decoys.tsv $d.aware.tsv
  for k in raw aware; do
    f=$d.results.tsv; [ $k = aware ] && f=$d.aware.tsv
    python3 $REPO/scripts/validate_run.py $f sheets/truth_$m.tsv score_i5/$m.$v.$k.refseq --by evidence level species_in_bfd genus_in_bfd phylum --title "Truth set $m $v $k" > /dev/null
  done
  X=""; [ $m = los ] && X="--excluded-rank species"; [ $m = log ] && X="--excluded-rank genus"
  B=sheets/truth_bfdlabel_$m.tsv; [ $m = indep ] && B=-
  python3 $REPO/scripts/review_disagreements.py $m.$v $d.aware.tsv sheets/truth_$m.tsv $B $REL/meta_full.tsv $GR/taxdump_20261030 $REPO/reports/label_audit score_i5/disagreements_draft.$m.$v.tsv $X --reviews $REPO/data/validation/truth_set_reviews.tsv
done
python3 $REPO/scripts/truth_set_release_table.py score_i5 $v score_i5/tables.md
grep -c . score_i5/tables.md
