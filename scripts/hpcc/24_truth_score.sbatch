#!/bin/bash -l
#SBATCH -J gmg_truth_score
#SBATCH -p short
#SBATCH -c 2
#SBATCH --mem=24G
#SBATCH -t 01:00:00
#SBATCH -o /bigdata/stajichlab/shared/projects/GmG_runs/logs/24_truth_score.%j.log
# #39: merge the shards of each run variant, score it against the RefSeq lineage (trusted truth) and against the index label of the
# same genome, with and without the phylum calls that only name a distant decoy (scripts/blank_decoy_calls.py), list every
# disagreement with its evidence (scripts/review_disagreements.py) and build the report tables (scripts/truth_set_tables.py).
# usage: REPO=<dir with scripts/ and reports_label_audit/> sbatch 24_truth_score.sbatch     (PREV=<earlier disagreements.tsv> for the prev_category column)
# needs pandas (validate_run.py); the committed numbers were produced with the same scripts on a workstation.
set -uo pipefail
source $HOME/projects/GmG/scripts/hpcc/env.sh
REPO=${REPO:-$HOME/projects/GmG}
AUDIT=${AUDIT:-$REPO/reports/label_audit}
REVIEWS=${REVIEWS:-$REPO/data/validation/truth_set_reviews.tsv}
VARIANTS=${VARIANTS:-"known novel known.nodecoy1 known.noaudit known.i3 novel.i3"}   # .i3 = index 2026.10.3 (production), the rest the first run on 2026.10.2 / 2026.10.1
T=$GR/truthset; cd $T
REL=$GR/release/2026.10.2
mkdir -p score
for v in $VARIANTS; do
  for m in loo los log indep; do
    d=runs/$m.$v
    [ -d $d ] || continue
    first=$(ls $d/[0-9]*/results.tsv | head -1)
    head -1 $first > $d.results.tsv
    for f in $d/[0-9]*/results.tsv; do tail -n +2 $f >> $d.results.tsv; done
    echo "== $m $v: $(($(wc -l < $d.results.tsv)-1)) rows"
    python3 $REPO/scripts/blank_decoy_calls.py $d.results.tsv $REL/meta_full.tsv $REPO/data/decoys/decoys.tsv $d.aware.tsv
    for k in raw aware; do
      f=$d.results.tsv; [ $k = aware ] && f=$d.aware.tsv
      python3 $REPO/scripts/validate_run.py $f sheets/truth_$m.tsv score/$m.$v.$k.refseq \
        --by evidence level species_in_bfd genus_in_bfd phylum --title "Truth set $m $v $k, RefSeq lineage as truth" > /dev/null
      [ $m = indep ] || python3 $REPO/scripts/validate_run.py $f sheets/truth_bfdlabel_$m.tsv score/$m.$v.$k.bfdlabel \
        --by evidence level species_in_bfd genus_in_bfd phylum --title "Truth set $m $v $k, index label as truth" > /dev/null
    done
    X=""; [ $m = los ] && X="--excluded-rank species"; [ $m = log ] && X="--excluded-rank genus"
    B=sheets/truth_bfdlabel_$m.tsv; [ $m = indep ] && B=-
    python3 $REPO/scripts/review_disagreements.py $m.$v $d.aware.tsv sheets/truth_$m.tsv $B $REL/meta_full.tsv $GR/taxdump_20261029 \
      $AUDIT score/disagreements_draft.$m.$v.tsv $X --reviews $REVIEWS
  done
done
python3 $REPO/scripts/truth_set_tables.py score score/tables.md score/disagreements.tsv ${PREV:-/dev/null}
