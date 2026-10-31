#!/bin/bash -l
#SBATCH -J gmg_score_i5
#SBATCH -p short
#SBATCH -c 4
#SBATCH --mem=32G
#SBATCH -t 01:00:00
#SBATCH -o /bigdata/stajichlab/shared/projects/GmG_runs/logs/score_i5.%j.log
# 2026.10.5 re-measurement, scoring of the leave-out bench (novelty curve, held-out calibration), the 1,257 post-freeze assemblies and the
# Coccidioides / RNA-seq read sets with the same scripts as the earlier reports. needs pandas.
set -uo pipefail
source $HOME/projects/GmG/scripts/hpcc/env.sh
REPO=$GR/precision/repo; O=$GR/precision/score; mkdir -p $O; cd $O
EX=$REPO/reports/label_audit/reference_exclusions.txt
python3 $REPO/scripts/analyze_bench.py $O/bench_i5 $GR/bench/i5/bench_*_s10000.tsv > /dev/null  # GMG_EXCLUDE below
GMG_EXCLUDE=$EX GMG_META=$GR/release/2026.10.5/meta_full.tsv python3 $REPO/scripts/analyze_bench.py $O/bench_i5 $GR/bench/i5/bench_*_s10000.tsv > /dev/null
python3 $REPO/scripts/novelty_curve.py $GR/bench/i5 $EX $O/novelty_curve
python3 $REPO/scripts/heldout_calibration.py $GR/bench/i5 $GR/release/2026.10.5/meta_full.tsv $EX $O/heldout_calibration 20 0.95
R=$GR/external_big/rel_2026.10.5; out=$O/post.results.tsv; first=$(ls $R/b*/results.tsv | head -1); head -1 $first > $out
for f in $R/b*/results.tsv; do tail -n +2 $f >> $out; done
echo "post-freeze: $(($(wc -l < $out)-1)) rows"
python3 $REPO/scripts/validate_run.py $out $REPO/data/validation/post_freeze.meta.tsv $O/post_freeze_balanced --by species_in_bfd genus_in_bfd level phylum --title "Post-freeze, 2026.10.5" > /dev/null
for N in 100000 500000 2000000; do
  python3 $REPO/scripts/validate_run.py $GR/reads/cocci_validation/n$N.i5/results.tsv $REPO/data/validation/cocci_reads.truth.tsv $O/cocci_reads_n$N --by --title "Cocci reads $N, 2026.10.5" > /dev/null
  python3 $REPO/scripts/validate_run.py $GR/reads/rna_validation/n$N.i5/results.tsv $GR/reads/rna_validation/rna.truth.tsv $O/rna_reads_n$N --by phylum --title "RNA reads $N, 2026.10.5" > /dev/null
done
ls $O
