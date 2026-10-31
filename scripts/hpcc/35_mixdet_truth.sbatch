#!/bin/bash -l
#SBATCH -J mixdet_truth
#SBATCH -p short
#SBATCH -c 4
#SBATCH --mem=12G
#SBATCH -t 02:00:00
#SBATCH -o /bigdata/stajichlab/shared/projects/GmG_runs/logs/35_mixdet_truth.%A_%a.log
# #31: pure single-organism assemblies as negatives: the 674-genome RefSeq truth set with the query's own genome (loo), species (los) or genus (log) hidden.
# usage: MODE=loo|los|log sbatch --array=0-13%8 35_mixdet_truth.sbatch   (shards of $GR/truthset/sheets)
set -uo pipefail
source $HOME/projects/GmG/scripts/hpcc/env.sh
R=$GR/mixdet; MG=${MG:-$R/bin/gmg}
MODE=${MODE:?}
INDEX=${INDEX:-$GR/release/2026.10.4/bfd-2026.10.4-s10000.gmgi}
S=$(printf "%03d" ${SLURM_ARRAY_TASK_ID})
SHEET=$GR/truthset/sheets/$MODE.$S.tsv
[ -f $SHEET ] || { echo "no shard $S"; exit 0; }
OUT=$R/out/truth_$MODE/$S
[ -s $OUT/results.tsv ] && exit 0
case $MODE in loo) L="--exclude-level self";; los) L="--exclude-level species";; log) L="--exclude-level genus";; *) L="";; esac
$MG --threads ${SLURM_CPUS_PER_TASK} run --sample-sheet $SHEET --index $INDEX -o $OUT --thresholds known $L --no-report -q
echo "shard $S: $(($(wc -l < $OUT/results.tsv)-1)) rows"
