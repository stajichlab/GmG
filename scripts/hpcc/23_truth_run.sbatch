#!/bin/bash -l
#SBATCH -J gmg_truth_run
#SBATCH -p short
#SBATCH -c 4
#SBATCH --mem=12G
#SBATCH -t 02:00:00
#SBATCH -o /bigdata/stajichlab/shared/projects/GmG_runs/logs/23_truth_run.%A_%a.log
# #39: classify one shard of the truth set with the query left out of the reference.
# usage: MODE=loo|los|log|indep TH=known|novel sbatch --array=0-13%8 23_truth_run.sbatch
#   loo/los/log: --exclude-level self/species/genus (sheet column `ref` = the index genome that is the query); indep: no exclusion
set -uo pipefail
source $HOME/projects/GmG/scripts/hpcc/env.sh
MODE=${MODE:?}; TH=${TH:-known}
INDEX=${INDEX:-$GR/release/2026.10.2/bfd-2026.10.2-s10000.gmgi}
T=$GR/truthset
S=$(printf "%03d" ${SLURM_ARRAY_TASK_ID})
SHEET=$T/sheets/$MODE.$S.tsv
[ -f $SHEET ] || { echo "no shard $S"; exit 0; }
OUT=$T/runs/$MODE.$TH/$S
[ -s $OUT/results.tsv ] && { echo "shard $S already done"; exit 0; }
case $MODE in loo) L="--exclude-level self";; los) L="--exclude-level species";; log) L="--exclude-level genus";; *) L="";; esac
$GMG --threads ${SLURM_CPUS_PER_TASK} run --sample-sheet $SHEET --index $INDEX -o $OUT --thresholds $TH $L --no-report -q
echo "shard $S: $(($(wc -l < $OUT/results.tsv)-1)) rows"
