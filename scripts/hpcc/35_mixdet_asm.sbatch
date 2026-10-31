#!/bin/bash -l
#SBATCH -J mixdet_asm
#SBATCH -p short
#SBATCH -c 8
#SBATCH --mem=32G
#SBATCH -t 02:00:00
#SBATCH -o /bigdata/stajichlab/shared/projects/GmG_runs/logs/35_mixdet_asm.%j.log
# #31: two-fungus assembly mixtures (scripts/mix_asm.py: 36 pairs in three classes, B at 0.1..50% of the bases, plus the parents alone) classified with
# the mixture-detector binary; the parents are hidden from the reference (--exclude-ids). Per-input output stays in GmG_runs/mixdet/out/asm.
set -uo pipefail
source $HOME/projects/GmG/scripts/hpcc/env.sh
R=$GR/mixdet; MG=${MG:-$R/bin/gmg}; S=${S:-$R/src/scripts}
INDEX=${INDEX:-$GR/release/2026.10.4/bfd-2026.10.4-s10000.gmgi}
TAG=${TAG:-asm}; W=$R/work/$TAG; mkdir -p $W $R/out/$TAG; cd $W
echo "binary $($MG --version); index $INDEX"
python3 $S/mix_asm.py $GR/truthset/truth.all.tsv $GR/truthset/genomes $W --per-class ${PER_CLASS:-12} || exit 1
$MG --threads ${SLURM_CPUS_PER_TASK} run --sample-sheet $W/sheet.tsv --index $INDEX -o $R/out/$TAG --thresholds known --exclude-ids $W/exclude_ids.txt --no-report -q --top 10
cp $W/truth.tsv $R/out/$TAG/truth.tsv
rm -rf $W/genomes
echo done
