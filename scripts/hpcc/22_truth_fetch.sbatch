#!/bin/bash -l
#SBATCH -J gmg_truth_dl
#SBATCH -p short
#SBATCH -c 2
#SBATCH --mem=8G
#SBATCH -t 02:00:00
#SBATCH -o /bigdata/stajichlab/shared/projects/GmG_runs/logs/22_truth_dl.%A_%a.log
# #39: download one batch (25 accessions, $GR/truthset/batches/batch_NNN) of the RefSeq truth-set genomes as $GR/truthset/genomes/ACC.fna.gz
# usage: sbatch --array=0-26%6 22_truth_fetch.sbatch
set -uo pipefail
source $HOME/projects/GmG/scripts/hpcc/env.sh
module load ncbi_datasets/18.30.1
B=$(printf "%03d" ${SLURM_ARRAY_TASK_ID})
T=$GR/truthset; W=$T/work/b$B; mkdir -p $T/genomes
[ -f $T/batches/batch_$B ] || { echo "no batch $B"; exit 0; }
rm -rf $W; mkdir -p $W
ok=0
for try in 1 2 3 4; do
  if datasets download genome accession --inputfile $T/batches/batch_$B --include genome --filename $W/d.zip --no-progressbar 2> $W/dl.err; then ok=1; break; fi
  sleep $((try*30))
done
[ $ok = 1 ] || { echo "download failed for batch $B"; cat $W/dl.err; exit 1; }
unzip -oq $W/d.zip -d $W/x
for f in $W/x/ncbi_dataset/data/GC*/*.fna; do
  a=$(basename $(dirname "$f")); gzip -c "$f" > $T/genomes/$a.fna.gz.part && mv $T/genomes/$a.fna.gz.part $T/genomes/$a.fna.gz
done
echo "batch $B: $(ls $T/genomes | wc -l) genomes so far; $(wc -l < $T/batches/batch_$B) requested"
rm -rf $W
