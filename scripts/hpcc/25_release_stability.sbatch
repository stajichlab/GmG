#!/bin/bash -l
#SBATCH -J gmg_relstab
#SBATCH -p short
#SBATCH -c 4
#SBATCH --mem=16G
#SBATCH -t 02:00:00
#SBATCH -o /bigdata/stajichlab/shared/projects/GmG_runs/logs/25_relstab.%A_%a.log
# #48: the 1,257 post-freeze assemblies (external_big batches) against every release index with the SAME binary: one download per
# batch, one `gmg run` per release. usage: sbatch --array=0-50%6 25_release_stability.sbatch   (VERS="2026.10 2026.10.1 ..." overrides)
set -uo pipefail
source $HOME/projects/GmG/scripts/hpcc/env.sh
module load ncbi_datasets/18.30.1
VERS=${VERS:-"2026.10 2026.10.1 2026.10.2 2026.10.3"}
B=$(printf "%03d" ${SLURM_ARRAY_TASK_ID})
R=$GR/external_big; W=$R/work/rs$B
todo=""
for v in $VERS; do [ -s $R/rel_$v/b$B/results.tsv ] || todo="$todo $v"; done
[ -n "$todo" ] || { echo "batch $B already done"; exit 0; }
rm -rf $W; mkdir -p $W/genomes
ok=0
for try in 1 2 3 4; do
  if datasets download genome accession --inputfile $R/batches/batch_$B --include genome --filename $W/d.zip --no-progressbar 2> $W/dl.err; then ok=1; break; fi
  sleep $((try*30))
done
[ $ok = 1 ] || { echo "download failed for batch $B"; cat $W/dl.err; exit 1; }
unzip -oq $W/d.zip -d $W/x
for f in $W/x/ncbi_dataset/data/GC*/*.fna; do mv "$f" $W/genomes/$(basename $(dirname "$f")).fna; done
for v in $todo; do
  $GMG --threads ${SLURM_CPUS_PER_TASK} -q run $W/genomes --index $GR/release/$v/bfd-$v-s10000.gmgi -o $R/rel_$v/b$B --no-report
done
rm -rf $W
