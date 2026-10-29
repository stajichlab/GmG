#!/bin/bash -l
#SBATCH -J gmg_extbig
#SBATCH -p short
#SBATCH -c 4
#SBATCH --mem=12G
#SBATCH -t 02:00:00
#SBATCH -o /bigdata/stajichlab/shared/projects/GmG_runs/logs/10_extbig.%A_%a.log
# Independent assemblies: download a batch of NCBI assemblies, classify them with `gmg run`, delete the FASTA.
# usage: sbatch --array=0-50%8 10_fetch_and_run.sbatch   (batches are $GR/external_big/batches/batch_NNN, 25 accessions each)
set -uo pipefail
source $HOME/projects/GmG/scripts/hpcc/env.sh
module load ncbi_datasets/18.30.1
B=$(printf "%03d" ${SLURM_ARRAY_TASK_ID})
R=$GR/external_big; W=$R/work/b$B; OUT=$R/${OUTNAME:-out}/b$B
INDEX=${INDEX:-$GR/index/bfd_s10000.v1.gmgi}   # OUTNAME/INDEX: rerun the same batches against another index (e.g. a decoy index)
[ -s $OUT/results.tsv ] && { echo "batch $B already done"; exit 0; }
rm -rf $W; mkdir -p $W/genomes $R/${OUTNAME:-out}
ok=0
for try in 1 2 3 4; do
  if datasets download genome accession --inputfile $R/batches/batch_$B --include genome --filename $W/d.zip --no-progressbar 2> $W/dl.err; then ok=1; break; fi
  sleep $((try*30))
done
[ $ok = 1 ] || { echo "download failed for batch $B"; cat $W/dl.err; exit 1; }
unzip -oq $W/d.zip -d $W/x
for f in $W/x/ncbi_dataset/data/GC*/*.fna; do mv "$f" $W/genomes/$(basename $(dirname "$f")).fna; done
echo "batch $B: $(ls $W/genomes | wc -l) genomes of $(wc -l < $R/batches/batch_$B) requested"
$GMG --threads ${SLURM_CPUS_PER_TASK} run $W/genomes --index $INDEX -o $OUT --no-report
rm -rf $W
