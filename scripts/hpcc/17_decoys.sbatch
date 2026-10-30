#!/bin/bash -l
#SBATCH -J gmg_decoys
#SBATCH -p short
#SBATCH -c 8
#SBATCH --mem=48G
#SBATCH -t 02:00:00
#SBATCH -o /bigdata/stajichlab/shared/projects/GmG_runs/logs/17_decoys.%j.log
# #79: download the non-fungal decoy genomes (data/decoys/decoys.tsv) and sketch them into the reference sketch store
# ($GR/sketch, scaled 1000, like the BFD genomes), and download the held-out negative genomes (data/decoys/heldout_negatives.txt)
# that are NOT decoys, for testing.
set -uo pipefail
source $HOME/projects/GmG/scripts/hpcc/env.sh
module load ncbi_datasets/18.30.1
R=$GR/decoys; mkdir -p $R/dl $R/held $R/work; cd $R/work
cut -f1 $GMG_REPO/data/decoys/decoys.tsv | tail -n +2 > decoy_acc.txt
grep -v '^#' $GMG_REPO/data/decoys/heldout_negatives.txt | cut -d' ' -f1 > held_acc.txt
for kind in decoy held; do
  dest=$R/dl; [ $kind = held ] && dest=$R/held
  rm -rf x_$kind d_$kind.zip
  for try in 1 2 3; do
    datasets download genome accession --inputfile ${kind}_acc.txt --include genome --filename d_$kind.zip --no-progressbar && break
    sleep 30
  done
  unzip -oq d_$kind.zip -d x_$kind
  for f in x_$kind/ncbi_dataset/data/GC*/*.fna; do mv "$f" $dest/$(basename $(dirname "$f")).fna; done
done
ls -la $R/dl $R/held
$GMG --threads ${SLURM_CPUS_PER_TASK} sketch -k 21 -s 1000 -o $GR/sketch $R/dl/*.fna
ls -la $GR/sketch/GCF_000001735.4.gsk $GR/sketch/GCA_000143045.1.gsk
cd $R && rm -rf work
