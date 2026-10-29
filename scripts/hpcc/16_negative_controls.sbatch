#!/bin/bash -l
#SBATCH -J gmg_neg
#SBATCH -p short
#SBATCH -c 4
#SBATCH --mem=16G
#SBATCH -t 01:00:00
#SBATCH -o /bigdata/stajichlab/shared/projects/GmG_runs/logs/16_neg.%j.log
# Negative controls (#26): real non-fungal genomes must not get a fungal call.
#   bacteria (E. coli K-12, B. subtilis 168), virus (phiX174), plant (A. thaliana), green alga (C. reinhardtii),
#   nematode (C. elegans), apicomplexan (P. falciparum 3D7). Human is not included (3 GB); see reports/validation/negative_controls.md.
set -uo pipefail
source $HOME/projects/GmG/scripts/hpcc/env.sh
module load ncbi_datasets/18.30.1
R=$GR/negative; W=$R/work; mkdir -p $W/g; cd $W
printf "%s\n" GCF_000005845.2 GCF_000009045.1 GCF_000819615.1 GCF_000001735.4 GCF_000002595.2 GCF_000002985.6 GCF_000002765.5 > acc.txt
datasets download genome accession --inputfile acc.txt --include genome --filename d.zip --no-progressbar || exit 1
unzip -oq d.zip -d x
for f in x/ncbi_dataset/data/GC*/*.fna; do mv "$f" g/$(basename $(dirname "$f")).fna; done
ls -la g
INDEX=${INDEX:-$GR/index/bfd_s10000.v1.gmgi}
$GMG --threads ${SLURM_CPUS_PER_TASK} run g --index $INDEX -o $R/out --report-format all
cp $R/out/results.tsv $R/results.tsv
rm -rf $W
