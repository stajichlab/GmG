#!/bin/bash -l
#SBATCH -J gmg_stress
#SBATCH -p short
#SBATCH -c 4
#SBATCH --mem=16G
#SBATCH -t 02:00:00
#SBATCH -o /bigdata/stajichlab/shared/projects/GmG_runs/logs/13_stress.%j.log
# Stress tests (#41): download 24 known assemblies, degrade them (scripts/stress_sim.py), classify with `gmg run`.
set -uo pipefail
source $HOME/projects/GmG/scripts/hpcc/env.sh
module load ncbi_datasets/18.30.1
R=$GR/stress; W=$R/work; mkdir -p $W/g; cd $W
datasets download genome accession --inputfile $GMG_REPO/data/validation/stress_hosts.txt --include genome --filename d.zip --no-progressbar || exit 1
unzip -oq d.zip -d x
for f in x/ncbi_dataset/data/GC*/*.fna; do mv "$f" g/$(basename $(dirname "$f")).fna; done
python3 $GMG_REPO/scripts/stress_sim.py g $GMG_REPO/data/validation/stress_hosts.meta.tsv $W/sim
$GMG --threads ${SLURM_CPUS_PER_TASK} run $W/sim/genomes --index $GR/index/bfd_s10000.v1.gmgi -o $R/out --no-report
cp $W/sim/truth.tsv $R/truth.tsv
rm -rf $W
