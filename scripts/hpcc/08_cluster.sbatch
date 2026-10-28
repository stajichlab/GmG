#!/bin/bash -l
#SBATCH -J gmg_cluster
#SBATCH -p short
#SBATCH -c 2
#SBATCH --mem=16G
#SBATCH -t 01:00:00
#SBATCH -o /bigdata/stajichlab/shared/projects/GmG_runs/logs/08_cluster.%j.log
# greedy ANI clustering of all BFD genomes (centroid priority = BFD rule: BUSCO, N50), at 95% (species-like) and 99%
set -uo pipefail
source /bigdata/stajichlab/shared/projects/GmG_runs/src/GmG/scripts/hpcc/env.sh
cd $GR/index
for A in 0.95 0.99; do
  $GMG cluster --index bfd_s10000.gmgi --sketches $GR/sketch --ani $A --order reps.order.txt -o clusters_ani${A}.tsv
done
ls -la clusters_*.tsv
