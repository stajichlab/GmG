#!/bin/bash -l
#SBATCH -J gmg_index
#SBATCH -p short
#SBATCH -c 4
#SBATCH --mem=64G
#SBATCH -t 02:00:00
#SBATCH -o /bigdata/stajichlab/shared/projects/GmG_runs/logs/03_index.%j.log
# builds indexes at several scaled values from the scaled=1000 sketches
set -euo pipefail
source $HOME/projects/GmG/scripts/hpcc/env.sh
echo "sketches: $(ls $GR/sketch | grep -c '\.gsk$')"
for S in 10000 2000 1000; do
  /usr/bin/time -v $GMG --threads 4 build --meta $GR/index/meta.tsv --sketches $GR/sketch --scaled $S -o $GR/index/bfd_s${S}.gmgi \
    2> $GR/logs/03_index.s${S}.err || true
  tail -4 $GR/logs/03_index.s${S}.err | head -3
  $GMG info $GR/index/bfd_s${S}.gmgi
done
