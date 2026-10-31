#!/bin/bash -l
#SBATCH -J gmg_parity
#SBATCH -p short
#SBATCH -c 8
#SBATCH --mem=16G
#SBATCH -t 01:00:00
#SBATCH -o /bigdata/stajichlab/shared/projects/GmG_runs/logs/28_parity.%j.log
# #49 on Linux: (1) hash sets of `gmg sketch` against Python sourmash 4.9 on the golden panel; (2) the golden-panel regression check
# with the deployed linux binary (compare with the macOS result: same calls, flags and ANI).
# usage: sbatch 28_parity_linux.sbatch [SCRIPT_DIR] [DATA_DIR]   (SCRIPT_DIR: scripts/ of a checkout; DATA_DIR: its data/golden)
set -uo pipefail
source $HOME/projects/GmG/scripts/hpcc/env.sh
SD=${1:-$GMG_REPO/scripts}
R=$GR/parity; mkdir -p $R/sk $R/sm
python3 $SD/golden_check.py --fetch-only --cache $R/cache || exit 1
$GMG --threads ${SLURM_CPUS_PER_TASK} -q sketch -k 21 -s 10000 -o $R/sk $R/cache/*.fna.gz
module load sourmash/4.9
for f in $R/cache/*.fna.gz; do
  i=$(basename $f .fna.gz)
  [ -s $R/sm/$i.sig ] || sourmash sketch dna -p k=21,scaled=10000 --name $i -o $R/sm/$i.sig $f 2> /dev/null
done
sourmash --version
python3 $SD/sourmash_parity.py $R/sk $R/sm $R/sourmash_parity.tsv
python3 $SD/golden_check.py --gmg $GMG --index $GR/release/2026.10.5/bfd-reps-2026.10.5-s20000.gmgi --cache $R/cache --threads ${SLURM_CPUS_PER_TASK}
echo "golden check exit: $?"
