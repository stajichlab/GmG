#!/bin/bash -l
#SBATCH -p short -N1 -c16 --mem=64G -t 4:00:00 -J f1genomes -o %x_%j.out
# F1 step 2: sketch the fusarioid BFD genomes, all-pairs sketch ANI, and gmg leave-out benchmarks over index resolution.
# usage: REPO=<checkout> W=<build dir> sbatch run_genomes.sbatch
# Writes $W/genomes (lists), $W/sketch (k21 s1000 .gsk, purgeable), $W/ani_s1000.tsv, $W/bench/*.tsv
set -euo pipefail
REPO=${REPO:?set REPO}; W=${W:?set W}
source "$REPO/scripts/hpcc/env.sh"   # GR, GMG, BFD_GENOMES
module load miniconda3
F1="$REPO/bench/F1_fusariumid"
T=$SLURM_CPUS_PER_TASK
mkdir -p "$W/sketch" "$W/bench"
[ -s "$W/genomes/paths.txt" ] || bash "$F1/make_lists.sh"
$GMG --threads $T -q sketch -k 21 -s 1000 --list "$W/genomes/paths.txt" -o "$W/sketch"
ls "$W/sketch" | wc -l
python3 -I "$F1/genome_ani.py" --scaled 1000 --threads $T -o "$W/ani_s1000.tsv" "$W"/sketch/*.gsk
for s in 1000 2000 10000 20000 50000; do
  for m in loo los; do
    $GMG --threads $T -q bench --index "$GR/index/bfd_s$s.gmgi" --sketches "$W/sketch" --ref-list "$GR/index/refs_2026.10.5.txt" \
      --mode $m -o "$W/bench/all_${m}_s$s.tsv"
    $GMG --threads $T -q bench --index "$GR/index/bfd_s$s.gmgi" --sketches "$W/sketch" --ref-list "$W/genomes/refs_fus.txt" \
      --mode $m -o "$W/bench/fus_${m}_s$s.tsv"
  done
done
ls -la "$W/bench"
