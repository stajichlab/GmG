#!/bin/bash
#SBATCH -p short -N1 -n1 -c24 --mem=48G -t 2:00:00 -J popsweep -o %x_%j.out
# Threshold sensitivity: re-run popsketch on cached sketches while varying one threshold at a time (no tree, no PCA).
# usage: sbatch --export=ALL,DS=cocci 35_threshold_sweep.sbatch   (DS = cocci | bd); summarise with 36_summarise_sweep.py
set -euo pipefail
P=/bigdata/stajichlab/jstajich/projects/GmG_popgen
BIN=$P/popsketch_build/target/release/popsketch
cd $P/popsketch_runs
mkdir -p sweep_$DS
run() {  # name, flags...
  n=$1; shift
  d=sweep_$DS/$n
  mkdir -p $d; [ -e $d/sketches ] || ln -s $P/popsketch_runs/${DS}_full/sketches $d/sketches
  $BIN run --sample-sheet ${DS}_samples.tsv --labels ${DS}_labels.tsv -o $d --threads 24 --resume --no-tree --no-pca --report-format md "$@" > $d/log.txt 2>&1 || echo "failed: $n"
}
run base
for v in 1 2 3 5 8;           do run minfloors_$v --min-floors $v; done
for v in 3 4 6 8 12;          do run sizez_$v --size-z $v; done
for v in 0.01 0.02 0.05 0.1;  do run sizefloor_$v --size-floor $v; done
for v in 3 4 8 12;            do run shoulderz_$v --shoulder-z $v; done
for v in 3 5 8 12;            do run mincov_$v --min-cov $v; done
