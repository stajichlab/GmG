#!/bin/bash
#SBATCH -p short -N1 -n1 -c24 --mem=48G -t 2:00:00 -J popsk -o %x_%j.out
# usage: sbatch --export=ALL,DS=cocci 28_popsketch_run.sbatch   (DS = cocci | bd)
set -euo pipefail
P=/bigdata/stajichlab/jstajich/projects/GmG_popgen
BIN=$P/popsketch_build/target/release/popsketch
IDX=/bigdata/stajichlab/shared/projects/GmG_runs/release/2026.10.5/bfd-2026.10.5-s10000.gmgi
cd $P/popsketch_runs
$BIN run --sample-sheet ${DS}_samples.tsv --labels ${DS}_labels.tsv -o ${DS}_full --index $IDX --threads 24 ${EXTRA:-}
