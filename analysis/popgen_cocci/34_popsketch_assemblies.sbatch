#!/bin/bash
#SBATCH -p short -N1 -n1 -c24 --mem=64G -t 2:00:00 -J popasm -o %x_%j.out
# usage: sbatch --export=ALL,DS=asm_cocci 34_popsketch_assemblies.sbatch   (DS = asm_cocci | asm_asp)
set -euo pipefail
P=/bigdata/stajichlab/jstajich/projects/GmG_popgen
cd $P/popsketch_runs
/usr/bin/time -v $P/popsketch_build/target/release/popsketch run --sample-sheet ${DS}_samples.tsv --labels ${DS}_labels.tsv -o ${DS}_full --threads 24 --max-samples 2000 ${EXTRA:-} 2> ${DS}_time.txt
tail -n 30 ${DS}_time.txt | grep -E "Elapsed|Maximum resident" || true
