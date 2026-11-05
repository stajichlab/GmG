#!/bin/bash
#SBATCH -p short -N1 -n1 -c24 --mem=96G -t 2:00:00 -J popbig -o %x_%j.out
# A larger set than the earlier validations: every 4th clean genome of the BFD set (about 2,500), no labels, no PCA.
set -euo pipefail
P=/bigdata/stajichlab/jstajich/projects/GmG_popgen
D=/bigdata/stajichlab/shared/projects/BFD/Fungi_BFD_runs/input_clean_genomes
cd $P/popsketch_runs
set +o pipefail  # head closes the pipe early
printf 'id\tpath\n' > asm_big_samples.tsv
ls $D | grep -E '\.fa\.gz$' | awk 'NR%4==0' | head -2500 | while read f; do printf '%s\t%s/%s\n' "${f%.fa.gz}" "$D" "$f"; done >> asm_big_samples.tsv
set -o pipefail
wc -l asm_big_samples.tsv
/usr/bin/time -v $P/popsketch_build/target/release/popsketch run --sample-sheet asm_big_samples.tsv -o asm_big_full --threads 24 --max-samples 3000 --no-pca --report-format md 2> asm_big_time.txt
grep -E "Elapsed|Maximum resident" asm_big_time.txt
