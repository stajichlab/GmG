#!/bin/bash -l
#SBATCH -J gmg_rnaex
#SBATCH -p short
#SBATCH -c 16
#SBATCH --mem=32G
#SBATCH -t 02:00:00
#SBATCH -o /bigdata/stajichlab/shared/projects/GmG_runs/logs/27_rna_excl.%j.log
# RNA-seq reads (#43), extending reports/validation/rna_reads_*.md: the same 345 species-labelled BFD RNA-seq sets, now also with the
# species / the genus of the sample hidden from the reference (`--exclude-level`, #85; the species comes from the `expected` column),
# so the species call is not a lookup of the same species. First 100k and 500k reads, default thresholds, newest index.
set -uo pipefail
source $HOME/projects/GmG/scripts/hpcc/env.sh
S=${READS_SCRIPTS:-$GMG_REPO/scripts}
REL=$GR/release/2026.10.3; [ -s $REL/bfd-2026.10.3-s10000.gmgi ] || REL=$GR/release/2026.10.2
INDEX=${INDEX:-$(ls $REL/bfd-20*-s10000.gmgi)}
V=$GR/reads/rna_validation; R=$GR/reads/rna_exclusion; mkdir -p $R; cd $R
echo "index $INDEX"
# sheet: id, path (the read file in the farm), expected (the species of the sample)
awk -F'\t' -v farm=$V/farm 'BEGIN{OFS="\t"; print "id","path","expected"} NR>1{print $1, farm"/"$8, $7}' $V/rna.truth.tsv > sheet.tsv
cut -f1-7 $V/rna.truth.tsv > truth.tsv
for N in 100000 500000; do
  for L in none species genus; do
    D=n${N}_$L
    if [ $L = none ]; then X=""; else X="--exclude-level $L"; fi
    [ -s $D/results.tsv ] || $GMG --threads ${SLURM_CPUS_PER_TASK} run --sample-sheet sheet.tsv --index $INDEX -o $D --no-report --max-reads $N $X
    $GMG validate $D/results.tsv truth.tsv -o v_$D --by phylum --title "RNA-seq reads, first $N reads, exclusion: $L"
  done
done
# inputs whose species/genus name matched nothing in the index (nothing was hidden: a caveat for the exclusion modes)
for D in n500000_species n500000_genus; do echo "$D: $(grep -c '"n_genomes":0' $D/results.jsonl) of $(wc -l < $D/results.jsonl) inputs hid no genome"; done > hidden_check.txt
cat hidden_check.txt
python3 $S/reads_summary.py all $(for N in 100000 500000; do for L in none species genus; do echo -n "n${N}_$L=v_n${N}_$L/validation.tsv "; done; done) > summary_all.md
cat summary_all.md
