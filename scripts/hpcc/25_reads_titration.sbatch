#!/bin/bash -l
#SBATCH -J gmg_titr
#SBATCH -p short
#SBATCH -c 16
#SBATCH --mem=32G
#SBATCH -t 02:00:00
#SBATCH -o /bigdata/stajichlab/shared/projects/GmG_runs/logs/25_titration.%j.log
# DNA reads titration (#42): real paired Illumina reads of ZyGoLife strains (lab label = file name), first N read pairs per strain,
# N = 3k ... 1.5M per file; scored at phylum..genus against the NCBI lineage of the label's genus. Read-only on the input folder.
set -uo pipefail
source $HOME/projects/GmG/scripts/hpcc/env.sh
S=${READS_SCRIPTS:-$GMG_REPO/scripts}   # scripts staged under GmG_runs/reads/scripts when the checkout is not on this branch
IN=/bigdata/stajichlab/shared/projects/ZyGoLife/LCG/Assembly/input
REL=$GR/release/2026.10.3; [ -s $REL/bfd-2026.10.3-s10000.gmgi ] || REL=$GR/release/2026.10.2
INDEX=${INDEX:-$(ls $REL/bfd-20*-s10000.gmgi)}
R=$GR/reads/titration; mkdir -p $R; cd $R
echo "index $INDEX"
[ -s sheet.tsv ] || python3 $S/make_reads_sheet.py $IN . ${PER_GENUS:-2} 5
python3 $S/name_truth.py labels.tsv $GR/taxdump_20261029 truth.tsv --index-meta $REL/meta_full.tsv
for N in 3000 15000 75000 300000 1500000; do
  [ -s n$N/results.tsv ] || $GMG --threads ${SLURM_CPUS_PER_TASK} run --sample-sheet sheet.tsv --index $INDEX -o n$N --no-report --max-reads $N
  $GMG validate n$N/results.tsv truth.tsv -o v$N --by genus_in_ref --title "ZyGoLife reads, first $N read pairs"
done
# R1 only against the merged pair at 75k reads: does the second read of a pair change the call?
awk -F'\t' 'NR==1 || $2 ~ /_R1\./' sheet.tsv > sheet_r1.tsv
[ -s n75000_r1/results.tsv ] || $GMG --threads ${SLURM_CPUS_PER_TASK} run --sample-sheet sheet_r1.tsv --index $INDEX -o n75000_r1 --no-report --max-reads 75000
$GMG validate n75000_r1/results.tsv truth.tsv -o v75000_r1 --by genus_in_ref --title "R1 only, 75000 reads"
python3 $S/reads_summary.py all $(for N in 3000 15000 75000 300000 1500000; do echo -n "n$N=v$N/validation.tsv "; done) n75000_r1=v75000_r1/validation.tsv > summary_all.md
for G in True False; do
  python3 $S/reads_summary.py genus_in_ref=$G $(for N in 3000 15000 75000 300000 1500000; do echo -n "n$N=v$N/validation.tsv "; done) > summary_genus_in_ref_$G.md
done
cat summary_all.md
