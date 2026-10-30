#!/bin/bash -l
#SBATCH -J rc_sim
#SBATCH -p short
#SBATCH -c 8
#SBATCH --mem=32G
#SBATCH -t 02:00:00
#SBATCH --array=0-2
#SBATCH -o /bigdata/stajichlab/shared/projects/GmG_runs/logs/33_covsim.%A_%a.log
# #32: accuracy of the coverage estimate and of the auto cutoff on simulated reads of KNOWN coverage at the production scaled (10000):
# 3 real assemblies (ZyGoLife aaftf/sort, small to mid size), coverage 0.1x ... 30x, substitution error rates 0.1 / 0.5 / 1%, 150 bp single-end
# (scripts/sim_reads.py). One array task = one genome. Output: $R/sim/<genome>/results_a{2,auto}.tsv and truth.tsv (id, coverage, error).
set -uo pipefail
source $HOME/projects/GmG/scripts/hpcc/env.sh
GMGBIN=${GMGBIN:-$GR/bin/gmg-readscal}
REPO=${REPO:-$GR/reads_cal_repo}
INDEX=${INDEX:-$GR/release/2026.10.4/bfd-2026.10.4-s10000.gmgi}
SRC=/bigdata/stajichlab/shared/projects/ZyGoLife/LCG/Annotation_2026/results/aaftf/sort
R=$GR/reads/calibration/sim; mkdir -p $R
# three assemblies of different size, taken at fixed ranks of the size-sorted listing (the large ones are skipped: 30x of 100 Mb is 3 Gb of reads)
f=$(ls -S $SRC | sed -n "$((200 + 150 * SLURM_ARRAY_TASK_ID))p")
g=$(basename $f .sorted.fasta.gz); W=$R/$g; mkdir -p $W/reads; cd $W
echo "genome $g ($f)"
[ -s genome.fa ] || zcat $SRC/$f > genome.fa
grep -v '>' genome.fa | tr -d '\n' | wc -c > size.txt
echo -e "id\tgenome\tcoverage\terror" > truth.tsv
for cov in 0.1 0.3 1 3 10 30; do
  for err in 0.001 0.005 0.01; do
    id=${g}_${cov}x_e${err}
    [ -s reads/$id.fastq.gz ] || python3 $REPO/scripts/sim_reads.py genome.fa reads/$id.fastq.gz $cov --err $err --seed 7
    echo -e "$id\t$g\t$cov\t$err" >> truth.tsv
  done
done
for ma in 1 2 auto; do
  $GMGBIN --threads ${SLURM_CPUS_PER_TASK} run reads --index $INDEX -o out_a$ma --no-report -q --kind reads --min-abund $ma
  cp out_a$ma/results.tsv results_a$ma.tsv
done
rm -f genome.fa
