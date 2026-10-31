#!/bin/bash -l
#SBATCH -J mixdet_reads
#SBATCH -p short
#SBATCH -c 8
#SBATCH --mem=32G
#SBATCH -t 02:00:00
#SBATCH -o /bigdata/stajichlab/shared/projects/GmG_runs/logs/35_mixdet_reads.%A_%a.log
# #31: the two-fungus read mixtures of #45 (plan from 33_mix_plan.sbatch; 32 pairs x 9 shares at 4x, 25 of them also at 20x, and 4 triples) classified with
# the mixture-detector binary, default read settings only. One array task takes every 4th pair of the plan (--array=0-3). Each variant is written as plain
# FASTQ (scripts/mix_reads.py), classified and deleted; output stays in GmG_runs/mixdet/out/reads/<pair>.
set -uo pipefail
source $HOME/projects/GmG/scripts/hpcc/env.sh
R=$GR/mixdet; MG=${MG:-$R/bin/gmg}; S=$GR/reads/mix/scripts
INDEX=${INDEX:-$GR/release/2026.10.4/bfd-2026.10.4-s10000.gmgi}
LOW=600000; HIGH=2700000
P=$GR/reads/mix/plan; W=$R/work/reads; OUT=$R/out/reads; mkdir -p $W $OUT
K=${SLURM_ARRAY_TASK_ID:-0}; NT=${NT:-4}
echo "binary $($MG --version); index $INDEX; task $K of $NT"
pair() {
  pid=$1; high=$2; a1=$3; a2=$4; b1=$5; b2=$6
  [ -s $OUT/$pid/results.jsonl ] && return 0
  d=$W/$pid; mkdir -p $d; : > $d/body.tsv
  for depth in low high; do
    [ $depth = high ] && { [ "$high" = 1 ] || continue; N=$HIGH; } || N=$LOW
    for pa in 100 99.9 99.5 99 95 90 80 50 0; do
      pb=$(python3 -c "print(round(100-$pa,4))")
      python3 $S/mix_reads.py $d ${pid}__${depth}__A$pa $N $pa,$pb $a1 $a2 $b1 $b2 >> $d/body.tsv || return 1
    done
  done
  (echo -e "id\tpath\texpected"; awk -F'\t' '{print $1"\t"$2"\t"}' $d/body.tsv) > $d/sheet.tsv
  $MG --threads 4 run --sample-sheet $d/sheet.tsv --index $INDEX -o $OUT/$pid --thresholds known --no-report -q --top 10 --max-reads 3000000
  rm -rf $d
}
triple() {
  tid=$1; shift
  [ -s $OUT/$tid/results.jsonl ] && return 0
  d=$W/$tid; mkdir -p $d
  python3 $S/mix_reads.py $d ${tid}__low__mix $LOW 60,30,10 "$@" > $d/body.tsv || return 1
  (echo -e "id\tpath\texpected"; awk -F'\t' '{print $1"\t"$2"\t"}' $d/body.tsv) > $d/sheet.tsv
  $MG --threads 4 run --sample-sheet $d/sheet.tsv --index $INDEX -o $OUT/$tid --thresholds known --no-report -q --top 10 --max-reads 3000000
  rm -rf $d
}
export -f pair triple; export S MG INDEX LOW HIGH OUT W
tail -n +2 $P/plan.tsv | awk -F'\t' -v k=$K -v n=$NT '(NR-1)%n==k {print $1, $3, $5, $6, $8, $9}' | xargs -P 2 -L 1 bash -c 'pair "$@"' _
tail -n +2 $P/triples.tsv | awk -F'\t' -v k=$K -v n=$NT '(NR-1)%n==k {print $1, $4, $5, $13, $14, $22, $23}' | xargs -P 2 -L 1 bash -c 'triple "$@"' _
echo done
