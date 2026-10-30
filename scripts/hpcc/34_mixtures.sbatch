#!/bin/bash -l
#SBATCH -J gmg_mix
#SBATCH -p short
#SBATCH -c 16
#SBATCH --mem=64G
#SBATCH -t 02:00:00
#SBATCH -o /bigdata/stajichlab/shared/projects/GmG_runs/logs/34_mixtures.%j.log
# #45 (reads half): two-fungus read mixtures from real read sets (plan from 33_mix_plan.sbatch) at six shares of parent A, pure controls, a second
# fungus at 0.5% and 0.1%, a high-depth subset, and 60:30:10 triples. Each variant is written as plain FASTQ (scripts/mix_reads.py), classified with the
# default read settings (min abundance 2) and with --min-abund 1, then deleted. Scored by scripts/mix_score.py.
set -uo pipefail
source $HOME/projects/GmG/scripts/hpcc/env.sh
S=${MIX_SCRIPTS:-$GR/reads/mix/scripts}
MG=$GR/reads/mix/bin/gmg; [ -x $MG ] || MG=$GMG
INDEX=${INDEX:-$GR/release/2026.10.4/bfd-2026.10.4-s10000.gmgi}
LOW=${LOW:-600000}; HIGH=${HIGH:-2700000}
R=$GR/reads/mix; P=$R/plan; W=$R/work; mkdir -p $W $R/out_default $R/out_abund1
echo "binary $($MG --version); index $INDEX; depths low=$LOW high=$HIGH pairs"
# a pair: write its variants, classify them in both abundance modes, delete them
pair() {
  pid=$1; high=$2; a1=$3; a2=$4; b1=$5; b2=$6
  d=$W/$pid; mkdir -p $d; : > $d/body.tsv
  for depth in low high; do
    [ $depth = high ] && { [ "$high" = 1 ] || continue; N=$HIGH; } || N=$LOW
    for pa in 100 99.9 99.5 99 95 90 80 50 0; do
      pb=$(python3 -c "print(round(100-$pa,4))")
      python3 $S/mix_reads.py $d ${pid}__${depth}__A$pa $N $pa,$pb $a1 $a2 $b1 $b2 >> $d/body.tsv || return 1
    done
  done
  (echo -e "id\tpath\texpected"; awk -F'\t' '{print $1"\t"$2"\t"}' $d/body.tsv) > $d/sheet.tsv
  $MG --threads 4 run --sample-sheet $d/sheet.tsv --index $INDEX -o $R/out_default/$pid --no-report -q --top 10 --max-reads 3000000
  $MG --threads 4 run --sample-sheet $d/sheet.tsv --index $INDEX -o $R/out_abund1/$pid --no-report -q --top 10 --max-reads 3000000 --min-abund 1
  rm -rf $d
}
triple() {
  tid=$1; shift
  d=$W/$tid; mkdir -p $d
  python3 $S/mix_reads.py $d ${tid}__low__mix $LOW 60,30,10 "$@" > $d/body.tsv || return 1
  (echo -e "id\tpath\texpected"; awk -F'\t' '{print $1"\t"$2"\t"}' $d/body.tsv) > $d/sheet.tsv
  $MG --threads 4 run --sample-sheet $d/sheet.tsv --index $INDEX -o $R/out_default/$tid --no-report -q --top 10 --max-reads 3000000
  $MG --threads 4 run --sample-sheet $d/sheet.tsv --index $INDEX -o $R/out_abund1/$tid --no-report -q --top 10 --max-reads 3000000 --min-abund 1
  rm -rf $d
}
export -f pair triple; export S MG INDEX LOW HIGH R W
# plan rows: pair_id distance high a_id a_r1 a_r2 b_id b_r1 b_r2 ...
tail -n +2 $P/plan.tsv | awk -F'\t' '{print $1, $3, $5, $6, $8, $9}' | xargs -P 4 -L 1 bash -c 'pair "$@"' _
tail -n +2 $P/triples.tsv | awk -F'\t' '{print $1, $4, $5, $13, $14, $22, $23}' | xargs -P 4 -L 1 bash -c 'triple "$@"' _
for M in default abund1; do
  for f in results hits; do
    awk 'FNR==1 && NR!=1 {next} {print}' $R/out_$M/*/$f.tsv > $R/${f}_$M.tsv
  done
done
wc -l $R/results_default.tsv $R/results_abund1.tsv
python3 $S/mix_score.py $P/plan.tsv $P/triples.tsv --run default=$R/results_default.tsv:$R/hits_default.tsv --run abund1=$R/results_abund1.tsv:$R/hits_abund1.tsv > $R/mix_tables.md
echo done
