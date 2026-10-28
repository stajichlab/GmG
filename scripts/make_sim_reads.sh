#!/bin/bash
# usage: make_sim_reads.sh GENOME_DIR ID_LIST OUT_DIR   (run in the popsketch conda env; needs target/release/gmg)
# For each id: simulate 150 bp reads at 0.5x, 1x, 5x, 20x, sketch with abundances (k=21 scaled=1000), delete the FASTQ.
set -eu
gd=$1; ids=$2; out=$3; mkdir -p $out/fq $out/sim
G=$(dirname $0)/../target/release/gmg
while read id; do
  for cov in 0.5 1 5 20; do
    [ -s $out/sim/${id}_${cov}x.gsk ] && continue
    python3 $(dirname $0)/sim_reads.py $gd/$id.fna $out/fq/${id}_${cov}x.fastq.gz $cov --seed 7 2>/dev/null
    $G sketch --abund -k 21 -s 1000 -o $out/sim $out/fq/${id}_${cov}x.fastq.gz 2>&1 | grep -E "FAIL|EMPTY" || true
    rm -f $out/fq/${id}_${cov}x.fastq.gz
  done
done < $ids
echo done > $out/sim.done
