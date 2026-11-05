#!/bin/bash
#SBATCH -p short -N1 -n1 -c4 --mem=8G -t 2:00:00 -J cocci_sk -o logs/sk_%A_%a.out
# usage: sbatch --array=1-N%40 01_sketch.sbatch   (task i = row i+1 of manifest.tsv, found==1 only)
set -euo pipefail
W=${W:-/bigdata/stajichlab/jstajich/projects/GmG_popgen/popgen_cocci}
GMG=${GMG:-/bigdata/stajichlab/shared/projects/GmG_runs/bin/gmg}
MAXREADS=${MAXREADS:-6000000}   # R1 reads; ~25-30x of a 28 Mb genome at 150 bp
cd "$W"
line=$(awk -F'\t' 'NR>1 && $6==1' manifest.tsv | sed -n "${SLURM_ARRAY_TASK_ID}p")
id=$(cut -f1 <<<"$line"); r1=$(cut -f4 <<<"$line")
# scaled=200 sketch (can be downsampled to 1000 later); abundance tracked, no filtering yet
mkdir -p sk200 sk1000
tmp=$(mktemp -d); ln -s "$r1" "$tmp/$id.fastq.gz"
$GMG sketch -q --threads 4 --abund -k 21 -s 200  --max-reads "$MAXREADS" -o sk200  "$tmp/$id.fastq.gz"
$GMG sketch -q --threads 4 --abund -k 21 -s 1000 --max-reads "$MAXREADS" -o sk1000 "$tmp/$id.fastq.gz"
rm -rf "$tmp"
