#!/bin/bash
#SBATCH -p stajichlab -N1 -n1 -c8 --mem=16G -t 3:00:00 -J classchg -o logs/%x_%A_%a.out
# #159 step 1: are the reads of a sample whose genome-size estimate changed (high-coverage fix, PR #158) from the species the sample is labelled as?
# One array task per row of classify/tasks.tsv: run, accession of the species' reference assembly (PopSeq catalog/genome_sizes.csv), label.
# Fetches the first N read pairs (default 4M), then reports (a) the fraction of reads that map to the species' reference assembly (minimap2 -x sr)
# and (b) the top gmg guess. A low mapped fraction supports "most reads are from another genome" (the explanation for the new 3-4 Mb estimates).
# usage: sbatch --array=1-4 10_classify_changed.sbatch      (results in classify/<run>/summary.tsv; temp files in $SCRATCH)
# CLASSDIR=<dir> (default classify) holds tasks.tsv and the results, so the same runs can be mapped to a second reference
# (hybrid check: map to both parents' references, see RESULTS_159_hybrids.md).
set -uo pipefail
module load sratoolkit/3.2.0 minimap2/2.30 samtools/1.22.1 2>/dev/null || true
P=/bigdata/stajichlab/jstajich/projects/GmG_popgen
GMG=/bigdata/stajichlab/shared/projects/GmG_runs/bin/gmg
INDEX=/bigdata/stajichlab/shared/projects/GmG_runs/index/bfd_s1000.gmgi
N=${N:-4000000}
cd $P/${CLASSDIR:-classify}
read -r run acc label < <(awk -F'\t' -v n=$SLURM_ARRAY_TASK_ID 'NR==n' tasks.tsv)
out=$run; mkdir -p $out
[[ -s $out/summary.tsv ]] && { echo "skip $run"; exit 0; }
work=$(mktemp -d "${SCRATCH:-/tmp}/classchg.$run.XXXXXX"); trap 'rm -rf "$work"' EXIT
curl -sSf -o $work/ref.zip "https://api.ncbi.nlm.nih.gov/datasets/v2/genome/accession/$acc/download?include_annotation_type=GENOME_FASTA" || { echo "FAIL ref $acc"; exit 1; }
unzip -q -o $work/ref.zip -d $work/ref && cat $work/ref/ncbi_dataset/data/$acc/*.fna > $work/ref.fa
fastq-dump -X $N --split-files --gzip -O $work/raw $run >$work/dump.log 2>&1 || { echo "FAIL dump $run"; tail -3 $work/dump.log; exit 1; }
r1=$work/raw/${run}_1.fastq.gz; r2=$work/raw/${run}_2.fastq.gz
minimap2 -t 8 -ax sr $work/ref.fa $r1 $r2 2>$work/mm.log | samtools flagstat - > $out/flagstat.txt
$GMG guess --index $INDEX --top 5 --min-abund auto --format json $r1 $r2 > $out/guess.json 2> $out/guess.err
mapped=$(awk '/ mapped \(/{sub(/\(/,"",$5); print $5; exit}' $out/flagstat.txt)
primary=$(awk '/ primary mapped \(/{print $5; exit}' $out/flagstat.txt)
printf 'run\tlabel\treference\tpairs_fetched\tmapped_pct\tprimary_mapped_pct\n%s\t%s\t%s\t%s\t%s\t%s\n' "$run" "$label" "$acc" "$N" "$mapped" "$primary" > $out/summary.tsv
cat $out/summary.tsv
