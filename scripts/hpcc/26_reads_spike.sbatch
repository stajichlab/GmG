#!/bin/bash -l
#SBATCH -J gmg_spike
#SBATCH -p short
#SBATCH -c 16
#SBATCH --mem=48G
#SBATCH -t 02:00:00
#SBATCH -o /bigdata/stajichlab/shared/projects/GmG_runs/logs/26_spike.%j.log
# Contaminated reads (#45): real fungal read pairs (8 Coccidioides with species truth, 12 ZyGoLife with genus truth) plus simulated bacterial
# (E. coli, B. subtilis) or human (chr20) read pairs, 1-50% of the pairs; classified with the default read settings (min abundance 2) and with
# --min-abund 1, scored against the host lineage. Needs the titration run (25_reads_titration.sbatch) for the ZyGoLife hosts.
set -uo pipefail
source $HOME/projects/GmG/scripts/hpcc/env.sh
module load ncbi_datasets/18.30.1
S=${READS_SCRIPTS:-$GMG_REPO/scripts}
REL=$GR/release/2026.10.3; [ -s $REL/bfd-2026.10.3-s10000.gmgi ] || REL=$GR/release/2026.10.2
INDEX=${INDEX:-$(ls $REL/bfd-20*-s10000.gmgi)}
T=$GR/reads/titration; R=$GR/reads/spike; W=$R/work; mkdir -p $W/donors; cd $R
echo "index $INDEX"
python3 $S/reads_spike_plan.py $GMG_REPO/data/validation/cocci_reads.truth.tsv $GR/reads/cocci_r1.list $T . ${N_COCCI:-8} ${N_ZYGO:-12}
cd $W/donors
if [ ! -s bacteria.fna ]; then
  printf 'GCF_000005845.2\nGCF_000009045.1\n' > bac.txt
  datasets download genome accession --inputfile bac.txt --include genome --filename b.zip --no-progressbar && unzip -oq b.zip -d bx && cat bx/ncbi_dataset/data/GC*/*.fna > bacteria.fna
  datasets download genome accession GCF_000001405.40 --chromosomes 20 --include genome --filename h.zip --no-progressbar && unzip -oq h.zip -d hx && cat hx/ncbi_dataset/data/GC*/*.fna > human.fna
fi
cd $R
# one host at a time: write its variants (plain FASTQ), classify them in both abundance modes, delete them
one() {
  id=$1; a=$2; b=$3
  d=$W/$id; mkdir -p $d
  python3 $S/reads_spike.py $id $a $b $W/donors/bacteria.fna $W/donors/human.fna $d > $d/sheet_body.tsv || return 1
  (echo -e "id\tpath\texpected"; awk -F'\t' '{print $1"\t"$2"\t"}' $d/sheet_body.tsv) > $d/sheet.tsv
  $GMG --threads 2 run --sample-sheet $d/sheet.tsv --index $INDEX -o $R/out_default/$id --no-report -q --max-reads 2000000
  $GMG --threads 2 run --sample-sheet $d/sheet.tsv --index $INDEX -o $R/out_abund1/$id --no-report -q --max-reads 2000000 --min-abund 1
  rm -rf $d
}
export -f one; export S W R GMG INDEX
mkdir -p $R/out_default $R/out_abund1
while IFS=$'\t' read -r id a b; do echo "$id $a $b"; done < hosts.tsv | xargs -P 8 -L 1 bash -c 'one "$0" "$1" "$2"'
for M in default abund1; do
  head -1 $(ls $R/out_$M/*/results.tsv | head -1) > results_$M.tsv
  for f in $R/out_$M/*/results.tsv; do tail -n +2 $f; done >> results_$M.tsv
  $GMG validate results_$M.tsv truth.tsv -o v_$M --by donor frac --title "Contaminated reads, min abundance $M"
done
cat v_default/validation.md | head -40
