#!/bin/bash -l
#SBATCH -J lr_real
#SBATCH -p short
#SBATCH -c 8
#SBATCH --mem=28G
#SBATCH -t 02:00:00
#SBATCH -o /bigdata/stajichlab/shared/projects/GmG_runs/logs/37_lr_real.%j.log
# Long reads (#44), part 2: the ~20 chosen public runs (data/validation/long_reads_runs.tsv, listed BEFORE the download). Steps: (1) harvest the first 200 Mb
# of every run by streaming (scripts/long_reads_harvest.py, total capped at 10 GB; a run whose first 200 Mb cannot be streamed is retried with
# `fastq-dump -X`, the SRA toolkit fallback that was pre-approved); (2) measure each run's real error by aligning 500 reads to the run's own assembly;
# (3) classify with the own genome hidden (`self`), its species or its genus hidden, min-abund 2, 1 and auto; (4) collect (scripts/long_reads_score.py).
# The raw reads are deleted at the end.
set -uo pipefail
source $HOME/projects/GmG/scripts/hpcc/env.sh
S=${LR_SCRIPTS:-$GR/long_reads/scripts}
GMGBIN=${GMGBIN:-$GR/long_reads/bin/gmg}   # private build of the branch (38_long_reads_build.sbatch)
INDEX=${INDEX:-$GR/release/2026.10.4/bfd-2026.10.4-s10000.gmgi}
RUNS=${RUNS:-$GR/long_reads/long_reads_runs.tsv}
W=$GR/long_reads/real; mkdir -p $W/raw; cd $W
echo "gmg: $($GMGBIN --version)  index: $INDEX"
python3 $S/long_reads_harvest.py $RUNS raw 200e6 3 10e9
# fallback for runs that failed to stream
module load sratoolkit 2>/dev/null
for f in raw/*.failed; do
  [ -e "$f" ] || continue
  run=$(basename $f .failed)
  ml=$(awk -F'\t' -v r=$run 'NR>1 && $1==r{print $18}' $RUNS); n=$(( 200000000 / ${ml:-10000} ))
  echo "fastq-dump -X $n for $run"
  fastq-dump -X $n --stdout $run 2>/dev/null | gzip -1 > raw/$run.fastq.gz && [ -s raw/$run.fastq.gz ] && rm $f && touch raw/$run.done
done
# (2) real error rate and the augmented run table
printf "run\tn_reads\tn_aligned\tmedian_len\tmean_len\terror_pct\tmismatch_pct\tindel_pct\n" > measured.tsv
tail -n +2 $RUNS | while IFS=$'\t' read -r run rest; do
  [ -e raw/$run.done ] || continue
  gid=$(awk -F'\t' -v r=$run '$1==r{print $11}' $RUNS)
  python3 $S/long_reads_measure.py $run raw/$run.fastq.gz $BFD_GENOMES/clean/$gid.purge.fasta.gz 500 >> measured.tsv
done
python3 - <<PY
import csv, os
runs = list(csv.DictReader(open("$RUNS"), delimiter="\t"))
m = {l[0]: l for l in list(csv.reader(open("measured.tsv"), delimiter="\t"))[1:]}
b = {l.split("\t")[0]: l.split("\t") for l in open("raw/bytes.tsv").read().splitlines()[1:]}
out = []
for r in runs:
    if r["run"] not in m: continue
    r = dict(r); r["measured_error_pct"] = m[r["run"]][5]; r["measured_mismatch_pct"] = m[r["run"]][6]; r["measured_indel_pct"] = m[r["run"]][7]
    r["base_count_used"] = b.get(r["run"], ["", "", "200000000"])[2]; r["bytes_read"] = b.get(r["run"], ["", ""])[1]
    out.append(r)
with open("runs_aug.tsv", "w") as f:
    w = csv.DictWriter(f, list(out[0]), delimiter="\t", lineterminator="\n"); w.writeheader(); w.writerows(out)
with open("sheet.tsv", "w") as f:
    f.write("id\tpath\texpected\tref\n")
    for r in out: f.write(f"{r['run']}\t{os.getcwd()}/raw/{r['run']}.fastq.gz\t{r['species']}\t{r['bfd_id']}\n")
PY
run() { tag=$1; shift; [ -s $tag/results.tsv ] || $GMGBIN --threads ${SLURM_CPUS_PER_TASK} run --sample-sheet sheet.tsv --index $INDEX -o $tag --no-report -q --kind reads "$@"; }
run self_a2 --exclude-level self --min-abund 2
run self_a1 --exclude-level self --min-abund 1
run self_aauto --exclude-level self --min-abund auto
run species_a2 --exclude-level species --min-abund 2
run genus_a2 --exclude-level genus --min-abund 2
python3 $S/long_reads_score.py collect-real runs_aug.tsv . queries_real.tsv $GR/release/2026.10.4/meta_full.tsv
awk -F'\t' 'NR>1{n+=$2}END{printf "total compressed bytes read: %.2f GB\n", n/1e9}' raw/bytes.tsv
[ -s queries_real.tsv ] && [ "${DELETE_RAW:-0}" = 1 ] && rm -rf raw && echo "raw reads deleted"
