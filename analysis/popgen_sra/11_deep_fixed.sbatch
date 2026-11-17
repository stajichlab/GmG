#!/bin/bash
#SBATCH -p stajichlab -N1 -n1 -c4 --mem=12G -t 4:00:00 -J deepfix -o logs/%x_%A_%a.out
# #159 step 2: the 6 Gb and whole-run sketches were deleted before the high-coverage fix (PR #158). Re-fetch deep runs (rows of depth/tasks.tsv: slug, run, base_count)
# and sketch the 6 Gb subsample and the whole run with the FIXED binary (FIXBIN), into depth/<slug>/s6g_fix/ and depth/<slug>/full_fix/.
# Compare with the _out_* (old binary) and s1g/s2g/s4g sketches using 09_validate_estimator.py or `popsketch check`.
# usage: FIXBIN=<fixed popsketch> sbatch --array=1,3,6,7,13,14 11_deep_fixed.sbatch      (temp files in $SCRATCH)
set -uo pipefail
module load sratoolkit/3.2.0 parallel-fastq-dump/0.6.7 2>/dev/null || true
P=/bigdata/stajichlab/jstajich/projects/GmG_popgen
BIN=${FIXBIN:-$P/popsketch_target_fix/release/popsketch}
cd $P
read -r slug run bases < <(awk -F'\t' -v n=$SLURM_ARRAY_TASK_ID 'NR==n' depth/tasks.tsv)
T=depth/$slug
[[ -s $T/full_fix/$run.json ]] && { echo "skip $run"; exit 0; }
echo "binary: $BIN ($($BIN --version 2>&1 | head -1)); run $run bases=$bases"
work=$(mktemp -d "${SCRATCH:-/tmp}/deepfix.$run.XXXXXX"); trap 'rm -rf "$work"' EXIT
parallel-fastq-dump --sra-id $run --threads 4 --split-files --gzip --tmpdir $work -O $work/raw >$work/dump.log 2>&1 || { echo "FAIL dump $run"; tail -5 $work/dump.log; exit 1; }
r1=$work/raw/${run}_1.fastq.gz; r2=$work/raw/${run}_2.fastq.gz
[[ -s $r1 && -s $r2 ]] || { echo "FAIL mates $run"; exit 1; }
f=$(awk -v t=6000000000 -v b=$bases 'BEGIN{printf "%.9f", (t>=b)?1:t/b}')
for m in 1 2; do
  zcat $work/raw/${run}_$m.fastq.gz | awk -v f=$f 'NR%4==1{i++; k=(int(i*f)!=int((i-1)*f))} k' | gzip > $work/s6_$m.fastq.gz
done
$BIN sketch --id $run -o $T/s6g_fix --max-bases 0 $work/s6_1.fastq.gz $work/s6_2.fastq.gz || { echo "FAIL sketch s6g $run"; exit 1; }
rm -f $work/s6_?.fastq.gz
$BIN sketch --id $run -o $T/full_fix --max-bases 0 $r1 $r2 || { echo "FAIL sketch full $run"; exit 1; }
echo "done $run"
