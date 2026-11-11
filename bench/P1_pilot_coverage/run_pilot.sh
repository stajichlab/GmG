#!/bin/bash
# P1 pilot: `popseq fetch --tier coverage --coverage 25` for one run of each of three species, on a
# laptop with no SRA toolkit (ENA is the only usable route). One run per species, the three in
# parallel, so the declared concurrency is 3. Everything here is purgeable: rm -rf work reads popseq.
#
#   bench/P1_pilot_coverage/run_pilot.sh        (after select_pilot.py has written plan.csv)
set -uo pipefail
HERE="$(cd "$(dirname "$0")" && pwd)"
export SCRATCH="$HERE/work/scratch"          # temp files go in $SCRATCH, never /tmp (PopSeq AGENTS section 5)
mkdir -p "$SCRATCH" "$HERE/reads"
RESULTS="$HERE/results.tsv"
printf "run\ttier\troute\tlayout\tbases_kept\treads_kept\tcoverage_target\tcoverage_reached\tseconds\tverified\toutcome\n" > "$RESULTS"
UV=$(command -v uv)
PROJ="$HERE/popseq"
# species taxid : run (the smallest run of each, from select_pilot.py)
for pair in 5478:ERR4669779 5007:ERR12524156 5476:SRR31272388; do
  tid=${pair%%:*}; run=${pair##*:}
  ( "$UV" run --project "$PROJ" popseq fetch "$run" --tier coverage --coverage 25 \
      --out "$HERE/reads/$tid" --plan "$HERE/plan.csv" --concurrency 3 --results "$RESULTS" \
      </dev/null > "$HERE/work/$run.out" 2> "$HERE/work/$run.err"; echo "$run exit $?" >> "$HERE/work/exits.txt" ) &
done
wait
echo "done"; cat "$HERE/work/exits.txt"
