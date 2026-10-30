#!/bin/bash
# Native performance numbers for reports/performance_*.md (#49): load time, peak memory and throughput.
# usage: perf_native.sh GMG_BIN FULL.gmgi REPS.gmgi GENOME_DIR     (GENOME_DIR: the golden panel cache, *.fna.gz)
# macOS: /usr/bin/time -l (bytes); Linux: /usr/bin/time -v (kB). Prints "what<TAB>seconds<TAB>peak_MB".
set -u
G=$1; FULL=$2; REPS=$3; D=$4
if /usr/bin/time -l true >/dev/null 2>&1; then T() { /usr/bin/time -l "$@" 2>&1 >/dev/null | awk '/maximum resident/ {printf "%.0f", $1/1048576}'; }; else T() { /usr/bin/time -v "$@" 2>&1 >/dev/null | awk '/Maximum resident/ {printf "%.0f", $6/1024}'; }; fi
now() { python3 -c 'import time; print(time.time())'; }
run() { local what=$1; shift; local t0 t1 mb; t0=$(now); mb=$(T "$@"); t1=$(now); printf '%s\t%.2f\t%s\n' "$what" "$(echo "$t1 - $t0" | bc -l)" "$mb"; }
bytes=$(cat $D/*.fna.gz | gzip -dc | wc -c)
echo "# panel: $(ls $D/*.fna.gz | wc -l) genomes, $((bytes / 1000000)) MB of sequence"
run "info, full index" $G info $FULL
run "info, reps index" $G info $REPS
for th in 1 8; do
  run "classify panel, full index, $th thread(s)" $G --threads $th -q run $D --index $FULL -o /tmp/perf_out_full --no-report
  run "classify panel, reps index, $th thread(s)" $G --threads $th -q run $D --index $REPS -o /tmp/perf_out_reps --no-report
done
run "sketch panel only (k=21, scaled 10000), 1 thread" $G --threads 1 -q sketch -k 21 -s 10000 -o /tmp/perf_sk $D/*.fna.gz
rm -rf /tmp/perf_out_full /tmp/perf_out_reps /tmp/perf_sk
