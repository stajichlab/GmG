#!/bin/bash
# Test for scripts/popsra_chunk.sh: every manifest row goes to exactly one chunk, in order, for many sizes.
set -euo pipefail
HERE="$(cd "$(dirname "$0")" && pwd)"
# shellcheck source=popsra_chunk.sh
. "$HERE/popsra_chunk.sh"
T=$(mktemp -d "${TMPDIR:-/tmp}/chunktest.XXXXXX"); trap 'rm -rf "$T"' EXIT
check() { local d=$1; shift; if "$@"; then echo "ok: $d"; else echo "FAIL: $d" >&2; exit 1; fi; }

for delim in "," $'\t'; do
  for n in 0 1 7 24 25 70; do
    {
      printf "run_accession${delim}study\n"
      for ((r = 1; r <= n; r++)); do printf "RUN%03d${delim}S%d\n" "$r" "$((r % 3))"; done
    } > "$T/m"
    for chunk in 1 4 7 25 100; do
      count=$(chunk_count "$T/m" "$chunk")
      : > "$T/all"
      for ((i = 0; i < count; i++)); do chunk_runs "$T/m" "$chunk" "$i" >> "$T/all"; done
      for ((r = 1; r <= n; r++)); do printf "RUN%03d\n" "$r"; done > "$T/want"
      check "delim '${delim//$'\t'/TAB}' rows $n chunk $chunk: every row once, in order" cmp -s "$T/all" "$T/want"
      check "  chunk_count is ceil(n/chunk)" test "$count" = "$(((n + chunk - 1) / chunk))"
      check "  a chunk past the end is empty" test -z "$(chunk_runs "$T/m" "$chunk" "$count")"
    done
  done
done
printf "study\trun_accession\nS1\tRUN001\n" > "$T/swapped"
check "the run_accession column is found wherever it is" test "$(chunk_runs "$T/swapped" 5 0)" = RUN001
printf "acc\nRUN001\n" > "$T/nocol"
check "a manifest with no run_accession column is refused" bash -c '. "$0"; ! chunk_runs "$1" 5 0 2>/dev/null' "$HERE/popsra_chunk.sh" "$T/nocol"
printf 'sample,geo,run_accession\nS1,"USA: California, Riverside",RUN001\nS2,"x, y, z",RUN002\n' > "$T/quoted"
check "a quoted comma does not shift the columns" test "$(chunk_runs "$T/quoted" 5 0 | tr '\n' ' ')" = "RUN001 RUN002 "
printf 'run_accession,study\nRUN001,S\n\nRUN002,S\n' > "$T/blank"
check "a blank row is not a run" test "$(chunk_count "$T/blank" 1)" = 2
printf 'run_accession,study\n' > "$T/headeronly"
check "a header-only manifest is zero chunks, not an error" test "$(chunk_count "$T/headeronly" 5)" = 0
echo "all chunk tests passed"
