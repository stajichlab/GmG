#!/bin/bash
# Test for scripts/popsra_sketch_hook.sh on popsketch's synthetic reads (no data, no network).
#   scripts/test_popsra_sketch_hook.sh [path/to/popsketch]
set -euo pipefail
HERE="$(cd "$(dirname "$0")" && pwd)"
export POPSKETCH=${1:-$HERE/../target/debug/popsketch}
[ -x "$POPSKETCH" ] || { echo "no popsketch binary at $POPSKETCH (cargo build -p popsketch)" >&2; exit 1; }
T=$(mktemp -d "${TMPDIR:-/tmp}/hooktest.XXXXXX"); trap 'rm -rf "$T"' EXIT
"$POPSKETCH" demo -o "$T/demo" >/dev/null 2>&1
export SKETCH_DIR="$T/sk" SKETCH_ARGS="--scaled 5" HOOK_LOG="$T/hook.log"
mapfile -t FILES < <(ls "$T/demo/reads"/*.fastq | head -8)
# check DESCRIPTION COMMAND...: a false condition fails the test loudly (a bare `[ ] && echo` would pass silently)
check() { local d=$1; shift; if "$@"; then echo "ok: $d"; else echo "FAIL: $d" >&2; exit 1; fi; }

# a paired call (two files) and a single-end call (one file), as popseq fetch makes them
check "paired call" "$HERE/popsra_sketch_hook.sh" PAIRED1 "${FILES[0]}" "${FILES[1]}"
for i in 2 3 4 5 6 7; do check "single-end call S$i" "$HERE/popsra_sketch_hook.sh" "S$i" "${FILES[$i]}"; done
for part in PAIRED1.json PAIRED1.gsk PAIRED1.half_a.gsk PAIRED1.half_b.gsk; do
  check "sketch file $part exists" test -s "$SKETCH_DIR/$part"
done
check "one log line per run (7)" test "$(tr -d ' ' < <(wc -l < "$HOOK_LOG"))" = 7
check "every log line has seconds and bytes" awk -F'\t' '$3 > 0 && $2 >= 0 {n++} END {exit n == 7 ? 0 : 1}' "$HOOK_LOG"

# the sketches are readable as a set: popsketch run --from-sketches takes the folder
check "popsketch run --from-sketches reads the folder" "$POPSKETCH" run --from-sketches "$SKETCH_DIR" -o "$T/out" --scaled 5 --max-bases 0 --flagged mark --no-pca --no-tree -q
check "and writes popsketch.json" test -s "$T/out/popsketch.json"

# failures are loud: no files, a missing file, and a missing SKETCH_DIR
check "no files is refused" bash -c '! "$0" X 2>/dev/null' "$HERE/popsra_sketch_hook.sh"
check "a missing file is refused" bash -c '! "$0" X "$1" 2>/dev/null' "$HERE/popsra_sketch_hook.sh" "$T/nope.fastq"
check "a missing SKETCH_DIR is refused" bash -c '! env -u SKETCH_DIR "$0" X "$1" 2>/dev/null' "$HERE/popsra_sketch_hook.sh" "${FILES[0]}"
echo "all hook tests passed"
