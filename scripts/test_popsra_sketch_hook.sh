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

# the sampling check: a genome the sketch matches passes; one it cannot cover fails with 75 and leaves no sketch behind
export CHECK_LOG="$T/check.tsv"
check "a matching genome size passes the sampling check" env POPSEQ_GENOME_MB=0.17 POPSEQ_COVERAGE=20 \
  "$HERE/popsra_sketch_hook.sh" CHK1 "${FILES[0]}"
check "the check line is logged with a header" test "$(head -1 "$CHECK_LOG" | cut -f1)" = id
rc=0; env POPSEQ_GENOME_MB=5 POPSEQ_COVERAGE=20 "$HERE/popsra_sketch_hook.sh" CHK2 "${FILES[0]}" 2>/dev/null || rc=$?
check "an undersampled sample exits 75" test "$rc" = 75
check "and its sketch files are removed" bash -c '! ls "$0"/CHK2.* >/dev/null 2>&1' "$SKETCH_DIR"
rc=0; env POPSEQ_GENOME_MB=0.01 POPSEQ_COVERAGE=20 "$HERE/popsra_sketch_hook.sh" CHK4 "${FILES[0]}" 2>/dev/null || rc=$?
check "a sample far larger than the genome is not assessed: exit 65" test "$rc" = 65
check "and its sketch files are removed too" bash -c '! ls "$0"/CHK4.* >/dev/null 2>&1' "$SKETCH_DIR"
rc=0; env POPSEQ_GENOME_MB=0 "$HERE/popsra_sketch_hook.sh" CHK5 "${FILES[0]}" 2>/dev/null || rc=$?
check "an error of the check itself (genome size 0) exits 2" test "$rc" = 2
check "and removes the sketch (it never passed the check)" bash -c '! ls "$0"/CHK5.* >/dev/null 2>&1' "$SKETCH_DIR"
check "and no blank line goes to the check log" bash -c '! grep -q "^$" "$0"' "$CHECK_LOG"
check "a sample that passed keeps its sketch" test -s "$SKETCH_DIR/CHK1.gsk"
check "POPSKETCH_CHECK=off skips the check" env POPSKETCH_CHECK=off POPSEQ_GENOME_MB=5 "$HERE/popsra_sketch_hook.sh" CHK3 "${FILES[0]}"

# hooks that finish together write the header once and every data line whole
export CHECK_LOG="$T/check_race.tsv"
for i in 0 1 2 3 4 5; do
  env POPSEQ_GENOME_MB=0.17 POPSEQ_COVERAGE=20 "$HERE/popsra_sketch_hook.sh" "RACE$i" "${FILES[$i]}" >/dev/null 2>&1 &
done
wait
check "six hooks at once: one header line" test "$(grep -c '^id	' "$CHECK_LOG")" = 1
check "and six data lines" test "$(grep -c '^RACE' "$CHECK_LOG")" = 6
check "the header is the first line" test "$(head -1 "$CHECK_LOG" | cut -f1)" = id
unset CHECK_LOG

# the species check (PopSeq#52). Sketching and the sampling check are real; `species-match` is a stub that returns the status in
# FAKE_STATUS (the verdict logic itself is tested in crates/popsketch/src/species_match.rs).
export REAL_POPSKETCH="$POPSKETCH"
cat > "$T/fake_popsketch" <<'STUB'
#!/bin/bash
if [ "$1" = species-match ]; then
  while [ $# -gt 0 ]; do case "$1" in --filed) F=$2; shift ;; -o) O=$2; shift ;; esac; shift; done
  [ "${FAKE_STATUS:-match}" = fail ] && { echo "boom" >&2; exit 3; }
  id=$(sed -n 2p "$F" | cut -f1); name=$(sed -n 2p "$F" | cut -f2)
  case "${FAKE_STATUS:-match}" in mismatch|host_dominated|no_fungal_call_suspect) x=yes ;; *) x=no ;; esac
  printf 'id\tfiled_name\tstatus\texclude_default\tdetected_genus\tdetected_species\tbest_shared_hashes\tbest_share\tpurity\treason\n%s\t%s\t%s\t%s\t\t\t10\t0.5\t0.1\treason text\n' "$id" "$name" "${FAKE_STATUS:-match}" "$x" > "$O"
  echo '{"tool":"fake"}' > "$O.meta.json"; exit 0
fi
exec "$REAL_POPSKETCH" "$@"
STUB
chmod +x "$T/fake_popsketch"
: > "$T/fake.gmgi"
printf 'id\tfiled\nSP1\tCryptococcus neoformans\nSP2\tCryptococcus neoformans\nSP3\tCryptococcus neoformans\nSP4\tCryptococcus neoformans\nSP6\tCryptococcus neoformans\nSP7\tCryptococcus neoformans\n' > "$T/filed.tsv"
sp() { env POPSKETCH="$T/fake_popsketch" SPECIES_LOG="$T/species.tsv" SPECIES_INDEX="$T/fake.gmgi" SPECIES_FILED="$T/filed.tsv" "$@"; }
check "a matching verdict keeps the run (exit 0)" sp FAKE_STATUS=match "$HERE/popsra_sketch_hook.sh" SP1 "${FILES[0]}"
check "and its sketch stays in the folder" test -s "$SKETCH_DIR/SP1.gsk"
check "the verdict row has a header and one data line" test "$(wc -l < "$T/species.tsv" | tr -d ' ')" = 2
check "the index and cutoffs are kept once, next to the log" test -s "$T/species.meta.json"
rc=0; sp FAKE_STATUS=mismatch "$HERE/popsra_sketch_hook.sh" SP2 "${FILES[0]}" 2>"$T/sp2.err" || rc=$?
check "an excluded verdict under refuse exits 66" test "$rc" = 66
check "the sketch is moved to refused/, not deleted" test -s "$SKETCH_DIR/refused/SP2.gsk"
check "and is no longer in the folder popsketch run reads" bash -c '! ls "$0"/SP2.* >/dev/null 2>&1' "$SKETCH_DIR"
check "the refusal says why" grep -q "species check refused the run (mismatch)" "$T/sp2.err"
check "the refused run still has its verdict row (exclude_default yes)" awk -F'\t' '$1 == "SP2" && $3 == "mismatch" && $4 == "yes" {n++} END {exit n == 1 ? 0 : 1}' "$T/species.tsv"
check "a set with a refused/ folder is still readable by popsketch run" "$REAL_POPSKETCH" run --from-sketches "$SKETCH_DIR" -o "$T/out2" --scaled 5 --max-bases 0 --flagged mark --no-pca --no-tree -q
rc=0; sp FAKE_STATUS=host_dominated SPECIES_MATCH=warn "$HERE/popsra_sketch_hook.sh" SP3 "${FILES[0]}" 2>"$T/sp3.err" || rc=$?
check "warn keeps an excluded run (exit 0)" test "$rc" = 0
check "and warns" grep -q "WARNING species check says host_dominated" "$T/sp3.err"
check "and the sketch stays" test -s "$SKETCH_DIR/SP3.gsk"
before=$(wc -l < "$T/species.tsv")
check "off skips the check and writes no row" sp SPECIES_MATCH=off FAKE_STATUS=mismatch "$HERE/popsra_sketch_hook.sh" SP4 "${FILES[0]}"
check "(no row added)" test "$(wc -l < "$T/species.tsv")" = "$before"
check "an unreadable index keeps the run and records unchecked" sp SPECIES_INDEX="$T/nope.gmgi" "$HERE/popsra_sketch_hook.sh" SP5 "${FILES[0]}"
check "(the row says unchecked with the reason)" awk -F'\t' '$1 == "SP5" && $3 == "unchecked" && $4 == "no" && $10 ~ /no readable index/ {n++} END {exit n == 1 ? 0 : 1}' "$T/species.tsv"
check "(and the sketch stays)" test -s "$SKETCH_DIR/SP5.gsk"
check "a run with no filed name is recorded as unchecked" sp "$HERE/popsra_sketch_hook.sh" SP8 "${FILES[0]}"
check "(with that reason)" awk -F'\t' '$1 == "SP8" && $3 == "unchecked" && $10 ~ /no filed species name/ {n++} END {exit n == 1 ? 0 : 1}' "$T/species.tsv"
check "species-match failing keeps the run and records unchecked" sp FAKE_STATUS=fail "$HERE/popsra_sketch_hook.sh" SP6 "${FILES[0]}"
check "(with the exit code and message)" awk -F'\t' '$1 == "SP6" && $3 == "unchecked" && $10 ~ /exit 3.*boom/ {n++} END {exit n == 1 ? 0 : 1}' "$T/species.tsv"
check "a bad SPECIES_MATCH value is refused (exit 2)" bash -c '! env SPECIES_MATCH=maybe "$0" X "$1" 2>/dev/null' "$HERE/popsra_sketch_hook.sh" "${FILES[0]}"
check "no SPECIES_LOG: no species check at all" env POPSKETCH="$T/fake_popsketch" FAKE_STATUS=mismatch SPECIES_INDEX="$T/fake.gmgi" "$HERE/popsra_sketch_hook.sh" SP9 "${FILES[0]}"
rm -f "$T/species_race.tsv"
for i in 0 1 2 3 4 5; do
  env POPSKETCH="$T/fake_popsketch" SPECIES_LOG="$T/species_race.tsv" SPECIES_INDEX="$T/fake.gmgi" SPECIES_FILED="$T/filed.tsv" \
    "$HERE/popsra_sketch_hook.sh" SP1 "${FILES[$i]}" >/dev/null 2>&1 &
done
wait
check "six species checks at once: one header line" test "$(grep -c '^id	' "$T/species_race.tsv")" = 1
check "and six whole data lines" test "$(grep -c '^SP1	' "$T/species_race.tsv")" = 6

# failures are loud: no files, a missing file, and a missing SKETCH_DIR
check "no files is refused" bash -c '! "$0" X 2>/dev/null' "$HERE/popsra_sketch_hook.sh"
check "a missing file is refused" bash -c '! "$0" X "$1" 2>/dev/null' "$HERE/popsra_sketch_hook.sh" "$T/nope.fastq"
check "a missing SKETCH_DIR is refused" bash -c '! env -u SKETCH_DIR "$0" X "$1" 2>/dev/null' "$HERE/popsra_sketch_hook.sh" "${FILES[0]}"
echo "all hook tests passed"
