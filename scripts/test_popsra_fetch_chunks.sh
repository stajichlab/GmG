#!/bin/bash
# Test for scripts/popsra_fetch_chunks.slurm with a fake popseq fetch: each run once with the right arguments, parallelism
# at the limit, failures do not fail the job, a bad manifest does.
set -uo pipefail
HERE="$(cd "$(dirname "$0")" && pwd)"
T=$(mktemp -d "${TMPDIR:-/tmp}/chunkjob.XXXXXX"); trap 'rm -rf "$T"' EXIT
check() { local d=$1; shift; if "$@"; then echo "ok: $d"; else echo "FAIL: $d" >&2; exit 1; fi; }

# a fake fetch: records its arguments, tracks how many run at once, exits 12 for a run called FAILxx
cat > "$T/fakefetch" <<'FAKE'
#!/bin/bash
run=$1
mkdir -p "$ACTIVE"; touch "$ACTIVE/$run"
n=$(ls "$ACTIVE" | wc -l | tr -d ' '); echo "$n" >> "$ACTIVE.max"
echo "$*" >> "$CALLS"
sleep 0.4
rm -f "$ACTIVE/$run"
case "$run" in FAIL*) exit 12;; esac
exit 0
FAKE
chmod +x "$T/fakefetch"
{ echo "run_accession,study"; for r in R01 R02 FAIL03 R04 R05 R06 R07 R08 R09 R10; do echo "$r,S"; done; } > "$T/manifest.csv"

export POPSEQ_ROOT=/x GMG_ROOT="$HERE/.." POPSKETCH=/x/popsketch WORK="$T/work" MANIFEST="$T/manifest.csv" SCRATCH="$T/scratch" \
       SKIP_MODULES=1 FETCH_CMD="$T/fakefetch" ACTIVE="$T/active" CALLS="$T/calls" CHUNK=5 PAR=2 CONCURRENCY=8 SLURM_CPUS_PER_TASK=8
for task in 0 1; do SLURM_ARRAY_TASK_ID=$task bash "$HERE/popsra_fetch_chunks.slurm" > "$T/job$task.out" 2>&1; check "chunk $task exits 0" test $? -eq 0; done

check "every run was fetched exactly once" test "$(cut -d' ' -f1 "$T/calls" | sort | uniq -d | wc -l | tr -d ' ')" = 0
check "all 10 runs were fetched" test "$(wc -l < "$T/calls" | tr -d ' ')" = 10
check "two at a time, no more" test "$(sort -n "$T/active.max" | tail -1)" = 2
check "the tier, coverage, sampling and pairs are passed" grep -q -E -- "--tier coverage --coverage 25 --sampling stride --out $T/work/reads --exec $HERE/../scripts/popsra_sketch_hook.sh --threads 4 --require-pairs --concurrency 8 --results $T/work/tasks/[01]/results.tsv --profile $T/work/tasks/[01]/profile.tsv" "$T/calls"
check "a failing run is in its task's exits.tsv with its code" grep -q "^FAIL03	12	" "$T/work/tasks/0/exits.tsv"
check "exits.tsv of the two tasks have a line per run" test "$(cat "$T/work/tasks/0/exits.tsv" "$T/work/tasks/1/exits.tsv" | wc -l | tr -d ' ')" = 10
check "no flat exits.tsv is written (tasks on different nodes must not share a file)" test ! -e "$T/work/exits.tsv"
check "popsra_merge.py joins them" python3 "$HERE/popsra_merge.py" "$T/work" > /dev/null
check "the merged exits.tsv has all 10 runs" test "$(wc -l < "$T/work/merged/exits.tsv" | tr -d ' ')" = 10
check "each run has its own log files" test -f "$T/work/logs/runs/R10.out" -a -f "$T/work/logs/runs/R10.err"
check "the threads are the cores divided by the parallelism" grep -q -- "--threads 4" "$T/calls"

# a chunk past the end does nothing and succeeds; a manifest with no run_accession column fails with 2
SLURM_ARRAY_TASK_ID=2 bash "$HERE/popsra_fetch_chunks.slurm" > "$T/job2.out" 2>&1
check "a chunk past the end exits 0 with nothing to do" grep -q "nothing to do" "$T/job2.out"
printf "acc\nR01\n" > "$T/bad.csv"
rc=0; MANIFEST="$T/bad.csv" SLURM_ARRAY_TASK_ID=0 bash "$HERE/popsra_fetch_chunks.slurm" > "$T/jobbad.out" 2>&1 || rc=$?
check "a bad manifest exits 2" test "$rc" = 2
rc=0; env -u CONCURRENCY SLURM_ARRAY_TASK_ID=0 bash "$HERE/popsra_fetch_chunks.slurm" > "$T/jobnc.out" 2>&1 || rc=$?
check "a missing CONCURRENCY is refused" test "$rc" != 0

# an accession is an argument, never shell code: a hostile one does not run anything (no "/" in it: it names a log file)
{ echo "run_accession,study"; echo 'X$(touch pwned);y,S'; } > "$T/evil.csv"
( cd "$T" && MANIFEST="$T/evil.csv" WORK="$T/workevil" SLURM_ARRAY_TASK_ID=0 bash "$HERE/popsra_fetch_chunks.slurm" > "$T/jobevil.out" 2>&1 )
check "a shell metacharacter in an accession is not executed" test ! -e "$T/pwned"
check "and the run was passed to fetch whole" grep -qF 'X$(touch pwned);y ' "$T/calls"

# a quoted comma in the manifest does not shift the columns: the chunk gets the accession, not a piece of another column
printf 'sample,geo,run_accession\nS1,"USA: California, Riverside",RQ01\n' > "$T/quoted.csv"
MANIFEST="$T/quoted.csv" WORK="$T/workq" SLURM_ARRAY_TASK_ID=0 bash "$HERE/popsra_fetch_chunks.slurm" > "$T/jobq.out" 2>&1
check "a quoted field with a comma gives the right accession" grep -q "^RQ01 " "$T/calls"

# every fetch exits 127 (the command is not there): that is a failed job, not ten failed runs
rc=0; FETCH_CMD="$T/does-not-exist" WORK="$T/work127" SLURM_ARRAY_TASK_ID=0 bash "$HERE/popsra_fetch_chunks.slurm" > "$T/job127.out" 2>&1 || rc=$?
check "all fetches exiting 127 fails the job" test "$rc" = 1
# but an earlier attempt's 127 lines do not count against a later, working attempt of the same task
FETCH_CMD="$T/fakefetch" WORK="$T/work127" SLURM_ARRAY_TASK_ID=0 bash "$HERE/popsra_fetch_chunks.slurm" > "$T/job127b.out" 2>&1
check "a later working attempt of the task succeeds" test "$?" = 0
echo "all chunk-job tests passed"
