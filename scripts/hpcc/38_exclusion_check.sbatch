#!/bin/bash -l
#SBATCH -J gmg_exclcheck
#SBATCH -p short
#SBATCH -c 8
#SBATCH --mem=24G
#SBATCH -t 02:00:00
#SBATCH -o /bigdata/stajichlab/shared/projects/GmG_runs/logs/38_exclcheck.%j.log
# Check of a release that leaves genomes out (2026.10.5: the 15 tier A/B genomes of the label hand-off report).
# usage: sbatch 38_exclusion_check.sbatch OLD NEW [N_LINES]
# N_LINES = how many last lines of reference_exclusions.txt are new in NEW (default 15). For each such genome, and up to 3 cluster mates of each
# (genomes that stay in both releases): `gmg run` on the BFD assembly against OLD and NEW (same binary). Expect: the excluded genome no longer
# finds itself (its top hit is another genome, ANI < 1), its cluster mates keep their calls.
set -uo pipefail
source $HOME/projects/GmG/scripts/hpcc/env.sh
OLD=${1:-2026.10.4}; NEW=${2:-2026.10.5}; N=${3:-15}
R=$GR/exclcheck_$NEW; mkdir -p $R/ex $R/mates; cd $R
O=$GR/release/$OLD/bfd-$OLD-s10000.gmgi; W=$GR/release/$NEW/bfd-$NEW-s10000.gmgi
tail -$N $GMG_REPO/reports/label_audit/reference_exclusions.txt > excluded_new.txt
python3 - <<PY
import csv
meta = list(csv.DictReader(open("$GR/release/$OLD/meta_full.tsv"), delimiter="\t"))
ex = set(open("excluded_new.txt").read().split())
by = {}
for r in meta:
    by.setdefault(r["cluster"], []).append(r["id"])
mates = []
for g in ex:
    c = next((r["cluster"] for r in meta if r["id"] == g), None)
    if c is None:
        continue
    mates += [m for m in by[c] if m not in ex][:3]
open("mates.txt", "w").write("\n".join(dict.fromkeys(mates)) + "\n")
PY
for i in $(cat excluded_new.txt); do ln -sf $BFD_GENOMES/$i.fa.gz ex/; done
for i in $(cat mates.txt); do ln -sf $BFD_GENOMES/$i.fa.gz mates/; done
for v in $OLD:$O $NEW:$W; do
  n=${v%%:*}; ix=${v#*:}
  for s in ex mates; do
    $GMG --threads ${SLURM_CPUS_PER_TASK} -q run $s --index $ix -o out_${s}_$n --no-report
  done
done
python3 - <<PY
import csv
def rd(p):
    return {r["query_id"]: r for r in csv.DictReader(open(p), delimiter="\t")}
rk = ["phylum", "class", "order", "family", "genus", "species"]
for s in ("ex", "mates"):
    a, b = rd(f"out_{s}_$OLD/results.tsv"), rd(f"out_{s}_$NEW/results.tsv")
    same = selfhit = 0
    print(f"== {s}: {len(a)} queries (old $OLD, new $NEW)")
    for q in a:
        ca = [a[q]["call_" + r] for r in rk]
        cb = [b[q]["call_" + r] for r in rk]
        same += ca == cb
        old_self = a[q]["top_hit"] == q
        new_self = b[q]["top_hit"] == q
        selfhit += new_self
        print(q[:36], "| old top", a[q]["top_hit"][:28], a[q]["ani"], "| new top", b[q]["top_hit"][:28], b[q]["ani"],
              "| calls same" if ca == cb else "| CALLS DIFFER: " + "/".join(ca) + " -> " + "/".join(cb), "| NEW SELF-MATCH" if new_self else "", "| old self" if old_self else "")
    print(f"{s}: calls identical {same} of {len(a)}; new top hit is the query itself: {selfhit}")
PY
