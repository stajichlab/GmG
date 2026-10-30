#!/bin/bash
# Reads validations (#42/#43): the same runs scored by call status, so the precision of confident (ok) calls is separate from low_confidence ones.
# Light (gmg validate on a few hundred rows): fine on the login node. usage: 29_reads_by_status.sh   (run where $GR/$GMG are defined)
set -uo pipefail
source $HOME/projects/GmG/scripts/hpcc/env.sh
S=${READS_SCRIPTS:-$GR/reads/scripts}
T=$GR/reads/titration; X=$GR/reads/rna_exclusion
for N in 3000 15000 75000 300000 1500000; do
  $GMG validate $T/n$N/results.tsv $T/truth.tsv -o $T/vs_$N --by status genus_in_ref --title "ZyGoLife reads, $N pairs, by status" > /dev/null 2>&1
done
for N in 100000 500000; do for L in none species genus; do
  $GMG validate $X/n${N}_$L/results.tsv $X/truth.tsv -o $X/vs_n${N}_$L --by status --title "RNA reads $N, exclusion $L, by status" > /dev/null 2>&1
done; done
for G in status=ok status=low_confidence; do
  echo "## titration $G"; python3 $S/reads_summary.py $G $(for N in 3000 15000 75000 300000 1500000; do echo -n "n$N=$T/vs_$N/validation.tsv "; done)
  echo "## RNA $G"; python3 $S/reads_summary.py $G $(for N in 100000 500000; do for L in none species genus; do echo -n "n${N}_$L=$X/vs_n${N}_$L/validation.tsv "; done; done)
done
