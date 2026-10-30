#!/bin/bash
# Tool comparison (#46): the reference set of every tool = the genomes of the GmG representatives index, without the decoys. The reps index of
# 2026.10.3 is the 2026.10.2 reps set minus the 57 ambiguous-label genomes (reports/label_audit/reference_exclusions.txt), so that is what is used
# when release/2026.10.3 does not exist yet. Writes $GR/compare/ref/{refids.txt,reflist.txt} (ids, FASTA paths). Light: run on the login node.
set -euo pipefail
source $HOME/projects/GmG/scripts/hpcc/env.sh
R=$GR/compare/ref; mkdir -p $R; cd $R
META=${META:-$GR/release/2026.10.2/meta_reps.tsv}
EXCL=${GMG_REPO:?}/reports/label_audit/reference_exclusions.txt
tail -n +2 $META | cut -f1 | grep -v -x -F -f $EXCL > refids.all.txt || true
while read id; do f=$BFD_GENOMES/$id.fa.gz; [ -s "$f" ] && echo "$id" || true; done < refids.all.txt > refids.txt
sed "s#^#$BFD_GENOMES/#; s#\$#.fa.gz#" refids.txt > reflist.txt
wc -l refids.all.txt refids.txt reflist.txt
