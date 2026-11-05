#!/bin/bash
# one sample: first NREADS reads of R1 dealt alternately into two files, sketched separately (scaled 200, abundance)
# usage: 20_splithalf_task.sh DATASET ID R1
set -uo pipefail
ds=$1; id=$2; r1=$3
P=/bigdata/stajichlab/jstajich/projects/GmG_popgen
GMG=/bigdata/stajichlab/shared/projects/GmG_runs/bin/gmg
NREADS=${NREADS:-8000000}
out=$P/splithalf/$ds/sk; mkdir -p "$out"
[ -s "$out/${id}_a.gsk" ] && [ -s "$out/${id}_b.gsk" ] && exit 0
tmp=$(mktemp -d "${TMPDIR:-/tmp}/sh.XXXXXX")
zcat "$r1" 2>/dev/null | head -n $((NREADS * 4)) | awk -v fa="$tmp/${id}_a.fastq" -v fb="$tmp/${id}_b.fastq" '{r=int((NR-1)/4)%2; print > (r ? fb : fa)}'
$GMG sketch -q --threads 2 --abund -k 21 -s 200 -o "$out" "$tmp/${id}_a.fastq" "$tmp/${id}_b.fastq"
rm -rf "$tmp"
