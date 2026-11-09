#!/bin/bash
# F1: list the BFD genomes of the fusarioid genera found in the FusariumID records.
# usage: REPO=<checkout> W=<build dir> bash make_lists.sh      (needs $W/fid/records.tsv from parse_fusariumid.py)
# writes $W/genomes/{meta.tsv,ids.txt,refs_fus.txt,paths.txt}
set -euo pipefail
REPO=${REPO:?set REPO}; W=${W:?set W}
source "$REPO/scripts/hpcc/env.sh"
mkdir -p "$W/genomes"
GEN=$(cut -f5 "$W/fid/records.tsv" | tail -n +2 | sort -u | paste -sd'|')
awk -F'\t' -v re="^($GEN)$" 'NR==1 || $7 ~ re' "$GR/index/meta.tsv" > "$W/genomes/meta.tsv"
tail -n +2 "$W/genomes/meta.tsv" | cut -f1 > "$W/genomes/ids.txt"
grep -F -x -f "$W/genomes/ids.txt" "$GR/index/refs_2026.10.5.txt" > "$W/genomes/refs_fus.txt"
: > "$W/genomes/paths.txt"
while read -r id; do f="$BFD_GENOMES/$id.fa.gz"; [ -s "$f" ] && echo "$f" >> "$W/genomes/paths.txt"; done < "$W/genomes/ids.txt"
echo "genomes: $(wc -l < "$W/genomes/ids.txt") ids, $(wc -l < "$W/genomes/paths.txt") files, $(wc -l < "$W/genomes/refs_fus.txt") in the reference set"
