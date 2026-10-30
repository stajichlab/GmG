#!/bin/bash -l
#SBATCH -J gmg_sratit
#SBATCH -p short
#SBATCH -c 16
#SBATCH --mem=32G
#SBATCH -t 02:00:00
#SBATCH -o /bigdata/stajichlab/shared/projects/GmG_runs/logs/32_sra_titration.%j.log
# Real SRA/ENA reads for #42, step 3: classify the harvested runs from the first 5k ... 500k read pairs with the run's own genome (`self`), its species or its
# genus hidden from the reference (`ref` column + --exclude-level, #85); paired against R1 only; score with `gmg validate`; depth tables; delete the raw reads.
set -uo pipefail
source $HOME/projects/GmG/scripts/hpcc/env.sh
S=${READS_SCRIPTS:-$GMG_REPO/scripts}
INDEX=${INDEX:-$GR/release/2026.10.4/bfd-2026.10.4-s10000.gmgi}
R=$GR/reads_sra; cd $R
echo "gmg: $($GMG --version)  index: $INDEX"
python3 - <<'PY'
import csv, os
# only runs whose two mates both reached the full 500,000 pairs (a short stream is never silently accepted)
full = {l.split("\t")[0] for l in open("raw/bytes.tsv") if l.split("\t")[0] != "run" and int(l.split("\t")[3]) >= 500000}
runs = [r for r in csv.DictReader(open("sra_runs.tsv"), delimiter="\t") if os.path.exists(f"raw/{r['run']}.done") and r["run"] in full]
with open("sheet.tsv", "w") as sh, open("sheet_r1.tsv", "w") as s1, open("truth.tsv", "w") as tr:
    sh.write("id\tpath\texpected\tref\n"); s1.write("id\tpath\texpected\tref\n")
    w = csv.writer(tr, delimiter="\t", lineterminator="\n")
    w.writerow(["id", "phylum", "class", "order", "family", "genus", "species", "center", "species_in_ref", "genus_in_ref"])
    for r in runs:
        p = os.path.abspath(f"raw/{r['run']}")
        sh.write(f"{r['run']}\t{p}_1.fastq.gz\t{r['species']}\t{r['bfd_id']}\n{r['run']}\t{p}_2.fastq.gz\t\t\n")
        s1.write(f"{r['run']}\t{p}_1.fastq.gz\t{r['species']}\t{r['bfd_id']}\n")
        w.writerow([r["run"], r["phylum"], r["cls"], r["order"], r["family"], r["genus"], r["species"], r["center"],
                    str(int(r["species_others_in_ref"]) > 0), str(int(r["genus_others_in_ref"]) > 0)])
print(len(runs), "runs with reads")
PY
SPECS=""
for M in self species genus; do
  for N in 5000 15000 50000 150000 500000; do
    D=${M}_pair_n$N
    [ -s $D/results.tsv ] || $GMG --threads ${SLURM_CPUS_PER_TASK} run --sample-sheet sheet.tsv --index $INDEX -o $D --no-report -q --max-reads $N --exclude-level $M
    $GMG validate $D/results.tsv truth.tsv -o v_$D --by center phylum species_in_ref --title "SRA reads, $M hidden, first $N pairs" > /dev/null
    SPECS="$SPECS $M:pair:$N=$D/results.tsv"
  done
  for N in 15000 150000 500000; do
    D=${M}_single_n$N
    [ -s $D/results.tsv ] || $GMG --threads ${SLURM_CPUS_PER_TASK} run --sample-sheet sheet_r1.tsv --index $INDEX -o $D --no-report -q --max-reads $N --exclude-level $M
    $GMG validate $D/results.tsv truth.tsv -o v_$D --by center phylum species_in_ref --title "SRA reads, $M hidden, R1 only, first $N reads" > /dev/null
    SPECS="$SPECS $M:single:$N=$D/results.tsv"
  done
done
python3 $S/sra_depth_table.py sra_runs.tsv truth.tsv depth_tables.md $SPECS
# the raw reads are scratch: delete them once everything is scored
if [ -s depth_tables.md ] && [ -s self_pair_n500000/results.tsv ]; then rm -rf raw; echo "raw reads deleted"; fi
