#!/bin/bash -l
#SBATCH -J gmg_sraplan
#SBATCH -p short
#SBATCH -c 2
#SBATCH --mem=8G
#SBATCH -t 00:40:00
#SBATCH -o /bigdata/stajichlab/shared/projects/GmG_runs/logs/30_sra_plan.%j.log
# Real SRA/ENA reads for #42, step 1: the run table (ENA portal API, read_run, metadata only), the NCBI BioSample metadata of the fungal assemblies, and the accession list.
set -uo pipefail
source $HOME/projects/GmG/scripts/hpcc/env.sh
S=${READS_SCRIPTS:-$GMG_REPO/scripts}
REL=$GR/release/2026.10.4
R=$GR/reads_sra; mkdir -p $R; cd $R
Q='tax_tree(4751) AND library_strategy="WGS" AND library_source="GENOMIC" AND instrument_platform="ILLUMINA" AND library_layout="PAIRED"'
F=run_accession,sample_accession,secondary_sample_accession,study_accession,scientific_name,tax_id,center_name,instrument_model,library_layout,library_strategy,library_source,read_count,base_count,fastq_ftp,fastq_bytes
if [ ! -s ena_runs.tsv ]; then
  for try in 1 2 3 4 5; do
    curl -sS -m 600 -G "https://www.ebi.ac.uk/ena/portal/api/search" --data-urlencode "result=read_run" --data-urlencode "query=$Q" \
      --data-urlencode "fields=$F" --data-urlencode "limit=0" --data-urlencode "format=tsv" -o ena_runs.tmp && head -c 12 ena_runs.tmp | grep -q run_accession && mv ena_runs.tmp ena_runs.tsv && break
    echo "ENA run table: attempt $try failed"; sleep $((try * 60))
  done
fi
[ -s ena_runs.tsv ] || { echo "ENA run table could not be fetched"; exit 1; }
echo "ENA runs: $(($(wc -l < ena_runs.tsv) - 1))"
[ -s assembly_summary.txt ] || curl -sS -m 600 -o assembly_summary.txt https://ftp.ncbi.nlm.nih.gov/genomes/genbank/fungi/assembly_summary.txt
python3 $S/sra_plan.py ena_runs.tsv assembly_summary.txt $REL/meta_full.tsv $GMG_REPO/reports/label_audit/reference_exclusions.txt sra_runs.tsv ${N:-100} 5 25 4
wc -l sra_runs.tsv
