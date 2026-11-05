#!/bin/bash
# dump the Bd SNP VCF as a dosage matrix (rows = biallelic SNPs, cols = samples; -1 = missing)
set -euo pipefail
module load bcftools/1.19
D=/bigdata/stajichlab/shared/projects/Population_Genomics/B_dendrobatidis
V=$D/vcf/Bd_2.All.SNP.combined_selected.vcf.gz
O=/bigdata/stajichlab/jstajich/projects/GmG_popgen/popgen_bd
bcftools query -l $V > $O/vcf_samples.txt
bcftools view -m2 -M2 -v snps -i 'F_MISSING<0.5' $V \
 | bcftools query -f '[%GT\t]\n' \
 | sed -e 's/0[\/|]0/0/g' -e 's/1[\/|]1/2/g' -e 's/[01][\/|][01]/1/g' -e 's/\.[\/|]\./-1/g' -e 's/\t$//' > $O/vcf_dosage.tsv
echo "rows: $(wc -l < $O/vcf_dosage.tsv)  samples: $(wc -l < $O/vcf_samples.txt)"
