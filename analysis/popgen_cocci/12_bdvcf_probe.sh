#!/bin/bash
# probe the Bd SNP VCF: samples, sites, metadata
module load bcftools/1.19
D=/bigdata/stajichlab/shared/projects/Population_Genomics/B_dendrobatidis
V=$D/vcf/Bd_2.All.SNP.combined_selected.vcf.gz
echo "samples: $(bcftools query -l $V | wc -l)"; bcftools query -l $V | head -8 | paste -sd' '
echo "sites: $(bcftools index -n $V)"
bcftools view -h $V | grep -E "^##(contig|FILTER)" | head -5
bcftools view -H $V | head -2 | cut -f1-12 | cut -c1-200
echo ---
head -3 $D/samples.csv; wc -l $D/samples.csv; head -3 $D/samples_metadata.csv | cut -c1-300; head -12 $D/population_sets.yaml
