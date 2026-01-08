#!/bin/bash

# Script: 12_variant_table_generation.sh
# Purpose: Generate a tab-delimited variant table for review

bcftools query \
-f '%CHROM\t%POS\t%REF\t%ALT\t%INFO/ANN\t%FORMAT/AF\t%INFO/DP\n' \
variants/SAMPLE.clinical.panel.vcf.gz \
> results/SAMPLE_variants.tsv
