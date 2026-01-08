#!/bin/bash

# Script: 11_clinical_variant_filtering.sh
# Purpose: Apply clinical-style filtering to annotated variants

# Filter variants by quality, depth, and allele frequency
bcftools view \
-i 'FILTER="PASS" && INFO/DP>=100 && FORMAT/AF>=0.05' \
variants/SAMPLE.annotated.vcf \
-Oz -o variants/SAMPLE.clinical.pass.vcf.gz

# Restrict variants to targeted gene panel regions
bcftools view \
-R panel/target_regions.bed \
variants/SAMPLE.clinical.pass.vcf.gz \
-Oz -o variants/SAMPLE.clinical.panel.vcf.gz

# Index final clinical VCF
bcftools index variants/SAMPLE.clinical.panel.vcf.gz
