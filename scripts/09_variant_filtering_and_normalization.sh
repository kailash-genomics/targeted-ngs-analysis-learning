#!/bin/bash

# Script: 09_variant_filtering_and_normalization.sh
# Purpose: Filter Mutect2 variant calls and normalize variants

# Filter Mutect2 calls
gatk FilterMutectCalls \
-R reference/human_reference.fa \
-V variants/SAMPLE.unfiltered.vcf.gz \
-O variants/SAMPLE.filtered.vcf.gz

# Normalize and left-align variants
bcftools norm \
-f reference/human_reference.fa \
-m -both \
variants/SAMPLE.filtered.vcf.gz \
-Oz -o variants/SAMPLE.filtered.norm.vcf.gz

# Index normalized VCF
bcftools index variants/SAMPLE.filtered.norm.vcf.gz
