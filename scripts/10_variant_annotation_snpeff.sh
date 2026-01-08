#!/bin/bash

# Script: 10_variant_annotation_snpeff.sh
# Purpose: Annotate filtered variants using snpEff

snpEff \
-canon \
-hgvs \
-stats snpEff_SAMPLE.html \
hg19 \
variants/SAMPLE.filtered.norm.vcf.gz \
> variants/SAMPLE.annotated.vcf
