#!/bin/bash

# Script: 08_mutect2_variant_calling.sh
# Purpose: Somatic variant calling using GATK Mutect2 (tumor-only mode)

gatk Mutect2 \
-R reference/human_reference.fa \
-I alignment/SAMPLE.dedup.bam \
-tumor SAMPLE \
-L panel/target_regions.bed \
--disable-read-filter MateOnSameContigOrNoMappedMateReadFilter \
-O variants/SAMPLE.unfiltered.vcf.gz
