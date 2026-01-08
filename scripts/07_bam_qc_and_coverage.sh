#!/bin/bash

# Script: 07_bam_qc_and_coverage.sh
# Purpose: Perform BAM quality control and coverage analysis

# Basic alignment statistics
samtools flagstat alignment/SAMPLE.dedup.bam > metrics/SAMPLE_flagstat.txt

# Chromosome-level read distribution
samtools idxstats alignment/SAMPLE.dedup.bam > metrics/SAMPLE_idxstats.txt

# BAM integrity check
samtools quickcheck alignment/SAMPLE.dedup.bam

# Coverage analysis over targeted regions
samtools depth \
-b panel/target_regions.bed \
alignment/SAMPLE.dedup.bam \
> metrics/SAMPLE_coverage.txt
