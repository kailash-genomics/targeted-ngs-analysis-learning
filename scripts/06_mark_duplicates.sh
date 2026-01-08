#!/bin/bash

# Script: 06_mark_duplicates.sh
# Purpose: Mark PCR duplicates in aligned BAM files

picard MarkDuplicates \
I=alignment/SAMPLE.RG.sorted.bam \
O=alignment/SAMPLE.dedup.bam \
M=metrics/SAMPLE_duplication_metrics.txt \
CREATE_INDEX=true
