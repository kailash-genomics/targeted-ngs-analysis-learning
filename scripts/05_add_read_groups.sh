#!/bin/bash

# Script: 05_add_read_groups.sh
# Purpose: Add or replace read group information required for GATK

picard AddOrReplaceReadGroups \
I=alignment/SAMPLE.sorted.bam \
O=alignment/SAMPLE.RG.sorted.bam \
RGID=SAMPLE \
RGLB=TargetedPanel \
RGPL=ILLUMINA \
RGPU=PU1 \
RGSM=SAMPLE

# Index BAM after adding read groups
samtools index alignment/SAMPLE.RG.sorted.bam
