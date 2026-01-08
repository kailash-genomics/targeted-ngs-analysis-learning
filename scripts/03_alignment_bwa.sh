#!/bin/bash

# Script: 03_alignment_bwa.sh
# Purpose: Align trimmed reads to the human reference genome using BWA-MEM

bwa mem -t 8 reference/human_reference.fa \
trimmed_fastq/SAMPLE_R1.trimmed.fastq.gz \
trimmed_fastq/SAMPLE_R2.trimmed.fastq.gz \
> alignment/SAMPLE.sam
