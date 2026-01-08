#!/bin/bash

# Script: 04_post_alignment_processing.sh
# Purpose: Convert SAM to BAM, sort and index aligned reads

# Convert SAM to BAM
samtools view -Sb alignment/SAMPLE.sam > alignment/SAMPLE.bam

# Sort BAM file
samtools sort alignment/SAMPLE.bam -o alignment/SAMPLE.sorted.bam

# Index sorted BAM
samtools index alignment/SAMPLE.sorted.bam
