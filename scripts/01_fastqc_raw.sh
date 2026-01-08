#!/bin/bash

# Script: 01_fastqc_raw.sh
# Purpose: Quality control of raw FASTQ files using FastQC

fastqc \
raw_fastq/SAMPLE_R1.fastq.gz \
raw_fastq/SAMPLE_R2.fastq.gz \
-o fastqc_raw/
