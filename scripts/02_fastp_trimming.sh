#!/bin/bash

# Script: 02_fastp_trimming.sh
# Purpose: Adapter trimming and quality filtering of raw FASTQ files

fastp \
-i raw_fastq/SAMPLE_R1.fastq.gz \
-I raw_fastq/SAMPLE_R2.fastq.gz \
-o trimmed_fastq/SAMPLE_R1.trimmed.fastq.gz \
-O trimmed_fastq/SAMPLE_R2.trimmed.fastq.gz \
--detect_adapter_for_pe \
--qualified_quality_phred 20 \
--length_required 50 \
--html fastqc_trimmed/SAMPLE_fastp.html \
--json fastqc_trimmed/SAMPLE_fastp.json
