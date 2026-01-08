# End-to-End Targeted NGS Panel Analysis (Learning Project)

## Overview
This repository contains a step-by-step implementation of a targeted DNA sequencing data analysis pipeline, starting from raw FASTQ files to annotated variants.

## Project Summary

This project demonstrates an end-to-end targeted clinical NGS data analysis workflow implemented as a documentation-first learning exercise. The pipeline follows SOP-style logic commonly used in somatic cancer testing laboratories, covering raw data quality control, read preprocessing, alignment, post-alignment processing, variant calling, annotation, clinical-style filtering, and final variant table generation. The project emphasizes workflow structure, quality control checkpoints, and responsible interpretation boundaries rather than automated clinical reporting.

## Disclaimer
This repository is intended for educational and skill-development purposes only.  
The pipeline and scripts demonstrated here are not validated for clinical diagnostics and should not be used for patient-related decision making.

## Dataset

This project is based on publicly available human targeted sequencing datasets used for training and learning purposes.

- Sequencing type: Targeted DNA sequencing (gene panel)
- Data format: Paired-end FASTQ
- Analysis type: Tumor-only somatic variant analysis
- Reference genome: hg19 / GRCh38
- Target regions: Cancer-related genes defined by a BED file

Raw sequencing data files are not included in this repository due to size and privacy considerations.

## Pipeline Overview
## Pipeline Workflow Summary

| Phase | Step | Description | Tools |
|------|------|------------|-------|
| Pre-processing | Raw data QC | Assessment of sequencing read quality | FastQC |
| Pre-processing | Read trimming | Adapter removal and quality filtering | fastp |
| Alignment | Read alignment | Mapping reads to human reference genome | BWA-MEM |
| BAM Processing | Post-alignment processing | SAM to BAM conversion, sorting, indexing | SAMtools |
| BAM Processing | Read groups | Add read group information | Picard |
| BAM Processing | Duplicate marking | Identify PCR duplicates | Picard |
| QC | BAM QC & coverage | Alignment metrics and target coverage | SAMtools |
| Variant Calling | Somatic calling | Tumor-only variant calling | GATK Mutect2 |
| Variant Processing | Variant filtering | Remove low-confidence variants | GATK, bcftools |
| Annotation | Functional annotation | Gene and effect annotation | snpEff |
| Clinical Filtering | Clinical thresholds | Apply DP, AF, PASS, panel BED filters | bcftools |
| Reporting | Final output | Generate variant table for review | bcftools |


The targeted NGS analysis pipeline implemented in this project follows a standard clinical-style workflow designed for tumor-only variant detection.

### Workflow Steps

1. **Raw Data Quality Control**
   - Assessment of sequencing read quality
   - Identification of low-quality bases and adapter contamination

2. **Read Trimming and Filtering**
   - Removal of adapters and low-quality bases
   - Retention of high-quality reads for downstream analysis

3. **Alignment to Reference Genome**
   - Mapping of filtered reads to the human reference genome (hg19 / GRCh38)

4. **Post-alignment Processing**
   - Conversion of SAM to BAM format
   - Sorting and indexing of aligned reads
   - Addition of read group information
   - Marking of PCR duplicates

5. **BAM Quality Control and Coverage Analysis**
   - Evaluation of alignment metrics
   - Assessment of coverage across targeted regions using BED files

6. **Variant Calling**
   - Identification of somatic single nucleotide variants (SNVs) and small indels
   - Tumor-only variant calling strategy

7. **Variant Filtering and Normalization**
   - Removal of low-confidence variant calls
   - Normalization and left-alignment of variants

8. **Variant Annotation**
   - Functional annotation of variants to assess potential biological impact

9. **Clinical-style Variant Filtering**
   - Filtering based on read depth, allele frequency, and quality metrics
   - Restriction of variants to target gene panel regions

10. **Variant Interpretation (Manual)**
    - Interpretation using public knowledge bases and literature
    - Classification into pathogenic, likely pathogenic, VUS, or benign categories
   
## Tools Used

The following open-source bioinformatics tools are used throughout the targeted NGS analysis pipeline:

- **FastQC** – Quality assessment of raw sequencing reads  
- **fastp** – Adapter trimming and quality filtering  
- **BWA-MEM** – Alignment of reads to the human reference genome  
- **SAMtools** – BAM file manipulation, sorting, indexing, and QC  
- **Picard** – Read group assignment and duplicate marking  
- **GATK (Mutect2)** – Somatic variant calling in tumor-only mode  
- **bcftools** – Variant filtering, normalization, and querying  
- **snpEff** – Functional annotation of genetic variants  
- **IGV** – Visualization of aligned reads and variants (manual review)

All tools listed are widely used in research and clinical genomics laboratories.

## Variant Interpretation

Variant interpretation is performed as a manual review step and is not fully automated in this project.

Filtered variants are reviewed using publicly available knowledge bases and literature to assess their potential clinical relevance. Variants are broadly categorized into:

- Pathogenic
- Likely Pathogenic
- Variant of Uncertain Significance (VUS)
- Benign / Likely Benign

Databases commonly referenced during interpretation include ClinVar, COSMIC, OncoKB, and relevant scientific literature.

## Limitations

This project is implemented as a learning and training exercise and has several limitations:

- Tumor-only analysis without a matched normal sample
- No copy number variation (CNV) or structural variant analysis
- Annotation limited to selected tools and databases
- Filtering thresholds are representative but may vary across laboratories
- The pipeline is not validated for clinical diagnostics

These limitations are acknowledged to ensure responsible and accurate representation of the project scope.



