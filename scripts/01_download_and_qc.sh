#!/bin/bash
# Download SRA data, convert to FASTQ, and run FastQC + MultiQC

# (You already did this part, so this is just for documentation)
# prefetch SRR12231697 SRR12231698 SRR12231699 SRR12231700 SRR12231701 SRR12231702
# fasterq-dump --split-files ...

mkdir -p qc/fastqc
fastqc data/raw/*.fastq.gz -o qc/fastqc -t 8
multiqc qc/fastqc -o qc/
