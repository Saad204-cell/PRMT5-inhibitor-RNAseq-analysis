#!/bin/bash
# Adapter and quality trimming with fastp

mkdir -p data/trimmed

for sample in SRR12231697 SRR12231698 SRR12231699 SRR12231700 SRR12231701 SRR12231702; do
    echo "Trimming $sample ..."
    fastp \
        -i data/raw/${sample}_1.fastq.gz \
        -I data/raw/${sample}_2.fastq.gz \
        -o data/trimmed/${sample}_1.trimmed.fastq.gz \
        -O data/trimmed/${sample}_2.trimmed.fastq.gz \
        --detect_adapter_for_pe \
        --thread 8 \
        --html qc/${sample}_fastp.html \
        --json qc/${sample}_fastp.json
done
