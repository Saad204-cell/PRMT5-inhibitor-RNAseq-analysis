#!/bin/bash
# Quantification with Salmon

for sample in SRR12231697 SRR12231698 SRR12231699 SRR12231700 SRR12231701 SRR12231702; do
    echo "Quantifying $sample ..."
    salmon quant \
        -i ref/gencode.vM36_salmon_index \
        -l A \
        -1 data/trimmed/${sample}_1.trimmed.fastq.gz \
        -2 data/trimmed/${sample}_2.trimmed.fastq.gz \
        -p 8 \
        --validateMappings \
        -o aligned/${sample}_salmon
done
