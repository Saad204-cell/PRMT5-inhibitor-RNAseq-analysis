# RNA-seq Analysis of PRMT5 Inhibitor (GSK3326595) in MYC-driven Mouse HCC

**Portfolio / Case Study Project**

Re-analysis of publicly available bulk RNA-seq data from a PRMT5 inhibitor study in a MYC-driven hepatocellular carcinoma (HCC) mouse model.

**Original Dataset:** [GSE154514](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE154514)  
**Original Paper:** Luo et al., *Hepatology* (2021) – “Myelocytomatosis-Protein Arginine N-Methyltransferase 5 Axis Defines the Tumorigenesis and Immune Response in Hepatocellular Carcinoma”

---

## Project Goal

Treat a real public RNA-seq dataset as a freelance client project and deliver a complete, reproducible analysis pipeline from raw data to biological interpretation.

---

## Experimental Design

| Group              | Samples                          | Treatment      |
|--------------------|----------------------------------|----------------|
| Vehicle (Control)  | SRR12231697, SRR12231698, SRR12231699 | Vehicle       |
| Treated            | SRR12231700, SRR12231701, SRR12231702 | GSK3326595    |

- Organism: *Mus musculus*
- Tissue: Liver tumors from MYC-driven HCC model
- Library: Paired-end total RNA-seq

---

## Analysis Pipeline

1. **Quality Control** – FastQC + MultiQC
2. **Trimming** – fastp (adapter & quality trimming)
3. **Quantification** – Salmon (selective alignment)
4. **Differential Expression** – tximport + DESeq2
5. **Visualization** – PCA, Volcano plot, Heatmap
6. **Pathway Enrichment** – g:Profiler (GO + KEGG)

---

## Key Results

- Clear separation between Vehicle and GSK3326595 samples in PCA
- 84 significantly differentially expressed genes (padj < 0.05)
- Pathway enrichment strongly pointed to **immune-related processes**, including:
  - Leukocyte migration
  - Th1/Th2 cell differentiation
  - IL-17 signaling
  - Cell adhesion
  - Inflammatory pathways

These findings are consistent with the original publication, which reported that PRMT5 inhibition enhances antitumor immunity (increased MHC class II expression and lymphocyte infiltration).

---

## Folder Structure
rnaseq_project/
├── data/
│   ├── raw/              # Original FASTQ files
│   └── trimmed/          # Cleaned FASTQ files
├── aligned/              # Salmon quantification output
├── counts/               # Sample metadata & count matrices
├── qc/                   # FastQC, MultiQC, fastp reports
├── results/              # DESeq2 results, plots, enrichment
├── scripts/              # Analysis scripts
└── README.md

---

## How to Reproduce

1. Download the six SRA runs listed above
2. Follow the pipeline steps described in the `scripts/` folder
3. All major outputs are saved in the `results/` directory

---

## Tools Used

- SRA Toolkit
- FastQC / MultiQC
- fastp
- Salmon
- R / Bioconductor (tximport, DESeq2, AnnotationDbi)
- gprofiler2

---

## Notes

This is an independent re-analysis performed for portfolio and learning purposes. It is not affiliated with the original authors. The goal was to demonstrate a complete, professional RNA-seq workflow on real public data.
