############################################################
# 04_deseq2_analysis.R
# Differential Expression Analysis of PRMT5 inhibitor 
# (GSK3326595) vs Vehicle in MYC-driven mouse HCC
############################################################

# Load libraries
library(tximport)
library(DESeq2)
library(tidyverse)
library(AnnotationDbi)
library(org.Mm.eg.db)
library(pheatmap)
library(ggplot2)
library(gprofiler2)

# Set working directory
setwd("~/rnaseq_project")

############################################################
# 1. Sample metadata
############################################################

sample_info <- data.frame(
  sample = c("SRR12231697", "SRR12231698", "SRR12231699",
             "SRR12231700", "SRR12231701", "SRR12231702"),
  condition = c("vehicle", "vehicle", "vehicle",
                "GSK3326595", "GSK3326595", "GSK3326595")
)

sample_info$condition <- factor(sample_info$condition, 
                                levels = c("vehicle", "GSK3326595"))

write.csv(sample_info, "counts/sample_info.csv", row.names = FALSE)

############################################################
# 2. Import Salmon quantifications
############################################################

files <- file.path("aligned", paste0(sample_info$sample, "_salmon"), "quant.sf")
names(files) <- sample_info$sample

# Check files exist
stopifnot(all(file.exists(files)))

# Import at transcript level
txi <- tximport(files, type = "salmon", txOut = TRUE)

############################################################
# 3. DESeq2 analysis
############################################################

dds <- DESeqDataSetFromTximport(txi,
                                colData = sample_info,
                                design = ~ condition)

# Pre-filter low counts
keep <- rowSums(counts(dds)) >= 10
dds <- dds[keep, ]

# Run DESeq2
dds <- DESeq(dds)

# Results
res <- results(dds, contrast = c("condition", "GSK3326595", "vehicle"))
res_ordered <- res[order(res$padj), ]

# Save full results
write.csv(as.data.frame(res_ordered), "results/deseq2_results_transcript.csv")

############################################################
# 4. Map to gene symbols
############################################################

transcripts <- rownames(res_ordered)
transcripts_clean <- gsub("\\..*", "", transcripts)

gene_symbols <- mapIds(org.Mm.eg.db,
                       keys = transcripts_clean,
                       column = "SYMBOL",
                       keytype = "ENSEMBLTRANS",
                       multiVals = "first")

res_df <- as.data.frame(res_ordered)
res_df$transcript <- rownames(res_df)
res_df$gene_symbol <- gene_symbols[transcripts_clean]

# Remove transcripts without gene symbol
res_df <- res_df[!is.na(res_df$gene_symbol), ]

# Keep best transcript per gene
res_gene <- res_df %>%
  group_by(gene_symbol) %>%
  slice_min(order_by = padj, n = 1, with_ties = FALSE) %>%
  ungroup() %>%
  as.data.frame()

rownames(res_gene) <- res_gene$gene_symbol

# Save gene-level results
write.csv(res_gene, "results/deseq2_results_gene_level.csv", row.names = FALSE)

res_gene_sig <- subset(res_gene, padj < 0.05)
write.csv(res_gene_sig, "results/deseq2_significant_genes.csv", row.names = FALSE)

cat("Number of significant genes (padj < 0.05):", nrow(res_gene_sig), "\n")

############################################################
# 5. PCA Plot
############################################################

vsd <- vst(dds, blind = FALSE)

pca_data <- plotPCA(vsd, intgroup = "condition", returnData = TRUE)
percentVar <- round(100 * attr(pca_data, "percentVar"))

ggplot(pca_data, aes(x = PC1, y = PC2, color = condition)) +
  geom_point(size = 4) +
  xlab(paste0("PC1: ", percentVar[1], "% variance")) +
  ylab(paste0("PC2: ", percentVar[2], "% variance")) +
  ggtitle("PCA - Vehicle vs GSK3326595") +
  theme_bw()

ggsave("results/PCA_plot.png", width = 7, height = 5)

############################################################
# 6. Volcano Plot
############################################################

res_df_plot <- as.data.frame(res)
res_df_plot$gene <- rownames(res_df_plot)
res_df_plot$significant <- ifelse(res_df_plot$padj < 0.05 & abs(res_df_plot$log2FoldChange) > 1, 
                                  "Significant", "Not significant")

ggplot(res_df_plot, aes(x = log2FoldChange, y = -log10(padj), color = significant)) +
  geom_point(alpha = 0.6) +
  scale_color_manual(values = c("grey70", "red3")) +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
  theme_bw() +
  ggtitle("Volcano Plot - GSK3326595 vs Vehicle")

ggsave("results/volcano_plot.png", width = 8, height = 6)

############################################################
# 7. Pathway Enrichment (g:Profiler)
############################################################

sig_genes <- res_gene_sig$gene_symbol

gost_res <- gost(query = sig_genes,
                 organism = "mmusculus",
                 significant = TRUE,
                 correction_method = "fdr")

# Clean and save results
gost_df <- gost_res$result %>%
  mutate(across(where(is.list), ~sapply(., function(x) paste(x, collapse = ";"))))

write.csv(gost_df, "results/gprofiler_enrichment.csv", row.names = FALSE)

# Simple bar plot of top pathways
top_terms <- gost_df %>%
  filter(significant == TRUE) %>%
  arrange(p_value) %>%
  head(15)

ggplot(top_terms, aes(x = reorder(term_name, -log10(p_value)), y = -log10(p_value))) +
  geom_col(fill = "steelblue") +
  coord_flip() +
  labs(title = "Top Enriched Pathways (g:Profiler)",
       x = "Pathway",
       y = "-log10(p-value)") +
  theme_bw()

ggsave("results/gprofiler_top_pathways.png", width = 10, height = 7)

############################################################
# Done
############################################################

cat("\nAnalysis complete. Results saved in the results/ folder.\n")
