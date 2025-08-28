# ============================================================
# DEG analysis PC vs TC5 (DESeq2 + export DEG table)
# Author: Puyam Tondonba Singh
# Date: 2025-08-18
# ============================================================

library(data.table)
library(DESeq2)
library(dplyr)

# Set working directory
file_path <- "E:/tony/HydroTC5_RNAseq/RSEM/MIX"
setwd(file_path)

# Load and prepare count matrix
cnt_files <- list.files(pattern = ".genes.results$", full.names = TRUE)
names(cnt_files) <- sub(".genes.results", "", basename(cnt_files))
cnt_mat <- rbindlist(lapply(cnt_files, fread, select = c("gene_id", "expected_count")), idcol = "sample")
cnt_mat <- dcast(cnt_mat, gene_id ~ sample, value.var = "expected_count")
cnt_mat <- as.data.frame(cnt_mat)
rownames(cnt_mat) <- cnt_mat$gene_id
cnt_mat <- cnt_mat[, -1]
cnt_mat <- round(as.matrix(cnt_mat))

# Filter lowly expressed genes
cnt_mat <- cnt_mat[rowSums(cnt_mat > 5) >= ncol(cnt_mat) / 2, ]

# Build sample metadata
sample_info <- data.table(ID = colnames(cnt_mat))
sample_info[, condition := ifelse(grepl("^PC", ID), "PC", "TC5")]
sample_info[, batch := sub(".*M(\\d+)$", "R\\1", ID)]
sample_info[, Group := condition]
sample_info <- as.data.frame(sample_info)
rownames(sample_info) <- sample_info$ID

# Run DESeq2
dds <- DESeqDataSetFromMatrix(countData = cnt_mat, colData = sample_info, design = ~ batch + condition)
dds <- DESeq(dds, test = "LRT", reduced = ~ batch)
res <- results(dds, contrast = c("condition", "TC5", "PC"))
res <- na.omit(as.data.frame(res))
res$Gene_ID <- rownames(res)

# Filter DEGs
deg_df <- res %>%
  filter(padj < 0.1, abs(log2FoldChange) > 0.5) %>%
  mutate(Regulation = case_when(
    log2FoldChange > 0.5 ~ "Up",
    log2FoldChange < -0.5 ~ "Down"
  ))

# Export DEG results to TSV
fwrite(deg_df, file = "TC5_vs_PC_DEG.tsv", sep = "\t")

message("DEG analysis finished. Results saved to TC5_vs_PC_DEG.tsv")


