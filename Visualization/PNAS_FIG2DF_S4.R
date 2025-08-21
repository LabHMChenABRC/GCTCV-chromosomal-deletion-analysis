# ============================================================
#DEG analysis PC vs TC5

# This visualizes up- and down-regulated genes with a barplot,
# conducts GO enrichment analysis on upregulated genes, and generates a heatmap
# of mean TPM values for genes located in a specified deletion region.

# ============================================================
# Author: Puyam Tondonba Singh
# Date: 2025-08-18
# ============================================================


library(data.table)
library(ggplot2)
library(dplyr)
library(clusterProfiler)
library(ComplexHeatmap)
library(circlize)

# Load DEG TSV
deg_df <- fread("E:/tony/HydroTC5_RNAseq/RSEM/MIX/TC5_vs_PC_DEG.tsv")

# DEG barplot
p <- ggplot(deg_df, aes(x = Regulation, fill = Regulation)) +
  geom_bar(width = 0.7, color = "black") +
  geom_text(stat = "count", aes(label = ..count..), 
            vjust = -0.5, size = 4, fontface = "bold") +
  scale_fill_manual(values = c("Up" = "#f78d8d", "Down" = "#2bbac2")) +
  labs(title = "Differential Expression: TC5 vs PC", x = "Regulation", y = "Gene Count") +
  theme_minimal(base_size = 10) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    axis.text.x = element_text(angle = 45, hjust = 1, size = 12),
    axis.text.y = element_text(size = 10),
    axis.title.x = element_text(size = 10, face = "bold"),
    axis.title.y = element_text(size = 10, face = "bold"),
    legend.position = "top",
    legend.title = element_blank(),
    legend.text = element_text(size = 10),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_line(color = "black"),
    axis.ticks = element_line(color = "black", size = 0.5)
  )
#Fig. S4B
ggsave("DEG_count_barplot.pdf", plot = p, width = 8, height = 12, units = "cm")



setwd("E:/tony/HydroTC5_RNAseq/RSEM/MIX")

file_path <- "E:/tony/VIGS/Transcriptome/Annotation_v4/Musa_acuminata_pahang_v4_go.txt"
# Check if file exists
file.exists(file_path)  # should return TRUE

# Read the file
gene2go <- read.delim(file_path, header = FALSE, stringsAsFactors = FALSE)
colnames(gene2go) <- c("GID", "GO")
gene2go$GID <- sub("\\.1$", "", gene2go$GID)
head(gene2go)


# Add default evidence code ("IEA" = Inferred from Electronic Annotation)
gene2go_df <- gene2go %>%
  mutate(EVIDENCE = "IEA")

# --------------------------------------
# Step 3: Build the OrgDb package
# --------------------------------------
library(AnnotationForge)


makeOrgPackage(
  go = gene2go_df,
  version = "0.1",
  maintainer = "tony <tonpuyam@gmail.com>",
  author = "Tony",
  outputDir = "E:/tony/HydroTC5_RNAseq/RSEM/MIX",
  tax_id = "4641",
  genus = "Musa",
  species = "acuminata",
  goTable = "go"
)

# --------

# Load custom OrgDb
library(org.Macuminata.eg.db)

# GO enrichment analysis
gene_list_clean <- deg_df %>% filter(Regulation == "Up") %>% pull(Gene_ID)
background_gene_list <- rownames(as.data.frame(cnt_mat))

ego_result <- enrichGO(
  gene          = gene_list_clean,
  universe      = background_gene_list,
  OrgDb         = org.Macuminata.eg.db,
  keyType       = "GID",
  ont           = "BP",
  pvalueCutoff  = 0.05,
  qvalueCutoff  = 0.05,
  pAdjustMethod = "BH"
)

ego_df <- as.data.frame(ego_result)
ego_df$GeneList <- strsplit(as.character(ego_df$geneID), "/")
custom_keep <- c("response to salicylic acid")

used <- rep(FALSE, nrow(ego_df))
filtered_terms <- list()

for (i in seq_len(nrow(ego_df))) {
  if (used[i]) next
  genes_i <- ego_df$GeneList[[i]]
  group <- i
  if (i < nrow(ego_df)) {
    for (j in (i+1):nrow(ego_df)) {
      if (used[j]) next
      genes_j <- ego_df$GeneList[[j]]
      similarity <- length(intersect(genes_i, genes_j)) / length(union(genes_i, genes_j))
      if (similarity >= 0.9) {
        group <- c(group, j)
        used[j] <- TRUE
      }
    }
  }
  group_terms <- ego_df$Description[group]
  keep_index <- if (any(group_terms %in% custom_keep)) {
    group[which(group_terms %in% custom_keep)[1]]
  } else {
    group[which.min(ego_df$p.adjust[group])]
  }
  filtered_terms[[length(filtered_terms) + 1]] <- ego_df[keep_index, ]
  used[i] <- TRUE
}

ego_df <- do.call(rbind, filtered_terms)

ego_df <- ego_df %>%
  filter(Count > 3) %>%
  mutate(
    GeneRatioNum = sapply(strsplit(as.character(GeneRatio), "/"), function(x) as.numeric(x[1]) / as.numeric(x[2])),
    BgRatioNum = sapply(strsplit(as.character(BgRatio), "/"), function(x) as.numeric(x[1]) / as.numeric(x[2])),
    log2_enrichment = log2(GeneRatioNum / BgRatioNum)
  )


#filter out similar GO terms manually
exclude_terms <- c(
  "regulation of cellular ketone metabolic process",
  "response to water deprivation",
  "defense response to bacterium",
  "response to cold",
  "regulation of response to stress"
)
filtered_df <- ego_df[!ego_df$Description %in% exclude_terms, ]

filtered_df <- filtered_df %>%
  mutate(
    GeneRatioPercent = GeneRatioNum * 100,
    log10_padj = -log10(p.adjust)
  )

TC5_upregulated <- ggplot(filtered_df, aes(GeneRatioPercent, y = reorder(Description, GeneRatioPercent))) +
  geom_point(aes(size = GeneRatioPercent, color = log10_padj)) +
  scale_color_gradientn(
    colours = c("#2c3793", "#a160bf", "#fc646f"),
    name = expression(-log[10]~"(adjusted p-value)"),
    limits = c(-0.2, 5.2),
    breaks = 0:5
  ) +
  scale_size_continuous(
    name = "Gene ratio (%)",
    range = c(4, 8),
    breaks = c(10, 13, 16),
    limits = c(10, 16)
  ) +
  theme_minimal(base_size = 10) +
  theme(
    axis.text.y = element_text(size = 10),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    panel.grid = element_blank(),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5),
    legend.position = "right"
  ) +
  labs(x = NULL, y = NULL, title = NULL)
#Fig. 2F
ggsave("TC5_upregulated_generatioPNAS1.pdf", plot = TC5_upregulated, paper = "A4", unit = "cm", height = 8, width = 16)



#####heatmap for deletion region

# ===== Define gene list =====
Deletion_region <- c("Macma4_04_g41380", "Macma4_04_g41390", "Macma4_04_g41400", 
                     "Macma4_04_g41410", "Macma4_04_g41420", "Macma4_04_g41430", "Macma4_04_g41440")

# ===== Load all files =====
all_files <- list.files(pattern = ".genes.results$")

# ===== Read TPMs =====
read_and_select <- function(filename) {
  fread(filename, sep = "\t")[, .(gene_id, TPM)][, TPM := as.integer(TPM)]
}
all_data.list <- lapply(all_files, read_and_select)
names(all_data.list) <- sub(".genes.results", "", basename(all_files))
all_data <- rbindlist(all_data.list, idcol = "sample")

# ===== Create wide-format TPM matrix =====
TPM.dt <- dcast(all_data, gene_id ~ sample, value.var = 'TPM')

# ===== Filter genes in the deletion region =====
TPM_deletion.dt <- TPM.dt %>% filter(gene_id %in% Deletion_region)


library(dplyr)


TPM_deletion.dt[, PC_mean := rowMeans(.SD, na.rm = TRUE), .SDcols = patterns("^PC_")]
TPM_deletion.dt[, TC5_mean := rowMeans(.SD, na.rm = TRUE), .SDcols = patterns("^TC5_")]

# ===== Create expression matrix with gene IDs as row names =====
mean_expr <- TPM_deletion.dt[, .(gene_id, PC_mean, TC5_mean)]
rownames(mean_expr) <- mean_expr$gene_id
mean_expr$gene_id <- NULL
heatmap_matrix <- as.matrix(mean_expr)

# ===== Handle NaN and Inf =====
heatmap_matrix[is.nan(heatmap_matrix)] <- 0
heatmap_matrix[is.infinite(heatmap_matrix)] <- max(heatmap_matrix[!is.infinite(heatmap_matrix)], na.rm = TRUE)

# ===== Plot Heatmap with gene names as row labels =====
#Fig. 2D
pdf("Deletion_Region_PC_TC5_Mean_Heatmap_with_Genes.pdf", width = 2, height = 3, paper = "a4")  # 6cm x 10cm in inches

Heatmap(
  heatmap_matrix,
  name = "Mean TPM",
  col = colorRamp2(c(0, max(heatmap_matrix)), c("#fae1e3", "#f20f0f")),
  cluster_rows = FALSE,
  cluster_columns = FALSE,
  row_names_gp = grid::gpar(fontsize = 10),     # gene names shown here with font size 10
  column_names_gp = grid::gpar(fontsize = 10)
)

dev.off()
shell.exec("Deletion_Region_PC_TC5_Mean_Heatmap_with_Genes.pdf")
