# ============================================================
# GO Enrichment 
# Author: Puyam Tondonba Singh
# Date: 2025-08-25
#
# This script performs GO enrichment analysis of upregulated genes 
# (from TC5 vs PC comparison) using clusterProfiler. It builds a 
# background gene set from expression counts, reduces redundancy 
# among GO terms, filters for biologically relevant categories, 
#

# ============================================================

# --- Libraries ---
library(data.table)
library(dplyr)
library(clusterProfiler)
library(org.Macuminata.eg.db)
library(ggplot2)

# --- Input paths ---
file_path <- "E:/Visualization/2F"
setwd(file_path)

# --- Build background gene list from expression matrix ---
cnt_files <- list.files(pattern = ".genes.results$", full.names = TRUE)
names(cnt_files) <- sub(".genes.results", "", basename(cnt_files))

cnt_mat <- rbindlist(
  lapply(cnt_files, fread, select = c("gene_id", "expected_count")),
  idcol = "sample"
) %>%
  dcast(gene_id ~ sample, value.var = "expected_count") %>%
  as.data.frame()

rownames(cnt_mat) <- cnt_mat$gene_id
cnt_mat <- round(as.matrix(cnt_mat[ , -1]))
background_gene_list <- rownames(cnt_mat)

# --- Load DEG file ---
deg_df <- fread("E:/Visualization/2F/TC5_vs_PC_DEG.tsv")

# --- Gene list for enrichment (only upregulated) ---
gene_list_clean <- deg_df %>%
  filter(Regulation == "Up") %>%
  pull(Gene_ID)

# --- Run GO enrichment ---
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
  
  for (j in (i+1):nrow(ego_df)) {
    if (used[j]) next
    genes_j <- ego_df$GeneList[[j]]
    similarity <- length(intersect(genes_i, genes_j)) / length(union(genes_i, genes_j))
    if (similarity >= 0.9) {
      group <- c(group, j)
      used[j] <- TRUE
    }
  }
  
  # keep custom terms if present, otherwise lowest padj
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

# --- Post-filtering ---
ego_df <- ego_df %>%
  filter(Count > 3) %>%
  mutate(
    GeneRatioNum   = sapply(strsplit(GeneRatio, "/"), function(x) as.numeric(x[1]) / as.numeric(x[2])),
    BgRatioNum     = sapply(strsplit(BgRatio, "/"), function(x) as.numeric(x[1]) / as.numeric(x[2])),
    log2_enrichment = log2(GeneRatioNum / BgRatioNum),
    GeneRatioPercent = GeneRatioNum * 100,
    log10_padj = -log10(p.adjust)
  )

# --- Manual exclusion of redundant terms ---
exclude_terms <- c(
  "regulation of cellular ketone metabolic process",
  "response to water deprivation",
  "defense response to bacterium",
  "response to cold",
  "regulation of response to stress"
)
filtered_df <- ego_df %>% filter(!Description %in% exclude_terms)

# --- Plot enrichment bubble plot ---
TC5_upregulated <- ggplot(filtered_df, aes(GeneRatioPercent, y = reorder(Description, GeneRatioPercent))) +
  geom_point(aes(size = GeneRatioPercent, color = log10_padj)) +
  scale_color_gradientn(
    colours = c("#2c3793", "#a160bf", "#fc646f"),
    name = expression(-log[10]~"(adjusted p-value)"),
    limits = c(0, 5),
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
    axis.text.y   = element_text(size = 10),
    axis.text.x   = element_blank(),
    axis.ticks.x  = element_blank(),
    panel.grid    = element_blank(),
    panel.border  = element_rect(color = "black", fill = NA, linewidth = 0.5),
    legend.position = "right"
  ) +
  labs(x = NULL, y = NULL, title = NULL)

ggsave("TC5_upregulated_generatioPNAS1.pdf",
       plot = TC5_upregulated,
       device = "pdf", height = 8, width = 16, units = "cm")
