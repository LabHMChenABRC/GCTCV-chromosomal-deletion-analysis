# ============================================================
# Collinearity Visualization (COLLINEARITY RESULTS FROM MCSCANX)
# ============================================================
# Author: Puyam Tondonba Singh
# Date: 2025-08-18
# ============================================================

# ==== Load Libraries ====
library(data.table)
library(ggplot2)
library(dplyr)
library(stringr)
library(cowplot)
library(patchwork)

# ============================================================
# Function to Parse Collinearity and Plot
# ============================================================
plot_collinearity <- function(gff_file, collinearity_file, xlab, ylab, colors, point_size = 0.5) {
  
  # ---- Load GFF ----
  gff <- fread(gff_file, header = FALSE)
  colnames(gff) <- c("Chr", "Gene", "Start", "End")
  gff <- gff %>% mutate(Mid = (Start + End) / 2)
  
  # ---- Parse Collinearity File ----
  lines <- readLines(collinearity_file)
  blocks <- data.frame(Block = integer(), Gene1 = character(), Gene2 = character(), stringsAsFactors = FALSE)
  block_id <- NA
  
  for (ln in lines) {
    # Detect new alignment block
    if (grepl("^## Alignment", ln)) {
      block_id <- as.integer(str_match(ln, "## Alignment\\s+(\\d+):")[, 2])
    }
    # Detect gene pairs
    else if (grepl("Musa|Macma", ln)) {
      parts <- str_match(ln, "^[^:]+:\\s+(\\S+)\\s+(\\S+)")
      if (!is.na(parts[1, 2]) && !is.na(parts[1, 3])) {
        blocks <- rbind(blocks, data.frame(
          Block = block_id,
          Gene1 = parts[1, 2],
          Gene2 = parts[1, 3],
          stringsAsFactors = FALSE
        ))
      }
    }
  }
  
  # ---- Merge with GFF ----
  blocks <- blocks %>%
    left_join(gff %>% select(Gene, Chr, Mid), by = c("Gene1" = "Gene")) %>%
    rename(Chr1 = Chr, Pos1 = Mid) %>%
    left_join(gff %>% select(Gene, Chr, Mid), by = c("Gene2" = "Gene")) %>%
    rename(Chr2 = Chr, Pos2 = Mid)
  
  # ---- Plot ----
  p <- ggplot(blocks, aes(x = Pos2/1e+6, y = Pos1/1e+6, color = as.factor(Block))) +
    geom_point(size = point_size) +
    theme_cowplot() +
    labs(x = xlab, y = ylab, color = "Block") +
    scale_color_manual(values = colors) +
    theme(
      panel.grid = element_blank(),
      axis.line = element_line(color = "black"),
      axis.text = element_text(size = 8),
      axis.title = element_text(size = 8),
      legend.position = "right"
    )
  
  return(p)
}

# ============================================================
# Generate Plots
# ============================================================

# Plot 1: BAKSHII vs Pahang
p1 <- plot_collinearity(
  gff_file = "E:/tony/WGS/horti_baxi_bw10000/baxixiao protein fasta/MCScanX-1.0.0/chr02.2/ban/r plot/ban_FM05.gff",
  collinearity_file = "E:/tony/WGS/horti_baxi_bw10000/baxixiao protein fasta/MCScanX-1.0.0/chr02.2/ban/r plot/ban_FM05.collinearity",
  xlab = "'DH Pahang' chr05 (Mb)",
  ylab = "Li et al., chr05 (Mb)",
  colors = c("#3529d9", "red", "green", "orange", "purple", "black")
)

# Plot 2: BXH2 vs Pahang
p2 <- plot_collinearity(
  gff_file = "E:/tony/WGS/horti_baxi_bw10000/baxixiao protein fasta/MCScanX-1.0.0/chr02.2/chr05.2/bxh2_FM.gff",
  collinearity_file = "E:/tony/WGS/horti_baxi_bw10000/baxixiao protein fasta/MCScanX-1.0.0/chr02.2/chr05.2/bxh2_FM.collinearity",
  xlab = "'DH Pahang' chr05 FM deletion region (Mb)",
  ylab = "Huang et al., 2023 chr05.2 FM deletion region (Mb)",
  colors = c("#3529d9", "red", "green", "orange", "purple")
)

# ============================================================
# Combine with Patchwork
# ============================================================
combined_plot <- p1 | p2

# Save
ggsave("Combined_chr05_collinearity.pdf", plot = combined_plot, width = 10, height = 8, units = "cm", paper = "a4")
