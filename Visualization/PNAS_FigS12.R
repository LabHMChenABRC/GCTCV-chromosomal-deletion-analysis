
# This script summarizes counts FOR  'Nongke No. 1' VS 'Baxi', and generates barplots and chromosome-level plots.


# ============================================================
# Author: Puyam Tondonba Singh
# Date: 2025-08-18
# ============================================================

library(data.table)
library(dplyr)
library(ggplot2)
library(scales)

setwd("E:/Visualization/S12") 

# --- LOAD DEG TABLES FROM TSV ---
deg_27H <- fread("DEG_results/DEG_27H.tsv")

deg_27H <- deg_27H %>%
  filter(grepl('chr',Chromosome)) %>%
  filter(padj < 0.1 & abs(log2FoldChange) > 0.5) %>%
  mutate(
    Start_Mbp = Start / 1e6,
    color = case_when(
      log2FoldChange > 0.5 ~ "Up-regulated",
      log2FoldChange < -0.5 ~ "Down-regulated"
    )
  )%>%
  mutate(
    color = factor(color,c("Up-regulated","Down-regulated"))
  )

highlight <- data.frame(Chromosome = "chr05", xmin = 0.956, xmax = 7.03)

DEG_plot_27H <- ggplot(deg_27H, aes(x = Start_Mbp, y = log2FoldChange, color = color)) +
  labs(color=NULL)+
  geom_rect(data = highlight, aes(xmin = xmin, xmax = xmax, ymin = -Inf, ymax = Inf), inherit.aes = FALSE, fill = "gold", alpha = 0.2) +
  geom_point(alpha = 0.8, size = 0.5) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +
  facet_wrap(~ Chromosome, scales = "free_x", ncol = 4) +
  scale_color_manual(values = c("Up-regulated" = "red", "Down-regulated" = "blue")) +
  scale_y_continuous(name = expression(Log[2]~Fold~Change), breaks = c(-5,0,5)) +
  scale_x_continuous(name = "'DH−Pahang' genomic position (Mbp)", breaks = pretty_breaks(5)) +
  theme_bw(base_size = 10) +
  theme(panel.grid = element_blank(), 
        strip.text = element_text(size = 10, face = "bold"), 
        axis.text.x = element_text(angle = 0, hjust = 0.5),
        legend.text = element_text(size = 10),
        legend.position = "top")

PDF_FILE="FigS12.DEG_along_chromosome_27H.pdf"
ggsave(PDF_FILE, DEG_plot_27H, width = 11, height = 6)
shell.exec(PDF_FILE)