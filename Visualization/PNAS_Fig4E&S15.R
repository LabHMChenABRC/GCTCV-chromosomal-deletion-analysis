
# This script summarizes counts FOR  'Nongke No. 1' VS 'Baxi', and generates barplots and chromosome-level plots.


# ============================================================
# Author: Puyam Tondonba Singh
# Date: 2025-08-18
# ============================================================


library(data.table)
library(dplyr)
library(tidyr)
library(ggplot2)
library(cowplot)
library(scales)

# --- HELPER FUNCTION ---
safe_count <- function(obj, key) {
  val <- obj[[key]]
  if (is.null(val) || is.na(val)) return(0)
  return(val)
}

# --- LOAD DEG TABLES FROM TSV ---
deg_tables <- list()
deg_counts <- list()

time_points <- c("27 H", "52 H")
for (tp in time_points) {
  tsv_file <- file.path("DEG_results", paste0("DEG_", gsub(" ", "", tp), ".tsv"))
  if (!file.exists(tsv_file)) next
  dt <- fread(tsv_file)
  deg_tables[[tp]] <- dt
  
  deg_all <- dt %>% filter(color %in% c("Upregulated", "Downregulated"))
  deg_region <- dt %>% filter(color %in% c("Upregulated", "Downregulated") & Chromosome == "chr05" & Start >= 956425 & End <= 7030674)
  
  deg_counts[[tp]] <- list(
    all = table(deg_all$color),
    region = table(deg_region$color)
  )
}

# --- SUMMARY TABLE ---
data <- data.frame(
  Condition = rep(c("27h", "51h"), times = 2),
  Region = rep(c("Whole genome", "FM deletion"), each = 2),
  Upregulated = c(
    safe_count(deg_counts[["27 H"]]$all, "Upregulated"),
    safe_count(deg_counts[["52 H"]]$all, "Upregulated"),
    safe_count(deg_counts[["27 H"]]$region, "Upregulated"),
    safe_count(deg_counts[["52 H"]]$region, "Upregulated")
  ),
  Downregulated = c(
    safe_count(deg_counts[["27 H"]]$all, "Downregulated"),
    safe_count(deg_counts[["52 H"]]$all, "Downregulated"),
    safe_count(deg_counts[["27 H"]]$region, "Downregulated"),
    safe_count(deg_counts[["52 H"]]$region, "Downregulated")
  )
)

data_long <- data %>%
  pivot_longer(cols = c("Upregulated", "Downregulated"), names_to = "Regulation", values_to = "Count") %>%
  mutate(
    Region = factor(Region, levels = c("Whole genome", "FM deletion")),
    Regulation = factor(Regulation, levels = c("Upregulated", "Downregulated"))
  )

# --- BARPLOT ---
p <- ggplot(data_long, aes(x = Condition, y = Count, fill = Regulation)) +
  geom_bar(stat = "identity") +
  facet_wrap(~Region, nrow = 1, scales = "free_y") +
  scale_fill_manual(values = c("Downregulated" = "#00CFC1", "Upregulated" = "#F8766D")) +
  theme_minimal(base_size = 12) +
  theme_cowplot()

ggsave("DEG_plot.pdf", plot = p, width = 12, height = 8, units = "cm")

# --- CHROMOSOME PLOT (27H) ---
deg_27H <- deg_tables[["27 H"]] %>%
  filter(padj < 0.1 & abs(log2FoldChange) > 0.5) %>%
  mutate(
    Start_Mbp = Start / 1e6,
    color = case_when(
      log2FoldChange > 0.5 ~ "Up-regulated",
      log2FoldChange < -0.5 ~ "Down-regulated"
    )
  )
highlight <- data.frame(Chromosome = "chr05", xmin = 0.956, xmax = 7.03)

DEG_plot_27H <- ggplot(deg_27H, aes(x = Start_Mbp, y = log2FoldChange, color = color)) +
  geom_rect(data = highlight, aes(xmin = xmin, xmax = xmax, ymin = -Inf, ymax = Inf), inherit.aes = FALSE, fill = "gold", alpha = 0.2) +
  geom_point(alpha = 0.8, size = 0.5) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +
  facet_wrap(~ Chromosome, scales = "free_x", ncol = 4) +
  scale_color_manual(values = c("Up-regulated" = "red", "Down-regulated" = "blue")) +
  scale_y_continuous(name = expression(Log[2]~Fold~Change), breaks = pretty_breaks(6)) +
  scale_x_continuous(name = "Genomic Position (Mbp)", breaks = pretty_breaks(5)) +
  labs(title = "Genomic Distribution of Differentially Expressed Genes - 27 H") +
  theme_bw(base_size = 10) +
  theme(panel.grid = element_blank(), strip.text = element_text(size = 10, face = "bold"), plot.title = element_text(hjust = 0.5, face = "bold"), axis.text.x = element_text(angle = 0, hjust = 0.5), legend.position = "bottom")

ggsave("DEG_along_chromosome_27H.pdf", DEG_plot_27H, width = 11, height = 6)


