# ============================================================
# ACA Phylogeny and Domain Visualization

# This script Visualize phylogenetic tree and protein domain structures of ACA/ECA genes 

# ============================================================
# Author: Puyam Tondonba Singh
# Date: 2025-08-18

# ----------------------------
# Load Libraries
# ----------------------------
library(ape)
library(ggtree)
library(ggplot2)
library(cowplot)
library(dplyr)
library(forcats)

# ----------------------------
# Set Paths
# ----------------------------
ref_dir <- "E:/tony/PHYLOGENY/ACA"
setwd(ref_dir)

tree_file <- "240617.nwk"
domain_file <- "ACA_domain.csv"

# ----------------------------
# Load Tree and Process Bootstrap
# ----------------------------
phylo_tree <- read.tree(tree_file)

# Convert tree to tibble for inspection (optional)
tree_tbl <- as_tibble(phylo_tree)

# Round bootstrap values
bs_values <- as.numeric(phylo_tree$node.label)
phylo_tree$node.label <- round(bs_values, 2)

# Initialize ggtree object
tree_plot <- ggtree(phylo_tree)

# ----------------------------
# Customize Tree Appearance
# ----------------------------
tree_plot <- tree_plot +
  coord_cartesian(xlim = c(-0.1, 5.5)) +
  geom_tiplab(size = 3, align = TRUE, color = "black") +
  geom_treescale(x = 0, y = 32.5, width = 0.2, color = "black",
                 fontsize = 3, linesize = 0.55, offset = 1) +
  geom_nodelab(nudge_y = 0.5, hjust = 1.4, size = 3) +
  scale_x_continuous(expand = expansion(mult = c(0, 0.34))) +
  theme(
    legend.position = "none",
    text = element_text(size = 10, color = "black"),
    panel.grid = element_blank(),
    axis.line = element_blank(),
    axis.text = element_blank(),
    axis.title = element_blank()
  )

# ----------------------------
# Load and Process Domain Data
# ----------------------------
domain_data <- read.csv(domain_file, header = TRUE, sep = ",")

# Clean Query names
domain_data$Query <- gsub("Q#\\d+ - >", "", domain_data$Query)

# Select relevant columns
domains <- domain_data %>%
  dplyr::select(Gene = Query, Start = From, End = To, Domain = Short.name)

# Convert positions to numeric
domains$Start <- as.numeric(domains$Start)
domains$End <- as.numeric(domains$End)

# ----------------------------
# Define Gene Order
# ----------------------------
gene_order <- c(
  "MaACA6", "MaACA3", "MaACA9", "MaACA7", "OsACA12", "OsACA5", "AtACA7", "AtACA2",
  "OsACA4", "AtACA1", "MaACA11", "MaACA8", "OsACA11", "OsACA10", "OsACA7", "OsACA1",
  "AtACA11", "AtACA4", "MaACA2", "MaACA1", "MaACA4", "MaACA12", "MaACA13", "OsACA6",
  "AtACA10", "AtACA8", "AtACA9", "MaACA10", "MaACA5", "OsACA9", "AtACA13", "AtACA12",
  "MaECA2", "MaECA7", "MaECA4", "MaECA1", "OsECA1", "AtECA04", "AtECA01", "MaECA6",
  "MaECA3", "AtECA02", "MaECA5", "AtECA03", "AmTrACA"
)

domains$Gene <- factor(domains$Gene, levels = rev(gene_order))

# Ensure proper ordering
domains <- domains[order(match(domains$Gene, levels(domains$Gene))), ]

# ----------------------------
# Domain Plot
# ----------------------------
domain_plot <- ggplot() +
  # Baseline for each gene
  geom_segment(
    data = data.frame(Gene = unique(domains$Gene), Start = 0, End = 1200),
    aes(x = Start, xend = End, y = Gene, yend = Gene),
    size = 0.2, color = "black"
  ) +
  # Domains colored by type
  geom_segment(
    data = domains,
    aes(x = Start, xend = End, y = Gene, yend = Gene, color = Domain),
    size = 3, lineend = "butt"
  ) +
  scale_color_manual(
    values = c(
      "#fca103", "#326e4a", "#5b63a6",
      "#855996", "#c7593e", "#7d6660",
      "#3395a6", "#a8407d", "#36b550"
    ),
    name = "Domain Type"
  ) +
  scale_x_continuous(expand = expansion(mult = c(0, 0))) +
  theme_minimal() +
  theme(
    legend.position = "right",
    text = element_text(size = 10, color = "black"),
    axis.title.y = element_blank(),
    axis.text.y = element_blank(),
    panel.grid = element_blank(),
    panel.border = element_blank()
  ) +
  labs(x = "Protein Position (aa)", y = NULL)

# ----------------------------
# Combine Tree and Domain Plots
# ----------------------------
combined_plot <- plot_grid(tree_plot, domain_plot, align = "vh", axis = "l")

# Display
print(combined_plot)

# ----------------------------
# Save Output
# ----------------------------
ggsave("ACA_combined_plot.pdf", combined_plot,
       device = "pdf", width = 30, height = 15, units = "cm", paper = "A4")
