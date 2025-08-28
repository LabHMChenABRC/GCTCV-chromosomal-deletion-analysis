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


# ----------------------------
# Set Paths
# ----------------------------
ref_dir <- "E:/Visualization/S5"
setwd(ref_dir)

tree_file <- "240617.nwk"
domain_file <- "ACA_domain.csv"

# ----------------------------
# Load Tree and Process Bootstrap
# ----------------------------
phylo_tree <- read.tree(tree_file)

# Round bootstrap values
bs_values <- as.numeric(phylo_tree$node.label)
phylo_tree$node.label <- round(bs_values, 2)
phylo_tree$tip.label <- sprintf(" %s",phylo_tree$tip.label )
# ----------------------------
# Customize Tree Appearance
# ----------------------------

tree_plot <- ggtree(phylo_tree,linewidth=0.25) +
  coord_cartesian(xlim = c(-0.1, 5.5)) +
  geom_tiplab(size = 3, align = TRUE, color = "black") +
  geom_treescale(x = 0, y = 35, width = 0.2, color = "black",
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
Accession_IDs=c(
  "pfam12515" = "Ca2 + ATPase N terminal\nautoinhibitory domain", #CaATP_NAI
  "cl38396"   = "Cation ATPase", #Cation_ATPase superfamily
  "pfam13246" = "Cation ATPase", #Cation_ATPase
  "pfam00690" = "N terminus Cation ATPase", # Cation_ATPase_N
  "pfam00689" = "C terminus Cation ATPase", # Cation_ATPase_C
  "cl37964"   = "E1-E2 ATPase", #E1-E2_ATPase superfamily
  "pfam00122" = "E1-E2 ATPase", #E1-E2_ATPase
  "pfam00702" = "haloacid dehalogenase-like\nhydrolase" # Hydrolase
)
domain_data <- domain_data[domain_data$Accession %in% names(Accession_IDs),]
domain_data$full_name<-Accession_IDs[domain_data$Accession]

# Select relevant columns
domains <- domain_data %>%
  dplyr::select(Gene = Query, Start = From, End = To, Domain = full_name)

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
domains$Domain<-factor(domains$Domain, levels = unique(Accession_IDs))
# Ensure proper ordering
domains <- domains[order(match(domains$Gene, levels(domains$Gene))), ]

# ----------------------------
# Domain Plot
# ----------------------------

colorset = setNames(
  c("#fca103", "#326e4a", "#c7593e","#855996","#3395a6","#a8407d"),
     # "#7d6660","#5b63a6",
   # "#3395a6", "#a8407d", "#36b550"),
  nm = levels(domains$Domain)
)

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
    values=colorset,
    name = "Domain"
  ) +
  scale_x_continuous(expand = expansion(mult = c(0, 0))) +
  theme_minimal() +
  theme(
    legend.title.position = 'top',
    legend.position = "top",
    text = element_text(size = 10, color = "black"),
    axis.title.y = element_blank(),
    axis.text.y = element_blank(),
    axis.line.x = element_line(linewidth = 0.25),
    axis.ticks.x = element_line(linewidth = 0.25),
    panel.grid = element_blank(),
    panel.border = element_blank()
  ) +
  labs(x = "Length (amino acids)", y = NULL)

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
