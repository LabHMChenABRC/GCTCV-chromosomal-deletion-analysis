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
protein_file <-'ACA_protein.fasta'

bs_values_show_cutoff <- 0.7
font.size=7
line.size=0.5/ggplot2::.pt*96/72

# ----------------------------
# Load Tree and Process Bootstrap
# ----------------------------
phylo_tree <- read.tree(tree_file)

bs_values <- as.numeric(phylo_tree$node.label)

# remove BS label <= cutoff
bs_values[bs_values  <= bs_values_show_cutoff] <-NA 

# Round bootstrap values
phylo_tree$node.label<- round(bs_values, 2)*100
phylo_tree$tip.label <- sprintf(" %s",phylo_tree$tip.label )

# ----------------------------
# Customize Tree Appearance
# ----------------------------
# tree scale y position
treescale_y_pos      <- length(phylo_tree$tip.label) - 1

tree_plot <- ggtree(phylo_tree,linewidth=line.size) +
  coord_cartesian(xlim = c(-0.1, 5.5)) +
  geom_tiplab(size = font.size/ggplot2::.pt, align = TRUE, color = "black",linesize = line.size,linetype = 2) +
  geom_treescale(x = 0, y = treescale_y_pos, width = 0.2, color = "black",
                 fontsize = font.size/ggplot2::.pt, linesize = line.size/0.5*0.75, offset = 0.5) +
  geom_nodelab(nudge_y = 0.5, hjust = 1.4, size = (font.size-1)/ggplot2::.pt) +
  scale_x_continuous(expand = expansion(mult = c(0, 0))) +
  coord_cartesian(clip = 'off')+
  theme(
    legend.position = "none",
    line = element_line(linewidth = line.size),
    text = element_text(size = font.size, color = "black"),
    panel.grid = element_blank(),
    axis.line = element_blank(),
    axis.text = element_blank(),
    axis.title = element_blank(),
    plot.margin = margin(r = 1.5,unit = 'lines')
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
  "pfam00702" = "haloacid dehalogenase-like\nhydrolase", # Hydrolase
  "cl37964"   = "E1-E2 ATPase", #E1-E2_ATPase superfamily
  "pfam00122" = "E1-E2 ATPase", #E1-E2_ATPase
  "pfam00689" = "C terminus Cation ATPase" # Cation_ATPase_C
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
# Baseline
# ----------------------------

protein.seq        <- read.FASTA(protein_file,type = 'AA')
names(protein.seq) <- sub("\\s+.+$","",names(protein.seq))
protein.length     <- sapply(protein.seq, length)
baseline.df        <- data.frame(Gene = factor(gene_order,levels = rev(gene_order)), Start = 0, End = protein.length[gene_order])


# ----------------------------
# Domain Plot
# ----------------------------

colorset = setNames(
  c("#fca103", "#326e4a", "#c7593e", "#a8407d", "#3395a6", "#855996"),
  nm = levels(domains$Domain)
)

domain_plot <- ggplot() +
  # Baseline for each protein
  geom_segment(
    data = baseline.df,
    aes(x = Start, xend = End, y = Gene, yend = Gene),
    size = 0.2, color = "black"
  ) +
  # Domains colored by type
  geom_segment(
    data = domains,
    aes(x = Start, xend = End, y = Gene, yend = Gene, color = Domain),
    size = (font.size*1.1)/ggplot2::.pt, lineend = "butt"
  ) +
  scale_color_manual(
    values=colorset,
    name = "Domain"
  ) +
  scale_x_continuous(expand = expansion(mult = c(0, NA)),breaks = seq(0,100,3)*100) +
  theme_minimal(base_size = font.size,base_line_size = line.size) +
  coord_cartesian(clip = 'off')+
  theme(
    legend.title.position = 'top',
    legend.justification.top = 'center',
    legend.position = "top",
    legend.title = element_text(hjust = 0.5),
    text = element_text(size = font.size, color = "black"),
    axis.title.y = element_blank(),
    axis.text.y = element_blank(),
    axis.text.x = element_text(size = font.size, color = "black"),
    axis.line.x = element_line(linewidth = line.size),
    axis.ticks.x = element_line(linewidth = line.size),
    panel.grid = element_blank(),
    panel.border = element_blank(),
    plot.margin = margin(l = 1,r=1,unit = 'lines')
    
  ) +
  labs(x = "Length (amino acids)", y = NULL)

# ----------------------------
# Combine Tree and Domain Plots
# ----------------------------
domain_legend <- get_plot_component(domain_plot+theme(legend.key.spacing.y = unit(0,'lines'),legend.key.height = unit(1,'lines')), pattern = 'guide-box-top')
combined_plot <- plot_grid(tree_plot, domain_plot+theme(legend.position = "none"),rel_widths = c(9,5), align = "vh", axis = "l")
combined_plot <- plot_grid(domain_legend,combined_plot,ncol = 1,rel_heights = c(0.15,1))
# Display
print(combined_plot)

# ----------------------------
# Save Output
# ----------------------------
pdf_file="FigS5.ACA_combined_plot.pdf"
ggsave(pdf_file, combined_plot,width = 14, height = 16.5, units = "cm", paper = "A4")
shell.exec(pdf_file)
