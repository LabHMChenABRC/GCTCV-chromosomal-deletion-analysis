# ============================================================
# Triple Synteny Plot with DH Pahang in Middle
# ============================================================

# --- 1. Load Libraries ---
library(data.table)
library(dplyr)
library(ggplot2)
library(cowplot)

# --- 2. Define Functions ---

#' Load and filter a GFF file based on a genomic region.
#'
#' @param path Path to the GFF file.
#' @param region A data.table with Chr, Start, and End for filtering.
#' @return A data.table of genes within the specified region.
load_gff <- function(path, region) {
  gff <- fread(path, header = FALSE)
  colnames(gff) <- c("Chr", "Gene", "Start", "End")
  gff[, Mid := (Start + End) / 2]
  setkey(gff, Chr, Start, End)
  # Filter genes that overlap with the defined region
  foverlaps(region, gff, type = 'any', nomatch = 0L)[, c("Chr", "Gene", "Start", "End", "Mid")]
}

#' Efficiently parse a collinearity file from MCScanX.
#'
#' @param path Path to the .collinearity file.
#' @param tag1 A tag for the first genome.
#' @param tag2 A tag for the second genome.
#' @return A data.table with parsed synteny blocks.
parse_collinearity <- function(path, tag1, tag2, prefix1) {
  # This function uses a vectorized data.table approach
  lines <- readLines(path)
  
  # Create a data.table for processing
  dt <- data.table(line = lines, block_group = cumsum(grepl("^## Alignment", lines)))
  
  # Extract block ID from header lines
  headers <- dt[grepl("^## Alignment", line), 
                .(block_group, Block = as.integer(sub(".*Alignment\\s+(\\d+):.*", "\\1", line)))]
  
  # Parse gene lines and extract gene pairs
  # NOTE: This assumes the gene columns in the file are in a consistent order
  # that matches the subsequent merge logic.
  gene_pairs <- dt[grepl("^\\s+\\d+-[^:]+:\\s+", line), ]
  gene_pairs[, c("Gene1", "Gene2") := tstrsplit(sub("^\\s+\\d+-[^:]+:\\s+", "", line), split="\\s+",keep = 1:2)]
  
  # The Gene1 should start with prefix1 
  if ( ! grepl(paste("^",prefix1,sep =""),gene_pairs$Gene1[1]) ){
    setnames(gene_pairs,c("Gene1","Gene2"),c("Gene2","Gene1"))
  }
  # Join to get the correct block ID and clean up
  result <- headers[gene_pairs, on = "block_group"
  ][!is.na(Gene1) & !is.na(Gene2),
  ][, Pair1 := tag1
  ][, Pair2 := tag2
  ][,.(block_group,Gene1,Gene2,Pair1,Pair2)]
  
  return(result)
}

# --- 3. Configuration ---
setwd("E:/Visualization/S11")

paths <- list(
  gff_ban = "ban_FM05.gff",
  gff_dh  = "DH_FM05.gff",
  gff_bxh = "chr05.2_FM05.gff",
  col_ban_dh = "ban_FM05.collinearity",
  col_bxh_dh = "bxh2_FM.collinearity"
)

plot_params <- list(
  font_size = 6,
  # Define line width in mm for ggplot2, converting from points (0.5pt)
  line_size = 0.5 / ggplot2::.pt / 96 * 72,
  colors = list(
    ban_chr = "#c6dbef", ban_gene = "darkblue",
    dh_chr  = "grey",   dh_gene  = "black",
    bxh_chr = "#c7e9c0", bxh_gene = "darkgreen"
  ),
  chr_offset  = 0.08,
  gene_offset = 0.05,
  tags=c(
    "bxh" = 'chr05.2',
    "dh"  = 'chr05',
    "ban" = 'Ban05'
  ),
  y_pos=c(
    "bxh"=3,
    "dh" =2,
    "ban"=1
  )
)
plot_params$link_offset = max(plot_params$chr_offset,plot_params$gene_offset)
plot_params$chr_levels  = plot_params$tags[names(sort(plot_params$y_pos,decreasing = TRUE))]

# --- 4. Data Loading and Processing ---

# Define the Feature of Interest (FM) deletion region
FM.DEL.dt <- data.table(
  Chr = c("chr05", "chr05.2", "Ban05"),
  Start = c(956425L, 830729L, 1067819L),
  End = c(7030674L, 6869408L, 7314107L)
)
setkey(FM.DEL.dt, Chr, Start, End)

# Load GFF files for the region
gff_ban <- load_gff(paths$gff_ban, FM.DEL.dt)
gff_dh  <- load_gff(paths$gff_dh, FM.DEL.dt)
gff_bxh <- load_gff(paths$gff_bxh, FM.DEL.dt)

# Parse collinearity files
blocks_ban_dh <- parse_collinearity(path = paths$col_ban_dh, tag1 = "ban", tag2 = "dh", prefix1 = "Musa")
blocks_bxh_dh <- parse_collinearity(path = paths$col_bxh_dh, tag1 = "bxh", tag2 = "dh", prefix1 = "BXH2")

# Merge blocks with gene coordinates to create synteny links
# Note the order of gff files matters and corresponds to Gene1 and Gene2 from parsing
links_ban_dh <- merge(blocks_ban_dh, gff_ban[, .(Gene, Pos1 = Mid)], by.x = "Gene1", by.y = "Gene") %>%
  merge(gff_dh[, .(Gene, Pos2 = Mid)], by.x = "Gene2", by.y = "Gene")

links_bxh_dh <- merge(blocks_bxh_dh, gff_bxh[, .(Gene, Pos1 = Mid)], by.x = "Gene1", by.y = "Gene") %>%
  merge(gff_dh[, .(Gene, Pos2 = Mid)], by.x = "Gene2", by.y = "Gene")


# --- 5. Calculate Statistics for Bar Chart ---

# Combine all genes within the DEL region to get a master list and counts
all_genes_in_del <- rbind(gff_dh, gff_bxh, gff_ban) %>%
  as.data.table() %>%
  foverlaps(FM.DEL.dt, type = 'any', nomatch = 0L)

gene_counts_in_del <- all_genes_in_del[, .(N = .N), by = Chr]

# Calculate collinear gene counts for each comparison
count_collinear_genes <- function(blocks, all_genes) {
  collinear_genes <- c(blocks$Gene1, blocks$Gene2) %>% unique()
  sum(collinear_genes %in% all_genes)
}

collinear_counts <- data.table(
  group = c("chr05.2", "Ban05"),
  collinear = c(
    count_collinear_genes(blocks_bxh_dh, all_genes_in_del$Gene),
    count_collinear_genes(blocks_ban_dh, all_genes_in_del$Gene)
  )
)

# Calculate total genes for each comparison pair
total_counts <- data.table(
  group = c("chr05.2", "Ban05"),
  total = c(
    gene_counts_in_del[Chr %in% c("chr05", "chr05.2"), sum(N)],
    gene_counts_in_del[Chr %in% c("chr05", "Ban05"), sum(N)]
  )
)

# Combine for plotting
collinear_result.dt <- merge(collinear_counts, total_counts, by = "group")
collinear_result.dt[, group := factor(group, levels = plot_params$chr_levels )]

# --- 6. Prepare Data for Plotting ---

# Create a metadata table to organize plot elements. This avoids repetition.
plot_meta <- data.table(
  id   = c("bxh", "dh", "ban"),
  genome_name = c("Cavendish", "DH Pahang", "Cavendish"),
  Chr      = c("chr05.2", "chr05", "Ban05"),
  gff_data = list(gff_bxh, gff_dh, gff_ban)
)[,y_pos:=plot_params$y_pos[id]][]



# Prepare chromosome ranges, gene models, and labels 
plot_data <- plot_meta[, rbindlist(gff_data), by = .(id, y_pos)]
plot_data[, `:=`(Start_Mb = Start / 1e6, End_Mb = End / 1e6)]

chr_ranges <- plot_data[, .(xmin = min(Start_Mb), xmax = max(End_Mb)), by = .(id, y_pos)
                        ][,id_fill:=paste(id,'chr',sep='_')][]


gene_models <- plot_data[, .(id,
                             y_pos,
                             xmin=Start_Mb,
                             xmax=End_Mb,
                             id_fill=paste(id, 'gene', sep = '_')
                             )]

# Prepare synteny links for plotting
links_combined.plot <- rbind(links_ban_dh,links_bxh_dh)
links_combined.plot[, `:=`(y1 = plot_params$y_pos[Pair1], 
                           y2 = plot_params$y_pos[Pair2],
                           Pos1 = Pos1 / 1e6,
                           Pos2 = Pos2 / 1e6) ]
links_combined.plot[,y1_offset:=ifelse(y1>y2,-1,+1)*plot_params$link_offset]
links_combined.plot[,y2_offset:=ifelse(y1>y2,+1,-1)*plot_params$link_offset]

# Prepare labels for y-axis and gene counts
plot_labels <- merge(FM.DEL.dt, gene_counts_in_del, by = "Chr", all.x = TRUE)
plot_labels <- merge(plot_labels, plot_meta[, .(Chr, y_pos, genome_name, id)], by = "Chr", all.x = TRUE)
plot_labels[, `:=`(
  y_label = sprintf("%s %s\n%.2f-%.2f Mb", genome_name, Chr, Start/1e6, End/1e6),
  count_label = sprintf("n=%d", N),
  xmax = End  / 1e6
)]
setorder(plot_labels, y_pos)

# --- 7. adjust position ---
chr_pos_shift<-chr_ranges[,setNames(xmin,id)]

chr_ranges[,`:=`(xmin_shift=xmin-chr_pos_shift[id],
                 xmax_shift=xmax-chr_pos_shift[id])]

gene_models[,`:=`(xmin_shift=xmin-chr_pos_shift[id],
                  xmax_shift=xmax-chr_pos_shift[id])]

links_combined.plot[,`:=`(Pos1_shift=Pos1-chr_pos_shift[Pair1],
                          Pos2_shift=Pos2-chr_pos_shift[Pair2])]

plot_labels[,`:=`(xmax_shift=xmax-chr_pos_shift[id])]

# --- 8. Generate Plots ---

## 8a. Synteny Plot
synteny_p <- ggplot() +
  # Chromosome backbones
  geom_rect(data = chr_ranges, aes(xmin = xmin_shift, 
                                   xmax = xmax_shift, 
                                   ymin = y_pos - plot_params$chr_offset, 
                                   ymax = y_pos + plot_params$chr_offset, fill = id_fill), color = NA) +
  # Gene models
  geom_rect(data = gene_models, aes(xmin = xmin_shift, 
                                    xmax = xmax_shift, 
                                    ymin = y_pos - plot_params$gene_offset, 
                                    ymax = y_pos + plot_params$gene_offset, fill = id_fill), color = NA) +
  # Synteny links
  geom_segment(data = links_combined.plot, aes(x    = Pos1_shift, 
                                               xend = Pos2_shift, 
                                               y    = y1 + y1_offset , 
                                               yend = y2 + y2_offset ),
               alpha = 0.5, color = "grey40", linewidth = plot_params$line_size/10) +
  # Gene count labels
  geom_text(data = plot_labels, aes(y = y_pos, 
                                    x = xmax_shift, 
                                    label = count_label),
            hjust = -0.15, size = plot_params$font_size, size.unit = 'pt') +
  # Aesthestics and Theming
  scale_fill_manual(
    values = plot_params$colors,
    aesthetics = "fill" # Apply to both fill and color aesthetics if needed
  ) +
  scale_y_continuous(breaks = plot_labels$y_pos, labels = plot_labels$y_label) +
  scale_x_continuous(expand = expansion(mult = c(0.01, 0.01))) +
  coord_cartesian(clip = 'off') +
  labs(x = NULL, y = NULL, title = 'FM deletion region') +
  theme_cowplot(font_size = plot_params$font_size, rel_small = 1, rel_tiny = 1, rel_large = 1) + 
  theme(
    plot.title = element_text(hjust = 0.5),
    axis.line = element_blank(),
    axis.ticks = element_blank(),
    axis.text.x = element_blank(),
    axis.text.y = element_text(face = "bold", vjust = 0.8),
    panel.grid = element_blank(),
    legend.position = "none",
    plot.margin = margin(t=1,r = 1.2,l=0.5,unit = 'lines')
    )

## 8b. Collinear Gene Count Bar Plot
collinear_cnt_p <- ggplot(collinear_result.dt) +
  geom_col(aes(y = collinear / total * 100, x = group, fill = group),
           width = 0.65, color = 'black', linewidth = plot_params$line_size) +
  geom_text(aes(y = collinear / total * 100, x = group, label = sprintf("%d/%d", collinear, total)),
            vjust = -0.5, size = plot_params$font_size, size.unit = 'pt') +
  geom_text(aes(y = 50, x = group, label = sprintf("%.0f%%", collinear / total * 100)),
            size = plot_params$font_size, size.unit = 'pt', color = "black") +
  scale_fill_manual(
    values = c("chr05.2" = plot_params$colors$bxh_chr, "Ban05" = plot_params$colors$ban_chr),
    breaks = c("chr05.2", "Ban05")
  ) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.1)), name = 'Collinear genes (%)') +
  labs(x = "vs 'DH Pahang'", y = 'Collinear genes (%)', title = 'FM deletion region') +
  theme_half_open(font_size = plot_params$font_size, line_size = plot_params$line_size, rel_small = 1, rel_tiny = 1, rel_large = 1) +
  theme(
    legend.position = 'none',
    plot.title = element_text(hjust = 0.5),
    plot.margin = margin(t=1,l=0.5,unit = 'lines')
  )

# --- 9. Combine and Save Final Plot ---

Merge_p <- plot_grid(
  synteny_p,
  collinear_cnt_p,
  rel_widths = c(17, 3),
  labels = c("A", "B"),
  label_size = 10
)

PDF_FILE <- "FigS11.synteny_plot.pdf"
ggsave(PDF_FILE, plot = Merge_p, width = 20, height = 5, units = "cm", paper = 'A4')
shell.exec(PDF_FILE)
