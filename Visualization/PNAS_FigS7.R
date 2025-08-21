# ============================================================
# Genome-wide Coverage Ratio Visualization mapped to DH Pahang
# ============================================================
# Author: Puyam Tondonba Singh
# Date: 2025-08-18
#
## This script processes BigWig coverage data from GCTCVs and Dwarf Cavendish,
## normalizes coverage by Z-score filtering, computes coverage ratios
## relative to a control sample, multiplies by 3 to represent triploidy,
## and generates chromosome-wise plots of genome-wide coverage ratios.
# ============================================================

# Load libraries
library(data.table)
library(rtracklayer)
library(ggplot2)
library(cowplot)

# ============================================================
# --- 1. Setup paths and data ---
bw_dir <- "E:/tony/WGS/Bigwig_all/10000bw"
bw_files <- list.files(bw_dir, pattern = "\\.bw$", full.names = TRUE)
names(bw_files) <- tools::file_path_sans_ext(basename(bw_files))
bw_files <- bw_files[!is.na(bw_files)]

control_id <- "PC.sorted.bam"
valid_chrs <- sprintf("chr%02d", 1:11)
binsize <- 100
Z_cutoff <- 1.64

# --- 2. Get reference lengths and genomic ranges ---
ref_dt <- as.data.table(as.data.frame(seqinfo(BigWigFile(bw_files[1]))), keep.rownames = "CHROM")
ref_len <- ref_dt[, .(CHROM, Length = seqlengths)]
ref_len <- ref_len[CHROM %in% valid_chrs]
ref_len[, start := 1]
genome_gr <- GRanges(ref_len[, .(seqnames = CHROM, start, end = Length)])

# --- 3. Fetch BigWig coverage data and process (main loop) ---
depth_dt_list <- list()
filtered_bins_list <- list()

for (file in bw_files) {
  sample_name <- tools::file_path_sans_ext(basename(file))
  
  rl <- import.bw(file, as = "RleList", which = genome_gr)
  rl <- rl[seqnames(genome_gr)]
  
  dt_chr_list <- list()
  for (chr in names(rl)) {
    rl_i <- rl[[chr]]
    nbins <- ceiling(length(rl_i) / binsize)
    views <- IRanges::Views(rl_i, start = seq(1, length(rl_i), by = binsize),
                            end = pmin(seq(1, length(rl_i), by = binsize) + binsize - 1, length(rl_i)))
    value_binned <- viewMeans(views)
    dt_chr_list[[chr]] <- data.table(Chr = chr, position = (seq_len(nbins) - 0.5) * binsize, value = value_binned)
  }
  
  dt <- rbindlist(dt_chr_list)
  
  # Z-score normalization
  dt[value > 0, Z := scale(value)]
  cov_mean <- dt[Z > -Z_cutoff & Z < Z_cutoff, mean(value, na.rm = TRUE)]
  
  # Apply Z-score filter and normalization
  dt[, value_filtered := fifelse(Z > -Z_cutoff & Z < Z_cutoff, value / cov_mean, NA_real_)]
  dt[, value_capped := pmin(value_filtered, 2)]
  dt[, filtered_out := is.na(value_filtered)]
  dt[, Sample := sample_name]
  
  # Store results
  depth_dt_list[[sample_name]] <- dt[, .(Sample, Chr, position, value = value_capped)]
  filtered_bins_list[[sample_name]] <- dt[filtered_out == TRUE, .(Sample, Chr, position)]
}

depth_dt <- rbindlist(depth_dt_list)
filtered_bins <- rbindlist(filtered_bins_list)

print("Summary of filtered bins:")
print(filtered_bins[, .N, by = Chr])

# --- 4. Bin coverage data ---
bin_size <- 300000
depth_binned <- depth_dt[value > 0]
depth_binned[, bin := cut(position, breaks = seq(0, max(position), by = bin_size), labels = FALSE)]
depth_binned <- depth_binned[, .(start = min(position), end = max(position), value = mean(value, na.rm = TRUE)),
                             by = .(Sample, Chr, bin)]

# --- 5. Make ratio matrix ---
wide_dt <- dcast(depth_binned, Chr + start + end + bin ~ Sample, value.var = "value")
sample_ids <- setdiff(names(wide_dt), c("Chr", "start", "end", "bin", control_id))
ratio_dt_list <- list()
for (id in sample_ids) {
  ratio_dt_list[[id]] <- wide_dt[, .(Chr, start, end, bin, Sample = id, ratio = get(id) / wide_dt[[control_id]])]
}
cov_ratio_melt <- rbindlist(ratio_dt_list)

# --- 6. Clean and prepare data for plotting ---
cov_ratio_melt[, ratio := ratio * 3]
cov_ratio_melt[, Sample := gsub("\\.sorted\\.bam$", "", Sample)]
cov_ratio_melt <- cov_ratio_melt[!Sample %in% c("119-3_1", "TC5-2")]
sample_order <- c("DwarfC", "TC5-1", "FM", "TC7", "217", "119-3")
cov_ratio_melt$Sample <- factor(cov_ratio_melt$Sample, levels = sample_order)
chromosomes <- sprintf("chr%02d", 1:11)

# --- 7. Generate plots and arrange ---
plot_matrix <- list()
plot_list_flat <- list()
plot_index <- 1

for (sample in sample_order) {
  for (chr in chromosomes) {
    
    dt <- cov_ratio_melt[Sample == sample & Chr == chr]
    chr_len <- ref_len[CHROM == chr, Length]
    
    highlight_region <- data.frame(
      xmin = 956425,
      xmax = 7030674,
      ymin = -Inf,
      ymax = Inf
    )
    
    p <- ggplot(dt, aes(x = (start + end) / 2, y = ratio, group = Sample)) +
      geom_hline(yintercept = 3, linetype = "dashed", color = "gray60", linewidth = 0.1)
    
    if (chr == "chr05") {
      p <- p + geom_rect(
        data = highlight_region,
        aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
        fill = "yellow", alpha = 0.3, inherit.aes = FALSE
      )
    }
    
    p <- p +
      geom_line(color = "black", linewidth = 0.2, alpha = 1) +
      geom_smooth(se = FALSE, span = 0.2, linewidth = 0.2, alpha = 0.6) +
      scale_y_continuous(limits = c(1, 5), breaks = seq(0, 5, 1)) +
      scale_x_continuous(
        breaks = seq(0, chr_len, by = 2e7),
        labels = function(x) x / 1e6,
        expand = c(0, 0),
        limits = c(0, chr_len)
      ) +
      theme_minimal(base_size = 6) +
      theme(
        axis.title = element_blank(),
        panel.grid.major = element_line(color = "gray90"),
        panel.grid.minor = element_blank(),
        axis.line = element_line(color = "black", linewidth = 0.2),
        legend.position = "none",
        strip.text = element_blank(),
        plot.margin = margin(1, 1, 1, 1)
      )
    
    # Adjust themes for grid alignment
    if (chr != chromosomes[1]) p <- p + theme(axis.text.y = element_blank(), axis.ticks.y = element_blank())
    if (sample != sample_order[length(sample_order)]) p <- p + theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())
    
    plot_list_flat[[plot_index]] <- p
    plot_index <- plot_index + 1
  }
}

aligned_plots <- align_plots(plotlist = plot_list_flat, align = "hv", axis = "tblr")
plots_grid <- plot_grid(plotlist = aligned_plots, ncol = length(chromosomes))

# --- 8. Add labels and save figure ---
chr_labels <- lapply(chromosomes, function(chr) ggdraw() + draw_label(chr, fontface = "bold", size = 8))
grid_with_chr <- plot_grid(plot_grid(plotlist = chr_labels, ncol = length(chromosomes)), plots_grid, ncol = 1, rel_heights = c(0.05, 1))

sample_labels <- lapply(sample_order, function(sample) ggdraw() + draw_label(sample, fontface = "bold", size = 8, angle = 0, hjust = 0))
grid_with_labels <- plot_grid(grid_with_chr, plot_grid(plotlist = sample_labels, ncol = 1), ncol = 2, rel_widths = c(1, 0.06))

Chromosomeall <- ggdraw(grid_with_labels) +
  draw_label("Coverage", x = 0.0, y = 0.5, angle = 90, vjust = 1, hjust = 0, size = 8) +
  draw_label("Genomic Position (Mbp)", x = 0.5, y = 0.0, vjust = 0, size = 8)

ggsave("figs5_pnas_clean.pdf", Chromosomeall, width = 18, height = 10, units = "cm", paper = "a4")