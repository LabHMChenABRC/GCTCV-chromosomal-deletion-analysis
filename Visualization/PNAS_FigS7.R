# ============================================================
# Genome-wide Coverage Ratio Visualization (Huang et al., 2023)
# ============================================================
# Author: Puyam Tondonba Singh
# Date: 2025-08-18
#
# This script processes BigWig coverage data from GCTCVs and Dwarf Cavendish,
# normalizes by Z-score filtering, computes ratios relative to a control,
# and visualizes coverage across chromosomes.
# ============================================================

# Load required libraries
library(data.table)
library(rtracklayer)
library(ggplot2)
library(cowplot)

# Set working directory and define input
bw_dir <- "E:/Visualization/S7"
setwd(bw_dir)

# Load and order BigWig files
bw_files_all <- list.files(bw_dir, pattern = "\\.bw$", full.names = TRUE)
names(bw_files_all) <- gsub(".bw", "", basename(bw_files_all))

desired_order <- c("DwarfC", "PC", "FM", "TC5", "TC7", "217", "119")
bw_files <- bw_files_all[desired_order]
names(bw_files)<-c("DC", "PC", "FM", "TC5", "TC7", "217", "119")

# Get genome info
ref_len_dt <- as.data.table(as.data.frame(seqinfo(BigWigFile(bw_files[1]))), keep.rownames = "CHROM")[grepl('^chr',CHROM), .(CHROM, Length = seqlengths)]
genome_gr <- GRanges(ref_len_dt[, .(seqnames = CHROM, start = 1, end = Length)])

# Process and normalize coverage
# position was divide by 100
raw_coverage_list <- lapply(bw_files, function(file) {
  message("Processing: ", basename(file))
  rl <- import(file,  which = genome_gr, as = "RleList")
  rl <- rl[runValue(seqnames(genome_gr))]
  
  dt <- rbindlist(lapply(names(rl), function(chr_name) {
    rl_i <- rl[[chr_name]]
    data.table(
      position = seq_len(sum(runLength(rl_i) / 100)),
      value = rep(runValue(rl_i), runLength(rl_i) / 100),
      Chr = chr_name
    )
  }))
  
  Z_cutoff <- 1.64
  dt[value > 0, Z := scale(value)]
  cov_mean <- dt[Z > -Z_cutoff & Z < Z_cutoff, mean(value, na.rm = TRUE)]
  
  dt[, value := ifelse(Z > -Z_cutoff & Z < Z_cutoff, value / cov_mean, 0)]
  dt[value > 2, value := 2]
  
  dt[, .(Chr, position, value)]
})

depth_dt <- rbindlist(raw_coverage_list,idcol = 'Sample')

# Bin coverage into windows
bin_size <- 10000
depth_dt[value == 0, value := NA]
depth_dt[, bin := cut(position,  breaks = 0:ceiling(max(position) / bin_size)*bin_size, labels = FALSE)]


depth_binned <- depth_dt[, .(
  start = min(position),
  end = max(position),
  value = mean(value, na.rm = TRUE)
), by = .(Sample, Chr, bin)]

# Compute coverage ratio to control
control_id <- "PC"
wide_dt    <- dcast(depth_binned, Chr + start + end + bin ~ Sample, value.var = "value")
sample_ids <- setdiff(names(wide_dt), c("Chr", "start", "end", "bin", control_id))

ratio_dt       <- wide_dt[, (sample_ids) := lapply(.SD, function(x) x / get(control_id)), .SDcols = sample_ids]
cov_ratio_melt <- melt(ratio_dt, 
                       id.vars = c("Chr", "start", "end", "bin"), 
                       measure.vars = sample_ids, 
                       variable.name = "Sample", 
                       value.name = "ratio")

# Tidy and transform
cov_ratio_melt[, ratio := ratio * 3]
desired_order_plot <- c("DC", "FM", "TC5", "TC7", "217", "119")
cov_ratio_melt$Sample <- factor(cov_ratio_melt$Sample, levels = desired_order_plot)

# --- Build per-chromosome bp->index conversion 
chr_lengths <- ref_len_dt[CHROM %in% unique(cov_ratio_melt$Chr)]
setnames(chr_lengths, c("CHROM","Length"), c("Chr","ChrLenBP"))

# max plotted unit per Chr (use end because bins cover full range)
chr_span_units <- depth_binned[, .(MaxUnit = max(end, na.rm = TRUE)), by = Chr]

# --- Define FM deletion in bp and convert to your x-units for chr05 ---
fm_chr <- "chr05"
fm_bp_min <- 956425
fm_bp_max <- 7030674
shade_df <- data.table(Chr=fm_chr,
                       xmin=fm_bp_min,
                       xmax=fm_bp_max,
                       ymin=-Inf,
                       ymax=Inf)

plot_chr_region<-ref_len_dt[,.(Chr=CHROM,
                               xmin=0,
                               xmax=Length)]

linesize=0.5
# position is multiple by 100
p <- ggplot(cov_ratio_melt, aes(x = (start + end)/2*100, y = ratio, group = Sample)) +
  geom_rect(
    data=plot_chr_region,
    aes(xmin = xmin, xmax = xmax),ymin = NA, ymax = NA,
    inherit.aes = FALSE
  )+
  geom_rect(
    data = shade_df,
    aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
    inherit.aes = FALSE,
    fill = "yellow", alpha = 0.18
  ) +
  geom_hline(yintercept = 3, linetype = "dashed", color = "gray60", linewidth = 0.3) +
  geom_line(color = "#3B528BFF", linewidth = linesize * 72 / 96, alpha = 0.9) +
  facet_grid(Sample ~ Chr, scales = "free_x", space = "free_x") +
  scale_y_continuous(limits = c(1, 5), breaks = seq(1, 5, 1)) +
  scale_x_continuous(
    breaks = seq(0, max(ref_len_dt$Length),2e7),
    labels = function(x) x / 1e6,
    expand = c(0, 0)
  ) +
  labs(x = "DH Pahang' genome (Mb)", y = "Normalized score") +
  theme_minimal(base_size = 10) +
  theme(
    strip.text.x = element_text(size = 10, face = "bold"),
    strip.text.y = element_text(size = 10, face = "bold", angle = 0),
    panel.grid.major = element_line(color = "gray90"),
    panel.grid.minor = element_blank(),
    axis.line = element_line(color = "black"),
    axis.ticks = element_line(color = "black"),
    axis.title = element_text(size = 10),
    axis.text = element_text(size = 8),
    panel.spacing = unit(0.6, "lines") 
  )

# Save and open PDF
PDF_FILE="FigS7.All_chr_copyratio_DHv4.pdf"
ggsave(PDF_FILE, p, width = 18, height = 10, units = "cm", paper = "a4")
shell.exec(PDF_FILE)
