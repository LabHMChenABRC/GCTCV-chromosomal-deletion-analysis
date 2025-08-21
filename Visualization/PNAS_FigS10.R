# ============================================================
# Genome-wide Coverage Ratio Visualization mapped to cavendish li et al., 2023
# ============================================================
# Author: Puyam Tondonba Singh
# Date: 2025-08-18
#
## This script processes BigWig coverage data from GCTCVs and Dwarf Cavendish,
## normalizes coverage by Z-score filtering, computes coverage ratios
## relative to a control sample
## and generates chromosome-wise plots of genome-wide coverage ratios.
# ============================================================

# Load required libraries
library(data.table)
library(rtracklayer)
library(ggplot2)
library(cowplot)


# Set your BigWig directory and define the control sample ID
bw_dir <- "E:/tony/WGS/nature_baxi_bw10000" # 
setwd(bw_dir)

# Load and order the BigWig files based on a desired sequence
bw_files_all <- list.files(bw_dir, pattern = "\\.bw$", full.names = TRUE)
names(bw_files_all) <- tools::file_path_sans_ext(basename(bw_files_all))

desired_order_files <- c("Dwarf C", "PC", "TC4", "TC5", "TC7", "217", "119")
bw_files <- bw_files_all[match(desired_order_files, names(bw_files_all))]
bw_files <- bw_files[!is.na(bw_files)]

#  Data Processing Pipeline ==========

# Get Genome Information ---
ref_len_dt <- as.data.table(as.data.frame(seqinfo(BigWigFile(bw_files[1]))), keep.rownames = "CHROM")[, .(CHROM, Length = seqlengths)]
genome_gr_all <- GRanges(ref_len_dt[, .(seqnames = CHROM, start = 1, end = Length)])

# Define valid chromosomes and subset the GRanges object
valid_chrs <- grep("Ban\\d+|Dh\\d+|Ze\\d+", as.character(seqnames(genome_gr_all)), value=TRUE)

genome_gr <- genome_gr_all[as.character(seqnames(genome_gr_all)) %in% valid_chrs]


#  Fetch and Normalize Coverage ---
raw_coverage_list <- lapply(bw_files, function(file) {
  message("  Processing: ", basename(file))
  rl <- import.bw(file, as = "RleList", which = genome_gr)
  
  # Convert RleList to a data.table with 100bp resolution
  dt <- rbindlist(lapply(names(rl), function(chr_name) {
    rl_i <- rl[[chr_name]]
    run_len <- runLength(rl_i)
    run_val <- runValue(rl_i)
    data.table(
      position = seq_len(sum(run_len / 100)),
      value = rep(run_val, run_len / 100),
      Chr = chr_name
    )
  }))
  
  # Normalize coverage based on Z-score filtering to find a stable mean
  Z_cutoff <- 1.64
  dt[value > 0, Z := scale(value)]
  cov_mean <- dt[Z > -Z_cutoff & Z < Z_cutoff, mean(value, na.rm = TRUE)]
  
  # Apply normalization and cap values at 2x the mean
  dt[, value := ifelse(Z > -Z_cutoff & Z < Z_cutoff, value / cov_mean, 0)]
  dt[value > 2, value := 2]
  
  # Add sample name and return selected columns
  dt[, Sample := tools::file_path_sans_ext(basename(file))]
  return(dt[, .(Sample, Chr, position, value)])
})

depth_dt <- rbindlist(raw_coverage_list)


# Bin Coverage into Larger Windows ---
bin_size <- 20000
depth_dt[value == 0, value := NA] # Exclude zero-coverage regions from mean calculation
depth_dt[, bin := cut(position, breaks = seq(0, max(position), by = bin_size), labels = FALSE)]
depth_binned <- depth_dt[, .(
  start = min(position),
  end = max(position),
  value = mean(value, na.rm = TRUE)
), by = .(Sample, Chr, bin)]

control_id <- "PC"

# Calculate Ratio to Control ---
wide_dt <- dcast(depth_binned, Chr + start + end + bin ~ Sample, value.var = "value")
sample_ids <- setdiff(names(wide_dt), c("Chr", "start", "end", "bin", control_id))
ratio_dt <- wide_dt[, (sample_ids) := lapply(.SD, function(x) x / get(control_id)), .SDcols = sample_ids]
cov_ratio_melt <- melt(ratio_dt,
                       id.vars = c("Chr", "start", "end", "bin"),
                       measure.vars = sample_ids,
                       variable.name = "Sample",
                       value.name = "ratio")


# Final Data Tidying and Preparation for Faceting ---
# Rename samples for cleaner labels
cov_ratio_melt[Sample == "Dwarf C", Sample := "DC"]

# Set sample order for plotting
desired_order_plot <- c("DC", "TC4", "TC5", "TC7", "217", "119")
cov_ratio_melt$Sample <- factor(cov_ratio_melt$Sample, levels = desired_order_plot)

# Split 'Chr' column into Chromosome and Haplotype for faceting
cov_ratio_melt[, c("Chromosome", "Haplotype") := tstrsplit(Chr, "\\.", fixed = TRUE)]
cov_ratio_melt[, Haplotype := paste("Haplotype", Haplotype)]

# ===== Plotting =====
fontsize <- 7
linesize <- 0.5
palette <- c("#3B528BFF", "#21908CFF", "#5DC863FF")  # Ze, Ban, Dh
chrms <- sprintf("Chr%02d", 1:11)
hap_prefix <- c("Ze", "Ban", "Dh")

# Store all chromosome plots
all_chr_plots <- list()

for (chr in chrms) {
  
  # Construct haplotype names: Ze05, Ban05, Dh05
  hap_chrs <- paste0(hap_prefix, sub("Chr", "", chr))
  hap_plots <- list()
  
  for (i in seq_along(hap_chrs)) {
    hap_chr <- hap_chrs[i]
    
    # Get chromosome length
    chr_len <- ref_len_dt$Length[ref_len_dt$CHROM == hap_chr]
    if (length(chr_len) != 1) {
      message("Skipping missing chromosome: ", hap_chr)
      next
    }
    
    # Subset coverage data
    dt_chr <- cov_ratio_melt[Chr == hap_chr]
    if (nrow(dt_chr) == 0) {
      message("No data for: ", hap_chr)
      next
    }
    
    # Generate ggplot for this haplotype
    p <- ggplot(dt_chr, aes(x = (start + end) / 2 * 100, y = ratio, group = Sample)) +
      geom_hline(yintercept = 1, linetype = "dashed", color = "gray60", linewidth = 0.3) +
      geom_line(color = "#3B528BFF", linewidth = linesize * 72 / 96, alpha = 0.9) +
      geom_smooth(se = FALSE, span = 0.2, color = palette[i], linewidth = 0.6, alpha = 0.6) +
      facet_wrap(~Sample, nrow = 1) +
      scale_y_continuous(limits = c(0, 2.2), breaks = seq(0, 2, 1)) +
      scale_x_continuous(
        breaks = seq(0, chr_len, by = 2e7),
        labels = function(x) x / 1e6,
        expand = c(0, 0),
        limits = c(0, chr_len)
      ) +
      labs(y = hap_chr, x = NULL) +
      theme_minimal(base_size = fontsize) +
      theme(
        strip.text = element_text(size = fontsize, face = "bold"),
        axis.title.x = element_blank(),
        axis.title.y = element_text(size = fontsize),
        panel.grid.major = element_line(color = "gray90"),
        panel.grid.minor = element_blank(),
        plot.margin = margin(5, 5, 5, 5),
        axis.line = element_line(color = "black", linewidth = 0.2)
      )
    
    hap_plots[[i]] <- p
  }
  
  if (length(hap_plots) == 0) next
  
  # Combine haplotypes vertically
  combined_hap <- plot_grid(plotlist = hap_plots, ncol = 1, align = "v")
  
  # Add x-axis label as a separate plot
  x_label <- ggplot() +
    theme_void() +
    annotate("text", x = 0.5, y = 0.5, label = "Genomic Position (Mb)", size = fontsize / ggplot2::.pt)
  
  all_chr_plots[[chr]] <- plot_grid(
    combined_hap, x_label, ncol = 1, rel_heights = c(1, 0.06),
    labels = chr, label_size = fontsize + 2
  )
}

# Combine all chromosomes into a single figure
p_all <- plot_grid(plotlist = all_chr_plots, ncol = 4)
ggsave("All_chr_copyratio.pdf", p_all, width = 20, height = 18, units= "cm", paper = 'a4')
getwd()
