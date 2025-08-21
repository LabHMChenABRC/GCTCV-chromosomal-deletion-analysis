# ============================================================
# Chromosome 5 Allele Ratio Visualization Whole genome


# This script generates a heatmap of variant allele ratios (VAR) 
# for FM deletion region on Chromosome 5 from whole-genome sequencing (WGS) data.
# It parses a VCF, filters high-quality SNPs, identifies heterozygous sites based on "Pei- Chiao" 
# calculates variant allele ratios for all samples, 
# and visualizes them as a heatmap.
# =========================================

# ============================================================
# Author: Puyam Tondonba Singh
# Date: 2025-08-18
# ============================================================

# ----------------------------
# Load Libraries
# ----------------------------
library(data.table)
library(dplyr)
library(ggplot2)
library(ggrastr)
library(patchwork)

read_vcf_table <- function(File, Samples) {
  rbindlist(lapply(File, function(x) {
    SampleName <- sub("(.+?)_.+", "\\1", basename(x))  # adjust if needed
    fread(x, header = TRUE, blank.lines.skip = TRUE) %>%
      setnames(sub("#", "", colnames(.))) %>%
      .[, Marker := SampleName] %>%
      .[, .SD, .SDcols = c("CHROM","POS","REF","ALT","FILTER","Marker", Samples)]
  }))[, .(Marker = paste0(unique(Marker), collapse = ",")),
       by = c("CHROM","POS","REF","ALT","FILTER", Samples)]
}

# ---- chr05 VCF only ----
library(data.table)

# ---- chr05 VCF only ----
VCF_File <- "E:/tony/WGS/wgs_vcf_bohan/snps_raw.chr05.HF.vcf"  # include extension
cultivars <- c("PC","TC5-2","TC4","TC7","GCTCV217","GCTCV119_L3")
allowed_genos <- c("0/1/1","0/0/1")

# Parse, keep biallelic SNPs, unique rows, sorted
VCF <- read_vcf_table(VCF_File, cultivars)[
  nchar(REF) == 1 & nchar(ALT) == 1
][
  , unique(.SD)
][
  order(CHROM, POS)
]

# ---- Quick per-sample GT extraction for filtering by PC genotype ----
genotype <- VCF[FILTER %in% c("PASS","SnpCluster")][
  , (cultivars) := lapply(.SD, function(col) tstrsplit(col, ":", fixed = TRUE, keep = 1L)[[1]]),
  .SDcols = cultivars
][
  order(CHROM, POS)
][
  , .SD, .SDcols = c("Marker","CHROM","POS","REF","ALT","FILTER", cultivars)
]

# Keep sites where PC is 0/1/1 or 0/0/1
VCF_filtered <- genotype[PC %in% allowed_genos, .(Marker, CHROM, POS)]

# --- Melt full fields for AD/DP parsing (limited to filtered sites) ---
vcf_long <- melt(
  VCF[FILTER %in% c("PASS","SnpCluster")],
  id.vars = c("Marker","CHROM","POS","REF","ALT","FILTER"),
  measure.vars = cultivars,
  variable.name = "Sample",
  value.name = "GenotypeField"
)[
  VCF_filtered, on = .(Marker, CHROM, POS)  # keep only loci where PC matched pattern
]

# Split GT:AD:DP:...
vcf_long[, c("GT","AD","DP","F4","Rest") := tstrsplit(GenotypeField, ":", fixed = TRUE)]
vcf_long[, c("R_AD","A1_AD") := tstrsplit(AD, ",", fixed = TRUE)]
vcf_long[, c("DP","R_AD","A1_AD") := lapply(.SD, as.numeric), .SDcols = c("DP","R_AD","A1_AD")]
vcf_long[, c("AD","F4","Rest") := NULL]

# --- Add DP_sum and allele ratio; check DP consistency ---
vcf_long[, DP_sum := R_AD + A1_AD]
vcf_long[, alleleratio := fifelse(!is.na(DP_sum) & DP_sum > 0, A1_AD / DP_sum, NA_real_)]
vcf_long[, correct := DP == DP_sum]

# --- Stats on DP mismatch ---
incorrect_count <- sum(!vcf_long$correct, na.rm = TRUE)
total_sites     <- nrow(vcf_long)
cat("Incorrect sites:", incorrect_count, "out of", total_sites, "\n")
cat("Percentage incorrect:", (incorrect_count / total_sites) * 100, "%\n")

# --- Filter for downstream analysis (coverage-supported) ---
GT <- vcf_long[DP_sum > 10]
GT[, alleleratio := A1_AD / DP_sum]

# Optional: how many CHROM:POS passed with PC in allowed_genos?
pc_pass_loci <- unique(genotype[PC %in% allowed_genos, .(CHROM, POS)])
cat("PC loci matching pattern (unique CHROM:POS):", nrow(pc_pass_loci), "\n")

GT_final <- GT[, .(CHROM, POS, Sample, alleleratio)]

# --- Cast to Wide Format ---
Cultivars_AR <- dcast(GT_final, CHROM + POS ~ Sample, value.var = "alleleratio") %>% na.omit()

# --- Semi-Join to Filtered Variants ---
Cultivars_AR <- semi_join(Cultivars_AR, VCF_filtered, by = "POS")

All_melt <- melt(Cultivars_AR, id.vars = c("CHROM", "POS"), variable.name = "Cultivar") %>% na.omit()
All_melt$Cultivar <- factor(All_melt$Cultivar, levels = c("PC","TC5-2",  "TC4", "TC7", "GCTCV217", "GCTCV119_L3"))

# --- Highlight Region for Plot ---
chr05_length <- VCF[CHROM == "chr05", max(POS, na.rm = TRUE)]
highlight_region <- data.frame(
  CHROM  = "chr05",
  Cultivar = c("TC4", "TC7", "GCTCV217", "GCTCV119_L3"),
  xmin   = c(956425, 1, 1, 1),
  xmax   = c(7030674, 13290586, 21500000, chr05_length),
  ymin = c(0, 0),
  ymax = c(1, 1)
)


highlight_region$Cultivar <- factor(highlight_region$Cultivar, levels = levels(All_melt$Cultivar))

library(ggrastr)

geom_plot <- ggplot(All_melt, aes(x = POS, y = value)) +
  
  geom_rect(
    data = highlight_region,
    aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
    inherit.aes = FALSE,
    fill = "#FFD6C1", alpha = 1
  ) +
  geom_point_rast(aes(colour = value), size = 0.2,stroke = NA, alpha = 0.1,raster.dpi = 300)+
  scale_x_continuous(
    breaks=seq(0,20*10^7,by=10*10^6),
    labels = seq(0,20*10^7,by=10*10^6)/10^6,
    expand = expansion(mult = c(0.01,0.01))
    )+
  scale_y_continuous(
    breaks = c(0, 0.33, 0.50, 0.67, 1),
    labels = c("", "0.33", "", "0.67", ""),
    expand = expansion(mult = c(0.1,0.1))
  ) +
  scale_color_gradientn(
    colors = c("navy"),
    breaks = c(0, 0.33, 0.50, 0.67, 1),
    labels = c("0", "0.33", "0.50", "0.67", "1"),
  ) +
  facet_grid(Cultivar ~ CHROM, scales = "free_y", space = "free_y") +
  labs(x = "Genomic Position (Mb)", y = "Allele Ratio") +
  theme_minimal(base_size = 8) +
  theme(
    panel.background = element_rect(fill = "grey95", colour = NA),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1, size = 6),
    axis.text.y = element_text(size = 6),
    axis.ticks = element_line(color = "black", size = 0.2),
    axis.line = element_line(color = "black", size = 0.3),
    axis.title = element_text(size = 8),
    strip.background = element_blank(),
    strip.text.y.right = element_text(size = 7, angle = 0),
    plot.title = element_text(size = 8, face = "bold", hjust = 0.5),
    legend.position = "none"
  )


ggsave("250815_chr5_geom.pdf", plot = geom_plot, width = 120, height = 80,units = "mm",paper='A4')
shell.exec("250815_chr5_geom.pdf")


# Lock facet order
lvl <- c("PC","TC5-2","TC4","TC7","GCTCV217","GCTCV119_L3")
All_melt$Cultivar <- factor(All_melt$Cultivar, levels = lvl)

# Per-cultivar highlight ranges on chr05
highlight_region <- data.frame(
  CHROM    = "chr05",
  Cultivar = c("TC4", "TC7", "GCTCV217", "GCTCV119_L3"),
  xmin     = c(956425, 1, 1, 1),
  xmax     = c(7030674, 13290586, 21500000, chr05_length),
  stringsAsFactors = FALSE
)

# Join and subset to the per-cultivar POS ranges
highlight_sub <- dplyr::inner_join(All_melt, highlight_region, by = c("CHROM","Cultivar")) %>%
  dplyr::filter(POS >= xmin, POS <= xmax)

# Reapply factor levels so join can't mess with order
highlight_sub$Cultivar <- factor(highlight_sub$Cultivar, levels = lvl)

# --- Final Density Plot (Highlighted Cultivars Only) ---
highlight_density <- ggplot(All_melt, aes(x = value)) +
  facet_grid(Cultivar ~ ., switch = "y") +  # respects 'lvl' order
  geom_density(fill = "navy", color = NA, alpha = 0.8) +
  geom_density(
    data = highlight_sub,
    aes(x = value),
    fill = "#FFD6C1", color = NA, alpha = 0.5
  ) +
  coord_flip() +
  scale_x_continuous(
    limits = c(0, 1),
    breaks = c(0, 0.33, 0.50, 0.67, 1),
    labels = c("", "", "", "", "")
  ) +
  scale_y_continuous(limits = c(0, 4)) +
  xlab("Heterozygote allele ratio") + ylab("Density") +
  theme_minimal(base_size = 8) +
  theme(
    panel.background = element_rect(fill = "grey95", colour = NA),
    panel.grid = element_blank(),
    axis.text.x  = element_text(size = 6),
    axis.ticks.x = element_line(color = "black", size = 0.2),
    axis.line.x  = element_line(color = "black", size = 0.3),
    axis.text.y  = element_blank(),
    axis.ticks.y = element_blank(),
    axis.line.y  = element_blank(),
    axis.title.y = element_blank(),
    axis.title.x = element_text(size = 8),
    strip.background = element_blank(),
    strip.placement = "outside",
    strip.text.y.left  = element_blank(),
    strip.text.y.right = element_text(size = 7, angle = 0),
    legend.position = "none"
  )


# --- 3. Combine using patchwork ---
final_plot <- geom_plot + highlight_density +
  plot_layout(ncol = 2, widths = c(4, 0.6), guides = "collect")


ggsave("250825_chr5_geom_density.pdf", plot = final_plot, width = 15, height = 10, dpi = 300, units = "cm", paper ="a4")

shell.exec("250825_chr5_geom_density.pdf")





















