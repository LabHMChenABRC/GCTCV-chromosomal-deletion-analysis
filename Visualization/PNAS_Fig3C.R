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
    fread(x, header = TRUE, blank.lines.skip = TRUE) %>%
      setnames(sub("#", "", colnames(.))) %>% 
      .[, .SD, .SDcols = c("CHROM","POS","REF","ALT","FILTER", Samples)]
  }))
}

# ---- chr05 VCF only ----
VCF_File <- "E:/Visualization/3D/snps_raw.chr05.HF.vcf.gz"  # include extension
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
  , .SD, .SDcols = c("CHROM","POS","REF","ALT","FILTER", cultivars)
]

# Keep sites where PC is 0/1/1 or 0/0/1
VCF_filtered <- genotype[PC %in% allowed_genos, .( CHROM, POS)]

# --- Melt full fields for AD/DP parsing (limited to filtered sites) ---
setkey(VCF,CHROM, POS)
vcf_long <- melt(
  VCF[VCF_filtered], # keep only loci where PC matched pattern
  id.vars = c("CHROM","POS","REF","ALT","FILTER"),
  measure.vars = cultivars,
  variable.name = "Sample",
  value.name = "GenotypeField"
)

# Split GT:AD:...
vcf_long[, c("GT","AD") := tstrsplit(GenotypeField, ":", fixed = TRUE,keep = 1:2)]
vcf_long[, c("R_AD","A1_AD") := tstrsplit(AD, ",", fixed = TRUE)]
vcf_long[, c("R_AD","A1_AD") := lapply(.SD, as.numeric), .SDcols = c("R_AD","A1_AD")]
vcf_long[, c("AD") := NULL]

# --- Add AD_sum and allele ratio ---
vcf_long[, AD_sum := R_AD + A1_AD]
vcf_long[, alleleratio := fifelse(!is.na(AD_sum) & AD_sum > 0, A1_AD / AD_sum, NA_real_)]

# --- Filter for downstream analysis (coverage-supported) ---
GT <- vcf_long[AD_sum > 10]

# Optional: how many CHROM:POS passed with PC in allowed_genos?
pc_pass_loci <- unique(genotype[PC %in% allowed_genos, .(CHROM, POS)])
cat("PC loci matching pattern (unique CHROM:POS):", nrow(pc_pass_loci), "\n")

GT_final <- GT[, .(CHROM, POS, Sample, alleleratio)]

# --- Cast to Wide Format ---
Cultivars_AR <- dcast(GT_final, CHROM + POS ~ Sample, value.var = "alleleratio") %>% na.omit()

# --- Semi-Join to Filtered Variants ---
Cultivars_AR <- semi_join(Cultivars_AR, VCF_filtered, by = "POS")

All_melt <- melt(Cultivars_AR, id.vars = c("CHROM", "POS"), variable.name = "Cultivar") %>% na.omit()
All_melt$Cultivar <- factor(All_melt$Cultivar, levels = c("PC","TC5-2",  "TC4", "TC7", "GCTCV217", "GCTCV119_L3"),labels = c("PC","TC5","FM","TC7","217","119"))

# --- Highlight Region for Plot ---
chr05_length <- 46513039
var_highlight_region <- data.frame(
  CHROM  = "chr05",
  Cultivar = c("FM", "TC7", "217", "119"),
  xmin   = c(956425, 1, 1, 1),
  xmax   = c(7030674, 13290586, 21500000, chr05_length),
  ymin = c(0, 0),
  ymax = c(1, 1)
)


var_highlight_region$Cultivar <- factor(var_highlight_region$Cultivar, levels = levels(All_melt$Cultivar))

library(ggrastr)

geom_plot <- ggplot(All_melt, aes(x = POS, y = value)) +
  
  geom_rect(
    data = var_highlight_region,
    aes(xmin = xmin, xmax = xmax, ymin = 0, ymax = 1),
    inherit.aes = FALSE,
    fill = "#f6c9cd", alpha = 1
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
  facet_grid(Cultivar ~ .) +
  labs(x = "Genomic Position (Mb)", y = "Allele Ratio") +
  theme_minimal(base_size = 8) +
  theme(
    panel.background = element_rect(fill = "grey95", colour = NA),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.text.x = element_text(angle = 0, hjust = 0.5, size = 6),
    axis.text.y = element_text(size = 6),
    axis.ticks = element_line(color = "black", size = 0.2),
    axis.line = element_line(color = "black", size = 0.3),
    axis.title = element_text(size = 8),
    strip.background = element_blank(),
    strip.text.y = element_text(size = 7, angle = 0),
    plot.title = element_text(size = 8, face = "bold", hjust = 0.5),
    legend.position = "none"
  )


# Lock facet order
lname<-c("PC","TC5","FM","TC7","217","119")

# Per-cultivar highlight ranges on chr05
highlight_region <- data.frame(
  CHROM    = "chr05",
  Cultivar = c("FM", "TC7", "217", "119"),
  xmin     = c(956425, 1, 1, 1),
  xmax     = c(7030674, 13290586, 21500000, chr05_length),
  stringsAsFactors = FALSE
)

# Join and subset to the per-cultivar POS ranges
highlight_sub <- dplyr::inner_join(All_melt, highlight_region, by = c("CHROM","Cultivar")) %>%
  dplyr::filter(POS >= xmin, POS <= xmax)

# Reapply factor levels so join can't mess with order
highlight_sub$Cultivar <- factor(highlight_sub$Cultivar, levels = lname)

# --- Final Density Plot (Highlighted Cultivars Only) ---
highlight_density <- ggplot(All_melt, aes(x = value)) +
  
  geom_density(fill = "navy", color = NA, alpha = 0.8) +
  geom_density(
    data = highlight_sub,
    aes(x = value),
    fill = "#e1171e", color = NA, alpha = 0.4
  ) +
  coord_flip() +
  facet_grid(Cultivar ~ .) +  
  scale_x_continuous(
    breaks = c(0, 0.33, 0.50, 0.67, 1),
    labels = c("", "", "", "", ""),
    expand = expansion(mult = c(0.1,0.1))
  ) +
  scale_y_continuous(limits = c(0, 4),expand = expansion(mult = c(0,NA))) +
  xlab(NULL) + ylab("Density") +
  theme_minimal(base_size = 8) +
  theme(
    panel.background = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.text.x = element_text(angle = 0, hjust = 0.5, size = 6),
    axis.text.y = element_blank(),
    axis.line.y = element_blank(),
    axis.ticks.x = element_line(color = "black", size = 0.2),
    axis.line.x = element_line(color = "black", size = 0.3),
    axis.title = element_text(size = 8),
    strip.background = element_blank(),
    strip.text.x = element_blank(),
    strip.text.y = element_text(size = 7,angle = 0),
    plot.title = element_text(size = 8, face = "bold", hjust = 0.5),
    legend.position = "none"
  )


# --- 3. Combine using patchwork ---


final_plot <- geom_plot+theme(strip.text.y.right = element_blank()) + 
  highlight_density +
  plot_layout(ncol = 2, widths = c(4, 0.6))

PDF_FILE="Fig3C.chr5_geom_density.pdf"
ggsave(PDF_FILE, plot = final_plot, width = 15, height = 10, dpi = 300, units = "cm", paper ="a4")
shell.exec(PDF_FILE)
