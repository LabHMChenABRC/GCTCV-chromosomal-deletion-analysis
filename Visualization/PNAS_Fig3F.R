# ==============================
# Heatmap of Variant Allele Ratios (VAR) from VCF 'Formosana Deletion region' Chromosome 5

# This script generates a heatmap of variant allele ratios (VAR) 
# for FM deletion region on Chromosome 5 from whole-genome sequencing (WGS) data.
# It parses a VCF, filters high-quality SNPs, identifies heterozygous sites based on "Pei- Chiao" 
# calculates variant allele ratios for all samples, 
# and visualizes them as a heatmap.
# ============================================================
# Author: Puyam Tondonba Singh
# Date: 2025-08-18

# ==============================

# --- Load libraries ---
library(data.table)
library(dplyr)
library(stringr)
library(ComplexHeatmap)
library(circlize)
library(grid)

# --- INPUT ---
vcf_path <- "E:/tony/WGS/vcf_chr5/genotypeChr5_filter.vcf/snps_raw.chr05.HF.vcf"

# --- Helper: Parse samples from VCF ---
read_vcf_table <- function(vcf_path) {
  vcf <- fread(vcf_path, skip = "#CHROM", header = TRUE, blank.lines.skip = TRUE)
  setnames(vcf, gsub("^#", "", names(vcf)))  # clean headers
  sample_cols <- setdiff(names(vcf), c("CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT"))
  message("Detected samples: ", paste(sample_cols, collapse = ", "))
  vcf <- vcf[, c("CHROM", "POS", "REF", "ALT", "FILTER", sample_cols), with = FALSE]
  return(list(vcf = vcf, samples = sample_cols))
}

parsed <- read_vcf_table(vcf_path)
VCF <- parsed$vcf
samples <- parsed$samples

# --- Filter SNPs: Single base, quality PASS or SnpCluster ---
VCF <- VCF[
  nchar(REF) == 1 & nchar(ALT) == 1 & FILTER %in% c("PASS", "SnpCluster")
]

# --- Melt into long format ---
vcf_long <- melt(
  VCF,
  id.vars = c("CHROM", "POS", "REF", "ALT", "FILTER"),
  measure.vars = samples,
  variable.name = "Sample",
  value.name = "GenotypeField"
)

# --- Extract GT, AD, DP ---
vcf_long[, c("GT", "AD", "DP") := tstrsplit(GenotypeField, ":", fixed = TRUE)[1:3]]
vcf_long[, c("R_AD", "A1_AD") := tstrsplit(AD, ",", fixed = TRUE)]
vcf_long[, c("DP", "R_AD", "A1_AD") := lapply(.SD, as.numeric), .SDcols = c("DP", "R_AD", "A1_AD")]

# --- Identify heterozygous sites in PC with DP > 20 ---
het_GT_patterns <- c("0/1/1", "0/0/1")
PC_hets <- vcf_long[Sample == "PC" & GT %in% het_GT_patterns & DP > 20, .(CHROM, POS)]

# --- Filter to those positions for all samples ---
vcf_long <- merge(vcf_long, PC_hets, by = c("CHROM", "POS"))
vcf_long[, DP_sum := R_AD + A1_AD]
# --- Apply per-sample DP filter (DP > 10) ---
vcf_long <- vcf_long[DP_sum > 10]

# --- Calculate Variant Allele Ratio ---
vcf_long[, alleleratio := A1_AD / DP]

# --- Output formatted VAR matrix ---
GT_final <- vcf_long[, .(CHROM, POS, Sample, alleleratio)]
GT_filt <- GT_final %>% filter(POS > 956425 & POS < 7030674)
GT_wide <- dcast(GT_filt, CHROM + POS ~ Sample, value.var = "alleleratio") %>% na.omit()

# --- Build matrix and transpose for heatmap ---
mat <- as.matrix(GT_wide[, -c("CHROM", "POS")])
rownames(mat) <- as.character(GT_wide$POS)
mat_transposed <- t(mat)

# --- Heatmap Settings ---
target_positions <- c("966452", "966451")
col_fun <- colorRamp2(c(0, 0.5, 1), c("#36428f", "white", "#a63232"))

# Selective column labels
custom_labels <- ifelse(colnames(mat_transposed) %in% target_positions,
                        format(as.numeric(colnames(mat_transposed)), big.mark = ","),
                        "")

# --- Plot ---
pdf("250725hotpotNEW.pdf", width = 9, height = 6, paper = "a4")

Heatmap(
  mat_transposed,
  col = col_fun,
  cluster_columns = TRUE,
  cluster_rows = FALSE,
  show_column_dend = FALSE,
  na_col = NA,
  column_names_side = "bottom",
  column_labels = custom_labels,
  column_names_gp = gpar(fontsize = 8),
  row_names_gp = gpar(fontsize = 8),
  show_row_names = FALSE,
  show_column_names = TRUE,
  heatmap_legend_param = list(
    at = c(0, 0.33, 0.5, 0.67, 1),
    labels = c("0 (LOH)", "0.33", "0.5", "0.67", "1 (LOH)"),
    title = "Variant allele ratio",
    legend_height = unit(4, "cm"),
    title_position = "leftcenter-rot"
  ),
  row_split = factor(
    rownames(mat_transposed),
    levels = c("PC", "TC5-2", "TC4", "TC7", "GCTCV217", "GCTCV119_L3")
  ),
  use_raster = TRUE
)

dev.off()

