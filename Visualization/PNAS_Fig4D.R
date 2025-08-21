# ============================================================
# Chromosome 5 Allele Ratio Visualization based on WGS

# This script visualizes allele ratios along Chromosome 5 
# from WGS-derived VCFs across multiple cultivars.
# It parses VCF files, filters heterozygous sites based on "Pei-Chiao", computes allele ratios
# visualize using point and density visualizations

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

# --- Helper: Parse VCF Function ---
read_vcf_table  <- function(File, Samples) {
  rbindlist(
    lapply(File, function(x) {
      SampleName <- sub("(.+?)_.+", "\\1", basename(x))
      fread(x, header = TRUE, blank.lines.skip = TRUE) %>%
        setnames(sub("#", "", colnames(.))) %>%
        .[, Marker := SampleName] %>%
        .[, .SD, .SDcols = c("CHROM", "POS", "REF", "ALT", "FILTER", "Marker", Samples)]
    })
  )[, .(Marker = paste0(unique(Marker), collapse = ",")), by = c("CHROM", "POS", "REF", "ALT", "FILTER", Samples)]
}

# --- Set Inputs ---
#setwd("YOUR/WORKING/DIRECTORY")

VCF_File <- "E:/tony/Marker_RNASeq/GCTCV/2025/filtered_variants_chr05.vcf/filtered_variants.vcf"

cultivars <- c(
  paste0("PC-", 1:4),
  paste0("TC5-", 1:4),
  paste0("TC7-", 1:2),
  paste0("TC4-", 1:3),
  paste0("GCTCV119-", 1:2),
  paste0("GCTCV217-", 1:2),
  paste0("BXCK27h", 1:3),
  paste0("NKCK27h", 1:3)
)

# --- Parse and Filter VCF ---
VCF <- read_vcf_table (VCF_File, cultivars)[
  nchar(REF) == 1 & nchar(ALT) == 1
][
  , unique(.SD)
][
  order(CHROM, POS)
]

genotype<-VCF[FILTER %in% c("PASS","SnpCluster")
][,eval(cultivars):=lapply(.SD,function(x) tstrsplit(x,":",fixed=TRUE,keep = 1L)[[1]]),.SDcols=cultivars
][order(CHROM,POS),
][,.SD,.SDcols=c("Marker","CHROM","POS","REF","ALT","FILTER",eval(cultivars))]


VCF_filtered <- genotype[
  (`PC-1` %in% c("0/1/1", "0/0/1")) &
    (`PC-2` %in% c("0/1/1", "0/0/1")) &
    (`PC-3` %in% c("0/1/1", "0/0/1")) &
    (`PC-4` %in% c("0/1/1", "0/0/1"))
]

# --- Melt Genotypes and Extract Fields ---
vcf_long <- melt(
  VCF[FILTER %in% c("PASS", "SnpCluster")],
  id.vars = c("Marker", "CHROM", "POS", "REF", "ALT", "FILTER"),
  measure.vars = cultivars,
  variable.name = "Sample",
  value.name = "GenotypeField"
)

vcf_long[, c("GT", "AD", "DP", "Field4", "Rest") := tstrsplit(GenotypeField, ":", fixed = TRUE)]
vcf_long[, c("R_AD", "A1_AD") := tstrsplit(AD, ",", fixed = TRUE)]
vcf_long[, c("DP", "R_AD", "A1_AD") := lapply(.SD, as.numeric), .SDcols = c("DP", "R_AD", "A1_AD")]
vcf_long[, c("AD", "Field4", "Rest") := NULL]

vcf_long[, DP_sum := R_AD + A1_AD]

# --- Filter and Compute Allele Ratio ---
GT <- vcf_long[DP_sum > 10]
# Check if DP equals the sum of R_AD and A1_AD

GT[, correct := DP == (R_AD + A1_AD)]
# Count how many are incorrect
incorrect_count <- sum(!GT$correct)
total_sites <- nrow(GT)



cat("Incorrect sites:", incorrect_count, "out of", total_sites, "\n")
cat("Percentage incorrect:", (incorrect_count / total_sites) * 100, "%\n")



GT[, alleleratio := A1_AD / (R_AD + A1_AD)]

GT_final <- GT[, .(CHROM, POS, Sample, alleleratio)]

# --- Plot Density (All Samples) ---
dENSITY_ALL <- ggplot(GT_final, aes(x = alleleratio)) +
  facet_wrap(~Sample, ncol = 4, scales = "free") +
  geom_density(fill = "black", color = "#e9ecef", alpha = 0.8) +
  theme_bw() +
  theme(
    panel.border = element_blank(), panel.grid = element_blank(),
    axis.line = element_line(colour = "black"),
    text = element_text(size = 10, color = "black"),
    strip.text = element_text(size = 10, color = "black"),
    panel.spacing = unit(1, "lines")
  ) +
  labs(x = "Heterozygote allele ratio", y = "Density") +
  scale_y_continuous(breaks = c(0, 1, 2, 2.5), labels = c("0", "1", "2", "2.5"))

ggsave("dENSITY_ALL.pdf", plot = dENSITY_ALL, width = 300, height = 200, dpi = 300, units = "mm")

# --- Cast to Wide Format ---
Cultivars_AR <- dcast(GT_final, CHROM + POS ~ Sample, value.var = "alleleratio") %>% na.omit()

# --- Semi-Join to Filtered Variants ---
Cultivars_AR <- semi_join(Cultivars_AR, VCF_filtered, by = "POS")

# --- Compute Mean Per Group ---
groups <- list(
  PC = paste0("PC-", 1:4),
  TC5 = paste0("TC5-", 1:4),
  TC7_1 = "TC7-1", TC7_2 = "TC7-2",
  TC4 = paste0("TC4-", 1:3),
  GCTCV119 = paste0("GCTCV119-", 1:2),
  GCTCV217 = paste0("GCTCV217-", 1:2),
  NKCK27h = paste0("NKCK27h", 1:3),
  BXCK27h = paste0("BXCK27h", 1:3)
)

meta_cols <- c("CHROM", "POS")
mean_expr <- lapply(groups, function(cols) rowMeans(Cultivars_AR[, ..cols], na.rm = TRUE))
Averaged_AR <- cbind(Cultivars_AR[, ..meta_cols], as.data.table(mean_expr))

# --- Final Table for Plotting ---
final_AR <- Averaged_AR[, .(CHROM, POS, PC, TC5, TC4, NKCK27h, BXCK27h)]
All_melt <- melt(final_AR, id.vars = c("CHROM", "POS"), variable.name = "Cultivar") %>% na.omit()
All_melt$Cultivar <- factor(All_melt$Cultivar, levels = c("PC", "BXCK27h","TC5",  "TC4", "NKCK27h"))

# --- Highlight Region for Plot ---
highlight_region <- data.frame(
  CHROM = c("chr05", "chr05"),
  Cultivar = c("TC4", "NKCK27h"),
  xmin = c(956425, 956425),
  xmax = c(7030674, 7030674),
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


ggsave("250727_chr5_geom.pdf", plot = geom_plot, width = 95, height = 100, dpi = 300, units = "mm")

# --- Final Density Plot (Highlighted Cultivars Only) ---
highlight_region <- filter(All_melt, POS >= 956425 & POS <= 7030674)

highlight_density <- ggplot(All_melt, aes(x = value)) +
  facet_grid(Cultivar ~ ., switch = "y") +  # Ensure vertical alignment
  geom_density(fill = "navy", color = NA, alpha = 0.8) +
  geom_density(
    data = filter(highlight_region, Cultivar %in% c("TC4", "NKCK27h")),
    aes(x = value),
    fill = "pink", color = NA, alpha = 0.5
  ) +
  coord_flip() +  # Reflect horizontally
  scale_x_continuous(
    limits = c(0, 1),
    breaks = c(0, 0.33, 0.50, 0.67, 1),
    labels = c("", "", "", "", "")
  ) +
  scale_y_continuous(limits = c(0, 4)) +  # Only once
  xlab("Heterozygote allele ratio") +
  ylab("Density") +
  theme_minimal(base_size = 8) +
  theme(
    panel.background = element_rect(fill = "grey95", colour = NA),
    panel.grid = element_blank(),
    
    axis.text.x = element_text(size = 6),       # restore x labels
    axis.ticks.x = element_line(color = "black", size = 0.2),
    axis.line.x = element_line(color = "black", size = 0.3),
    
    axis.text.y = element_blank(),             
    axis.ticks.y = element_blank(),
    axis.line.y = element_blank(),
    
    axis.title.y = element_blank(),
    axis.title.x = element_text(size = 8),
    
    # Strip layout
    strip.background = element_blank(),
    strip.placement = "outside",
    strip.text.y.left = element_blank(),
    strip.text.y.right = element_text(size = 7, angle = 0),
    
    legend.position = "none"
  )


# --- 3. Combine using patchwork ---
final_plot <- geom_plot + highlight_density +
  plot_layout(ncol = 2, widths = c(4, 0.6), guides = "collect")


ggsave("250815_chr5_geom_density.pdf", plot = final_plot, width = 10, height = 7, dpi = 300, units = "cm", paper ="a4")
getwd()





















