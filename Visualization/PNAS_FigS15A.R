# ============================================================
# Chromosome 5 Synteny Visualization (genome alignment)

# This script visualizes synteny relationships on Chromosome 5
# between DH Pahang, Cavendish, and BXH assembly.


# ============================================================
# Author: Puyam Tondonba Singh
# Date: 2025-08-18



# ==== Load Libraries ====
library(data.table)
library(ggplot2)
library(ggforce)

setwd("E:/tony/WGS/horti_baxi_bw10000/baxixiao protein fasta/linkview")

# ==== Load Alignment File ====
align <- fread("alignment.txt", header = FALSE)
setnames(align, c("RefChr", "RefStart", "RefEnd", "QryChr", "QryStart", "QryEnd"))

# ==== Load Karyotype File ====
karyo <- fread("karyotype.txt", header = FALSE)
setnames(karyo, "V1")
karyo[, c("Chr", "Start", "End") := tstrsplit(V1, ":", fixed = TRUE)]
karyo[, `:=`(Start = as.integer(Start), End = as.integer(End))]
karyo[, Length := End - Start + 1]
karyo[, LengthMb := Length / 1e6]

# ==== Add Pretty Labels (BXH_chr05.1 REMOVED) ====
label_map <- c(
  "Pahang_chr05"    = "DH Pahang chr05",
  "Cavendish_chr05" = "Li et al. 2023 — Ban05 FM del",
  "BXH_chr05.2"     = "Huang et al. 2023 — chr05.2 FM del"
)
karyo[, Label := label_map[Chr]]

# Keep only chromosomes we are plotting
karyo <- karyo[Chr %in% names(label_map)]

# ==== Convert Alignment Positions to Mb ====
align[, `:=`(
  RefStartMb = RefStart / 1e6,
  RefEndMb   = RefEnd / 1e6,
  RefMidMb   = (RefStart + RefEnd) / 2 / 1e6,
  QryStartMb = QryStart / 1e6,
  QryEndMb   = QryEnd / 1e6,
  QryMidMb   = (QryStart + QryEnd) / 2 / 1e6
)]

# Drop any rows involving BXH_chr05.1 from align, and keep only target chroms
align <- align[!(RefChr %in% "BXH_chr05.1" | QryChr %in% "BXH_chr05.1")]
align <- align[RefChr %in% names(label_map) & QryChr %in% names(label_map)]

# ==== Assign Y Positions ====
y_positions <- data.table(
  Chr = names(label_map),
  Y = c(0, -1, 1.3)  # DH = 0, Cavendish = -1, chr05.2 = 1.3
)

# ==== Add X Offsets for Horizontal Shift ====
x_offsets <- data.table(
  Chr = names(label_map),
  X_offset = c(0, 0, 0)  # shift BXH_chr05.2 for spacing
)

# ==== Merge Y and X offsets with karyotype ====
bar_data <- merge(karyo, y_positions, by = "Chr", all.x = TRUE)
bar_data <- merge(bar_data, x_offsets, by = "Chr", all.x = TRUE)
bar_data[, `:=`(
  xstart = 0,
  xend = LengthMb,
  xstart_shifted = 0 + X_offset,
  xend_shifted = LengthMb + X_offset
)]

# ==== Merge Y and X offsets with alignment ====
align <- merge(align, y_positions, by.x = "RefChr", by.y = "Chr", all.x = TRUE)
setnames(align, "Y", "Y_Ref")
align <- merge(align, y_positions, by.x = "QryChr", by.y = "Chr", all.x = TRUE)
setnames(align, "Y", "Y_Qry")
align <- merge(align, x_offsets, by.x = "RefChr", by.y = "Chr", all.x = TRUE)
setnames(align, "X_offset", "X_RefOffset")
align <- merge(align, x_offsets, by.x = "QryChr", by.y = "Chr", all.x = TRUE)
setnames(align, "X_offset", "X_QryOffset")

# ==== Prepare Bezier Curves ====
bezier_data <- align[, .(
  x = c(RefMidMb + X_RefOffset, 
        RefMidMb + X_RefOffset + 1.5 * sign(QryMidMb - RefMidMb), 
        QryMidMb + X_QryOffset - 1.5 * sign(QryMidMb - RefMidMb), 
        QryMidMb + X_QryOffset),
  y = c(Y_Ref, 
        Y_Ref + 0.6 * sign(Y_Qry - Y_Ref), 
        Y_Qry - 0.6 * sign(Y_Qry - Y_Ref), 
        Y_Qry),
  group = .I,
  QryChr = QryChr
), by = .I]

# ==== Highlight Region on DH Pahang ====
highlight_bar <- data.table(
  Chr = "Pahang_chr05",
  xstart = 956425 / 1e6,
  xend   = 7030674 / 1e6,
  fill   = "highlight"
)
highlight_bar <- merge(highlight_bar, bar_data[, .(Chr, Y, X_offset)], by = "Chr")
highlight_bar[, `:=`(
  xstart = xstart + X_offset,
  xend = xend + X_offset,
  ymin = Y - 0.05,
  ymax = Y + 0.05
)]

# ==== Plot ====
a <- ggplot() +
  # Bezier synteny links
  geom_bezier(
    data = bezier_data,
    aes(x = x, y = y, group = group, color = QryChr),
    size = 0.1, alpha = 0.6, lineend = "round") +
  # Chromosome bars
  geom_segment(
    data = bar_data,
    aes(x = xstart_shifted, xend = xend_shifted, y = Y, yend = Y, color = Chr),
    size = 2
  ) +
  # Bar caps
  geom_rect(
    data = bar_data,
    aes(xmin = xstart_shifted, xmax = xend_shifted, ymin = Y - 0.03, ymax = Y + 0.03),
    fill = "white", color = "black", size = 0.6
  ) +
  # Chromosome labels
  geom_text(
    data = bar_data,
    aes(x = xend_shifted + 0.3, y = Y, label = Label),
    hjust = 0, size = 4.5
  ) +
  # X-axis with 2 Mb scale (labels won’t show with theme_void, but kept for scale bar sanity)
  scale_x_continuous(name = NULL, breaks = seq(0, 50, 2), expand = c(0, 0)) +
  # Y-axis hidden
  scale_y_continuous(name = NULL, breaks = NULL, expand = expansion(0.2)) +
  # Custom colors (only remaining chroms)
  scale_color_manual(values = c(
    "Pahang_chr05"    = "#EBCB8B",
    "Cavendish_chr05" = "#B48EAD",
    "BXH_chr05.2"     = "#4C566A"
  )) +
  # Highlighted DH region
  geom_rect(
    data = highlight_bar,
    aes(xmin = xstart, xmax = xend, ymin = ymin, ymax = ymax),
    fill = "yellow", color = NA, alpha = 0.5
  ) +
  # Coordinate settings
  coord_cartesian(xlim = c(min(bar_data$xstart_shifted) - 1, max(bar_data$xend_shifted) + 3)) +
  # Theme and annotations
  theme_void() +
  theme(
    plot.margin = margin(10, 40, 10, 10),
    legend.position = "none"
  ) +
  # Add scale bar
  annotate("segment", x = 1, xend = 3, y = 2.5, yend = 2.5, size = 0.8, color = "black") +
  annotate("text", x = 2, y = 2.65, label = "2 Mb", size = 3.5)

# ==== Save ====
ggsave("chr05gene_synteny1.pdf", plot = a, width = 12, height = 15, units = "cm", paper ="a4")
