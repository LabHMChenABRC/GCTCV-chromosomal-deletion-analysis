# ============================================================
# Triple Synteny Plot with DH Pahang in Middle
# ============================================================

library(data.table)
library(dplyr)
library(ggplot2)
library(cowplot)
library(stringr)


# --- Input paths ---
gff_path_ban <- "E:/Visualization/S14/ban_FM05.gff"
gff_path_dh  <- "E:/Visualization/S14/DH_FM05.gff"
gff_path_bxh <- "E:/Visualization/S14/chr05.2_FM05.gff"

collinearity_ban_dh <- "E:/Visualization/S14/ban_FM05.collinearity"
collinearity_bxh_dh <- "E:/Visualization/S14/bxh2_FM.collinearity"

# --- Load GFF files ---
load_gff <- function(path) {
  gff <- fread(path, header = FALSE)
  colnames(gff) <- c("Chr", "Gene", "Start", "End")
  gff %>% mutate(Mid = (Start + End)/2)
}

gff_ban <- load_gff(gff_path_ban)
gff_dh  <- load_gff(gff_path_dh)
gff_bxh <- load_gff(gff_path_bxh)

# --- Parse collinearity file with consistent orientation ---
parse_collinearity <- function(path, tagA, tagB) {
  lines <- readLines(path)
  blocks <- data.frame(Block = integer(), Gene1 = character(), Gene2 = character(), stringsAsFactors = FALSE)
  block_id <- NA
  for (ln in lines) {
    if (grepl("^## Alignment", ln)) {
      block_id <- as.integer(str_match(ln, "## Alignment\\s+(\\d+):")[,2])
    } else if (grepl("Musa|Macma", ln)) {
      parts <- str_match(ln, "^[^:]+:\\s+(\\S+)\\s+(\\S+)")
      if (!is.na(parts[1,2]) && !is.na(parts[1,3])) {
        # Force consistent orientation: Gene1 = tagA genome, Gene2 = tagB genome
        if (tagA == "dh") {  # DH should always be second
          g1 <- parts[1,3]
          g2 <- parts[1,2]
        } else {
          g1 <- parts[1,2]
          g2 <- parts[1,3]
        }
        blocks <- rbind(blocks, data.frame(
          Block = block_id,
          Gene1 = g1,
          Gene2 = g2,
          stringsAsFactors = FALSE
        ))
      }
    }
  }
  blocks$Pair <- paste(tagA, tagB, sep = "_")
  blocks
}

blocks_ban_dh <- parse_collinearity(collinearity_ban_dh, "ban", "dh")
blocks_bxh_dh  <- parse_collinearity(collinearity_bxh_dh, "bxh", "dh")

# --- Merge with gene coordinates ---
merge_blocks <- function(blocks, gff1, gff2, tagA, tagB) {
  blocks <- merge(blocks, gff1[, .(Gene, ChrA = Chr, PosA = Mid, StartA = Start, EndA = End)],
                  by.x = "Gene1", by.y = "Gene", all.x = TRUE)
  blocks <- merge(blocks, gff2[, .(Gene, ChrB = Chr, PosB = Mid, StartB = Start, EndB = End)],
                  by.x = "Gene2", by.y = "Gene", all.x = TRUE)
  blocks$TagA <- tagA
  blocks$TagB <- tagB
  blocks
}

links_ban_dh <- merge_blocks(blocks_ban_dh, gff_ban, gff_dh, "ban", "dh")
links_bxh_dh <- merge_blocks(blocks_bxh_dh, gff_dh, gff_bxh, "dh", "bxh")


# ============================================================
# Prepare plotting data
# ============================================================

# Chromosome ranges
chr_ban <- gff_ban %>% summarise(xmin = min(Start)/1e6, xmax = max(End)/1e6) %>% mutate(y = 3)
chr_dh  <- gff_dh  %>% summarise(xmin = min(Start)/1e6, xmax = max(End)/1e6) %>% mutate(y = 2)
chr_bxh <- gff_bxh %>% summarise(xmin = min(Start)/1e6, xmax = max(End)/1e6) %>% mutate(y = 1)

# Gene models
genes_ban <- gff_ban %>% mutate(y = 3, Start = Start/1e6, End = End/1e6)
genes_dh  <- gff_dh  %>% mutate(y = 2, Start = Start/1e6, End = End/1e6)
genes_bxh <- gff_bxh %>% mutate(y = 1, Start = Start/1e6, End = End/1e6)

# Links
links_ban_dh <- links_ban_dh %>% mutate(y1 = 3, y2 = 2, PosA = PosA/1e6, PosB = PosB/1e6)
links_bxh_dh <- links_bxh_dh %>% mutate(y1 = 1, y2 = 2, PosA = PosA/1e6, PosB = PosB/1e6)

# ============================================================
# Plot
# ============================================================

p <- ggplot() +
  # Synteny links
  geom_segment(data = links_ban_dh,
               aes(x = PosA, xend = PosB, y = y1, yend = y2),
               alpha = 0.5, color = "grey40", size = 0.05) +
  geom_segment(data = links_bxh_dh,
               aes(x = PosB, xend = PosA, y = y1, yend = y2),
               alpha = 0.5, color = "grey40", size = 0.05) +
  geom_rect(data = chr_ban, aes(xmin = xmin, xmax = xmax, ymin = y-0.08, ymax = y+0.08),
            fill = "#c6dbef", color = NA) +
  geom_rect(data = chr_dh, aes(xmin = xmin, xmax = xmax, ymin = y-0.08, ymax = y+0.08),
            fill = "grey", color = NA) +
  geom_rect(data = chr_bxh, aes(xmin = xmin, xmax = xmax, ymin = y-0.08, ymax = y+0.08),
            fill = "#c7e9c0", color = NA) +
  
  # Gene models
  geom_rect(data = genes_ban, aes(xmin = Start, xmax = End, ymin = y-0.05, ymax = y+0.05), fill = "darkblue") +
  geom_rect(data = genes_dh, aes(xmin = Start, xmax = End, ymin = y-0.05, ymax = y+0.05), fill = "black") +
  geom_rect(data = genes_bxh, aes(xmin = Start, xmax = End, ymin = y-0.05, ymax = y+0.05), fill = "darkgreen") +
  
  
  scale_y_continuous(breaks = c(1,2,3),
                     labels = c("Cavendish chr0.2", "'DH Pahang' chr05", "Cavendish Ban05")) +
  labs(x = "Genomic position (Mb)", y = "") +
  theme_cowplot() +
  theme(
    axis.text.y = element_text(size = 10, face = "bold"),
    axis.title.x = element_text(size = 10),
    panel.grid = element_blank(),
    legend.position = "none"
  )

ggsave("synteny_triple_DHmiddle.pdf", plot = p, width = 20, height = 6, units = "cm", paper ="a4")

shell.exec("synteny_triple_DHmiddle.pdf")




