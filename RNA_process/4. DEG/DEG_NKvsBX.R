# ============================================================
# DEG analysis 'Nongke No. 1' VS 'Baxi'
#
# This script performs DEG analysis for 'Nongke No. 1' VS 'Baxi'
# We load the GTF annotation here to map each gene ID to its genomic
# coordinates (chromosome, start, end, strand). This is necessary because:
#   - We want to filter DEGs located specifically in the chr05 deletion region
#     (region_start–region_end).
#   - We use this to generate chromosome-level plots where DEGs are placed along
#     their genomic positions.

# Author: Puyam Tondonba Singh
# Date: 2025-08-18
# ============================================================

library(data.table)
library(DESeq2)
library(dplyr)

# --- CONFIGURATION ---
setwd("E:/tony/NK1 RNA seq/RSEM")   # project root
gtf_file <- "Musa_acuminata_pahang_v4.CorrectedMaACA8.gtf"
chr_target <- "chr05"
region_start <- 956425
region_end <- 7030674
time_points <- list("27 H" = "27", "52 H" = "51")

# --- IMPORT GTF ---
gtf_data <- import(gtf_file)
gene_data <- gtf_data[gtf_data$type == "transcript", ]
gene_info <- data.frame(
  Gene_ID = mcols(gene_data)$gene_id,
  Chromosome = as.character(seqnames(gene_data)),
  Start = start(gene_data),
  End = end(gene_data),
  Strand = as.character(strand(gene_data))
)

# --- FUNCTIONS ---
read_counts <- function(file_path, prefix) {
  NK_files <- sprintf("%s/NKCK%sh%d.genes.results", file_path, prefix, 1:3)
  BX_files <- sprintf("%s/BXCK%sh%d.genes.results", file_path, prefix, 1:3)
  files <- c(NK_files, BX_files)
  
  missing <- files[!file.exists(files)]
  if (length(missing) > 0) {
    warning("Missing files: ", paste(missing, collapse = ", "))
    return(NULL)
  }
  
  read_and_select <- function(file) {
    fread(file, sep = "\t")[, .(gene_id, expected_count)][, expected_count := as.integer(expected_count)]
  }
  
  data_list <- lapply(files, read_and_select)
  names(data_list) <- sub(".genes.results", "", basename(files))
  all_data <- rbindlist(data_list, idcol = "sample")
  
  dt <- dcast(all_data, gene_id ~ sample, value.var = "expected_count")
  mat <- as.matrix(dt[,-1])
  rownames(mat) <- dt$gene_id
  mat <- mat[rowSums(mat > 5) >= ncol(mat)/2, ]
  return(mat)
}

build_sample_info <- function(count_matrix) {
  sample_info <- data.frame(sampleName = colnames(count_matrix))
  sample_info$batch <- factor(sub(".*h", "", sample_info$sampleName))
  sample_info$condition <- factor(sub("h.*", "", sample_info$sampleName))
  rownames(sample_info) <- sample_info$sampleName
  return(sample_info)
}

run_deseq_lrt <- function(mat, sample_info) {
  dds <- DESeqDataSetFromMatrix(countData = round(mat), colData = sample_info, design = ~ batch + condition)
  dds <- DESeq(dds, test = "LRT", reduced = ~batch)
  res <- results(dds)
  res <- as.data.frame(res)
  res$Gene_ID <- rownames(res)
  return(res)
}

classify_deg <- function(df) {
  df %>% mutate(
    DEG_status = ifelse(padj < 0.1, "DEG", "Not DEG"),
    color = case_when(
      DEG_status == "DEG" & log2FoldChange > 0.5 ~ "Upregulated",
      DEG_status == "DEG" & log2FoldChange < -0.5 ~ "Downregulated",
      TRUE ~ "Unchanged"
    )
  )
}

filter_by_region <- function(df) {
  df %>% filter(Chromosome == chr_target & Start >= region_start & End <= region_end)
}

safe_count <- function(obj, key) {
  val <- obj[[key]]
  if (is.null(val) || is.na(val)) return(0)
  return(val)
}

# --- OUTPUT DIRECTORY ---
output_dir <- file.path(getwd(), "DEG_results")
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

# --- ANALYSIS LOOP ---
deg_counts <- list()
deg_tables <- list()

for (tp_label in names(time_points)) {
  tp_prefix <- time_points[[tp_label]]
  file_path <- file.path(getwd(), tp_label)
  
  cat("Processing:", tp_label, "\n")
  if (!dir.exists(file_path)) {
    warning("Directory not found:", file_path)
    next
  }
  
  mat <- read_counts(file_path, tp_prefix)
  if (is.null(mat)) next
  
  sample_info <- build_sample_info(mat)
  res <- run_deseq_lrt(mat, sample_info)
  res <- classify_deg(res)
  combined <- merge(res, gene_info, by = "Gene_ID", all.x = TRUE)
  
  # save DEG table as TSV
  out_file <- file.path(output_dir, paste0("DEG_", gsub(" ", "", tp_label), ".tsv"))
  fwrite(combined, file = out_file, sep = "\t")
  cat("Saved:", out_file, "\n")
  
  deg_tables[[tp_label]] <- combined
  
  deg_all <- combined %>% filter(color %in% c("Upregulated", "Downregulated"))
  deg_region <- filter_by_region(combined) %>% filter(color %in% c("Upregulated", "Downregulated"))
  
  deg_counts[[tp_label]] <- list(
    all = table(deg_all$color),
    region = table(deg_region$color)
  )
}



