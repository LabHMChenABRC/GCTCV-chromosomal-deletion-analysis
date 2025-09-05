# ============================================================
# Build Custom OrgDb Package for Musa acuminata (Banana)
# ============================================================
# Author: Puyam Tondonba Singh
# Date: 2025-08-25
#
# This script creates a custom organism annotation package
# (OrgDb) using gene-to-GO mappings for Musa acuminata.
#
# Required: AnnotationForge, GO.db

if (!requireNamespace("AnnotationForge", quietly = TRUE)) {
  BiocManager::install("AnnotationForge")
}

library(data.table)
library(dplyr)
library(GO.db)
library(AnnotationForge)

# --- GO annotation file ---
file_path <- "E:/tony/VIGS/Transcriptome/Annotation_v4/Musa_acuminata_pahang_v4_go.txt"

# Read GO mapping
gene2go <- read.delim(file_path, header = FALSE, stringsAsFactors = FALSE)
colnames(gene2go) <- c("GID", "GO")
gene2go$GID <- sub("\\.1$", "", gene2go$GID)

# Add evidence column
gene2go_df <- gene2go %>%
  mutate(EVIDENCE = "IEA")

# --- Build OrgDb source package ---
makeOrgPackage(
  go         = gene2go_df,
  version    = "0.1",
  maintainer = "Tony <tonpuyam@gmail.com>",
  author     = "Tony",
  outputDir  = "E:/tony/HydroTC5_RNAseq/RSEM/MIX",
  tax_id     = "4641",   # NCBI taxonomy ID for Musa acuminata
  genus      = "Musa",
  species    = "acuminata",
  goTable    = "go"
)


# --- Install the OrgDb you just built ---
install.packages("E:/tony/HydroTC5_RNAseq/RSEM/MIX/org.Macuminata.eg.db", 
                 repos = NULL, type = "source")

# --- Test loading ---
library(org.Macuminata.eg.db)

