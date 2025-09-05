# ============================================================
# Heatmap of Deletion Region Genes (PC vs TC5)
# Author: Puyam Tondonba Singh
# Date: 2025-08-25
#
# This script generates a heatmap of gene expression (TPM) 
# for genes located in a defined deletion region. 
#

# ======================================================
library(data.table)
library(ComplexHeatmap)
library(circlize)

file_path <- "E:/Visualization/2D"
setwd(file_path)

# ===== Define gene list =====
TC5_DEL_Gene <- c("MaACA7"         ="Macma4_04_g41380", 
                  "Uncharacterized"="Macma4_04_g41390", 
                  "MaAPT1"         ="Macma4_04_g41400", 
                  "MaLBD"          ="Macma4_04_g41410", 
                  "MaP4H13"        ="Macma4_04_g41420", 
                  "Uncharacterized"="Macma4_04_g41430", 
                  "Uncharacterized"="Macma4_04_g41440")
ID2Name<-setNames(names(TC5_DEL_Gene),nm = TC5_DEL_Gene)

# ===== Load all files =====
all_files <- list.files(pattern = ".genes.results$")

read_and_select <- function(filename) {
  dt <- fread(filename, sep = "\t")[, .(gene_id, TPM)]
  dt[, TPM := as.numeric(TPM)]
  sample_name <- sub(".genes.results", "", basename(filename))
  setnames(dt, "TPM", sample_name)
  return(dt)
}

all_data.list <- lapply(all_files, read_and_select)
TPM.dt <- Reduce(function(x, y) merge(x, y, by = "gene_id", all = TRUE), all_data.list)

# ===== Subset genes expression =====
TPM_deletion.dt <- TPM.dt[gene_id %in% TC5_DEL_Gene]
TPM_deletion.dt[,gene_name:=ID2Name[gene_id]]
# ===== Check sample columns =====
print(colnames(TPM_deletion.dt))

# ===== Calculate mean expression per group =====
TPM_deletion.dt[, PC := rowMeans(.SD, na.rm = TRUE), .SDcols = patterns("^PC", cols = colnames(TPM_deletion.dt))]
TPM_deletion.dt[, TC5 := rowMeans(.SD, na.rm = TRUE), .SDcols = patterns("^TC5", cols = colnames(TPM_deletion.dt))]

# ===== Create expression matrix with gene IDs as rownames =====
mean_expr <- as.matrix(TPM_deletion.dt[, .(gene_name, PC, TC5)],rownames = 'gene_name')


# ===== Plot Heatmap =====

HM.p<-Heatmap(
  mean_expr,
  name = "Expression(TPM)",
  col = colorRamp2(c(0, max(mean_expr)), c("#fae1e3", "#f20f0f")),
  cluster_rows = FALSE,
  cluster_columns = FALSE,
  show_row_names = TRUE,
  row_names_side = "right",
  row_names_gp = grid::gpar(fontsize = 10, fontface = "bold.italic"),
  show_column_names = TRUE,
  column_names_gp = grid::gpar(fontsize = 10),
  column_names_centered = TRUE,
  column_names_rot= 0,
  rect_gp = gpar(col = "black", lwd = 0.5),
  heatmap_legend_param = list(direction = "horizontal",
                              title_gp = gpar(fontsize = 10, fontface = "plain"),
                              at=c(80,40,0),
                              legend_width=unit('3','cm'))
  
)
HM.p<-draw(HM.p,heatmap_legend_side = "top")

PDF_FILE="Fig2D.Deletion_Region_PC_TC5_Mean_Heatmap_with_Genes.pdf"
pdf(PDF_FILE, width = 1.8, height = 2.6,paper='A4')
print(HM.p)
dev.off()
shell.exec(PDF_FILE)
