
# This script summarizes counts FOR  'Nongke No. 1' VS 'Baxi', and generates barplots and chromosome-level plots.


# ============================================================
# Author: Puyam Tondonba Singh
# Date: 2025-08-18
# ============================================================


library(data.table)
library(ggplot2)
library(cowplot)


setwd("E:/Visualization/4E") 

# --- LOAD DEG TABLES FROM TSV ---
DEG.files <-list.files(path = 'DEG_results',pattern = 'DEG.*.tsv',full.names = TRUE)
names(DEG.files)<-ifelse(grepl('27',basename(DEG.files)),'27 H', '52 H')

DEG.dt <- lapply(DEG.files,fread)|>rbindlist(idcol = 'Condition')
DEG.dt <- DEG.dt[grepl('chr',Chromosome)][color %in%  c("Upregulated", "Downregulated")]

DEG.DEL.dt<-DEG.dt[Chromosome == "chr05" & Start >= 956425 & End <= 7030674]
Count.data<-rbind(
  DEG.dt[,.(Region='Whole genome',Count=.N),by=.(Condition,Regulation=color)],
  DEG.DEL.dt[,.(Region='NK/FM deletion',Count=.N),by=.(Condition,Regulation=color)]
)

Count.data[,Region    := factor(Region, levels = c("Whole genome", "NK/FM deletion"))]
Count.data[,Regulation:= factor(Regulation, levels = c("Upregulated", "Downregulated"),labels = c("Up-regulated", "Down-regulated"))]

# --- BARPLOT ---
p <- ggplot(Count.data, aes(x = Condition, y = Count, fill = Regulation)) +
  geom_bar(stat = "identity") +
  labs(y='Number of DEGs',x='NK vs BX',fill=NULL)+
  facet_wrap(~Region, nrow = 1, scales = "free_y") +
  scale_fill_manual(values = c("Down-regulated" = "#00CFC1", "Up-regulated" = "#F8766D")) +
  theme_minimal(base_size = 12) +
  scale_y_continuous(expand = expansion(mult = c(0,NA)))+
  theme_cowplot()+
  theme(strip.background = element_blank(),
        legend.position = 'top')

PDF_FILE="Fig4E.DEG_plot.pdf"
ggsave(PDF_FILE, plot = p, width = 12, height = 8, units = "cm")
shell.exec(PDF_FILE)