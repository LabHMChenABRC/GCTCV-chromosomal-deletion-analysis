# ============================================================
# Chromosome 5 bin count gene density

# This script visualizes compare of gene density in FM chr05 deletion regions across three genome assemblies.Gene density is shown in 50 kb bins.

# ============================================================
# Author: BoHan Hou 
# Date: 2025-08-18



library(data.table)
library(GenomicRanges)
library(rtracklayer) # import.gff3()


# Function: Generate fixed-width bins and count genes
count_genes_per_bin <- function(target_region, genes, bin_size = 50000,as.datatable=FALSE) {
  # target_region: GRanges of length 1 (seqname, start, end)
  # genes: GRanges object with gene features
  # bin_size: fixed bin size in bp (last bin may be shorter)
  
  if (!is(target_region, "GRanges")) {
    stop("target_region must be a GRanges object")
  }
  if (length(target_region) != 1) {
    stop("target_region must have exactly 1 range")
  }
  
  # Generate fixed-width bins
  bins_rel <- tileGenome(
    seqlengths = setNames(width(target_region), as.character(seqnames(target_region))),
    tilewidth = bin_size,
    cut.last.tile.in.chrom = TRUE
  )
  
  # Shift bins to start at target_region@start
  start_bins <- start(bins_rel) + start(target_region) - 1L
  end_bins   <- pmin(end(bins_rel) + start(target_region) - 1, end(target_region))
  bins       <- GRanges(seqnames = seqnames(bins_rel),IRanges(start=start_bins,end=end_bins))
  # Keep bins only within target_region
  bins <- pintersect(bins, target_region)
  
  # Count overlaps
  mcols(bins)$gene_count <- countOverlaps(bins, genes)
  if(as.datatable){
    return(as.data.table(bins))
  }else{
    return(bins)
  }
  
}

Li_Baxi_gff_file="R:/NGS-Reference/Musa_acuminata_BaxiJiao/X_Li_2023/Cavendish_gene.gff3"
Hu_Baxi_H2_gff_file="R:/NGS-Reference/Musa_acuminata_BaxiJiao/HR_Huang2023/Baxijiao.assembly_haplotype2.gff3"
DH_gff_file="R:/NGS-Reference/Musa_acuminata_DH-Pahang_v4.3/Musa_acuminata_pahang_v4.20240506.gff3"

DH_gff.gr <-import.gff3(DH_gff_file)
Li_Baxi_gff.gr <-import.gff3(Li_Baxi_gff_file)
Hu_Baxi_H2_gff.gr <-import.gff3(Hu_Baxi_H2_gff_file)

# Gene position
DH_gene.gr         <-DH_gff.gr[mcols(DH_gff.gr)$type=='gene']
Li_Baxi_Ban_gene.gr<-Li_Baxi_gff.gr[mcols(Li_Baxi_gff.gr)$type=='gene' & grepl('Ban',seqnames(Li_Baxi_gff.gr))]
Hu_Baxi_H2_gene.gr <-Hu_Baxi_H2_gff.gr[mcols(Hu_Baxi_H2_gff.gr)$type=='gene']

# target region
DEL_Region.gr<-data.table(
  seqnames=c("chr05","chr05.1","chr05.2",'Ban05'),
  start =c(956425, ,830729, 1067819),
  end   =c(7030674,6869408,7314107)
)|>GRanges()
names(DEL_Region.gr)<-c("chr05","chr05.2","Ban05")

Dh.bincount.dt    <-count_genes_per_bin(target_region = DEL_Region.gr['chr05'],genes = DH_gene.gr,bin_size = 50*10^3,as.datatable=TRUE) 
Hu.H2.bincount.dt <-count_genes_per_bin(target_region = DEL_Region.gr['chr05.2'],genes = Hu_Baxi_H2_gene.gr,bin_size = 50*10^3,as.datatable=TRUE) 
Li.Ban.bincount.dt<-count_genes_per_bin(target_region = DEL_Region.gr['Ban05'],genes = Li_Baxi_Ban_gene.gr,bin_size = 50*10^3,as.datatable=TRUE) 

all.dt <- rbind(Dh.bincount.dt, Hu.H2.bincount.dt, Li.Ban.bincount.dt
                )[,seqnames:=factor(seqnames,
                                    levels=names(DEL_Region.gr),
                                    labels=c("chr05 (DH pahang)","Chr05.2 (Huang 2023)","Ban05 (Li 2023)"))]


ggplot2_linewidth_base=ggplot2::.pt * 72 / 96 # 72/96 in pixel
# Convert 1 point to ggplot2 linewidth units:
# linewidth = 1 / ggplot2::.pt * 72 / 96

# Step 1: 1 / ggplot2::.pt
#   - Converts 1 point to millimeters (ggplot2's internal unit)
#   - 1 pt = 1/72 inch, ggplot2::.pt handles the pt-to-mm conversion

# Step 2: * 72 / 96  
#   - Adjusts for screen vs print DPI difference
#   - 72 = standard print DPI, 96 = standard screen DPI
#   - Ensures consistent visual size across devices

# Result: linewidth value that renders as 1 point thick visually


library(ggplot2)
Gene_BinCount.p<-ggplot(all.dt, aes(x=start, y=gene_count, fill=seqnames)) +
  geom_col(position = position_dodge(width = 2,preserve ='single'),color=NA) +
  facet_wrap(~seqnames, scales="free_x", ncol=1) +
  labs(x="Position within deletion region (Mb)", y="Number of genes per 50kb bin") +
  scale_y_continuous(expand = expansion(mult = c(0,0.05)),breaks = c(0,15))+
  scale_x_continuous(breaks = seq(0, max(all.dt$start), 2000000),labels = seq(0, max(all.dt$start), 2000000)/10^6,expand = expansion(mult = c(0.01,0.01)))+
  cowplot::theme_half_open(line_size = 0.5/ggplot2_linewidth_base, font_size = 6,rel_small = 1,rel_tiny = 1,rel_large = 1)+
  theme(legend.position = 'none')


fill_colorset=setNames(c("#F8766D","#00D552", "#619CFF"),
                       c("chr05 (DH pahang)","Chr05.2 (Huang 2023)","Ban05 (Li 2023)"))
# change colors
# Gene_BinCount.p<-Gene_BinCount.p+scale_fill_manual(values = fill_colorset)



Gene_BinCount.p
Gene_BinCount.pdf="R:/Qsync/Banana/FOC_RNASeq/MS_ACAs/Figure_PNAS/Gene_BinCount_withinFM-DEL-chr05.pdf"
ggsave(Gene_BinCount.pdf,Gene_BinCount.p,width = 6,height = 5,units = 'cm',paper="A4")
shell.exec(Gene_BinCount.pdf)
