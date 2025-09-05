# RNAseq_process

This folder contains the RNA-seq analysis pipeline used to study gene expression changes in banana somaclonal variants with resistance to Fusarium wilt.  
It includes trimming, alignment, quantification, variant calling, and differential expression analysis.

---

## 📁 1. preprocess

Preprocessing of raw RNA-seq reads using **Trimmomatic**.
- Script: `RNA_preprocess.sh`
- Input: Tab-delimited file (`GCTCV.tsv`) listing sample IDs, input FASTQ files, and paths.

### TSV File Format
The TSV file contains the following columns (tab-delimited):

| ID   | SizeMin | FileName1          | FileName2        | Path |
|------|---------|--------------------|------------------|------|
| PC-1 | 100     | Sample_R1.fastq.gz | Sample_R2.fastq.gz | /path/to/files/ |

Example snippet from `GCTCV.tsv`:

- **Tasks**:  
  - Remove adapters  
  - Quality trimming  
  - Drop short reads

- Script: `RNA_preprocess.sh`  

---

## 📁 2. mapping and quantification

Alignment of trimmed reads to the reference genome using **STAR**.

### 📁 3.1. Generate STAR index
- Script: `STAR_index.sh`

### 📁 3.2. STAR align for GATK calling
- **Tasks**:   
  - 2-pass alignment mode  
  - Output BAM sorted by coordinate
  - Index BAM with samtools

- Script: `STAR_align.sh`
  
### 📁 3.3. STAR align + RSEM quantification
- **Tasks**:  
  - Prepare RSEM reference from STAR index  
  - 2-pass alignment mode  
  - Output BAM sorted by coordinate
  - Estimate gene and transcript expression levels  
  
- Script: 
  - `RSEM_ref.sh`  
  - `STAR_RSEM_quant.sh`  
               
---
## 📁 3. variant_calling

Call variants from RNA-seq BAM files using **GATK**.

### 📁 3.1.PreparationGATK

Prepare genome and interval files.

- Script:
  - `01_split_intervals.sh`
  - `02_gatk_dict.sh`
  
Preprocess BAMs for RNA-seq variant calling:

- Script: 
  - `03_mark_duplicates.sh`
  - `04_split_n_cigar_reads.sh`


### 📁 3.2. gHC (GATK HaplotypeCaller)

Call variants in GVCF mode.

- Script: `05_process_vcf05.sh`

Merges per-interval GVCF outputs into a single file.

- Script: `06_merge_intervals.sh`

Combines GVCFs from all samples for joint genotyping

- Script: `07_merge_gvcfs.sh`

### 📁 3.3. GT (GenotypeGVCFs)

Joint genotyping across GVCFs.

- Script: `08_jointgenotyping.sh`


### 📁 3.4. HF (Hard Filtering)

Apply hard filters to raw variant calls.

- Script: `09_filter_variants.sh`
---

## 📁 4. DEG_analysis

Differential expression analysis using **DESeq2** in R.

- **PC vs TC5 Comparison**:  
  - Import and merge RSEM count tables (`*.genes.results`)
  - Filter lowly expressed genes  
  - Build metadata (condition: PC vs TC5; batch effect correction)  
  - Run **DESeq2**  
  - Extract DEGs (|log2FC| > 0.5, FDR < 0.1) 
  
- Script: `DEG_analysis.R`  
- Output: DEG table `TC5_vs_PC_DEG.tsv` 
  
- **Nongke No.1 vs Baxi Comparison**:
  - Load GTF annotation to fetch **Genomic coordinates of genes**  
  - Import time-series RNA-seq counts (27 h & 52 h)  
  - Build sample metadata (condition, batch)  
  - Run **DESeq2** 
  - Annotate DEGs with genomic position  
  - Extract DEGs (|log2FC| > 0.5, FDR < 0.1)  
  - Save results per time point 

- Script: `DEG_analysis.R`  
- Output: DEG tables in `DEG_results/` directory (e.g., `DEG_27h.tsv`, `DEG_52h.tsv`)

---

## 📁 5. Build Musa acuminata OrgDb.

Custom Musa acuminata OrgDb built with AnnotationForge

- **Build OrgDb**:  
  - Input: GO annotation file Musa_acuminata_pahang_v4_go.txt
  - Build package using makeOrgPackage()
  - Install the package locally (org.Macuminata.eg.db)
  
- Script: `musaOrgDb.R`  
  

## 🔧 Software Requirements

| Tool                     | Version                     |
|--------------------------|-----------------------------|
| **Trimmomatic**          | version 0.36                |
| **STAR**                 | version 2.7.0f              |
| **RSEM**                 | version 1.3.0               |
| **samtools**             | version 1.19.2              |
| **picard**               | version 2.21.2              |
| **GATK**                 | version 4.1.4.1             |
| **R**                    | version 4.3.2               |
| **Bioconductor**         | version 3.18                |
| **DESeq2**               | version 1.42 (R-Package)    |
| **AnnotationForge**      | version 1.44 (R-Package)    |
---

## 🧬 Reference Genome Sources

| Description                         | Source Link                                                   |
|-------------------------------------|---------------------------------------------------------------|
| *DH-Pahang* reference genome v4     | [banana-genome-hub](https://banana-genome-hub.southgreen.fr/) |

## MaACA8
The gene model from DH-Pahang genome v2 was mapped to the genomic coordinates of genome v4 and is available as: MaACA8.DHv4.gtf

---

## 👤 Author

Puyam Tondonba Singh
*This pipeline was developed and maintained for integrative RNA-seq analysis of banana Fusarium wilt resistance.*

---
