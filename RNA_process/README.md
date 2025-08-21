# RNAseq_process

This folder contains the RNA-seq analysis pipeline used to study gene expression changes in banana somaclonal variants with resistance to Fusarium wilt.  
It includes trimming, alignment, quantification, variant calling, and differential expression analysis.

---

## 📁 1.preprocess

Preprocessing of raw RNA-seq reads using **Trimmomatic**.

This step uses **Trimmomatic v0.36** to trim raw FASTQ reads based on adapters and quality filters.  
The script loads fasta file from a **tab-delimited file** (`GCTCV.tsv`) that lists the sample IDs, input FASTQ files, and their paths.

### TSV File Format
The TSV file contains the following columns (tab-delimited):

| ID   | SizeMin | FileName1         | FileName2       | Path |
|------|---------|------------------|------------------|------|
| PC-1 | 100     | Sample_R1.fastq.gz | Sample_R2.fastq.gz | /path/to/files/ |

Example snippet from `GCTCV.tsv`:

- **Tasks**:  
  - Remove adapters  
  - Quality trimming  
  - Drop short reads

- **Script**:  
  - `RNA_preprocess.sh`  

---

## 📁 2.mapping and quantification

Alignment of trimmed reads to the reference genome using **STAR**.

### 2.1. Generate STAR and RSEM index
- **Scripts**: 
`RSEM_index.sh`
`STAR_index.sh`

### 2.2. STAR align for GATK calling
- **Tasks**:   
  - 2-pass alignment mode  
  - Output BAM sorted by coordinate
  - Index BAM with samtools

- **Script**:  
  - `STAR_align_gatk.sh`
  
### 2.3. STAR align + RSEM quantification
- **Tasks**:  
  - Prepare RSEM reference from STAR index  
  - 2-pass alignment mode  
  - Output BAM sorted by coordinate
  - Estimate gene and transcript expression levels  
  
- **Script**:
   `STAR_RSEM_quant.sh`  
               
---
## 📁 3.variant_calling

Call variants from RNA-seq BAM files using **GATK**.

### 📁 3.1.PreparationGATK

Prepare genome dict and interval files.

- **Scripts**:
  `01_gatk_dict.sh`
  `02_split_intervals.sh`
  
  
Preprocess BAMs for RNA-seq variant calling:

- **Scripts**:
  `03_mark_duplicates.sh`
  `04_split_n_cigar_reads.sh`


### 📁 3.2.gHC (GATK HaplotypeCaller)

Call variants in GVCF mode.

- **Script**:  
  - `05_process_vcf05.sh`

Merges per-interval GVCF outputs into a single file.

- **Script**:  
  - `06_merge_intervals.sh`

Combines GVCFs from all samples for joint genotyping

- **Scripts**:  
  - `07_merge_gvcfs.sh`


### 📁 3.3.GT (GenotypeGVCFs)

Joint genotyping across GVCFs.

- **Scripts**:  
  - `08_jointgenotyping.sh`


### 📁 3.4.HF (Hard Filtering)

Apply hard filters to raw variant calls.

- **Scripts**:  
  - `09_filter_variants.sh`
---

## 📁 4.DEG_analysis

Differential expression analysis using **DESeq2** in R.

- **PC vs TC5 Comparison**:  
  - Import and merge RSEM count tables (`*.genes.results`)
  - Filter lowly expressed genes  
  - Build metadata (condition: PC vs TC5; batch effect correction)  
  - Run **DESeq2**  
  - Extract DEGs (|log2FC| > 0.5, FDR < 0.1) 
  
- **Script**:  
  - `DEG_analysis.R`
  **Output**
  - DEG table: `TC5_vs_PC_DEG.tsv` 
  
Nongke No.1 vs Baxi Comparison  
  - Load GTF annotation to map **Gene IDs → Genomic coordinates**  
  - Import time-series RNA-seq counts (27 h & 52 h)  
  - Build sample metadata (condition, batch)  
  - Run **DESeq2** 
  - Annotate DEGs with genomic position  
  - Extract DEGs (|log2FC| > 0.5, FDR < 0.1)  
  - Save results per time point 

- **Script**:  
  - `DEG_analysis.R` (Design matrix: `sample_metadata.tsv`)
  **Output**
  - DEG tables in `DEG_results/` directory (e.g., `DEG_27h.tsv`, `DEG_52h.tsv`)

---

## 🔧 Software Requirements

| Tool            | Version                     |
|-----------------|-----------------------------|
| **Trimmomatic** | version 0.36                |
| **STAR**        | version 2.7.0f              |
| **RSEM**        | version 1.3.0               |
| **samtools**    | version 1.19.2              |
| **picard**      | version 2.21.2              |
| **GATK**        | version 4.1.4.1             |
| **DESeq2**      | version 1.24.0 (R-Package)  |

---

## 🧬 Reference Genome Sources

| Description                                              | Source Link                                                                 |
|----------------------------------------------------------|-----------------------------------------------------------------------------|
| *DH-Pahang* reference genome v4                          | [banana-genome-hub](https://banana-genome-hub.southgreen.fr/)               |

## MaACA8
The gene model from DH-Pahang genome v2 was mapped to the genomic coordinates of genome v4 and is available as: MaACA8.DHv4.gtf

---

## 👤 Author

Puyam Tondonba Singh
*This pipeline was developed and maintained by Puyam Tondonba Singh for integrative RNA-seq analysis of banana Fusarium wilt resistance.*

---
