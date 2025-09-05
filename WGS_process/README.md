# WGS_process

This folder contains the step-by-step whole-genome sequencing (WGS) analysis pipeline used in the study of **chromosomal deletions in banana somaclonal variants associated with Fusarium wilt resistance**.

The pipeline includes preprocessing, alignment, coverage calculation, and variant calling. All steps are organized sequentially and numerically.

---
## 🚀 Pipeline Overview

This pipeline processes raw whole-genome sequencing data to produce analysis-ready, filtered variant calls. The general data flow is as follows:

`Raw FASTQ` -> **1. Preprocessing** -> `Cleaned FASTQ` -> **2. Mapping** -> `Aligned CRAM` -> **3. Coverage** & **4. GATK** -> `BigWig` & `Filtered VCF`

---
## ▶️ How to Run

This pipeline is designed to be run on a High-Performance Computing (HPC) cluster using the SLURM workload manager.

1.  **Configuration**: Before running, you may need to edit the `*.sbatch` scripts to set correct paths for your reference genomes, input data, and output directories.
2.  **Execution**: Submit the `*.sbatch` script for each step in numerical order. The sbatch script handles job submission and passes the correct parameters to its corresponding `*.one.sh` core script.

    ```bash
    cd 1.preprocess
    sbatch WGS_preprocess.sbatch
    # ...and so on for the subsequent steps.
    ```

## 📁 1.preprocess

Prepare raw FASTQ files before alignment.

- **Tasks**:  
  - Remove duplicated reads  
  - Quality trimming
  - Drop short trimmed reads

- **Input**: Raw paired-end FASTQ files.
- **Output**: Quality-trimmed, deduplicated FASTQ files.

- **Scripts**:  
  - Submit script: `WGS_preprocess.sbatch`  
  - Core script: `WGS_preprocess.one.sh`
  - sample list file: `sample.list` per sample per line

---

## 📁 2.mapping

Alignment of sequencing reads to the reference genome.

### 📁 2.1.mk_index

Build reference genome index files.

- **For full reference genomes**:  
  - `build_index_[Genome Code].sbatch`

- **For haplotype-specific reference copies**:  
  - `build_index_[Genome Code].split.sbatch`

---

### 📁 2.2.Mapping_merged_copies

This step aligns reads to a primary reference genome assembly.

- **Input**: Cleaned FASTQ files from Step 1.
- **Output**: Coordinate-sorted CRAM files (`.cram`) with index (`.crai`).
Align reads to full reference genomes.

- **Scripts**:  
  - Submit script: `mapping.[Genome Code].sbatch`  
  - Core script: `mapping.one.sh`  
  - Read group file: `WGS.RG.tsv`

---

### 📁 2.3.Mapping_individual_copies

This step is for aligning reads to haplotype-specific assemblies, if available.

- **Input**: Cleaned FASTQ files from Step 1.
- **Output**: Coordinate-sorted CRAM files (`.cram`) with index (`.crai`).
Align reads to haplotype-specific reference genomes.

- **Scripts**:  
  - Submit script: `mapping.[Genome Code].sbatch`  
  - Core script: `mapping.one.sh`  
  - Read group file: `WGS.RG.tsv`

---

## 📁 3. Coverage Calculation

- **Input**: CRAM files from Step 2.
- **Output**: Normalized coverage tracks in BigWig (`.bw`) format.

Estimate sequencing depth based on CRAM files from `2.2.Mapping_merged_copies`.

- **Step 1: Calculate genome size (for normalization)**  
  - Script: `1.calculate_genome_size.sh`  
  - Tool: `faCount` (UCSC tool)

- **Step 2: Compute coverage using deepTools**  
  - Submit script: `2.coverage.[Genome Code].sbatch`  
  - Core script: `coverage.one.sh`

---

## 📁 4.GATK

- **Input**: CRAM files from Step 2.
- **Output**: Hard-filtered SNP and INDEL VCF files.

GATK-based variant calling pipeline.

### 📁 4.1.PreparationGATK

Prepare genome and interval files.

- **Scripts**:  
  - Create `.dict`: `1.Genome.dict.sbatch`  
  - Create `.interval_list`: `2.mk.interval.sbatch`

---

### 📁 4.2.gHC (GATK HaplotypeCaller)

Call variants in GVCF mode.

- **Scripts**:  
  - Submit script: `HC.[Genome Code].sbatch`  
  - Core script: `gHC.one.sh`

---

### 📁 4.3.GT (GenotypeGVCFs)

Joint genotyping across GVCFs.

- **Scripts**:  
  - Submit script: `GT.[Genome Code].sbatch`  
  - Core script: `GT.one.sh`  
  - BED file for filtering high-variant regions: `DHv4.XL.bed`

---

### 📁 4.4.HF (Hard Filtering)

Apply hard filters to raw variant calls.

- **Scripts**:  
  - Submit script: `HF.[Genome Code].sbatch`  
  - Core script: `HF.one.sh`

---

## 🧬 Reference Genome Sources

| Genome Code    | Description                                              | Source Link                                                                                         |
|----------------|----------------------------------------------------------|-----------------------------------------------------------------------------------------------------|
| **DHv4**       | *DH-Pahang* reference genome v4 (Belser etal., 2021)     | [banana-genome-hub](https://banana-genome-hub.southgreen.fr/)                                       |
| **Baxi-Huang** | *Musa Cavendish Baxijiao* (Huang et al., 2023)           | [banana-genome-hub](https://banana-genome-hub.southgreen.fr/)                                       |
| **Baxi-Li**    | *Musa Cavendish Baxijiao* (Li et al., 2023)              | [figshare](https://figshare.com/articles/dataset/The_assembly_and_annotation_of_Cavendish/24113352) |

---

## 🔧 Software Requirements

| Tool          | Version        |
|---------------|----------------|
| **BBTools**   | 39.26          |
| **bwa-mem2**  | 2.1.1          |
| **htslib**    | 1.19.2         |
| **samtools**  | 1.19.2         |
| **bcftools**  | 1.19.2         |
| **GATK**      | 4.1.4.1        |
| **Picard**    | 2.21.2         |
| **bedtools**  | 2.31.1         |
| **deepTools** | 3.4.3          |
| **faCount**   | (UCSC Tools)   |

---

## 👤 Author

BoHan Hou  
*This pipeline was developed and maintained by BoHan Hou for the study of chromosomal deletions linked to Fusarium wilt resistance in banana.*

---
## 📌 Notes

- `[REF-name]` refers to the reference genome used (e.g., `DHv4`, `Baxi-Huang`, `Baxi-Li`)
- Scripts with a `.sbatch` extension are submission scripts for a SLURM-based HPC environment. They wrap the corresponding `.one.sh` scripts, which contain the core command-line logic.
- This pipeline supports both full reference genomes and haplotype-resolved assemblies.
- The `WGS.RG.tsv` file is a **tab-separated** file that provides input paths and read group information for the mapping script. It must not have a header and should contain the following **four** columns:
  1.  **Sample ID**: A unique identifier for the sample (e.g., `PC`).
  2.  **R1 Path**: The full path to the forward-strand (R1) FASTQ file.
  3.  **R2 Path**: The full path to the reverse-strand (R2) FASTQ file.
  4.  **Read Group String**: The complete read group string that will be passed to `bwa-mem2`. The string itself must contain literal tabs (`\t`) as required by the SAM/BAM specification.
---
