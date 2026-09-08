# Bash Pipeline for QIIME2 and PICRUSt2

Complete workflow for processing 16S rRNA amplicon sequencing data from raw FASTQ files to functional predictions.

---

## Overview

This pipeline automates the entire 16S data processing workflow:

1. **Manifest Creation** - Generate QIIME2 manifest file from FASTQ files
2. **Data Import** - Import sequences into QIIME2
3. **Denoising** - DADA2 or Deblur (user choice)
4. **Taxonomy Classification** - SILVA database (full-length or region-specific)
5. **Phylogenetic Tree** - Build tree for diversity analyses
6. **Functional Prediction** - PICRUSt2 (optional)

---

## Scripts

- `qiime2_picrust2_pipeline.sh` - Main pipeline script
- `submit_pipeline.slurm` - SLURM submission script for clusters

---

## Requirements

### Software
- **Conda** - Environment management
- **QIIME2** (2023.9 or later)
- **PICRUSt2** (optional, for functional prediction)

### Databases
- **SILVA** (v138.2) - Automatically downloaded if not found

---

## Installation

```bash
# Clone repository
git clone https://github.com/tiagofmota/myMicrobiomeHelpers.git
cd myMicrobiomeHelpers/bash

# Make scripts executable
chmod +x qiime2_picrust2_pipeline.sh
chmod +x submit_pipeline.slurm

# Set up conda environments
conda create -n qiime2-amplicon-2025.7 -c conda-forge -c bioconda -c qiime2 qiime2
conda create -n picrust2 -c bioconda -c conda-forge picrust2  # optional

```

## Usage
### Required Arguments
Argument	Description	Example

```bash
--work-dir	Working directory (absolute path)	/home/user/project
--seqs-dir	Directory containing FASTQ files	/home/user/data/reads
--forward-pattern	Pattern for forward reads	_R1, _1.fq.gz
--reverse-pattern	Pattern for reverse reads	_R2, _2.fq.gz
--n-char	Characters to trim for sample ID	5
--o-manifest	Output manifest path (no extension)	results/manifest
--o-qiime	Output directory for QIIME2 results	results/qiime2
Optional Arguments
Argument	Description	Default
--d-denoise	Denoising method: deblur or dada2	deblur
--p-picrust	Run PICRUSt2: true or false	false
--threads	Number of CPU threads	4
--qiime-env	Conda environment for QIIME2	qiime2-amplicon-2025.7
--picrust-env	Conda environment for PICRUSt2	picrust2
Primer Options (Region-Specific Classification)
Argument	Description
--forward-primer	Forward primer sequence (e.g., CCTACGGGNGGCWGCAG)
--reverse-primer	Reverse primer sequence (e.g., GGACTACNVGGGTWTCTAAT)

```
Important: Both primers must be provided together. If omitted, the full-length classifier will be used.

## Examples
Example 1: Basic Pipeline (Full-Length Classifier)

```bash
./qiime2_picrust2_pipeline.sh \
    --work-dir /project/16S_data \
    --seqs-dir /project/raw_data \
    --forward-pattern _R1_001.fastq.gz \
    --reverse-pattern _R2_001.fastq.gz \
    --n-char 5 \
    --o-manifest results/manifest \
    --o-qiime results/qiime2 \
    --threads 8
```

Example 2: V3-V4 Region (DADA2 + Region-Specific Classifier)
Region is specified by primers provided in arguments --forward-primer and --reverse-primer

```bash
./qiime2_picrust2_pipeline.sh \
    --work-dir /project/16S_data \
    --seqs-dir /project/raw_data \
    --forward-pattern _R1 \
    --reverse-pattern _R2 \
    --n-char 5 \
    --o-manifest results/manifest \
    --o-qiime results/qiime2 \
    --d-denoise dada2 \
    --forward-primer CCTACGGGNGGCWGCAG \
    --reverse-primer GGACTACNVGGGTWTCTAAT \
    --threads 8
```

Example 3: Full Workflow (DADA2 + PICRUSt2)
bash

```bash
./qiime2_picrust2_pipeline.sh \
    --work-dir /project/16S_data \
    --seqs-dir /project/raw_data \
    --forward-pattern _R1 \
    --reverse-pattern _R2 \
    --n-char 5 \
    --o-manifest results/manifest \
    --o-qiime results/qiime2 \
    --d-denoise dada2 \
    --p-picrust true \
    --threads 8
```
