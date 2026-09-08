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
