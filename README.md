# myMicrobiomeHelpers

**A collection of R functions and bash scripts for reproducible 16S rRNA microbiome data analysis**

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.XXXXXXX.svg)](https://doi.org/10.5281/zenodo.XXXXXXX)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![R](https://img.shields.io/badge/R-4.0+-blue.svg)](https://www.r-project.org/)
[![QIIME2](https://img.shields.io/badge/QIIME2-2023.9+-orange.svg)](https://qiime2.org/)

---

## Overview

This repository contains scripts I've developed for my microbiome research. It bridges the gap between raw 16S rRNA sequencing data processing and downstream statistical analysis and visualization.

### What's Inside

| Component | Description |
|-----------|-------------|
| **R Functions** | Helper functions for diversity analysis, visualization, and statistical testing |
| **Bash Pipeline** | Complete QIIME2 + PICRUSt2 workflow for 16S data processing |
| **SLURM Scripts** | Job submission scripts for cluster environments |

---

## Quick Navigation

- [**R Functions**](./R/README.md) - Import QIIME2 data, calculate diversity, create visualizations
- [**Bash Pipeline**](./bash/README.md) - Process raw FASTQ files through QIIME2 and PICRUSt2
- [**Installation**](#installation) - Set up dependencies
- [**Usage Examples**](#usage-examples) - Get started quickly
- [**Citation**](#citation) - How to cite this work

---

## Installation

### R Functions

```r
# Option 1: Source directly from GitHub
source("https://raw.githubusercontent.com/tiagofmota/myMicrobiomeHelpers/main/R/microbiome_helpers.R")

# Option 2: Clone the repository
git clone https://github.com/tiagofmota/myMicrobiomeHelpers.git
cd myMicrobiomeHelpers
# Then source the R file
R
> source("R/microbiome_helpers.R")
