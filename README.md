# myMicrobiomeHelpers

[![R-project](https://shields.io)](https://r-project.org)
[![License: MIT](https://shields.io)](https://opensource.org)

A minimalist R package designed to streamline, automate, and scale **microbiome analyses** and **bootstrapping networks**. This tool acts as an optimized, wrapper pipeline for integrating `phyloseq` datasets with downstream analysis. This repository contains custom R utility functions for microbiome analysis workflows of 16S rRNA sequencing data.

## 🚀 Key Features

* **Composition analysis:** Uses heatmaps to visualize the differences in abundance of the most abundant bacteria between groups at any given taxonomic level.
* **Differential abundance analysis:** For a given taxonomic level of interest, a differential abundance analysis using LinDA is performed.
* **Pairwise PERMANOVA:** A pairwise PERMANOVA between groups for the beta diversity using a distance matrix.
* **Automated Bootstrapping:** Resamples sample matrices with replacement to generate robust, reproducible edge stability profiles.
* **Smart Resource Management:** Automatically detects system architecture limits and optimizes multi-core processing assignments (defaulting safely to 80% CPU capacity).
* **Fault-Tolerant Execution:** Built-in checkpointing skips previously computed experimental groups to safeguard against pipeline interruptions.
* **Defensive Programming:** Implements strict data structure and type validation checks up front to prevent deep-loop execution failures.

## 🛠️ Installation

You can install this portfolio utility directly from GitHub using the `remotes` package:

```R
# Install remotes if you haven't already
if (!requireNamespace("remotes", quietly = TRUE)) {
    install.packages("remotes")
}

# Install this package
remotes::install_github("tiagofmota/myMicrobiomeHelpers")
```

## 💻 Quick Start Example

```R
# Run composition pipeline on Genus breakdowns
res <- compos(
  phylo            = my_phyloseq_obj,                       # Phyloseq object
  level            = "Genus",                               # Taxonomic level as written in colnames(my_phyloseq_obj@tax_table)
  group            = "Treatment",                           # Can be meta column name or custom group vector
  transform.method = "CLR",                                 # Transformation method. So far, VST, ILR and CLR are accepted
  transform.offset = 0.5,                                   # Transformation offset to be used in case transform.method is either "CLR" or "ILR"
  group.order      = c("Control", "Low_Dose", "High_Dose")  # A vector with groups to set order of factors. The first will be for comparisons with all other like a negative control group
)

# Extract structured outputs
transformed_data <- res@transformDF
ComplexHeatmap::draw(res@heatmap, heatmap_legend_side = "bottom")

# Run differential abundance analysis using automated CPU tracking
analysis_output <- difab(
  phylo        = my_phyloseq_obj,                    # Phyloseq object
  level        = "Genus",                            # Taxonomic level as written in colnames(my_phyloseq_obj@tax_table)
  formula      = "~ TreatmentGroup + (1|SubjectID)", # Mixed-effect models benefit from multi-core parsing
  var          = "Treatment",                        # Name of Variable (colname) of interest
  groups       = c("Control", "Treated"),            # A vector with groups to set order of factors. The first will be for comparisons with all other like a negative control group
  da.alpha     = 0.05,                               # Alpha threshold for adjusted p-value
  FCtreshold   = 1.0,                                # Fold Change threshold
  cpus         = NULL                                # Leave NULL to automatically manage 80% system resource load
)

# Extract structured results
print(analysis_output@daPlot)

# Export ASV table to biom format for picrust2 using any phyloseq object
write_biom_csv(
  ps,    # Phyloseq object
  file   # Character with name for the biom file to be generated. It can be a full directory with filename in the end
)

# Performs post-hoc pairwise PERMANOVA tests directly on a pre-calculated distance/dissimilarity matrix (`dist` or symmetric matrix format) across multi-level group vectors. It automatically partitions distance subsets, applies multi-test corrections, and outputs multivariate homogeneity group dispersion values side-by-side.

## Generate your custom distance matrix (e.g., PhiLR Euclidean distance)
# (Or use any other metric like phyloseq::distance(ps, method = "wunifrac"))
philr_coordinates <- philr(t(otu_table(my_phyloseq_obj)), tree, part, groups)
philr_dist_matrix <- dist(philr_coordinates, method = "euclidean")

## Extract matching metadata grouping vector
group_factors <- as.character(sample_data(my_phyloseq_obj)\$Treatment)

## Execute post-hoc comparisons using the pre-calculated distance matrix
pairwise_results <- pairwise.adonis(
  x          = as.matrix(dist), # Accepts dist object or full square distance matrix
  factors    = group_factors,   # Grouping variable
  p.adjust.m = "bonferroni"     # p-adjust method as from stats::p.adjust()
)

# View results table containing adjusted p-values and beta-dispersion p-values
print(pairwise_results)

# Run the bootstrapping network pipeline across a named list of phyloseq objects
boot_network(
  phy_list   = my_phyloseq_list, # Named list (e.g., list(Control = ps1, Treated = ps2))
  nboot      = 100,              # Number of bootstrap iterations
  outputDir  = "./network_res",  # Directory to save .rds outputs and checkpoints
  keep_taxa  = 10,               # Minimum abundance threshold to keep an ASV
  cpus       = NULL              # Leave NULL to automatically allocate 80% of system cores
)
```

## 📋 Package Dependencies

This package integrates multiple core ecology and data science frameworks. Since these dependencies span CRAN, Bioconductor, and GitHub, please ensure they are installed using the commands below:

| Package | Source / Installation Command | Used In Function(s) |
| :--- | :--- | :--- |
| `parallel` | *Built-in (Base R — No installation required)* | `boot_network` |
| `phyloseq` | `BiocManager::install("phyloseq")` | `boot_network`, `compos`, `difab`, `write_biom_csv` |
| `ggClusterNet` | `remotes::install_github("taowenmicro/ggClusterNet")` | `boot_network` |
| `mixOmics` | `BiocManager::install("mixOmics")` | `compos`, `difab` |
| `microbiome` | `BiocManager::install("microbiome")` | `compos`, `difab` |
| `ComplexHeatmap` | `BiocManager::install("ComplexHeatmap")` | `compos`  |
| `MicrobiomeStat` | `remotes::install_github("cafferychen777/MicrobiomeStat")` | `difab` |
| `vegan` | `install.packages("vegan")` | `pairwise.adonis` |
| `DESeq2` | `BiocManager::install("DESeq2")` | `compos`  |
| `SpiecEasi` | `remotes::install_github("zdk123/SpiecEasi")` | `boot_network`*If SparCC is used* |
| `ggpubr` | `install.packages("ggpubr")` | `compos`, `difab` |
| `tidyverse` | `install.packages("tidyverse")` | `write_biom_csv` |

## ⚠️ Disclaimer & Support

This repository serves as a personal research utility and a professional portfolio piece. 
* **Support Status:** This software is provided **"as-is"** without active technical support. GitHub Issues and pull requests are disabled.
* **Modifications:** You are welcome to **fork** this repository to adapt, extend, or fix the code for your own research constraints.

## 📄 Citation

If you use this workflow or codebase to support your academic research, please cite it as follows:

> **Mota, Tiago Feitosa.** (2026). *myMicrobiomeHelpers: An R package for microbiome and network bootstrapping workflow automation*. GitHub repository: `https://github.com`

