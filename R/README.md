# R Functions for Microbiome Analysis

This folder contains R helper functions for analyzing microbiome data from QIIME2 and PICRUSt2 outputs.

---

## File Contents

### `micro_functions.R`
This file contains all helper functions I have built to improve my microbiome analysis workflow. Key features and usage are described bellow.

#### 🚀 Key Features

* **Composition analysis:** Uses heatmaps to visualize the differences in abundance of the most abundant bacteria between groups at any given taxonomic level.
* **Differential abundance analysis:** For a given taxonomic level of interest, a differential abundance analysis using LinDA is performed.
* **Pairwise PERMANOVA:** A pairwise PERMANOVA between groups for the beta diversity using a distance matrix.
* **Automated Bootstrapping:** Resamples sample matrices with replacement to generate robust, reproducible edge stability profiles.
* **Smart Resource Management:** Automatically detects system architecture limits and optimizes multi-core processing assignments (defaulting safely to 20% CPU capacity).
* **Fault-Tolerant Execution:** Built-in checkpointing skips previously computed experimental groups to safeguard against pipeline interruptions.
* **Defensive Programming:** Implements strict data structure and type validation checks up front to prevent deep-loop execution failures.

### `micro_analysis.R`
This file contains the whole microbiome analysis script I have written during the past years and it includes the import of .qza data from the Qiime2 output files to make a phyloseq object until co-occurrence network analysis using SparCC. The script contains the following steps:

| Analysis | Brief description | Packages involved |
| :--- | :--- | :--- |
| `Composition` | `CLR transformed counts are summarized by groups and plotted in heatmaps at any taxonomic level` | `phyloseq`; `microbiome`; `ComplexHeatmap` |
| `Alpha diversity` | `Shannon, Inverse Simpson and Faith's phylogenetic diversity metrics are computed and plotted in violin boxplots` | `picante`; `phyloseq`; `ggplot2`; `ggpubr` |
| `Beta diversity` | `PhiLR transformation is applied and using a using the Euclidian distance a PCoA scatter plot is produced. A pairwise PERMANOVA is used to statistically assess differences between groups` | `philr`; `mixOmics`; `vegan`; `ggplot2` |
| `Differential abundance` | `LinDA is used in order to assess deferentially abundant bacteria and plotted in Diverging Bar Plot` | `MicrobiomeStat`; `ggplot2` |
| `Functional prediction` | `The output file from PICRUSt2's functional prediction is used and LinDA is used to assess deferentially abundant KOs which are applied in an enrichment analysis` | `MicrobiomeProfiler`; `ggplot2`; `clusterProfiler` |
| `Co-occurrence network` | `Bootstrapped co-occurrence networks are built and used for plotting as well as statistical comparison between groups using the global network metrics` | `igraph`; `SpiecEasi`; `ggClusterNet`; `parallel` |

## 🛠️ Source all functions
```R
# Option 1: Source directly from GitHub
source("https://raw.githubusercontent.com/tiagofmota/myMicrobiomeHelpers/main/R/microbiome_helpers.R")

# Option 2: Clone the repository
git clone https://github.com/tiagofmota/myMicrobiomeHelpers.git
cd myMicrobiomeHelpers
# Then source the R file
R
> source("R/microbiome_helpers.R")
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

## 📋 Dependencies

These scripts integrate multiple core ecology and data science frameworks. Since these dependencies span CRAN, Bioconductor, and GitHub, please ensure they are installed using the commands below:

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
| `qiime2R` | `remotes::install_github("jbisanz/qiime2R")` |  |
| `picante` | `install.packages("picante")` |  |
| `philr` | `BiocManager::install("philr")` |  |
| `MicrobiomeProfiler` | `BiocManager::install("MicrobiomeProfiler")` |  |
| `clusterProfiler` | `BiocManager::install("clusterProfiler")` |  |
| `ggplot2` | `install.packages("ggplot2")` |  |
| `igraph` | `install.packages("igraph")` |  |
| `ape` | `install.packages('ape')` |  |
| `decontam` | `BiocManager::install("decontam")` |  |

## ⚠️ Disclaimer & Support

This repository serves as a personal research utility and a professional portfolio piece. 
* **Support Status:** This software is provided **"as-is"** without active technical support. GitHub Issues and pull requests are disabled.
* **Modifications:** You are welcome to **fork** this repository to adapt, extend, or fix the code for your own research constraints.

## 📄 Citation

If you use this workflow or codebase to support your academic research, please cite it as follows:

Mota, T.F. (2026). QIIME2 + PICRUSt2 Pipeline for 16S Analysis [Computer software]. GitHub. https://github.com/tiagofmota/myMicrobiomeHelpers

Also cite the tools used:

    QIIME 2 ------ Reproducible, interactive, scalable and extensible microbiome data science using QIIME 2. Bolyen et al., 2019, Nature Biotechnology; https://doi.org/10.1038/s41587-019-0209-9
    PICRUSt2 ------ PICRUSt2 for prediction of metagenome functions. Douglas et al., 2020, Nature Biotechnology; https://doi.org/10.1038/s41587-020-0548-6
    DADA2 ------ DADA2: High-resolution sample inference from Illumina amplicon data. Callahan et al., 2016, Nature Methods; https://doi.org/10.1038/nmeth.3869
    Deblur ------ Deblur Rapidly Resolves Single-Nucleotide Community Sequence Patterns. Amir et al., 2017, Novel Systems Biology Techniques; https://doi.org/10.1128/msystems.00191-16
    SILVA database ------ SILVA in 2026: a global core biodata resource for rRNA within the DSMZ digital diversity. Chuvochina et al., 2026, Nucleic Acids Research; https://doi.org/10.1093/nar/gkaf1247
    phyloseq: phyloseq ------ An R Package for Reproducible Interactive Analysis and Graphics of Microbiome Census Data, McMurdie PJ & Holmes S, 2013, PLOS ONE; https://doi.org/10.1371/journal.pone.0061217
    ComplexHeatmap ------ Complex heatmaps reveal patterns and correlations in multidimensional genomic data, Gu et al., 2016, Bioinformatics; https://doi.org/10.1093/bioinformatics/btw313
    ape ------ ape 5.0: an environment for modern phylogenetics and evolutionary analyses in R, Paradis & Schliep, 2019, Bioinformatics; https://doi.org/10.1093/bioinformatics/bty633
    picante ------ Picante: R tools for integrating phylogenies and ecology, Kembel et al., 2010, Bioinformatics; https://doi.org/10.1093/bioinformatics/btq166
    philr ------ A phylogenetic transform enhances analysis of compositional microbiota data, 2017, eLife; https://doi.org/10.7554/eLife.21887
    DESeq2 ------ Moderated estimation of fold change and dispersion for RNA-seq data with DESeq2, Love et al., 2014, Genome Biology; https://doi.org/10.1186/s13059-014-0550-8
    LinDA ------ LinDA: linear models for differential abundance analysis of microbiome compositional data. Zhou et al., 2022, Genome Biology; https://doi.org/10.1186/s13059-022-02655-5
    clusterProfiler ------ Thirteen years of clusterProfiler, Yu, 2024, The Innovation; https://doi.org/10.1016/j.xinn.2024.100722
    decontam ------ Simple statistical identification and removal of contaminant sequences in marker-gene and metagenomics data, Davis et al., 2018, Microbiome; https://doi.org/10.1186/s40168-018-0605-2

Please cite packages according to citation()

