# Load packages ####
{
  library(phyloseq)
  library(microbiome)
  library(ggpubr) 
  library(vegan)
  library(qiime2R)
  library(mixOmics)
  library(MicrobiomeStat)
  library(philr)
  library(ape) 
  library(ComplexHeatmap)
  library(pROC)
  library(ConQuR)
  library(doParallel) 
  library(ggsankey)
  library(ggvenn)
  library(igraph) 
  library(ggClusterNet)
  library(caret)
  library(ggrepel)
  library(dplyr)
  library(tibble)
  library(picante)
  library(decontam)
  library(MicrobiomeProfiler)
  source("R/micro_functions.R")
}

# Load Qiime2 data ####
# Set directory variables for Qimme2 files
setwd("~/work_dir")

# Pliefunction .get_colorase with full directory. No "~/"
dir_features <- "qiime_outputs/table-v3v4.qza"
dir_taxonomy <- "qiime_outputs/taxonomy-v3v4.qza"
dir_metadata <- "qiime_outputs/metadata_file2.tsv" # Must be .tsv without quote marks
dir_tree <- "qiime_outputs/fasttree-tree-deblur.qza"

phylo <- qza_to_phyloseq(features = dir_features
                         , taxonomy = dir_taxonomy
                         , metadata = dir_metadata
                         , tree = dir_tree
                         )

# Initial processing ####
# generates the tax ID to use in plots if needed
{
  taxic <- as.data.frame(phylo@tax_table)
  otu.df <- abundances(phylo)
  
  # make a dataframe for OTU information.
  otu.df <- as.data.frame(otu.df)
  
  # check the rows and columns
  # head(otu.df)
  
  # Add the OTU ids from OTU table into the taxa table at the end.
  taxic$ASV <- row.names(otu.df)
  
  # You can see that we now have extra taxonomy levels.
  colnames(taxic)
  
  # convert it into a matrix.
  taxmat <- as.matrix(taxic)
  
  # convert into phyloseq compatible file.
  new.tax <- tax_table(taxmat)
  
  # incroporate into phyloseq Object
  tax_table(phylo) <- new.tax
}

# If metadata cleaning is need it can be done directly in phylo@sam_data
# Here a some examples for dataset filtering/cleaning

# Filter some samples from the phyloseq object
phylo <- prune_samples(!grepl("NTC", phylo@sam_data$ID), phylo)
phylo <- prune_samples(names(which(sample_sums(phylo) >= 1000)), phylo)

# Filter some taxa from the phyloseq object
phylo <- prune_taxa(taxa_sums(phylo) >= 1, phylo)
phylo <- subset_taxa(phylo, !ASV %in% taxa_names(subset_taxa(phylo, Genus == "Incertae_Sedis")))

# Creating some variables from other with cleaner names
phylo@sam_data$Infection.status <- ifelse(phylo@sam_data$Infection.status == "Bac_+ve", "Positive",
                                          ifelse(phylo@sam_data$Infection.status == "Neg_cntrl", "Control",
                                                 "Negative"))

phylo@sam_data$life_stage <- ifelse(grepl("nymph", phylo@sam_data$lab.id), "Nymph", 
                                    ifelse(grepl("male|female", phylo@sam_data$lab.id), "Adult", 
                                           "is_neg"))

# Detect and remove contaminants from ntc samples
sample_data(phylo)$is.neg <- sample_data(phylo)$Infection.status == "NTC"
contamdf.prev <- isContaminant(phylo, method="prevalence", neg="is.neg")
table(contamdf.prev$contaminant)

# Remove ASV assigned as contaminant and remove the NTC samples from phyloseq object
phylo <- prune_taxa(!contamdf.prev$contaminant, phylo)
phylo <- prune_samples(sample_data(phylo)$Infection.status != "Control", phylo)

# Prune samples and taxa to remove possibly empty ones
# Keep phylo before making any other changes in it. Work further with ps2
ps2 = prune_taxa(taxa_sums(phylo) >= 1, phylo)
ps2 = prune_samples(names(which(sample_sums(ps2) >= 1)), ps2)
ps2 <- subset_taxa(ps2, !ASV %in% taxa_names(subset_taxa(ps2, Genus == "Incertae_Sedis")))

# Check how many samples and taxa passed the contaminant removal
nsamples(ps2)
ntaxa(ps2)

# Compositional analysis ####
# Filtering ASVs to improve visualization. remove the low frequent bacteria
pseq2 <- core(ps2, 
              detection = 10, 
              prevalence = .01, # Increase this parameter to have less bacteria in the heatmaps
              include.lowest = T
)
ntaxa(pseq2)

# Run the compos function in a loop to get the heatmap for all taxonomic levels from phylum to genus
for (i in colnames(pseq2@tax_table)[2:6]) {
  assign(paste0(i, ".clr"), compos(pseq2, level = i, 
                                   group = as.character(pseq2@sam_data$Infection.status), 
                                   transform.method = "CLR", 
                                   transform.offset = .5, 
                                   group.order = c("Negative", "Positive"),
                                   ht_subtittle = i,
                                   color.range = c(-1,0,1)))
}

# Plot the heatmaps
draw(Genus.clr@heatmap, merge_legend = TRUE, heatmap_legend_side = "bottom")
draw(Family.clr@heatmap, merge_legend = TRUE, heatmap_legend_side = "bottom")
draw(Order.clr@heatmap, merge_legend = TRUE, heatmap_legend_side = "bottom")
draw(Class.clr@heatmap, merge_legend = TRUE, heatmap_legend_side = "bottom")
draw(Phylum.clr@heatmap, merge_legend = TRUE, heatmap_legend_side = "bottom")

# Set color blind friendly palette for downstream analysis
cblind_palette8 <- c("#0072B2", "#D55E00", "#999999", "#E69F00",
                     "#56B4E9", "#009E73", "#F0E442", "#CC79A7")

# Alpha diversity ####
## Rarefaction and plotting ####
# Parameters to bootstrap the alpha diversity metrics
# This is necessary if there is big difference in library depth among samples
n_bootstrap <- 1000
rarefaction_depth <- min(sample_sums(ps2))  # or whatever depth you choose

# Initialize storage: bootstrap iterations x samples
shannon_matrix <- as.data.frame(matrix(ncol = nsamples(ps2)))
colnames(shannon_matrix) <- rownames(ps2@sam_data)
inverse_matrix <- as.data.frame(matrix(ncol = nsamples(ps2)))
colnames(inverse_matrix) <- rownames(ps2@sam_data)
pd_matrix <- as.data.frame(matrix(ncol = nsamples(ps2)))
colnames(pd_matrix) <- rownames(ps2@sam_data)

# Make an empty df to store the summary of all bootstraps, i.e., how many times there was a 
# statistically significant difference (Mann-Whitney) and which group was the median higher
alpha_summary <- as.data.frame(matrix(ncol = 6))
colnames(alpha_summary) <- c("shannon_p", "shannon_dir", "invsimpson_p", "invsimpson_dir", "pd_p", "pd_dir")

# Bootstrap loop
set.seed(123)  # for reproducibility
for (i in 1:n_bootstrap) {
   # Rarefy each sample to the target depth. In this case the lowest depth was chosen
  ps2_rare <- rarefy_even_depth(ps2, sample.size = min(sample_sums(ps2)))
  ps2_rare <- prune_taxa(taxa_sums(ps2_rare) > 10, ps2_rare)
  
  # Calculate alpha diversity metrics per sample
  shannon_matrix[i,] <- estimate_richness(ps2_rare, 
                                           measures = "Shannon")$Shannon
  inverse_matrix[i,] <- estimate_richness(ps2_rare, 
                                           measures = "InvSimpson")$InvSimpson
  pd_matrix[i,] <- picante::pd(t(ps2_rare@otu_table), ps2_rare@phy_tree, F)$PD
  
  temp_df <- data.frame(vc = ps2_rare@sam_data$vc,
                        cbind(t(shannon_matrix[i,]),
                        t(inverse_matrix[i,]),
                        t(pd_matrix[i,]))
                        )
  colnames(temp_df) <- c("Infection.status", "shannon", "invsimpson", "pd")
  
  row2add <- data.frame(shannon_p = ifelse(wilcox.test(shannon~Infection.status, data = temp_df)$p.value < 0.05, 
                                           TRUE, 
                                           FALSE),
                        shannon_dir = ifelse(median(temp_df$shannon[temp_df$Infection.status != "Positive"], na.rm = T) > median(temp_df$shannon[temp_df$Infection.status == "Positive"], na.rm = T), 
                                             "Lower in Positive", 
                                             "Higher in Positive"),
                        invsimpson_p = ifelse(wilcox.test(invsimpson~Infection.status, data = temp_df)$p.value < 0.05, 
                                              TRUE, 
                                              FALSE),
                        invsimpson_dir = ifelse(median(temp_df$invsimpson[temp_df$Infection.status != "Positive"], na.rm = T) > median(temp_df$invsimpson[temp_df$Infection.status == "Positive"], na.rm = T), 
                                                "Lower in Positive", 
                                                "Higher in Positive"),
                        pd_p = ifelse(wilcox.test(pd~Infection.status, data = temp_df)$p.value < 0.05, 
                                      TRUE, 
                                      FALSE),
                        pd_dir = ifelse(median(temp_df$pd[temp_df$Infection.status != "Positive"], na.rm = T) > median(temp_df$pd[temp_df$Infection.status == "Positive"], na.rm = T), 
                                        "Lower in Positive", 
                                        "Higher in Positive"))
  
  alpha_summary <- rbind(alpha_summary, row2add)
  
}

rm(temp_df, ps2_rare)

alpha_summary <- alpha_summary[-1,]
# Collapse to median per sample and join the vector competence column
alpha_boot_mean <- data.frame(Shannon = apply(t(shannon_matrix), 1, mean, na.rm = TRUE),
                                InverseSimpson = apply(t(inverse_matrix), 1, mean, na.rm = TRUE),
                                pd = apply(t(pd_matrix), 1, mean, na.rm = TRUE),
                                Treatment = sample_data(ps2)$Treatment)

alpha_boot_median <- data.frame(Shannon = apply(t(shannon_matrix), 1, median, na.rm = TRUE),
                              InverseSimpson = apply(t(inverse_matrix), 1, median, na.rm = TRUE),
                              pd = apply(t(pd_matrix), 1, median, na.rm = TRUE),
                              Treatment = sample_data(ps2)$Treatment)

# Plot a boxplot with jitter points for each metric using the median of each sample throughout the 100 bootstraps
# Shannon
alpha.sh <- ggplot(alpha_boot_median,
                   aes(x=Infection.status, 
                       y=Shannon, 
                       fill = Infection.status)) +
  geom_violin(width = .35, 
              alpha = .3) +
  geom_boxplot(width = .08, 
               alpha = .6,
               outlier.shape = NA) +
  geom_point(position = position_jitterdodge(jitter.width = .25),
             shape = 21, 
             color = "black",
             size = 2,
             alpha = .6) +
  scale_fill_manual(values = cblind_palette8[-3]) +
  theme_classic2() +
  ylab("Shannon") +
  geom_pwc(hide.ns = T,
           p.adjust.method = "none",
           bracket.nudge.y = -.05) +
  theme(legend.position = "none",
        axis.title.x = element_blank(),
        axis.text = element_text(size = 12),
        axis.title = element_text(size = 12))

# Inverse Simpson
alpha.is <- ggplot(alpha_boot_median,
                   aes(x=Infection.status, 
                       y=InverseSimpson, 
                       fill = Infection.status)) +
  geom_violin(width = .35, 
              alpha = .3) +
  geom_boxplot(width = .08, 
               alpha = .6,
               outlier.shape = NA) +
  geom_point(position = position_jitterdodge(jitter.width = .25),
             shape = 21, 
             color = "black",
             size = 2,
             alpha = .6) +
  scale_fill_manual(values = cblind_palette8[-3]) +
  theme_classic2() +
  ylab("Inverse Simpson") +
  geom_pwc(hide.ns = T,
           p.adjust.method = "none",
           bracket.nudge.y = -.05) +
  theme(legend.position = "none",
        axis.title.x = element_blank(),
        axis.text = element_text(size = 12),
        axis.title = element_text(size = 12))
# PD
alpha.pd <- ggplot(alpha_boot_median,
                   aes(x=Infection.status, 
                       y=pd, 
                       fill = Infection.status)) +
  geom_violin(width = .35, 
              alpha = .3) +
  geom_boxplot(width = .08, 
               alpha = .6,
               outlier.shape = NA) +
  geom_point(position = position_jitterdodge(jitter.width = .25),
             shape = 21, 
             color = "black",
             size = 2,
             alpha = .6) +
  scale_fill_manual(values = cblind_palette8[-3]) +
  theme_classic2() +
  ylab("Faith's phylogenetic diversity") +
  geom_pwc(hide.ns = T,
           p.adjust.method = "none",
           bracket.nudge.y = -.05) +
  theme(legend.position = "none",
        axis.title.x = element_blank(),
        axis.text = element_text(size = 12),
        axis.title = element_text(size = 12))

# Plot both boxplots side by side 
ggarrange(alpha.sh,
          alpha.is,
          alpha.pd,
          cols = 1
         )

# Beta diversity ####
{
  # Positive vs Negative
  ps2_beta <- ps2
  ps2_beta@otu_table <- ps2_beta@otu_table + .5
  ps2_philr <- philr(ps2_beta, part.weights='enorm.x.gm.counts', ilr.weights='blw.sqrt')
  philr_dist <- dist(ps2_philr, method="euclidean")
  
  philr.pcoa <- pcoa(philr_dist, correction = "cailliez")
  philr.pcoa.df <- as.data.frame(philr.pcoa$vectors[,1:2])
  colnames(philr.pcoa.df) <- paste0("PC", 1:2)
  philr.pcoa.df$Infection.status <- ps2_beta@sam_data$Infection.status
  
  # Beta diversity by timepoints 
  ggplot(data = philr.pcoa.df, 
         aes(x=PC1, y=PC2, 
             color=Infection.status)) +
    geom_point(size=3) + 
    stat_ellipse(linetype = 15, linewidth = 1.5) +
    theme_bw() +
    #facet_grid(~ dev_stage) +
    labs(title = "Euclidean distances with PhILR") +
    scale_color_manual(values= cblind_palette8[-3], name=NULL) +
    xlab(paste('PCo1 ', round(philr.pcoa$values$Relative_eig[1]*100, digits = 2), "%", sep = "")) +
    ylab(paste('PCo2 ', round(philr.pcoa$values$Relative_eig[2]*100, digits = 2), "%", sep = "")) +
    theme(aspect.ratio = 1)
  
  # PERMANOVA
  adonis.beta <- pairwise.adonis(as.matrix(philr_dist), 
                                 philr.pcoa.df$Infection.status, 
                                 p.adjust.m = "BH",
                                 sim.method = "euclidean")
  
  # ANOVA Dispersion test
  disp_philr <- vegan::permutest(vegan::betadisper(philr_dist, 
                                                   philr.pcoa.df$Infection.status),
                                 pairwise = T)
  
  adonis.beta$beta_disp <- unname(disp_philr$pairwise$permuted)
  
}

# DA analysis ####
# Prepare aggregated taxa according to taxonomic level of interest
# Change taxa names to wanted level
## Genus ####
ps2_genus <- aggregate_taxa(ps2, "Genus")
ps2_genus <- subset_taxa(ps2, Genus != "Unknown")
ps2_genus@sam_data$Infection.status <- factor(ps2_genus@sam_data$Infection.status, levels = c("Negative", "Positive"))

# LinDA 
linda.obj.genus <- linda(phyloseq.obj = ps2_genus, 
                         formula = "~ Infection.status + Feeding.status + other_variable2control", # If the exp design is longitudinal, you can use + (1|ID) in the formula 
                         alpha = 0.05, 
                         feature.dat.type = "count", 
                         p.adj.method = "BH")

# Prepare data frame for graphical visualization
# Genus
tick.da.genus <- linda.obj.genus$output$Infection.statusPositive
tick.da.genus$bacs <- row.names(tick.da.genus)
tick.da.filt <- subset(tick.da.genus[abs(tick.da.genus$log2FoldChange) > 1,], 
                       reject == TRUE)

# Plot bar plot of log2 fold change and bacteria which were statistically 
# significant for any of the groups
# Genus
ggplot(tick.da.filt, 
       aes(x=reorder(bacs, log2FoldChange, sort), 
           y=log2FoldChange, 
           fill=Infection.status)) +
  geom_col(aes(fill = factor(Infection.status)),
           position = position_dodge2(preserve = "single")) +
  scale_fill_manual(values=c("#D55E00", "#0072B2"), name = NULL, breaks = c("Negative", "Positive")) +
  theme_bw() +
  ylab(expression(Log[2]*" Fold change")) +
  theme(axis.title.y = element_blank(),
        panel.grid.major.y = element_blank(),
        strip.text = element_text(size = 11),
        axis.text = element_text(size = 11),
        legend.text = element_text(size = 11),
        plot.title = element_text(size = 12, hjust = 0.5)) +
  guides(fill=guide_legend(ncol=1)) +
  ggtitle("Differently abundant genera") +
  coord_flip()

# Functional analysis ####
{
  # Load kegg brite map for further characterization of pathways
  kegg_brite_map <- read.csv("~/Microbiome_metaanalysis/picrust1_KO_BRITE_map.tsv",
                             sep = "\t",
                             header = T)
  
  # Prepare mapping data frame to map Brite hierarchy levels 1 and 2 for each pathway
  pathway_list <- unlist(strsplit(kegg_brite_map[, "metadata_KEGG_Pathways"], "\\|"))
  pathway_list <- pathway_list[!duplicated(pathway_list)]
  
  path_map <- data.frame()
  
  for(pathway in 1:length(pathway_list)) {
    
    level1 <- strsplit(pathway_list[pathway], ";")[[1]][1]
    level2 <- strsplit(pathway_list[pathway], ";")[[1]][2]
    level3 <- strsplit(pathway_list[pathway], ";")[[1]][3]
    
    path_map <- rbind(path_map, 
                      c(level1, level2, level3))
    
  }
  
  colnames(path_map) <- c("Level1", "Level2", "Level3")
  
  # Manualy update some pathways without levels 1 and 2 information
  path_map <- rbind(path_map, 
                    c("Metabolism", "Glycan biosynthesis and metabolism", "Teichoic acid biosynthesis"),
                    c("Metabolism", "Global and overview maps", "Biosynthesis of amino acids"),
                    c("Metabolism", "Global and overview maps", "Biosynthesis of nucleotide sugars"),
                    c("Human Diseases", "Drug resistance: antimicrobial", "Vancomycin resistance"),
                    c("Metabolism", "Global and overview maps", "Biosynthesis of cofactors"),
                    c("Cellular Processes", "Cellular community - prokaryotes", "Quorum sensing"),
                    c("Metabolism", "Global and overview maps", "Carbon metabolism"),
                    c("Metabolism", "Metabolism of other amino acids", "D-Amino acid metabolism"),
                    c("Human Diseases", "Drug resistance: antimicrobial", "Cationic antimicrobial peptide (CAMP) resistance"),
                    c("Metabolism", "Glycan biosynthesis and metabolism", "Exopolysaccharide biosynthesis"),
                    c("Metabolism", "Amino acid metabolism", "Arginine biosynthesis"),
                    c("Cellular Processes", "Cellular community - prokaryotes", "Biofilm formation - Escherichia coli"),
                    c("Cellular Processes", "Cellular community - prokaryotes", "Biofilm formation - Pseudomonas aeruginosa"),
                    c("Cellular Processes", "Cellular community - prokaryotes", "Biofilm formation - Vibrio cholerae"),
                    c("Metabolism", "Metabolism of terpenoids and polyketides", "Pinene, camphor and geraniol degradation"),
                    c("Metabolism", "Global and overview maps", "Degradation of aromatic compounds"),
                    c("Metabolism", "Global and overview maps", "2-Oxocarboxylic acid metabolism"),
                    c("Human Diseases", "Neurodegenerative disease", "Parkinson disease"),
                    c("Human Diseases", "Cardiovascular disease", "Diabetic cardiomyopathy"),
                    c("Human Diseases", "Cancer: overview", "Chemical carcinogenesis - reactive oxygen species"),
                    c("Organismal Systems", "Environmental adaptation", "Thermogenesis"),
                    c("Human Diseases", "Neurodegenerative disease", "Prion disease"),
                    c("Human Diseases", "Neurodegenerative disease", "Huntington disease"),
                    c("Human Diseases", "Neurodegenerative disease", "Pathways of neurodegeneration - multiple diseases"),
                    c("Human Diseases", "Neurodegenerative disease", "Amyotrophic lateral sclerosis"),
                    c("Human Diseases", "Neurodegenerative disease", "Alzheimer disease"),
                    c("Organismal Systems", "Nervous system", "Retrograde endocannabinoid signaling"),
                    c("Human Diseases", "Endocrine and metabolic disease", "Non-alcoholic fatty liver disease"),
                    c("Human Diseases", "Endocrine and metabolic disease", "Insulin resistance"),
                    c("Metabolism", "Xenobiotics biodegradation and metabolism", "Steroid degradation"),
                    c("Metabolism", "Glycan biosynthesis and metabolism", "Lipoarabinomannan (LAM) biosynthesis"),
                    c("Metabolism", "Metabolism of terpenoids and polyketides", "Nonribosomal peptide structures"),
                    c("Metabolism", "Biosynthesis of other secondary metabolites", "Monobactam biosynthesis"),
                    c("Metabolism", "Lipid metabolism", "Fatty acid degradation"),
                    c("Metabolism", "Metabolism of cofactors and vitamins", "Porphyrin metabolism"),
                    c("Metabolism", "Metabolism of terpenoids and polyketides", "Biosynthesis of vancomycin group antibiotics"),
                    c("Metabolism", "Biosynthesis of other secondary metabolites", "Staurosporine biosynthesis"))
  
  # Load PICURSt2 results
  picrust2 <- read.delim("~/Microbiome_metaanalysis/qiime_outputs_new/picrust2_out_pipeline/pathways_out/path_abun_unstrat_descripKO.tsv", 
                         row.names = 1)[,-1]
  
  # Differential abundance for microbiome functional profile
  linda.fun <- linda(feature.dat = picrust2,
                     meta.dat = meta(ps2_cor),
                     formula = "~ vc + sex + origin + washing + (1|tick_species)",
                     feature.dat.type = "count")
  
  ko_comp <- linda.fun$output$vcCompetent[linda.fun$output$vcCompetent$pvalue < .05 & 
                                            abs(linda.fun$output$vcCompetent$log2FoldChange) > 1.4,]
  
  ko_comp.df <- data.frame(Ko = rownames(ko_comp),
                           Log2FC = ko_comp$log2FoldChange,
                           pvalue = ko_comp$pvalue,
                           padj = ko_comp$padj,
                           group = ifelse(ko_comp$log2FoldChange > 0, "Up regulated in Competent", "Down regulated in Competent")
  )
  
  ko_comp.df$Ko <- gsub("ko:", "", ko_comp.df$Ko)
  
  openxlsx::write.xlsx(ko_comp.df, 
                       file = "Results/Differentially_abundant_KEGGIDs.xlsx")
  
  # Enrichment
  ## List of compared clusters
  ck_comp <- clusterProfiler::compareCluster(Ko~group,
                                             data=ko_comp.df,
                                             pAdjustMethod = "BH",
                                             qvalueCutoff = 0.2,
                                             pvalueCutoff = 0.05,
                                             fun='enrichKO')
  
  # Make dotplot of enriched pathways using enrichplot package function
  path_data <- enrichplot::dotplot(ck_comp,
                                   showCategory= Inf,
                                   font.size = 5) + 
    ggplot2::theme(aspect.ratio = 1.8)
  
  # Retrieve data in a data frame to make a custom dotplot with ggplot2
  pathways_data <- data.frame(Pathway = as.character(path_data$data$Description),
                              GeneRatio = as.numeric(sapply(path_data$data$GeneRatio,function(x) eval(parse(text=x)))),
                              BH = path_data$data$p.adjust,
                              value = path_data$data$Cluster,
                              Associated_GeneID = path_data$data$geneID
  )
  
  # Order the pathways data dataframe
  pathways_data <- pathways_data[order(pathways_data$GeneRatio, decreasing = TRUE), ]
  
  # Map Brite hierarchy levels 1 and 2 for each enriched pathway
  for(lev3 in unique(pathways_data$Pathway)){
    pathways_data$level1[pathways_data$Pathway == lev3] <- path_map$Level1[path_map$Level3 == lev3][2]
    pathways_data$level2[pathways_data$Pathway == lev3] <- path_map$Level2[path_map$Level3 == lev3][2]
  }
  
  pathways_data$value <- gsub("\\n\\([0-9][0-9][0-9]\\)$","", 
                              gsub("\\n\\([0-9][0-9]\\)$", "", 
                                   pathways_data$value))
  
  pathways_data$value <- ifelse(pathways_data$value == "Up regulated in Competent", 
                                "Up", 
                                "Down")
  
  # Enriched Pathway dotplot. Use factor to sort the x axis according to time points. 
  # Value in parenthesis should be changed according to results.
  pathway_plot <- ggplot(pathways_data, aes(x = value,
                                            y = Pathway,
                                            fill = level1,
                                            color = level1)) +
    geom_point(alpha = 2,
               shape = 21, 
               size = 4) +
    # facet_wrap(~factor(comparison, levels=c("Dead-HIV",
    #                                        "Survivor-HIV")),
    #            scales = "free_x") +
    theme_bw() +
    theme(axis.title.y = element_blank(),
          axis.text.x = element_text(size = 14, angle = 45, hjust = 1),
          axis.text.y = element_text(size = 14),
          strip.text = element_text(size = 16),
          panel.grid.minor = element_blank(),
          legend.title = element_blank(),
          legend.position = "top",
          legend.justification='right',
          legend.text = element_text(size = 14)
    ) + 
    guides(fill = guide_legend(nrow = 2)) +
    labs(x = "", y = "")+
    scale_fill_manual(values = ggsci::pal_lancet(palette = c("lanonc"), alpha = 1)(5)) +
    scale_color_manual(values = ggsci::pal_lancet(palette = c("lanonc"), alpha = 1)(5))
  
  pathway_plot
  
}  
# Co-occurence network ####
### Building and plotting the network ####
# Change between the Competent and Non-competent groups

ASVs_2keep <- taxa_names(core(ps2, 
                         detection = 10, 
                         prevalence = .01, 
                         include.lowest = T
))

ps2_genus <- aggregate_taxa(prune_taxa(taxa_names(ps2) %in% ASVs_2keep, 
                                       ps2), 
                            "Genus")
ps2_genus <- subset_taxa(ps2_genus, Genus != "Unknown")
ps2_genus@sam_data$Infection.status <- as.factor(ps2_genus@sam_data$Infection.status)

ntaxa_min <- 1

ps2_genus_nymph <- prune_samples(ps2_genus@sam_data$life_stage == "Nymph", ps2_genus)

ps2_genus_pos_fed <- subset_samples(ps2_genus_nymph, blood_status == "Fed" & 
                                      Infection.status == "Positive")
ps2_genus_pos_fed <- prune_taxa(taxa_sums(ps2_genus_pos_fed) >= ntaxa_min, ps2_genus_pos_fed)
ps2_genus_pos_fed <- prune_samples(names(which(sample_sums(ps2_genus_pos_fed) >= 1)), ps2_genus_pos_fed)
ntaxa(ps2_genus_pos_fed)

ps2_genus_neg_fed <- subset_samples(ps2_genus_nymph, blood_status == "Fed" & 
                                      Infection.status != "Positive")
ps2_genus_neg_fed <- prune_taxa(taxa_sums(ps2_genus_neg_fed) >= ntaxa_min, ps2_genus_neg_fed)
ps2_genus_neg_fed <- prune_samples(names(which(sample_sums(ps2_genus_neg_fed) >= 1)), ps2_genus_neg_fed)
ntaxa(ps2_genus_neg_fed)

ps2_genus_pos_unfed <- subset_samples(ps2_genus_nymph, blood_status != "Fed" & 
                                      Infection.status == "Positive")
ps2_genus_pos_unfed <- prune_taxa(taxa_sums(ps2_genus_pos_unfed) >= ntaxa_min, ps2_genus_pos_unfed)
ps2_genus_pos_unfed <- prune_samples(names(which(sample_sums(ps2_genus_pos_unfed) >= 1)), ps2_genus_pos_unfed)
ntaxa(ps2_genus_pos_unfed)

ps2_genus_neg_unfed <- subset_samples(ps2_genus_nymph, blood_status != "Fed" & 
                                      Infection.status != "Positive")
ps2_genus_neg_unfed <- prune_taxa(taxa_sums(ps2_genus_neg_unfed) >= ntaxa_min, ps2_genus_neg_unfed)
ps2_genus_neg_unfed <- prune_samples(names(which(sample_sums(ps2_genus_neg_unfed) >= 1)), ps2_genus_neg_unfed)
ntaxa(ps2_genus_neg_unfed)

ps2_list <- list("Pos_Fed" = ps2_genus_pos_fed,
                 "Neg_Fed" = ps2_genus_neg_fed,
                 "Pos_Unfed" = ps2_genus_pos_unfed,
                 "Neg_Unfed" = ps2_genus_neg_unfed)

NoSleepR::with_nosleep({
  boot_network(ps2_list, 
               nboot = 100,
               outputDir = "boot_outputs",
               keep_taxa = 10
  )
})


# Key taxa ####
source("C:/Users/tiago/OneDrive/Microbiome_metaanalysis/HCIC_keyNodes2.R")

# start_time <- Sys.time()
competent_boot <- bootstrap_network_hcic(
  physeq = ps2_genus_comp,
  group_var = "vc",
  group_label = "Competent",
  n_bootstrap = 100,
  cor_threshold = 0.4,
  min_reads = 10,
  min_prevalence = 0.1,
  n_cpus = 15
)
# end_time <- Sys.time()
# 
# time_comp <- end_time - start_time
saveRDS(competent_boot, "C:/Users/tiago/OneDrive/Microbiome_metaanalysis/Net_boot_results/competent_boot.rds")

# 
# start_time <- Sys.time()
noncompetent_boot <- bootstrap_network_hcic(
  physeq = ps2_genus_ncomp,
  group_var = "vc",
  group_label = "Non_competent",
  n_bootstrap = 100,
  cor_threshold = 0.4,
  min_reads = 10,
  min_prevalence = 0.1,
  n_cpus = 15
)
# end_time <- Sys.time()
# 
# time_ncomp <- end_time - start_time

saveRDS(noncompetent_boot, "C:/Users/tiago/OneDrive/Microbiome_metaanalysis/Net_boot_results/noncompetent_boot.rds")

# After bootstrapping in background, RDSs need to be read to include in the current environment
competent_boot <- readRDS("~/Microbiome_metaanalysis/Net_boot_results/competent_boot.rds")
noncompetent_boot <- readRDS("~/Microbiome_metaanalysis/Net_boot_results/noncompetent_boot.rds")
competent_bootMetrics <- readRDS("~/Microbiome_metaanalysis/Net_boot_results/competent_boot_metrics.rds")
noncompetent_bootMetrics <- readRDS("~/Microbiome_metaanalysis/Net_boot_results/noncompetent_boot_metrics.rds")

# Compare groups
comparison <- compare_hcic_groups(competent_boot, 
                                  noncompetent_boot, 
                                  metric = "HCIC")

# View taxa that change keystone status
comparison[comparison$role_change == TRUE, ]

# Inspect results
print(competent_boot)

#

intersect_bac <- intersect(colnames(competent_boot$bootstrap_data$HCIC), 
                           colnames(noncompetent_boot$bootstrap_data$HCIC))

comp_bac_pars <- lapply(1:length(competent_boot$stable_taxa), function(i) {
  result <- sapply(competent_boot$bootstrap_data, function(mat) mat[, i])
  result <- as.matrix(result)
  colnames(result) <- names(competent_boot$bootstrap_data)
  rownames(result) <- paste0("Boot_comp", seq(1, nrow(result), 1))
  
  return(result)
})

# Name list with bac names
names(comp_bac_pars) <- colnames(competent_boot$bootstrap_data$HCIC)

# Impute the median of each variable for PCA
for(j in which(names(comp_bac_pars) %in% intersect_bac)) {
  for (u in 1:ncol(comp_bac_pars[[j]])) {
    comp_bac_pars[[j]][which(is.na(comp_bac_pars[[j]][,u])),u] <- median(comp_bac_pars[[j]][,u], na.rm = T)
  }
}

# Run a PCA for each bac and make a matrix of PC1s of each bac and bootstrap
for (bac.index in which(names(comp_bac_pars) %in% intersect_bac)) {
  bac_name <- names(comp_bac_pars)[bac.index]
  
  pca_bac <- mixOmics::pca(comp_bac_pars[[bac.index]])
  
  if(!exists("mat_comp")){
    mat_comp <- as.data.frame(pca_bac$x[,1]) 
  }else{
    col2add <- as.data.frame(pca_bac$x[,1]) 
    mat_comp <- cbind(mat_comp, col2add)
  }
  
}

names(mat_comp) <- names(comp_bac_pars[which(names(comp_bac_pars) %in% intersect_bac)])

ncomp_bac_pars <- lapply(1:length(noncompetent_boot$stable_taxa), function(i) {
  result <- sapply(noncompetent_boot$bootstrap_data, function(mat) mat[, i])
  result <- as.matrix(result)
  colnames(result) <- names(noncompetent_boot$bootstrap_data)
  rownames(result) <- paste0("Boot_ncomp", seq(1, nrow(result), 1))
  
  return(result)
})

# Name list with bac names
names(ncomp_bac_pars) <- colnames(noncompetent_boot$bootstrap_data$HCIC)

# Impute the median of each variable for PCA
for(j in which(names(ncomp_bac_pars) %in% intersect_bac)) {
  for (u in 1:ncol(ncomp_bac_pars[[j]])) {
    ncomp_bac_pars[[j]][which(is.na(ncomp_bac_pars[[j]][,u])),u] <- median(ncomp_bac_pars[[j]][,u], na.rm = T)
  }
}

# Run a PCA for each bac and make a matrix of PC1s of each bac and bootstrap
for (bac.index in which(names(ncomp_bac_pars) %in% intersect_bac)) {
  bac_name <- names(ncomp_bac_pars)[bac.index]
  
  pca_bac <- mixOmics::pca(ncomp_bac_pars[[bac.index]])
  
  if(!exists("mat_ncomp")){
    mat_ncomp <- as.data.frame(pca_bac$x[,1]) 
  }else{
    col2add <- as.data.frame(pca_bac$x[,1]) 
    mat_ncomp <- cbind(mat_ncomp, col2add)
  }
  
}

names(mat_ncomp) <- names(ncomp_bac_pars[which(names(ncomp_bac_pars) %in% intersect_bac)])

mat_keystone <- as.data.frame(cbind(t(mat_ncomp),t(mat_comp)))
meta_keystone <- data.frame("Sample" = colnames(mat_keystone),
                            "Class" = c(rep("Non Competent", 100), rep("Competent", 100)))

mdp_comp_vs_ncomp <- mdp(mat_keystone, meta_keystone, "Non Competent")

meta_keystone <- data.frame("Sample" = colnames(mat_keystone),
                            "Class" = c(rep("Competent", 100), rep("Non Competent", 100)))

mdp_ncomp_vs_comp <- mdp(mat_keystone, meta_keystone, "Non Competent")

abundance_comp <- list()
abundance_ncomp <- list()
prevalence_comp <- list()
prevalence_ncomp <- list()

for(chunk in 1:10){
  min_index <- (chunk*10)-9
  max_index <- chunk*10
  abundance_comp[[chunk]] <- colMeans(competent_boot$bootstrap_data$abundance[min_index:max_index,], na.rm = T)
  abundance_comp[[chunk]] <- colMeans(noncompetent_bootMetrics$bootstrap_data$abundance[min_index:max_index,], na.rm = T)
  prevalence_comp[[chunk]] <- colMeans(competent_boot$bootstrap_data$prevalence[min_index:max_index,], na.rm = T)
  prevalence_ncomp[[chunk]] <- colMeans(noncompetent_bootMetrics$bootstrap_data$prevalence[min_index:max_index,], na.rm = T)
}


