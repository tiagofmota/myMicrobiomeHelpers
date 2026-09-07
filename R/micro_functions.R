# Composition analysis. It returns a data frame with CLR transformed abundances ####
# aggregated by group and a heatmap from CLR z score abundances
compos <- function(phylo, 
                   level = "Genus", 
                   group, 
                   ht_subtittle = NULL, 
                   transform.method = "VST", 
                   transform.offset = .5, 
                   group.order,
                   color.range = c(min(scale(OTU.df.clr[,-ncol(OTU.df.clr)], center = F)), 0, 
                                   max(scale(OTU.df.clr[,-ncol(OTU.df.clr)], center = F)))){
  
  level <- match.arg(level, choices = c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species"))
  stopifnot("group should be a character vector" = is.vector(group, mode = "character"))
  transform.method <- match.arg(transform.method, choices = c("none", "CLR", "ILR", "VST"))
  stopifnot("transform.offset should be a number higher than 0" = isTRUE(is.numeric(transform.offset) & transform.offset > 0))
  stopifnot("group.order should be a character vector" = is.vector(group.order, mode = "character"))
  stopifnot("group.order is not included in group vector" = all(group.order %in% group))
  stopifnot("ht_subtittle should be a character vector to name the heatmap subtittle" = is.vector(group.order, mode = "character"))
  
  if(level == "Species"){
    pseq.fam <- aggregate_taxa(phylo, level)
    pseq.fam <- subset_taxa(pseq.fam, Species != "Unknown")
  }else if(level == "Genus"){
    pseq.fam <- aggregate_taxa(phylo, level)
    pseq.fam <- subset_taxa(pseq.fam, Genus != "Unknown")
  }else if(level == "Family"){
    pseq.fam <- aggregate_taxa(phylo, level)
    pseq.fam <- subset_taxa(pseq.fam, Family != "Unknown")
  }else if(level == "Order"){
    pseq.fam <- aggregate_taxa(phylo, level)
    pseq.fam <- subset_taxa(pseq.fam, Order != "Unknown")
  }else if(level == "Class"){
    pseq.fam <- aggregate_taxa(phylo, level)
    pseq.fam <- subset_taxa(pseq.fam, Class != "Unknown")
  }else if(level == "Phylum"){
    pseq.fam <- aggregate_taxa(phylo, level)
    pseq.fam <- subset_taxa(pseq.fam, Phylum != "Unknown")
  }else if(level == "Kingdom"){
    pseq.fam <- aggregate_taxa(phylo, level)
    pseq.fam <- subset_taxa(pseq.fam, Kingdom != "Unknown")
  }
  
  OTU.df <- data.frame(groups = group, t(pseq.fam@otu_table))
  
  if(transform.method == "VST") {
    
    for (i in colnames(OTU.df[,-1])) {
      OTU.df[,i] <- as.integer(OTU.df[,i])
    }
    
    OTU.df.vst <- DESeq2::varianceStabilizingTransformation(as.matrix(OTU.df[,-1]+1))
    OTU.df.vst <- data.frame(groups = OTU.df$groups, OTU.df.vst) |>
      group_by(groups) |> 
      summarise_all(median) |>
      as.data.frame()
    
    OTU.df <- OTU.df.vst
  }else if(transform.method == "CLR"){
    OTU.df <- OTU.df |>
      group_by(groups) |> 
      summarise_all(mean) |>
      as.data.frame()  
    
    OTU.df.clr <- logratio.transfo(OTU.df[,-1], 
                                   logratio = "CLR", 
                                   offset = transform.offset)
    class(OTU.df.clr) <- "matrix"
    OTU.df <- data.frame(groups = OTU.df$groups, OTU.df.clr)
  }else if(transform.method == "ILR"){
    OTU.df <- OTU.df |>
      group_by(groups) |> 
      summarise_all(median) |>
      as.data.frame()  
    
    OTU.df.ilr <- logratio.transfo(OTU.df[,-1], 
                                   logratio = "ILR", 
                                   offset = transform.offset)
    class(OTU.df.ilr) <- "matrix"
    OTU.df <- data.frame(groups = OTU.df$groups, OTU.df.ilr)
  }else if(transform.method == "none"){
    OTU.df <- OTU.df |>
      group_by(groups) |> 
      summarise_all(median) |>
      as.data.frame()  
  }
  
  ht <- Heatmap(t(scale(OTU.df[,-1], center = F)),
                name = paste0(transform.method, " abundances\n        Z-score"),
                col= circlize::colorRamp2(color.range,
                                          c("darkblue","darkgrey","yellow")),
                top_annotation = columnAnnotation(Groups = anno_text(OTU.df$groups, 
                                                                     just = "center", 
                                                                     rot = 0,
                                                                     location = .5,
                                                                     gp = gpar(border = "darkgrey", 
                                                                               lwd = 2,
                                                                               fill = "grey",
                                                                               col = "darkred"),
                                                                     height = max_text_height(group)*2)),
                heatmap_legend_param = list(direction = "horizontal", legend_width = unit(3, "cm")),
                column_split = factor(OTU.df$groups, levels = group.order),
                row_names_gp = grid::gpar(fontsize = 10),
                cluster_rows = T,
                row_labels = as.expression(lapply(colnames(OTU.df[,-1]), function(a) bquote(italic(.(a))))),
                cluster_columns = F,
                row_names_side = "right",
                show_column_names = T,
                row_title = NULL,
                column_title = ht_subtittle,
                # dsitances availiables: 
                # "euclidean" "maximum"   "manhattan" "canberra"  "binary"   
                # "minkowski" "pearson" "spearman" "kendall"
                clustering_distance_columns = "euclidean")
  
  setClass("Microb.composition", slots = c(transformDF="data.frame", 
                                           heatmap="Heatmap"))
  myObj <- new("Microb.composition", transformDF = OTU.df, heatmap = ht)
  
  return(myObj)
}

# Function for whole DA analysis microbiome ####
difab <- function(phylo, level = "Genus", formula, da.alpha = .05, da.dat.type = "count",
                  p.adj.method = "BH", groups, var, FCtreshold = 1){
  
  stopifnot("phylo should be a phyloseq object" = isTRUE(class(phylo)[1] == "phyloseq"))
  level <- match.arg(level, choices = c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species"))
  stopifnot("groups should be a two character vector, c(control,condition)" = is.vector(groups, mode = "character"))
  stopifnot("formula should be a character similar to a regression formula as ~condition1+characteristic2adjust+(1|ID)" = is.character(formula))
  stopifnot("da.alpha should be a number between 0 and 1" = isTRUE(is.numeric(da.alpha) & da.alpha > 0 & da.alpha <= 1))
  da.dat.type <- match.arg(da.dat.type, choices = c("count", "proportion"))
  p.adj.method <- match.arg(p.adj.method, choices = c("holm", "hochberg", "hommel", "bonferroni", "BH", "BY", "fdr", "none"))
  stopifnot("var should be a character of variable of interest colname the same way as in formula" = is.character(var))
  stopifnot("FCtreshold should be a number higher than 0" = isTRUE(is.numeric(FCtreshold) & FCtreshold > 0))
  
  if(level == "Species"){
    pseq.fam <- aggregate_taxa(phylo, level)
    pseq.fam <- subset_taxa(pseq.fam, Species != "Unknown")
    level.plural <- "Species"
  }else if(level == "Genus"){
    pseq.fam <- aggregate_taxa(phylo, level)
    pseq.fam <- subset_taxa(pseq.fam, Genus != "Unknown")
    level.plural <- "Genera"
  }else if(level == "Family"){
    pseq.fam <- aggregate_taxa(phylo, level)
    pseq.fam <- subset_taxa(pseq.fam, Family != "Unknown")
    level.plural <- "Families"
  }else if(level == "Order"){
    pseq.fam <- aggregate_taxa(phylo, level)
    pseq.fam <- subset_taxa(pseq.fam, Order != "Unknown")
    level.plural <- "Orders"
  }else if(level == "Class"){
    pseq.fam <- aggregate_taxa(phylo, level)
    pseq.fam <- subset_taxa(pseq.fam, Class != "Unknown")
    level.plural <- "Classes"
  }else if(level == "Phylum"){
    pseq.fam <- aggregate_taxa(phylo, level)
    pseq.fam <- subset_taxa(pseq.fam, Phylum != "Unknown")
    level.plural <- "Phyla"
  }else if(level == "Kingdom"){
    pseq.fam <- aggregate_taxa(phylo, level)
    pseq.fam <- subset_taxa(pseq.fam, Kingdom != "Unknown")
    level.plural <- "Kingdoms"
  }
  
  pseq.fam = prune_samples(names(which(sample_sums(pseq.fam) > 0)), pseq.fam)
  
  linda.obj <- linda(phyloseq.obj = pseq.fam, formula = formula, 
                     alpha = da.alpha, feature.dat.type = da.dat.type, 
                     p.adj.method = p.adj.method)
  
  da.output <- paste0(var, groups[2])
  
  da.df <- as.data.frame(linda.obj$output[da.output])
  colnames(da.df) <- c("baseMean", "log2FoldChange", "lfcSE", "stat",
                       "pvalue", "padj", "reject", "df")
  da.df$bacs <- row.names(da.df)
  da.df$varname[da.df$log2FoldChange < 0] <- groups[1]
  da.df$varname[da.df$log2FoldChange > 0] <- groups[2]
  da.df <- as.data.frame(da.df)
  
  da.plot <- ggplot(subset(da.df[da.df$log2FoldChange < -FCtreshold | 
                                   da.df$log2FoldChange > FCtreshold,], reject == TRUE), 
                    aes(x=reorder(bacs, log2FoldChange, sort), y=log2FoldChange, fill=varname)) +
    geom_col(aes(fill = factor(varname)),
             position = position_dodge2(preserve = "single")) +
    scale_fill_manual(values=c("#0072B2", "#D55E00"), name = NULL, breaks = c(groups[1], groups[2])) +
    theme_bw() +
    ylab(expression(Log[2]*" Fold change")) +
    theme(axis.title.y = element_blank(),
          panel.grid.major.y = element_blank(),
          strip.text = element_text(size = 11),
          axis.text = element_text(size = 11),
          legend.text = element_text(size = 11),
          plot.title = element_text(size = 12, hjust = 0.5),
          axis.text.y = element_text(face = "italic")) +
    guides(fill=guide_legend(ncol=1)) +
    ggtitle(paste("Differently abundant", level.plural)) +
    coord_flip()
  
  setOldClass("gg")
  setClass("Microb.da", slots = c(daResult="data.frame",
                                  daPlot="gg",
                                  linda.obj = "list"))
  
  myObj <- new("Microb.da", daResult = da.df, daPlot = da.plot, linda.obj = linda.obj)
  
  return(myObj)
  
}

# Export ASV table to biom for picrust2 ####
write_biom_csv <- function(ps, file, sep = "; ") {
  phyloseq::otu_table(ps) %>%
    as.data.frame() %>%
    rownames_to_column("#OTU ID") %>%
    left_join(phyloseq::tax_table(ps) %>% 
                as.data.frame() %>%
                rownames_to_column("#OTU ID") %>% 
                tidyr::unite("taxonomy", !`#OTU ID`, sep = sep)) -> phyloseq_biom
  
  readr::write_tsv(phyloseq_biom, file = file)
}

## Pairwise Adonis - PERMANOVA ####

pairwise.adonis = function(x,factors, sim.function = 'vegdist', sim.method = 'bray', p.adjust.m ='bonferroni')
{
  #library(vegan)
  
  co = combn(unique(sort(as.character(factors))),2)
  pairs = c()
  F.Model =c()
  R2 = c()
  p.value = c()
  
  
  for(elem in 1:ncol(co)){
    if(sim.function == 'daisy'){
      library(cluster); x1 = daisy(x[factors %in% c(co[1,elem],co[2,elem]),],metric=sim.method)
    } else{x1 = vegan::vegdist(x[factors %in% c(co[1,elem],co[2,elem]),],method=sim.method)}
    
    ad = vegan::adonis2(x1 ~ factors[factors %in% c(co[1,elem],co[2,elem])] );
    pairs = c(pairs,paste(co[1,elem],'vs',co[2,elem]));
    F.Model =c(F.Model,ad$F[1]);
    R2 = c(R2,ad$R2[1]);
    p.value = c(p.value,ad$`Pr(>F)`[1])
  }
  p.adjusted = p.adjust(p.value,method=p.adjust.m)
  sig = c(rep('',length(p.adjusted)))
  sig[p.adjusted <= 0.05] <-'.'
  sig[p.adjusted <= 0.01] <-'*'
  sig[p.adjusted <= 0.001] <-'**'
  sig[p.adjusted <= 0.0001] <-'***'
  
  pairw.res = data.frame(pairs,F.Model,R2,p.value,p.adjusted,sig)
  print("Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1")
  return(pairw.res)
  
}

## p-value matrix p adjustment ####
mat_padjust <- function(x, method = "BH", ...) {
  stopifnot(class(x)[1] == "matrix")
  x[upper.tri(x)] <- p.adjust(x[upper.tri(x)], method = method)
  x[lower.tri(x)] <- p.adjust(x[lower.tri(x)], method = method)
  x[is.na(x)] <- 1
  return(x)
}

# Bootstrap networks using ####
boot_network <- function(phy_list, nboot, outputDir, keep_taxa, cpus = NULL, ...){
  
  if (missing(phy_list)) {
    stop("Error: 'phy_list' is missing.\n",
         "  -> Expected: A named list containing one or more 'phyloseq' objects (one per experimental group) where names of elements represents the group names.")
  }
  
  if (missing(nboot)) {
    stop("Error: 'nboot' is missing.\n",
         "  -> Expected: An integer specifying the number of bootstrap iterations (e.g., nboot = 100).")
  }
  
  if (missing(outputDir)) {
    stop("Error: 'outputDir' is missing.\n",
         "  -> Expected: A character string indicating the directory path where results and checkpoints should be saved.")
  }
  
  if (missing(keep_taxa)) {
    stop("Error: 'keep_taxa' is missing.\n",
         "  -> Expected: A numeric threshold for filtering ASVs. Taxa with total counts below this value within a bootstrap iteration will be discarded.")
  }
  
  available_cores <- parallel::detectCores()
  
  if (is.null(cpus)) {
    # If cpus is left blank, automatically assign 80% of system capacity
    recommended_cpus <- floor(available_cores * 0.8)
    if (recommended_cpus < 1) recommended_cpus <- 1
    
    message(sprintf("Notice: 'cpus' parameter not specified. Automatically utilizing 80%% of available cores (%d/%d).", 
                    recommended_cpus, available_cores))
    cpus <- recommended_cpus
    
  } else {
    # Validate manually entered cpus input
    if (!is.numeric(cpus) || cpus %% 1 != 0 || cpus <= 0) {
      stop("Error in 'cpus': Must be a positive integer representing CPU cores. Received value: ", cpus)
    }
    
    # Check if manually entered cpus input exceeds machine limits
    if (cpus > available_cores) {
      recommended_cpus <- floor(available_cores * 0.8)
      if (recommended_cpus < 1) recommended_cpus <- 1
      
      warning(sprintf(
        "Requested cpus (%d) exceeds available system cores (%d).\n  -> Automatically capping allocation to 80%% capacity: using %d cpus instead.",
        cpus, available_cores, recommended_cpus
      ), immediate. = TRUE)
      
      cpus <- recommended_cpus
    }
  }
  
  # Validate phy_list (Must be a list and MUST be named)
  if (!is.list(phy_list)) {
    stop("Error in 'phy_list': The provided object is not a list.")
  }
  if (is.null(names(phy_list)) || any(names(phy_list) == "")) {
    stop("Error in 'phy_list': The list must be named to identify experimental groups (e.g., list(GroupA = physeq1, GroupB = physeq2)).")
  }
  
  # Validate nboot (Must be numeric, a whole integer, and positive)
  if (!is.numeric(nboot) || nboot %% 1 != 0 || nboot <= 0) {
    stop("Error in 'nboot': Must be a positive integer (e.g., 100). Received value: ", nboot)
  }
  
  # Validate keep_taxa (Must be a single positive numeric value or zero)
  if (!is.numeric(keep_taxa) || length(keep_taxa) != 1 || keep_taxa < 0) {
    stop("Error in 'keep_taxa': Must be a single numeric value greater than or equal to zero.")
  }
  
  # Validate outputDir (Must be a valid single character string)
  if (!is.character(outputDir) || length(outputDir) != 1 || outputDir == "") {
    stop("Error in 'outputDir': Must be a valid character string pointing to a directory path.")
  }
  
  corMicro_args <- list(...)
  
  if (!"method" %in% names(corMicro_args)) corMicro_args$method <- "sparcc"
  if (!"r.threshold" %in% names(corMicro_args)) corMicro_args$r.threshold <- 0.4
  if (!"R" %in% names(corMicro_args)) corMicro_args$R <- 100
  if (!"ncpus" %in% names(corMicro_args)) corMicro_args$ncpus <- cpus
  
  n_ASVtables <- length(phy_list)
  group_names <- names(phy_list)
  check <- "__checkpoints"
  dir.create(file.path(outputDir), showWarnings = FALSE)
  dir.create(file.path(outputDir, check), showWarnings = FALSE)
  
  cat("===Starting network bootstraps===", 
      paste0(">>> Number of ASV tables: ", n_ASVtables),
      paste0(">>> Number of bootstraps: ", nboot),
      paste0(">>> Number of CPUs to be used: ", cpus),
      paste0(">>> Group names identified: ", paste(group_names, collapse = ", ")),
      sep = "\n")
  
  for(group in group_names){
    if(file.exists(file.path(outputDir, check, group))){
      cat("Group", group, 
          "has already been done previously. Jumping to next:", 
          group_names[which(group_names %in% group) + 1])
      next
    }else{
      # Extract OTU table and metadata
      phy_obj <- phy_list[[which(group_names %in% group)]]
      otu_mat <- as(otu_table(phy_obj), "matrix")
      if (taxa_are_rows(phy_obj)) {
        otu_mat <- t(otu_mat)
      }
      
      metadata <- as(sample_data(phy_obj), "data.frame")
      
      original_taxa <- colnames(otu_mat)
      n_samples_orig <- nrow(otu_mat)
      
      # Progress tracking
      pb <- txtProgressBar(min = 1, max = nboot, style = 3)
      
      # Bootstrap ASV tables to make 100 bootstrapped correlation matrices with sparcc
      occor_ref = NULL
      corR_ref = list()
      count_ref = list()
      for(i in 1:nboot){
        
        # Resample samples with replacement
        set.seed(i)
        boot_idx <- sample(1:n_samples_orig, size = n_samples_orig, replace = TRUE)
        
        # Get bootstrapped OTU table
        otu_boot <- otu_mat[boot_idx, , drop = FALSE]
        rownames(otu_boot) <- paste0("BootSample_", 1:n_samples_orig, "_", i)
        
        taxa2keep <- colSums(otu_boot) >= keep_taxa
        
        # Filter OTU table
        otu_boot_filtered <- otu_boot[, taxa2keep, drop = FALSE]
        
        # Create phyloseq object for this bootstrap
        metadata_boot <- metadata[boot_idx, , drop = FALSE]
        rownames(metadata_boot) <- rownames(otu_boot_filtered)
        
        physeq_boot <- phyloseq(
          otu_table(otu_boot_filtered, taxa_are_rows = FALSE),
          sample_data(metadata_boot)
        )
        
        # 2. Combine the dynamic phyloseq object with the captured ellipsis args
        # The primary position argument for corMicro is 'ps'
        final_cor_args <- c(list(ps = physeq_boot), corMicro_args)
        
        # 3. Call corMicro dynamically using do.call
        occor_ref <- do.call(ggClusterNet::corMicro, final_cor_args)
        
        r_thresh <- corMicro_args$r.threshold
        
        corR_ref[[i]] = ifelse(abs(occor_ref[[1]]) > r_thresh & occor_ref[[4]] < 0.05, occor_ref[[1]], 0 )
        count_ref[[i]] = otu_boot_filtered
        
        setTxtProgressBar(pb, i)
      }
      
      close(pb)
      
      cat("===Exporting list of networks and boostrapped ASV tables for group:", group, "===")
      saveRDS(corR_ref, file.path(outputDir, paste0("corR_", group, "_ref.rds")))
      saveRDS(count_ref, file.path(outputDir, paste0("count_", group, "_ref.rds")))
      
      if(file.exists(file.path(outputDir, paste0("corR_", group, "_ref.rds")))){
        cat("===List of networks successfully exported for group:", group, "===")
        file.create(file.path(outputDir, check, paste0(group, "_Nets_OK.txt")))
      }else{
        cat(paste0("===Some error to export list of networks for group: ", group, "==="),
            "===Please check error message and repeat from this group===",
            sep = "\n")
        quit(status = 1)
      }
    }
  }
}

Add microbiome helper function
