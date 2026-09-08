# Composition analysis. It returns a data frame with CLR transformed abundances aggregated by group and a heatmap from CLR z score abundances
# Define S4 return container
setClass("Microb.composition", 
         slots = c(transformDF = "data.frame", heatmap = "ANY"))

compos <- function(phylo, 
                   level = "Genus", 
                   group, 
                   ht_subtittle = NULL, 
                   transform.method = "VST", 
                   transform.offset = .5, 
                   group.order = NULL,
                   color.range = NULL) {
  
  # Input validation and defensive programming
  if (missing(phylo)) stop("Error: 'phylo' dataset is missing.")
  if (missing(group)) stop("Error: 'group' matching vector or sample metadata column is missing.")
  
  # Clean parameter matching
  level <- match.arg(level, choices = c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species"))
  transform.method <- match.arg(transform.method, choices = c("none", "CLR", "ILR", "VST"))
  
  # Extract or validate group mappings
  if (is.character(group) && length(group) == 1) {
    # If a column name string is passed, extract from phyloseq sample data
    meta_df <- as(phyloseq::sample_data(phylo), "data.frame")
    if (!group %in% colnames(meta_df)) {
      stop(sprintf("Error: Sample column '%s' not found inside phyloseq metadata.", group))
    }
    group_vec <- as.character(meta_df[[group]])
  } else if (is.vector(group) && (is.character(group) || is.factor(group))) {
    group_vec <- as.character(group)
    if (length(group_vec) != phyloseq::nsamples(phylo)) {
      stop("Error: Provided 'group' vector length does not match sample count in phyloseq object.")
    }
  } else {
    stop("Error: 'group' must be a character string metadata name or an explicit matching group vector.")
  }
  
  if (!is.numeric(transform.offset) || length(transform.offset) != 1 || transform.offset <= 0) {
    stop("Error: 'transform.offset' parameter must be a positive numeric value.")
  }
  
  # Standardize unique groups
  unique_groups <- unique(group_vec)
  if (is.null(group.order)) {
    group.order <- sort(unique_groups)
  } else {
    if (!is.vector(group.order, mode = "character") || !all(group.order %in% unique_groups)) {
      stop("Error: 'group.order' values must perfectly match the values contained inside the group identifier matrix.")
    }
  }
  
  if (!is.null(ht_subtittle) && (!is.character(ht_subtittle) || length(ht_subtittle) != 1)) {
    stop("Error: 'ht_subtittle' must be a single descriptive character string.")
  }

  # Unified taxonomy aggregation and cleaning
  pseq.fam <- microbiome::aggregate_taxa(phylo, level)
  
  # Filter out "Unknown" classifications
  filter_expr <- sprintf("%s != 'Unknown'", level)
  pseq.fam <- phyloseq::subset_taxa(pseq.fam, eval(parse(text = filter_expr)))
  
  # Transformations and groups summarization
  otu_mat <- as(phyloseq::otu_table(pseq.fam), "matrix")
  if (phyloseq::taxa_are_rows(pseq.fam)) {
    otu_mat <- t(otu_mat)
  }
  
  OTU.df <- data.frame(groups = group_vec, otu_mat, check.names = FALSE)
  
  if (transform.method == "VST") {
    # Coerce to integer values safely for DESeq2 matrix compatibility
    numeric_cols <- colnames(OTU.df)[-1]
    for (col in numeric_cols) {
      OTU.df[[col]] <- as.integer(OTU.df[[col]])
    }
    
    # Calculate regularized transformations
    vst_matrix <- DESeq2::varianceStabilizingTransformation(as.matrix(OTU.df[, -1] + 1))
    OTU.df <- data.frame(groups = OTU.df$groups, vst_matrix, check.names = FALSE) |>
      dplyr::group_by(groups) |> 
      dplyr::summarise_all(mean) |>
      as.data.frame()
    
  } else if (transform.method %in% c("CLR", "ILR")) {
    # Aggregate data profiles first
    OTU.df <- OTU.df |>
      dplyr::group_by(groups) |> 
      dplyr::summarise_all(mean) |>
      as.data.frame()  
    
    # Run targeted logratio transformation from mixOmics package
    transformed_mat <- mixOmics::logratio.transfo(as.matrix(OTU.df[, -1]), 
                                                   logratio = transform.method, 
                                                   offset = transform.offset)
    transformed_mat <- as.matrix(transformed_mat)
    OTU.df <- data.frame(groups = OTU.df$groups, transformed_mat, check.names = FALSE)
    
  } else if (transform.method == "none") {
    OTU.df <- OTU.df |>
      dplyr::group_by(groups) |> 
      dplyr::summarise_all(dplyr::median) |>
      as.data.frame()  
  }
  
  # Color range configuration and heatmap generation
  scaled_matrix <- t(scale(OTU.df[, -1], center = FALSE))
  
  # Dynamically configure default color boundaries if not specified by the user
  if (is.null(color.range)) {
    color.range <- c(min(scaled_matrix, na.rm = TRUE), 0, max(scaled_matrix, na.rm = TRUE))
  }
  
  ht <- ComplexHeatmap::Heatmap(
    matrix = scaled_matrix,
    name = paste0(transform.method, " abundances\nZ-score"),
    col = circlize::colorRamp2(color.range, c("darkblue", "darkgrey", "yellow")),
    top_annotation = ComplexHeatmap::columnAnnotation(
      Groups = ComplexHeatmap::anno_text(
        OTU.df$groups, just = "center", rot = 0, location = 0.5,
        gp = grid::gpar(border = "darkgrey", lwd = 2, fill = "grey", col = "darkred"),
        height = grid::max_text_height(OTU.df$groups) * 2
      )
    ),
    heatmap_legend_param = list(direction = "horizontal", legend_width = grid::unit(3, "cm")),
    column_split = factor(OTU.df$groups, levels = group.order),
    row_names_gp = grid::gpar(fontsize = 10),
    cluster_rows = TRUE,
    row_labels = as.expression(lapply(rownames(scaled_matrix), function(a) bquote(italic(.(a))))),
    cluster_columns = FALSE,
    row_names_side = "right",
    show_column_names = TRUE,
    row_title = NULL,
    column_title = ht_subtittle,
    clustering_distance_columns = "euclidean"
  )
  
  # Build S4 response payload object
  myObj <- new("Microb.composition", transformDF = OTU.df, heatmap = ht)
  return(myObj)
}

# Function for differential abundance analysis of the microbiome using the LinDA method
setOldClass("gg")
setClass("Microb.da", 
         slots = c(daResult  = "data.frame",
                   daPlot    = "gg",
                   linda.obj = "list"))

difab <- function(phylo, 
                  level = "Genus", 
                  formula, 
                  da.alpha = .05, 
                  da.dat.type = "count",
                  p.adj.method = "BH", 
                  groups, 
                  var, 
                  FCtreshold = 1,
                  cpus = NULL) {
  
  # Input validation and defensive programming
  if (missing(phylo)) stop("Error: 'phylo' dataset parameter is missing.")
  if (!inherits(phylo, "phyloseq")) stop("Error: 'phylo' must be a valid phyloseq object.")
  if (missing(groups) || !is.character(groups) || length(groups) != 2) {
    stop("Error: 'groups' must be a character vector of length 2 containing exactly: c('control', 'condition').")
  }
  if (missing(formula)) stop("Error: Modeling 'formula' string or expression is missing.")
  if (missing(var) || !is.character(var) || length(var) != 1) {
    stop("Error: 'var' must be a single character string matching the main variable column name in your metadata.")
  }

  level <- match.arg(level, choices = c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species"))
  da.dat.type <- match.arg(da.dat.type, choices = c("count", "proportion"))
  p.adj.method <- match.arg(p.adj.method, choices = c("holm", "hochberg", "hommel", "bonferroni", "BH", "BY", "fdr", "none"))
  
  if (!is.numeric(da.alpha) || length(da.alpha) != 1 || da.alpha <= 0 || da.alpha > 1) {
    stop("Error: 'da.alpha' must be a numeric threshold between 0 and 1.")
  }
  if (!is.numeric(FCtreshold) || length(FCtreshold) != 1 || FCtreshold < 0) {
    stop("Error: 'FCtreshold' must be a positive numeric value.")
  }

  # Dynamic CPU cores configurations (80% Safety Cap)
  available_cores <- parallel::detectCores()
  
  if (is.null(cpus)) {
    recommended_cpus <- floor(available_cores * 0.8)
    if (recommended_cpus < 1) recommended_cpus <- 1
    message(sprintf("Notice: 'cpus' parameter not specified. Automatically utilizing 80%% of available cores (%d/%d) for LinDA calculations.", 
                    recommended_cpus, available_cores))
    cpus <- recommended_cpus
  } else {
    if (!is.numeric(cpus) || cpus %% 1 != 0 || cpus <= 0) {
      stop("Error in 'cpus': Must be a positive integer representing CPU cores. Received value: ", cpus)
    }
    if (cpus > available_cores) {
      recommended_cpus <- floor(available_cores * 0.8)
      if (recommended_cpus < 1) recommended_cpus <- 1
      warning(sprintf("Requested cpus (%d) exceeds available system cores (%d).\n  -> Automatically capping allocation to 80%% capacity: using %d cpus instead.",
                      cpus, available_cores, recommended_cpus), immediate. = TRUE)
      cpus <- recommended_cpus
    }
  }

  # Map plural tags for clean plot titling
  plural_mappings <- c(Kingdom = "Kingdoms", Phylum = "Phyla", Class = "Classes", 
                       Order = "Orders", Family = "Families", Genus = "Genera", Species = "Species")
  level.plural <- plural_mappings[level]

  # Taxonomy aggregation and zero-omission filtering
  pseq.fam <- microbiome::aggregate_taxa(phylo, level)
  
  filter_expr <- sprintf("%s != 'Unknown'", level)
  pseq.fam <- phyloseq::subset_taxa(pseq.fam, eval(parse(text = filter_expr)))
  pseq.fam <- phyloseq::prune_samples(phyloseq::sample_sums(pseq.fam) > 0, pseq.fam)
  
  # LinDA regression models and data parsing
  formal_formula <- if (is.character(formula)) as.formula(formula) else formula

  # n.cores is passed dynamically here to MicrobiomeStat::linda
  linda.obj <- MicrobiomeStat::linda(
    phyloseq.obj     = pseq.fam, 
    formula          = formal_formula, 
    alpha            = da.alpha, 
    feature.dat.type = da.dat.type, 
    p.adj.method     = p.adj.method,
    n.cores          = cpus
  )
  
  da.output <- paste0(var, groups[2])
  
  if (!da.output %in% names(linda.obj$output)) {
    stop(sprintf("Error: Coefficient output matrix named '%s' was not found in LinDA results. Check your variable naming.", da.output))
  }
  
  da.df <- as.data.frame(linda.obj$output[[da.output]])
  colnames(da.df) <- c("baseMean", "log2FoldChange", "lfcSE", "stat", "pvalue", "padj", "reject", "df")
  
  da.df$bacs <- row.names(da.df)
  da.df$varname <- ifelse(da.df$log2FoldChange < 0, groups[1], groups[2])
  
  # Plotting pipeline
  plot_data <- subset(da.df, reject == TRUE & (log2FoldChange < -FCtreshold | log2FoldChange > FCtreshold))
  
  da.plot <- ggplot2::ggplot(plot_data, ggplot2::aes(x = reorder(bacs, log2FoldChange), y = log2FoldChange, fill = varname)) +
    ggplot2::geom_col(position = ggplot2::position_dodge2(preserve = "single")) +
    ggplot2::scale_fill_manual(values = c("#0072B2", "#D55E00"), name = NULL, breaks = groups) +
    ggplot2::theme_bw() +
    ggplot2::labs(y = expression(Log[2]*" Fold change"), x = NULL, title = paste("Differently abundant", level.plural)) +
    ggplot2::theme(
      panel.grid.major.y = ggplot2::element_blank(),
      strip.text         = ggplot2::element_text(size = 11),
      axis.text          = ggplot2::element_text(size = 11),
      legend.text        = ggplot2::element_text(size = 11),
      plot.title         = ggplot2::element_text(size = 12, hjust = 0.5),
      axis.text.y        = ggplot2::element_text(face = "italic")
    ) +
    ggplot2::guides(fill = ggplot2::guide_legend(ncol = 1)) +
    ggplot2::coord_flip()
  
  myObj <- new("Microb.da", daResult = da.df, daPlot = da.plot, linda.obj = linda.obj)
  return(myObj)
}

# Export ASV table to biom for picrust2 using any phyloseq object
write_biom_csv <- function(ps, file) {
  phyloseq::otu_table(ps) %>%
    as.data.frame() %>%
    rownames_to_column("#OTU ID") %>%
    left_join(phyloseq::tax_table(ps) %>% 
                as.data.frame() %>%
                rownames_to_column("#OTU ID") %>% 
                tidyr::unite("taxonomy", !`#OTU ID`, sep = "; ")) -> phyloseq_biom
  
  readr::write_tsv(phyloseq_biom, file = file)
}

## Pairwise Adonis - PERMANOVA
pairwise.adonis <- function(x, 
                            factors, 
                            p.adjust.m = "bonferroni") {
  
  # Input validation and defensive programming
  if (missing(x)) stop("Error: Input distance matrix 'x' is missing.")
  if (missing(factors)) stop("Error: Grouping vector 'factors' is missing.")
  
  p.adjust.m <- match.arg(p.adjust.m, choices = c("holm", "hochberg", "hommel", "bonferroni", "BH", "BY", "fdr", "none"))
  
  # Ensure the input is converted/treated as a full matrix for row/col indexing
  full_dist_mat <- as.matrix(x)
  factors_vec   <- as.character(factors)
  
  if (nrow(full_dist_mat) != length(factors_vec)) {
    stop(sprintf("Error: Distance matrix dimensions (%dx%d) do not match 'factors' vector length (%d).", 
                 nrow(full_dist_mat), ncol(full_dist_mat), length(factors_vec)))
  }
  
  if (nrow(full_dist_mat) != ncol(full_dist_mat)) {
    stop("Error: Input 'x' must be a square distance matrix or a valid 'dist' object.")
  }
  
  unique_groups <- unique(sort(factors_vec))
  if (length(unique_groups) < 2) {
    stop("Error: 'factors' must contain at least 2 unique groups to perform pairwise comparisons.")
  }
  
  # Combination generation and memory pre-allocation
  co <- combn(unique_groups, 2)
  n_comparisons <- ncol(co)
  
  pairs_vec <- vector("character", n_comparisons)
  F_Model   <- vector("numeric", n_comparisons)
  R2        <- vector("numeric", n_comparisons)
  p_value   <- vector("numeric", n_comparisons)
  
  # Pairwise distance matrix subsetting and adonis loop
  for (elem in 1:n_comparisons) {
    group1 <- co[1, elem]
    group2 <- co[2, elem]
    
    # Identify which samples belong to the current pair
    keep_indices <- factors_vec %in= c(group1, group2)
    
    # Subset BOTH rows and columns of the distance matrix to isolate the pair
    sub_dist_mat <- full_dist_mat[keep_indices, keep_indices, drop = FALSE]
    sub_factors  <- factors_vec[keep_indices]
    
    # Convert back to a 'dist' object as required by adonis2
    sub_dist_obj <- as.dist(sub_dist_mat)
    
    # Execute adonis2 directly on the subsetted distance structure
    ad <- vegan::adonis2(sub_dist_obj ~ sub_factors)
    
    pairs_vec[elem] <- paste(group1, "vs", group2)
    F_Model[elem]   <- ad$F[1]
    R2[elem]        <- ad$R2[1]
    p_value[elem]   <- ad$`Pr(>F)`[1]
  }
  
  # Multivariate homogeneity of groups dispersions on the full distance object
  full_dist_obj <- as.dist(full_dist_mat)
  mod_disp      <- vegan::betadisper(full_dist_obj, group = factors_vec)
  disp_perm     <- vegan::permutest(mod_disp, pairwise = TRUE)
  
  # Extract raw pairwise dispersion results and map names
  disp_p_matrix <- disp_perm$pairwise$permuted
  disp_lookup   <- setNames(unname(disp_p_matrix), gsub("-", " vs ", rownames(disp_p_matrix)))
  
  # Safely match labels back into our pre-allocated vector
  beta_disp_p <- vector("numeric", n_comparisons)
  for (j in seq_along(pairs_vec)) {
    alt_pair <- paste(co[2, j], "vs", co[1, j])
    
    if (pairs_vec[j] %in% names(disp_lookup)) {
      beta_disp_p[j] <- disp_lookup[pairs_vec[j]]
    } else if (alt_pair %in% names(disp_lookup)) {
      beta_disp_p[j] <- disp_lookup[alt_pair]
    } else {
      beta_disp_p[j] <- NA
    }
  }

  # P-value adjustment and export
  p.adjusted <- p.adjust(p_value, method = p.adjust.m)
  
  sig <- rep("", length(p.adjusted))
  sig[p.adjusted <= 0.05]   <- "."
  sig[p.adjusted <= 0.01]   <- "*"
  sig[p.adjusted <= 0.001]  <- "**"
  sig[p.adjusted <= 0.0001] <- "***"
  
  pairw.res <- data.frame(
    pairs       = pairs_vec, 
    F.Model     = F_Model, 
    R2          = R2, 
    p.value     = p_value, 
    p.adjusted  = p.adjusted, 
    sig         = sig,
    beta_disp_p = beta_disp_p,
    stringsAsFactors = FALSE
  )
  
  cat("Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1\n")
  
  return(pairw.res)
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
    if(file.exists(file.path(outputDir, check, paste0(group, "_Nets_OK.txt")))){
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
