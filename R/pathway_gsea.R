#' Gene Set Enrichment Analysis for PICRUSt2 output
#'
#' This function performs Gene Set Enrichment Analysis (GSEA) on PICRUSt2 predicted functional data
#' to identify enriched pathways between different conditions.
#'
#' @param abundance A data frame containing KO/EC/MetaCyc abundance data, with features as rows and samples as columns
#' @param metadata A data frame containing sample metadata
#' @param group A character string specifying the column name in metadata that contains the grouping variable
#' @param pathway_type A character string specifying the pathway type: "KEGG", "MetaCyc", or "GO"
#' @param method A character string specifying the GSEA method: "fgsea", "GSEA", or "clusterProfiler"
#' @param rank_method A character string specifying the ranking statistic: "signal2noise", "t_test", "log2_ratio", or "diff_abundance"
#' @param nperm An integer specifying the number of permutations
#' @param min_size An integer specifying the minimum gene set size
#' @param max_size An integer specifying the maximum gene set size
#' @param p.adjust A character string specifying the p-value adjustment method
#' @param seed An integer specifying the random seed for reproducibility
#'
#' @return A data frame containing GSEA results
#' @export
#'
#' @examples
#' \dontrun{
#' # Load example data
#' data(ko_abundance)
#' data(metadata)
#' 
#' # Prepare abundance data
#' abundance_data <- as.data.frame(ko_abundance)
#' rownames(abundance_data) <- abundance_data[, "#NAME"]
#' abundance_data <- abundance_data[, -1]
#' 
#' # Run GSEA analysis
#' gsea_results <- pathway_gsea(
#'   abundance = abundance_data,
#'   metadata = metadata,
#'   group = "Environment",
#'   pathway_type = "KEGG",
#'   method = "fgsea"
#' )
#' 
#' # Visualize results
#' visualize_gsea(gsea_results, plot_type = "enrichment_plot", n_pathways = 10)
#' }
pathway_gsea <- function(abundance, 
                        metadata, 
                        group, 
                        pathway_type = "KEGG", 
                        method = "fgsea", 
                        rank_method = "signal2noise", 
                        nperm = 1000, 
                        min_size = 10, 
                        max_size = 500, 
                        p.adjust = "BH",
                        seed = 42) {
  
  # Input validation
  if (!is.data.frame(abundance) && !is.matrix(abundance)) {
    stop("'abundance' must be a data frame or matrix")
  }
  
  if (!is.data.frame(metadata)) {
    stop("'metadata' must be a data frame")
  }
  
  if (!group %in% colnames(metadata)) {
    stop(paste("Group variable", group, "not found in metadata"))
  }
  
  if (!pathway_type %in% c("KEGG", "MetaCyc", "GO")) {
    stop("pathway_type must be one of 'KEGG', 'MetaCyc', or 'GO'")
  }
  
  if (!method %in% c("fgsea", "GSEA", "clusterProfiler")) {
    stop("method must be one of 'fgsea', 'GSEA', or 'clusterProfiler'")
  }
  
  if (!rank_method %in% c("signal2noise", "t_test", "log2_ratio", "diff_abundance")) {
    stop("rank_method must be one of 'signal2noise', 't_test', 'log2_ratio', or 'diff_abundance'")
  }
  
  # Check if required packages are installed
  required_packages <- list(
    "fgsea" = "fgsea",
    "GSEA" = "clusterProfiler",
    "clusterProfiler" = "clusterProfiler"
  )
  
  pkg_name <- required_packages[[method]]
  if (!requireNamespace(pkg_name, quietly = TRUE)) {
    stop(paste("Package", pkg_name, "is required for method", method, 
               ". Please install it using BiocManager::install('", pkg_name, "').", sep = ""))
  }
  
  # Set seed for reproducibility
  set.seed(seed)
  
  # Prepare data
  # Ensure abundance is a matrix with samples as columns
  abundance_mat <- as.matrix(abundance)
  
  # Extract group information
  Group <- factor(metadata[[group]])
  names(Group) <- rownames(metadata)
  
  # Check if sample names match between abundance and metadata
  if (!all(colnames(abundance_mat) %in% names(Group))) {
    stop("Sample names in abundance data do not match sample names in metadata")
  }
  
  # Subset abundance to include only samples in metadata
  abundance_mat <- abundance_mat[, names(Group)]
  
  # Prepare gene sets
  gene_sets <- prepare_gene_sets(pathway_type)
  
  # Calculate ranking metric
  ranked_list <- calculate_rank_metric(abundance_mat, metadata, group, rank_method)
  
  # Run GSEA based on selected method
  if (method == "fgsea") {
    if (!requireNamespace("fgsea", quietly = TRUE)) {
      stop("Package 'fgsea' is required. Please install it using BiocManager::install('fgsea').")
    }
    
    results <- run_fgsea(ranked_list, gene_sets, nperm, min_size, max_size)
    
  } else if (method == "GSEA" || method == "clusterProfiler") {
    if (!requireNamespace("clusterProfiler", quietly = TRUE)) {
      stop("Package 'clusterProfiler' is required. Please install it using BiocManager::install('clusterProfiler').")
    }
    
    # Convert ranked list to format required by clusterProfiler
    gene_list <- sort(ranked_list, decreasing = TRUE)
    
    # Run GSEA using clusterProfiler
    gsea_result <- clusterProfiler::GSEA(
      geneList = gene_list,
      TERM2GENE = data.frame(
        term = rep(names(gene_sets), sapply(gene_sets, length)),
        gene = unlist(gene_sets)
      ),
      minGSSize = min_size,
      maxGSSize = max_size,
      pvalueCutoff = 1,
      pAdjustMethod = p.adjust,
      nPermSimple = nperm,
      seed = seed
    )
    
    # Convert results to data frame
    if (length(gsea_result) > 0) {
      results <- as.data.frame(gsea_result)
      
      # Rename columns to match fgsea output
      results <- data.frame(
        pathway_id = results$ID,
        pathway_name = results$Description,
        size = results$setSize,
        ES = results$enrichmentScore,
        NES = results$NES,
        pvalue = results$pvalue,
        p.adjust = results$p.adjust,
        leading_edge = results$core_enrichment,
        stringsAsFactors = FALSE
      )
    } else {
      results <- data.frame(
        pathway_id = character(),
        pathway_name = character(),
        size = integer(),
        ES = numeric(),
        NES = numeric(),
        pvalue = numeric(),
        p.adjust = numeric(),
        leading_edge = character(),
        stringsAsFactors = FALSE
      )
    }
  }
  
  # Add method information
  results$method <- method
  
  return(results)
}

#' Prepare gene sets for GSEA
#'
#' @param pathway_type A character string specifying the pathway type: "KEGG", "MetaCyc", or "GO"
#' @param organism A character string specifying the organism (only relevant for KEGG and GO)
#'
#' @return A list of pathway gene sets
#' @keywords internal
prepare_gene_sets <- function(pathway_type = "KEGG", organism = "ko") {
  
  if (pathway_type == "KEGG") {
    # Load KEGG pathway to KO mapping
    if (!exists("ko_to_kegg_reference")) {
      data("ko_to_kegg_reference", package = "ggpicrust2", envir = environment())
    }
    
    # Convert to list format required for GSEA
    ko_to_kegg_reference <- as.data.frame(ko_to_kegg_reference)
    
    # Create a list where each element is a pathway containing KO IDs
    gene_sets <- list()
    
    for (i in 1:nrow(ko_to_kegg_reference)) {
      pathway_id <- ko_to_kegg_reference[i, 1]
      ko_ids <- as.character(ko_to_kegg_reference[i, -1])
      ko_ids <- ko_ids[!is.na(ko_ids) & ko_ids != ""]
      
      if (length(ko_ids) > 0) {
        gene_sets[[pathway_id]] <- ko_ids
      }
    }
    
  } else if (pathway_type == "MetaCyc") {
    # For MetaCyc, we need to create a mapping from MetaCyc pathways to EC numbers
    # This would require additional reference data
    # For now, we'll return a placeholder
    gene_sets <- list()
    warning("MetaCyc pathway gene sets not yet implemented")
    
  } else if (pathway_type == "GO") {
    # For GO, we would need to map KO IDs to GO terms
    # This would require additional reference data
    # For now, we'll return a placeholder
    gene_sets <- list()
    warning("GO pathway gene sets not yet implemented")
  }
  
  return(gene_sets)
}

#' Calculate rank metric for GSEA
#'
#' @param abundance A matrix of abundance data
#' @param metadata A data frame of metadata
#' @param group A character string specifying the grouping variable
#' @param method A character string specifying the ranking method
#'
#' @return A named vector of ranking statistics
#' @keywords internal
calculate_rank_metric <- function(abundance, 
                                 metadata, 
                                 group, 
                                 method = "signal2noise") {
  
  # Extract group information
  Group <- factor(metadata[[group]])
  names(Group) <- rownames(metadata)
  
  # Ensure abundance is a matrix with samples as columns
  abundance <- as.matrix(abundance)
  
  # Subset abundance to include only samples in metadata
  abundance <- abundance[, names(Group)]
  
  # Get unique group levels
  levels <- levels(Group)
  
  if (length(levels) != 2) {
    stop("GSEA currently only supports two-group comparisons")
  }
  
  # Split samples by group
  group1_samples <- names(Group)[Group == levels[1]]
  group2_samples <- names(Group)[Group == levels[2]]
  
  # Calculate ranking statistic based on method
  if (method == "signal2noise") {
    # Signal-to-noise ratio: (mean1 - mean2) / (sd1 + sd2)
    mean1 <- rowMeans(abundance[, group1_samples, drop = FALSE])
    mean2 <- rowMeans(abundance[, group2_samples, drop = FALSE])
    sd1 <- apply(abundance[, group1_samples, drop = FALSE], 1, stats::sd)
    sd2 <- apply(abundance[, group2_samples, drop = FALSE], 1, stats::sd)
    
    # Handle zero standard deviations
    sd1[sd1 == 0] <- 0.00001
    sd2[sd2 == 0] <- 0.00001
    
    metric <- (mean1 - mean2) / (sd1 + sd2)
    
  } else if (method == "t_test") {
    # t-test statistic
    metric <- numeric(nrow(abundance))
    names(metric) <- rownames(abundance)
    
    for (i in 1:nrow(abundance)) {
      t_test <- stats::t.test(abundance[i, group1_samples], abundance[i, group2_samples])
      metric[i] <- t_test$statistic
    }
    
  } else if (method == "log2_ratio") {
    # Log2 fold change
    mean1 <- rowMeans(abundance[, group1_samples, drop = FALSE])
    mean2 <- rowMeans(abundance[, group2_samples, drop = FALSE])
    
    # Add small constant to avoid division by zero
    mean1[mean1 == 0] <- 0.00001
    mean2[mean2 == 0] <- 0.00001
    
    metric <- log2(mean1 / mean2)
    
  } else if (method == "diff_abundance") {
    # Simple difference in abundance
    mean1 <- rowMeans(abundance[, group1_samples, drop = FALSE])
    mean2 <- rowMeans(abundance[, group2_samples, drop = FALSE])
    metric <- mean1 - mean2
  }
  
  return(metric)
}

#' Run fast GSEA implementation
#'
#' @param ranked_list A named vector of ranking statistics
#' @param gene_sets A list of pathway gene sets
#' @param nperm An integer specifying the number of permutations
#' @param min_size An integer specifying the minimum gene set size
#' @param max_size An integer specifying the maximum gene set size
#'
#' @return A data frame of fgsea results
#' @keywords internal
run_fgsea <- function(ranked_list, 
                     gene_sets, 
                     nperm = 1000, 
                     min_size = 10, 
                     max_size = 500) {
  
  # Check if fgsea package is available
  if (!requireNamespace("fgsea", quietly = TRUE)) {
    stop("Package 'fgsea' is required. Please install it using BiocManager::install('fgsea').")
  }
  
  # Run fgsea
  fgsea_result <- fgsea::fgsea(
    pathways = gene_sets,
    stats = ranked_list,
    minSize = min_size,
    maxSize = max_size,
    nperm = nperm
  )
  
  # Convert to data frame
  results <- as.data.frame(fgsea_result)
  
  # Rename columns for consistency
  results <- data.frame(
    pathway_id = results$pathway,
    pathway_name = results$pathway,  # Will be updated later with annotation
    size = results$size,
    ES = results$ES,
    NES = results$NES,
    pvalue = results$pval,
    p.adjust = results$padj,
    leading_edge = sapply(results$leadingEdge, function(x) paste(x, collapse = ";")),
    stringsAsFactors = FALSE
  )
  
  return(results)
}
