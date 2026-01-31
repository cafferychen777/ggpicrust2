#' Validate group sizes for statistical analysis
#'
#' This function checks group sizes and balance to ensure statistical reliability.
#' Follows Linus principle: fail fast with clear reasons, don't hide problems.
#'
#' @param group_vector A factor vector with group assignments
#' @param group_name A character string with the group variable name for error messages
#' @keywords internal
validate_group_sizes <- function(group_vector, group_name) {
  group_counts <- table(group_vector)
  group_names <- names(group_counts)
  n_groups <- length(group_counts)
  
  # Absolute minimum: need at least 2 samples per group for any comparison
  min_samples_per_group <- 2
  insufficient_groups <- group_counts < min_samples_per_group
  
  if (any(insufficient_groups)) {
    insufficient_names <- paste(group_names[insufficient_groups], ":", group_counts[insufficient_groups], collapse = "; ")
    stop(sprintf("Groups with <2 samples detected in '%s' (%s). Statistical comparison impossible. Need at least 2 samples per group.", 
                 group_name, insufficient_names))
  }
  
  # Statistical power warning: less than 3 samples per group
  small_groups <- group_counts < 3
  if (any(small_groups)) {
    small_names <- paste(group_names[small_groups], ":", group_counts[small_groups], collapse = "; ")
    warning(sprintf("Small group sizes in '%s' (%s). Statistical power severely limited. Recommend n>=3 per group for reliable results.",
                    group_name, small_names), call. = FALSE)
  }
  
  # Group imbalance warning (only for two-group comparisons)
  if (n_groups == 2) {
    imbalance_ratio <- max(group_counts) / min(group_counts)
    if (imbalance_ratio > 3) {
      warning(sprintf("Severe group imbalance in '%s' (ratio %.1f:1). Results may be biased. Consider balancing groups or using robust methods.",
                      group_name, imbalance_ratio), call. = FALSE)
    }
  }
  
  # Multi-group warning
  if (n_groups > 2) {
    warning(sprintf("Multiple groups detected in '%s' (%d groups). GSEA will use pairwise comparisons. Consider running separate two-group analyses.",
                    group_name, n_groups), call. = FALSE)
  }
}

#' Gene Set Enrichment Analysis for PICRUSt2 output
#'
#' This function performs Gene Set Enrichment Analysis (GSEA) on PICRUSt2 predicted functional data
#' to identify enriched pathways between different conditions.
#'
#' @param abundance A data frame containing gene/enzyme abundance data, with features as rows and samples as columns.
#'   For KEGG analysis: features should be KO IDs (e.g., K00001).
#'   For MetaCyc analysis: features should be EC numbers (e.g., EC:1.1.1.1 or 1.1.1.1), NOT pathway IDs.
#'   For GO analysis: features should be KO IDs that will be mapped to GO terms.
#'   NOTE: This function requires gene-level data, not pathway-level abundances.
#'   For pathway abundance analysis, use \code{\link{pathway_daa}} instead
#' @param metadata A data frame containing sample metadata
#' @param group A character string specifying the column name in metadata that contains the grouping variable
#' @param pathway_type A character string specifying the pathway type: "KEGG", "MetaCyc", or "GO"
#' @param method A character string specifying the GSEA method:
#'   \itemize{
#'     \item \code{"camera"}: Competitive gene set test using limma's camera function (recommended).
#'           Accounts for inter-gene correlations and provides more reliable p-values.
#'     \item \code{"fry"}: Fast approximation to rotation gene set testing using limma's fry function.
#'           Self-contained test that is computationally efficient.
#'     \item \code{"fgsea"}: Fast preranked GSEA implementation. Note: preranked methods may produce
#'           unreliable p-values due to not accounting for inter-gene correlations (Wu et al., 2012).
#'     \item \code{"GSEA"} or \code{"clusterProfiler"}: clusterProfiler's GSEA implementation.
#'   }
#' @param covariates A character vector specifying column names in metadata to use as covariates
#'   for adjustment. Only used when method is "camera" or "fry". Default is NULL (no covariates).
#'   Example: \code{covariates = c("age", "sex", "BMI")}
#' @param contrast For multi-group comparisons with "camera" or "fry" methods, specify the contrast
#'   to test. Can be a character string naming a group level, or a numeric vector of contrast weights.
#'   Default is NULL (automatic: compares second group to first).
#' @param inter.gene.cor Numeric value specifying the inter-gene correlation for camera method.
#'   Default is 0.01. Use NA to estimate correlation from data for each gene set.
#' @param rank_method A character string specifying the ranking statistic for preranked methods
#'   (fgsea, GSEA, clusterProfiler): "signal2noise", "t_test", "log2_ratio", or "diff_abundance"
#' @param nperm An integer specifying the number of permutations (for fgsea/clusterProfiler methods)
#' @param min_size An integer specifying the minimum gene set size
#' @param max_size An integer specifying the maximum gene set size
#' @param p.adjust A character string specifying the p-value adjustment method
#' @param seed An integer specifying the random seed for reproducibility
#' @param go_category A character string specifying GO category: "all" (default),
#'   "BP" (Biological Process), "MF" (Molecular Function), or "CC" (Cellular Component).
#'   Note: coverage varies by category depending on the reference data.
#' @param organism A character string specifying the organism for KEGG analysis (default: "ko" for KEGG Orthology)
#'
#' @return A data frame containing GSEA results with columns:
#'   \itemize{
#'     \item \code{pathway_id}: Pathway identifier
#'     \item \code{pathway_name}: Pathway name/description
#'     \item \code{size}: Number of genes in the pathway
#'     \item \code{direction}: Direction of enrichment ("Up" or "Down", for camera/fry)
#'     \item \code{pvalue}: Raw p-value
#'     \item \code{p.adjust}: Adjusted p-value (FDR)
#'     \item \code{method}: The method used for analysis
#'   }
#'   For fgsea/clusterProfiler methods, additional columns include ES, NES, and leading_edge.
#'
#' @details
#' \strong{Method Selection:}
#'
#' The \code{camera} method (default) is recommended for most analyses because:
#' \itemize{
#'   \item It accounts for inter-gene correlations, providing more accurate p-values
#'   \item It supports covariate adjustment through the design matrix
#'   \item It performs a competitive test (genes in set vs. genes not in set)
#' }
#'
#' The \code{fry} method is a fast alternative that:
#' \itemize{
#'   \item Performs a self-contained test (are genes in the set differentially expressed?)
#'   \item Is computationally very efficient for large numbers of gene sets
#'   \item Also supports covariate adjustment
#' }
#'
#' The preranked methods (\code{fgsea}, \code{GSEA}) are included for compatibility but
#' users should be aware that Wu et al. (2012) demonstrated these can produce "spectacularly
#' wrong p-values" even with low inter-gene correlations.
#'
#' \strong{Covariate Adjustment:}
#'
#' When using \code{method = "camera"} or \code{method = "fry"}, you can adjust for
#' confounding variables by specifying them in the \code{covariates} parameter.
#' This is particularly important in microbiome studies where factors like age, sex,
#' BMI, and batch effects can influence results.
#'
#' @references
#' Wu, D., & Smyth, G. K. (2012). Camera: a competitive gene set test accounting for
#' inter-gene correlation. Nucleic Acids Research, 40(17), e133.
#'
#' Wu, D., Lim, E., Vaillant, F., Asselin-Labat, M. L., Visvader, J. E., & Smyth, G. K. (2010).
#' ROAST: rotation gene set tests for complex microarray experiments. Bioinformatics, 26(17), 2176-2182.
#'
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
#' # Method 1: Using camera (recommended) - accounts for inter-gene correlations
#' gsea_results <- pathway_gsea(
#'   abundance = abundance_data,
#'   metadata = metadata,
#'   group = "Environment",
#'   pathway_type = "KEGG",
#'   method = "camera"
#' )
#'
#' # Method 2: Using camera with covariate adjustment
#' gsea_results_adj <- pathway_gsea(
#'   abundance = abundance_data,
#'   metadata = metadata,
#'   group = "Disease",
#'   covariates = c("age", "sex"),
#'   pathway_type = "KEGG",
#'   method = "camera"
#' )
#'
#' # Method 3: Using fry for fast self-contained testing
#' gsea_results_fry <- pathway_gsea(
#'   abundance = abundance_data,
#'   metadata = metadata,
#'   group = "Environment",
#'   pathway_type = "KEGG",
#'   method = "fry"
#' )
#'
#' # Method 4: Using fgsea (preranked, less reliable p-values)
#' gsea_results_fgsea <- pathway_gsea(
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
                        method = "camera",
                        covariates = NULL,
                        contrast = NULL,
                        inter.gene.cor = 0.01,
                        rank_method = "signal2noise",
                        nperm = 1000,
                        min_size = 5,
                        max_size = 500,
                        p.adjust = "BH",
                        seed = 42,
                        go_category = "all",
                        organism = "ko") {
  
  # Input validation using unified functions
  validate_abundance(abundance, min_samples = 4)
  validate_metadata(metadata)
  validate_group(metadata, group, min_groups = 2)
  
  if (!pathway_type %in% c("KEGG", "MetaCyc", "GO")) {
    stop("pathway_type must be one of 'KEGG', 'MetaCyc', or 'GO'")
  }

  valid_methods <- c("camera", "fry", "fgsea", "GSEA", "clusterProfiler")
  if (!method %in% valid_methods) {
    stop("method must be one of: ", paste(valid_methods, collapse = ", "))
  }

  # Validate rank_method only for preranked methods
  if (method %in% c("fgsea", "GSEA", "clusterProfiler")) {
    if (!rank_method %in% c("signal2noise", "t_test", "log2_ratio", "diff_abundance")) {
      stop("rank_method must be one of 'signal2noise', 't_test', 'log2_ratio', or 'diff_abundance'")
    }
  }

  # Validate covariates
 if (!is.null(covariates)) {
    if (method %in% c("fgsea", "GSEA", "clusterProfiler")) {
      warning("Covariates are only supported for 'camera' and 'fry' methods. ",
              "Covariates will be ignored for method '", method, "'.",
              call. = FALSE)
    } else {
      # Check if all covariates exist in metadata
      missing_covs <- setdiff(covariates, colnames(metadata))
      if (length(missing_covs) > 0) {
        stop("Covariates not found in metadata: ", paste(missing_covs, collapse = ", "))
      }
    }
  }

  # Validate inter.gene.cor for camera
  if (method == "camera" && !is.null(inter.gene.cor) && !is.na(inter.gene.cor)) {
    if (!is.numeric(inter.gene.cor) || inter.gene.cor < 0 || inter.gene.cor > 1) {
      stop("inter.gene.cor must be a numeric value between 0 and 1, or NA")
    }
  }

  # Validate GO category parameter when pathway_type is "GO"
  if (pathway_type == "GO") {
    valid_go_categories <- c("BP", "MF", "CC", "all")
    if (!go_category %in% valid_go_categories) {
      stop(sprintf("Invalid go_category '%s'. Must be one of: %s",
                   go_category, paste(valid_go_categories, collapse = ", ")))
    }
  }

  # Check if required package is installed
  method_packages <- list(
    "camera" = "limma",
    "fry" = "limma",
    "fgsea" = "fgsea",
    "GSEA" = "clusterProfiler",
    "clusterProfiler" = "clusterProfiler"
  )
  require_package(method_packages[[method]], purpose = paste("GSEA method", method))

  # Warning for preranked methods about p-value reliability
  if (method %in% c("fgsea", "GSEA", "clusterProfiler")) {
    message("Note: Preranked GSEA methods (fgsea, clusterProfiler) do not account for ",
            "inter-gene correlations, which may lead to unreliable p-values ",
            "(Wu et al., 2012). Consider using method='camera' or method='fry' for ",
            "more reliable statistical inference.")
  }
  
  # Set seed for reproducibility
  set.seed(seed)
  
  # Prepare data
  # Handle #NAME column commonly found in PICRUSt2 output
  if (ncol(abundance) > 0 && colnames(abundance)[1] == "#NAME") {
    # Convert tibble to data.frame if necessary and set proper rownames
    abundance <- as.data.frame(abundance)
    rownames(abundance) <- abundance[, 1]
    abundance <- abundance[, -1]
  }
  
  # Ensure abundance is a matrix with samples as columns
  abundance_mat <- as.matrix(abundance)
  
  # Input validation for MetaCyc pathway data
  if (pathway_type == "MetaCyc") {
    # Check if features look like MetaCyc pathway IDs instead of EC numbers
    feature_names <- rownames(abundance_mat)
    if (length(feature_names) > 0) {
      # Check first 10 features for pathway patterns
      check_features <- feature_names[1:min(10, length(feature_names))]
      if (any(grepl("-PWY$|-PWY[0-9]+$", check_features))) {
        warning(
          "Input appears to contain MetaCyc pathway IDs (e.g., ending in '-PWY'). ",
          "pathway_gsea() requires gene/enzyme level data (EC numbers). ",
          "For pathway-level abundance analysis, use pathway_daa() instead. ",
          "To run GSEA, use EC abundance data from PICRUSt2 output ",
          "(EC_metagenome_out/pred_metagenome_unstrat.tsv). ",
          "See ?pathway_gsea for details.", 
          call. = FALSE, immediate. = TRUE
        )
      }
    }
  }

  # Align samples between abundance and metadata using unified function
  aligned <- align_samples(abundance_mat, metadata)
  abundance_mat <- as.matrix(aligned$abundance)
  metadata <- aligned$metadata

  # Enforce minimum sample requirement
  if (aligned$n_samples < 4) {
    stop("Insufficient samples (", aligned$n_samples,
         "). Need at least 4 samples for statistical analysis.")
  }

  # Extract group information
  Group <- factor(metadata[[group]])
  names(Group) <- colnames(abundance_mat)

  # Validate group balance for statistical reliability
  group_counts <- table(Group)
  if (any(group_counts < 2)) {
    stop("Each group must have at least 2 samples. Current group sizes: ",
         paste(names(group_counts), "=", group_counts, collapse = ", "))
  }
  
  # Prepare gene sets
  gene_sets <- prepare_gene_sets(pathway_type, organism = organism, go_category = go_category)

  # Run analysis based on selected method
  if (method %in% c("camera", "fry")) {
    # =========================================================================
    # Limma-based methods (camera/fry) with covariate support
    # =========================================================================
    results <- run_limma_gsea(
      abundance_mat = abundance_mat,
      metadata = metadata,
      group = group,
      covariates = covariates,
      contrast = contrast,
      gene_sets = gene_sets,
      method = method,
      inter.gene.cor = inter.gene.cor,
      min_size = min_size,
      max_size = max_size,
      p.adjust.method = p.adjust
    )

  } else if (method == "fgsea") {
    # =========================================================================
    # Preranked methods (fgsea)
    # =========================================================================

    # Calculate ranking metric for preranked methods
    ranked_list <- calculate_rank_metric(abundance_mat, metadata, group, method = rank_method)
    results <- run_fgsea(ranked_list, gene_sets, nperm, min_size, max_size)

  } else if (method == "GSEA" || method == "clusterProfiler") {
    # =========================================================================
    # clusterProfiler GSEA
    # =========================================================================

    # Calculate ranking metric for preranked methods
    ranked_list <- calculate_rank_metric(abundance_mat, metadata, group, method = rank_method)

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
      results <- create_empty_gsea_result(method, full = TRUE)
    }
  }

  # Add method information to populated results
  # (empty results already have method from factory)
  if (nrow(results) > 0) {
    results$method <- method
  }
  
  return(results)
}

#' Prepare gene sets for GSEA
#'
#' @param pathway_type A character string specifying the pathway type: "KEGG", "MetaCyc", or "GO"
#' @param organism A character string specifying the organism (only relevant for KEGG and GO)
#' @param go_category A character string specifying the GO category: "BP" (Biological Process), "MF" (Molecular Function), "CC" (Cellular Component), or "all"
#'
#' @return A list of pathway gene sets
#' @export
prepare_gene_sets <- function(pathway_type = "KEGG", organism = "ko", go_category = "all") {

  # Validate pathway_type
  valid_types <- c("KEGG", "MetaCyc", "GO")
  if (!pathway_type %in% valid_types) {
    stop("pathway_type must be one of: ", paste(valid_types, collapse = ", "))
  }

  if (pathway_type == "KEGG") {
    # Load KEGG reference using unified loader
    ko_to_kegg_reference <- load_reference_data("ko_to_kegg")

    # Create gene sets: list where each element is a pathway containing KO IDs
    gene_sets <- split(ko_to_kegg_reference$ko_id, ko_to_kegg_reference$pathway_id)

  } else if (pathway_type == "MetaCyc") {
    # Load MetaCyc to EC mapping using unified loader
    metacyc_to_ec_reference <- load_reference_data("metacyc_to_ec")

    # Create gene sets list from MetaCyc pathways
    gene_sets <- list()
    for (i in seq_len(nrow(metacyc_to_ec_reference))) {
      pathway_id <- metacyc_to_ec_reference[i, "pathway"]
      ec_string <- as.character(metacyc_to_ec_reference[i, "ec_numbers"])

      # Skip pathways with no EC mappings
      if (is.na(ec_string) || ec_string == "" || ec_string == "NA") {
        next
      }

      # Split EC numbers by semicolon and clean
      ec_numbers <- strsplit(ec_string, ";")[[1]]
      ec_numbers <- trimws(ec_numbers)
      ec_numbers <- ec_numbers[ec_numbers != ""]

      if (length(ec_numbers) > 0) {
        # Add EC: prefix if not present for consistency
        ec_numbers <- ifelse(grepl("^EC:", ec_numbers), ec_numbers, paste0("EC:", ec_numbers))
        gene_sets[[pathway_id]] <- ec_numbers
      }
    }

  } else if (pathway_type == "GO") {
    # Load GO reference using unified loader
    go_reference <- load_reference_data("ko_to_go")
    
    # Validate and filter by GO category if specified
    valid_go_categories <- c("BP", "MF", "CC", "all")
    if (!is.null(go_category) && !go_category %in% valid_go_categories) {
      stop(sprintf("Invalid go_category '%s'. Must be one of: %s", 
                   go_category, paste(valid_go_categories, collapse = ", ")))
    }
    
    if (!is.null(go_category) && go_category != "all" && go_category %in% c("BP", "MF", "CC")) {
      if ("category" %in% colnames(go_reference)) {
        go_reference <- go_reference[go_reference$category == go_category, ]
      }
    }

    if (nrow(go_reference) == 0) {
      full_ref <- load_reference_data("ko_to_go")
      available <- if ("category" %in% colnames(full_ref)) {
        paste(unique(full_ref$category), collapse = ", ")
      } else {
        "unknown"
      }
      stop(sprintf(
        "No GO terms found for category '%s'. Available categories in reference data: %s. ",
        go_category, available),
        "Try go_category = \"all\" or go_category = \"MF\".",
        call. = FALSE)
    }

    # Create gene sets list for each GO term
    gene_sets <- list()

    # Use seq_len to handle empty data frame correctly (1:0 returns c(1,0), not empty)
    for (i in seq_len(nrow(go_reference))) {
      go_id <- go_reference$go_id[i]

      # Get KO members for this GO term
      if ("ko_members" %in% colnames(go_reference)) {
        ko_string <- as.character(go_reference[i, "ko_members"])
        if (!is.na(ko_string) && nchar(ko_string) > 0) {
          ko_ids <- unlist(strsplit(ko_string, ";"))
          ko_ids <- ko_ids[!is.na(ko_ids) & ko_ids != ""]

          if (length(ko_ids) > 0) {
            gene_sets[[go_id]] <- ko_ids
          }
        }
      } else {
        # Fallback: assume other columns contain KO IDs
        ko_ids <- as.character(go_reference[i, -1])
        ko_ids <- ko_ids[!is.na(ko_ids) & ko_ids != ""]

        if (length(ko_ids) > 0) {
          gene_sets[[go_id]] <- ko_ids
        }
      }
    }
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
  common_samples <- intersect(colnames(abundance), names(Group))
  abundance <- abundance[, common_samples, drop = FALSE]
  Group <- Group[common_samples]

  # Handle problematic values - replace with 0 (undetected features)
  if (any(is.infinite(abundance))) {
    warning("Abundance data contains Inf values. Replacing with 0.", call. = FALSE)
    abundance[is.infinite(abundance)] <- 0
  }
  if (any(is.na(abundance))) {
    warning("Abundance data contains NA values. Replacing with 0.", call. = FALSE)
    abundance[is.na(abundance)] <- 0
  }
  if (any(abundance < 0)) {
    warning("Abundance data contains negative values. Replacing with 0.", call. = FALSE)
    abundance[abundance < 0] <- 0
  }

  # Group size validation already done above in main function
  
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
    
    # Handle zero standard deviations using GSEA-style minimum
    # GSEA uses: sigma_min = 0.2 * |mu|, where mu=0 is adjusted to mu=1
    # This ensures the floor is proportional to the data scale
    adjusted_mean1 <- ifelse(mean1 == 0, 1, abs(mean1))
    adjusted_mean2 <- ifelse(mean2 == 0, 1, abs(mean2))
    sd1 <- pmax(sd1, 0.2 * adjusted_mean1)
    sd2 <- pmax(sd2, 0.2 * adjusted_mean2)
    
    metric <- (mean1 - mean2) / (sd1 + sd2)
    # Ensure names are preserved - critical for fgsea
    names(metric) <- rownames(abundance)
    
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

    # Calculate pseudocount using unified function
    # Note: GSEA uses mean1/mean2 direction (group1/group2)
    pseudocount <- calculate_pseudocount(c(mean1, mean2))
    mean1[mean1 == 0] <- pseudocount
    mean2[mean2 == 0] <- pseudocount

    metric <- log2(mean1 / mean2)
    # Ensure names are preserved - critical for fgsea
    names(metric) <- rownames(abundance)
    
  } else if (method == "diff_abundance") {
    # Simple difference in abundance
    mean1 <- rowMeans(abundance[, group1_samples, drop = FALSE])
    mean2 <- rowMeans(abundance[, group2_samples, drop = FALSE])
    metric <- mean1 - mean2
    # Ensure names are preserved - critical for fgsea
    names(metric) <- rownames(abundance)
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
  # Run fgsea
  fgsea_result <- fgsea::fgsea(
    pathways = gene_sets,
    stats = ranked_list,
    minSize = min_size,
    maxSize = max_size,
    nperm = nperm
  )
  
  # Convert to data frame and handle empty results
  results <- as.data.frame(fgsea_result)
  
  # Check if we have any results
  if (nrow(results) > 0) {
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
  } else {
    results <- create_empty_gsea_result(full = TRUE)
  }

  return(results)
}

#' Run limma-based gene set analysis (camera/fry)
#'
#' This internal function implements limma's camera and fry methods for gene set
#' enrichment analysis with support for covariates.
#'
#' @param abundance_mat A matrix of abundance data with features as rows and samples as columns
#' @param metadata A data frame containing sample metadata
#' @param group A character string specifying the grouping variable column name
#' @param covariates A character vector of covariate column names (optional)
#' @param contrast Contrast specification for multi-group comparisons (optional)
#' @param gene_sets A named list of gene sets (pathway -> gene IDs)
#' @param method Either "camera" or "fry"
#' @param inter.gene.cor Inter-gene correlation for camera (default 0.01)
#' @param min_size Minimum gene set size
#' @param max_size Maximum gene set size
#' @param p.adjust.method P-value adjustment method
#'
#' @return A data frame containing gene set analysis results
#' @keywords internal
run_limma_gsea <- function(abundance_mat,
                           metadata,
                           group,
                           covariates = NULL,
                           contrast = NULL,
                           gene_sets,
                           method = "camera",
                           inter.gene.cor = 0.01,
                           min_size = 5,
                           max_size = 500,
                           p.adjust.method = "BH") {
  # Align samples using unified function
  aligned <- align_samples(abundance_mat, metadata, verbose = FALSE)
  abundance_mat <- as.matrix(aligned$abundance)
  metadata <- aligned$metadata

  if (aligned$n_samples < 4) {
    stop("Insufficient matching samples between abundance data and metadata")
  }

  # Build design matrix
  design <- build_design_matrix(metadata, group, covariates)

  # Determine which coefficient to test
  # For a design with intercept + group, the group coefficient is typically column 2
  # For a design without intercept (0 + group), we need to construct a contrast
  group_levels <- levels(factor(metadata[[group]]))

  if (is.null(contrast)) {
    # Default: test the group effect (second column for intercept model)
    # or construct a contrast for no-intercept model
    if (ncol(design) >= 2 && grepl(group, colnames(design)[2])) {
      # Intercept model: test second coefficient
      contrast_coef <- 2
    } else {
      # Create contrast: second group vs first group
      contrast_coef <- 2
    }
  } else if (is.numeric(contrast)) {
    contrast_coef <- contrast
  } else if (is.character(contrast)) {
    # Find the column matching the contrast
    contrast_col <- grep(contrast, colnames(design))
    if (length(contrast_col) == 0) {
      stop("Contrast '", contrast, "' not found in design matrix columns: ",
           paste(colnames(design), collapse = ", "))
    }
    contrast_coef <- contrast_col[1]
  } else {
    contrast_coef <- 2
  }

  # Convert gene sets to index format for limma
  # gene_sets is a named list: pathway_id -> vector of gene IDs
  # We need to convert to row indices in abundance_mat
  feature_names <- rownames(abundance_mat)
  gene_set_indices <- lapply(gene_sets, function(genes) {
    which(feature_names %in% genes)
  })

  # Filter gene sets by size
  gene_set_sizes <- sapply(gene_set_indices, length)
  valid_sets <- gene_set_sizes >= min_size & gene_set_sizes <= max_size

  if (sum(valid_sets) == 0) {
    warning("No gene sets passed the size filter (min_size=", min_size,
            ", max_size=", max_size, "). Returning empty results.")
    return(create_empty_gsea_result(method))
  }

  gene_set_indices <- gene_set_indices[valid_sets]
  message(sprintf("Testing %d gene sets (filtered from %d by size constraints)",
                  length(gene_set_indices), length(gene_sets)))

  # Handle problematic values before voom transformation
  # voom/limma do not accept NA, Inf, or negative values in count data

  # Handle Inf/-Inf values first
  if (any(is.infinite(abundance_mat))) {
    inf_count <- sum(is.infinite(abundance_mat))
    warning(sprintf("Abundance data contains %d Inf/-Inf values. ",
                    inf_count),
            "Replacing with 0.",
            call. = FALSE)
    abundance_mat[is.infinite(abundance_mat)] <- 0
  }

  # Handle NA values
  if (any(is.na(abundance_mat))) {
    na_count <- sum(is.na(abundance_mat))
    na_features <- sum(rowSums(is.na(abundance_mat)) > 0)
    warning(sprintf("Abundance data contains %d NA values in %d features. ",
                    na_count, na_features),
            "Replacing NA with 0 (assuming undetected features).",
            call. = FALSE)
    abundance_mat[is.na(abundance_mat)] <- 0
  }

  # Handle negative values (voom expects count-like data)
  if (any(abundance_mat < 0)) {
    neg_count <- sum(abundance_mat < 0)
    warning(sprintf("Abundance data contains %d negative values. ",
                    neg_count),
            "Replacing with 0 (voom expects non-negative count data).",
            call. = FALSE)
    abundance_mat[abundance_mat < 0] <- 0
  }

  # Apply voom transformation for count data
  # First, add a small pseudocount to avoid log(0)
  abundance_mat <- abundance_mat + 0.5

  # Use limma-voom for count data transformation
  # This estimates the mean-variance relationship and computes precision weights
  tryCatch({
    v <- limma::voom(abundance_mat, design, plot = FALSE)
  }, error = function(e) {
    # Fallback: use log2 transformation if voom fails
    warning("voom transformation failed, using log2 transformation instead: ", e$message)
    v <<- list(
      E = log2(abundance_mat),
      weights = matrix(1, nrow = nrow(abundance_mat), ncol = ncol(abundance_mat))
    )
    class(v) <<- "EList"
  })

  # Run the selected method
  if (method == "camera") {
    # Camera: competitive gene set test
    # Accounts for inter-gene correlations
    results <- limma::camera(
      y = v,
      index = gene_set_indices,
      design = design,
      contrast = contrast_coef,
      inter.gene.cor = inter.gene.cor
    )

  } else if (method == "fry") {
    # Fry: fast rotation gene set test (self-contained)
    results <- limma::fry(
      y = v,
      index = gene_set_indices,
      design = design,
      contrast = contrast_coef
    )
  }

  # Format results
  if (nrow(results) > 0) {
    # Add pathway IDs from rownames
    results$pathway_id <- rownames(results)

    # Standardize column names
    formatted_results <- data.frame(
      pathway_id = results$pathway_id,
      pathway_name = results$pathway_id,  # Will be annotated later if mapping available
      size = results$NGenes,
      direction = results$Direction,
      pvalue = results$PValue,
      p.adjust = stats::p.adjust(results$PValue, method = p.adjust.method),
      method = method,
      stringsAsFactors = FALSE
    )

    # Add FDR from results if available (camera/fry compute their own FDR)
    if ("FDR" %in% colnames(results)) {
      formatted_results$FDR_original <- results$FDR
    }

    # Add synthetic NES for visualization compatibility
    # Camera/fry don't produce NES, so we create a synthetic measure based on
    # direction and -log10(pvalue) for visualization purposes
    formatted_results$NES <- ifelse(
      formatted_results$direction == "Up",
      -log10(pmax(formatted_results$pvalue, 1e-300)),  # Positive for Up
      log10(pmax(formatted_results$pvalue, 1e-300))     # Negative for Down
    )

    # Add a leading_edge placeholder for visualization compatibility
    formatted_results$leading_edge <- ""

    # Sort by p-value
    formatted_results <- formatted_results[order(formatted_results$pvalue), ]
    rownames(formatted_results) <- NULL

  } else {
    formatted_results <- create_empty_gsea_result(method)
  }

  return(formatted_results)
}

#' Build design matrix for limma analysis
#'
#' Creates a design matrix incorporating the group variable and optional covariates.
#'
#' @param metadata A data frame containing sample metadata
#' @param group A character string specifying the grouping variable column name
#' @param covariates A character vector of covariate column names (optional)
#'
#' @return A design matrix suitable for limma
#' @keywords internal
build_design_matrix <- function(metadata, group, covariates = NULL) {

  # Ensure group is a factor
  metadata[[group]] <- factor(metadata[[group]])

  # Build formula
  if (is.null(covariates) || length(covariates) == 0) {
    # Simple formula: ~ group
    formula_str <- paste("~", group)
  } else {
    # Formula with covariates: ~ group + cov1 + cov2 + ...
    # Ensure covariates are properly formatted
    for (cov in covariates) {
      if (is.character(metadata[[cov]])) {
        metadata[[cov]] <- factor(metadata[[cov]])
      }
    }
    formula_str <- paste("~", group, "+", paste(covariates, collapse = " + "))
  }

  # Create design matrix
  design <- stats::model.matrix(stats::as.formula(formula_str), data = metadata)

  # Clean up column names for readability
  colnames(design) <- gsub("\\(Intercept\\)", "Intercept", colnames(design))

  return(design)
}
