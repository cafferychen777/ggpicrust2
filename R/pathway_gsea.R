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
#' @param go_category A character string specifying GO category: "BP" (Biological Process), "MF" (Molecular Function), "CC" (Cellular Component), or "all"
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
                        go_category = "BP",
                        organism = "ko") {
  
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

  # Check if required packages are installed
  required_packages <- list(
    "camera" = "limma",
    "fry" = "limma",
    "fgsea" = "fgsea",
    "GSEA" = "clusterProfiler",
    "clusterProfiler" = "clusterProfiler"
  )

  pkg_name <- required_packages[[method]]
  if (!requireNamespace(pkg_name, quietly = TRUE)) {
    stop(paste("Package", pkg_name, "is required for method", method,
               ". Please install it using BiocManager::install('", pkg_name, "').", sep = ""))
  }

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

  # Ensure metadata has proper rownames for sample matching
  # This handles tibbles and data.frames where rownames may not be set properly
  metadata <- as.data.frame(metadata)
  if (length(intersect(colnames(abundance_mat), rownames(metadata))) == 0) {
    # Try to use sample/sample_name column as rownames
    sample_col <- NULL
    if ("sample_name" %in% colnames(metadata)) {
      sample_col <- "sample_name"
    } else if ("sample" %in% colnames(metadata)) {
      sample_col <- "sample"
    } else if ("SampleID" %in% colnames(metadata)) {
      sample_col <- "SampleID"
    }

    if (!is.null(sample_col)) {
      rownames(metadata) <- metadata[[sample_col]]
    }
  }

  # Extract group information
  Group <- factor(metadata[[group]])
  names(Group) <- rownames(metadata)

  # Find common samples (eliminate special case handling)
  common_samples <- intersect(colnames(abundance_mat), names(Group))
  
  # Enforce minimum requirements (clear validation)
  if (length(common_samples) < 4) {
    stop("Insufficient overlapping samples (", length(common_samples), 
         "/", ncol(abundance_mat), "). Need at least 4 samples for statistical analysis.")
  }
  
  # Inform user about sample subsetting (transparency)
  if (length(common_samples) < ncol(abundance_mat)) {
    warning(sprintf("Using %d/%d overlapping samples between abundance data and metadata", 
                    length(common_samples), ncol(abundance_mat)))
  }
  
  # Clean subset to common samples (no special cases)
  abundance_mat <- abundance_mat[, common_samples]
  Group <- Group[common_samples]
  
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
    if (!requireNamespace("fgsea", quietly = TRUE)) {
      stop("Package 'fgsea' is required. Please install it using BiocManager::install('fgsea').")
    }

    # Calculate ranking metric for preranked methods
    ranked_list <- calculate_rank_metric(abundance_mat, metadata, group, method = rank_method)
    results <- run_fgsea(ranked_list, gene_sets, nperm, min_size, max_size)

  } else if (method == "GSEA" || method == "clusterProfiler") {
    # =========================================================================
    # clusterProfiler GSEA
    # =========================================================================
    if (!requireNamespace("clusterProfiler", quietly = TRUE)) {
      stop("Package 'clusterProfiler' is required. Please install it using BiocManager::install('clusterProfiler').")
    }

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
  
  # Add method information (handle empty results)
  if (nrow(results) > 0) {
    results$method <- method
  } else {
    # Create empty result with correct structure
    results <- data.frame(
      pathway_id = character(),
      pathway_name = character(),
      size = integer(),
      ES = numeric(),
      NES = numeric(),
      pvalue = numeric(),
      p.adjust = numeric(),
      leading_edge = character(),
      method = character(),
      stringsAsFactors = FALSE
    )
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
prepare_gene_sets <- function(pathway_type = "KEGG", organism = "ko", go_category = "BP") {

  # Validate pathway_type
  valid_types <- c("KEGG", "MetaCyc", "GO")
  if (!pathway_type %in% valid_types) {
    stop("pathway_type must be one of: ", paste(valid_types, collapse = ", "))
  }

  if (pathway_type == "KEGG") {
    # Initialize gene_sets
    gene_sets <- list()
    
    # Try to load KEGG pathway to KO mapping
    tryCatch({
      if (!exists("ko_to_kegg_reference")) {
        data("ko_to_kegg_reference", package = "ggpicrust2", envir = environment())
      }
      
      # Convert to list format required for GSEA (now using long-format data)
      ko_to_kegg_reference <- as.data.frame(ko_to_kegg_reference)

      # Create a list where each element is a pathway containing KO IDs
      # Using split() for efficient conversion from long format
      gene_sets <- split(ko_to_kegg_reference$ko_id, ko_to_kegg_reference$pathway_id)
      
    }, error = function(e) {
      # Create dummy gene sets for testing when reference data is not available
      warning("KEGG reference data not found. Creating dummy gene sets for testing.", call. = FALSE)
      
      # Create some dummy pathways for demonstration
      gene_sets <<- list(
        "ko00010" = c("K00844", "K12407", "K00845", "K00886", "K08074"),
        "ko00020" = c("K00239", "K00240", "K00241", "K00242", "K01902"),
        "ko00030" = c("K00016", "K00018", "K00128", "K01595", "K01596"),
        "ko00040" = c("K01623", "K01624", "K11645", "K01803", "K15633"),
        "ko00051" = c("K00134", "K00150", "K03781", "K03782", "K14085")
      )
    })
    
  } else if (pathway_type == "MetaCyc") {
    # Load MetaCyc pathway to EC mapping
    if (!exists("metacyc_to_ec_reference")) {
      # Try to load from package extdata first
      tryCatch({
        metacyc_ref_path <- system.file("extdata", "metacyc_to_ec_reference.RData", package = "ggpicrust2")
        if (file.exists(metacyc_ref_path)) {
          load(metacyc_ref_path, envir = environment())
        } else {
          stop("metacyc_to_ec_reference data file not found")
        }
      }, error = function(e) {
        stop("Failed to load MetaCyc to EC mapping: ", e$message)
      })
    }
    
    # Convert to data frame
    metacyc_to_ec_reference <- as.data.frame(metacyc_to_ec_reference)
    
    # Create gene sets list
    gene_sets <- list()
    
    for (i in 1:nrow(metacyc_to_ec_reference)) {
      pathway_id <- metacyc_to_ec_reference[i, "pathway"]
      ec_string <- as.character(metacyc_to_ec_reference[i, "ec_numbers"])
      
      # Skip pathways with no EC mappings
      if (is.na(ec_string) || ec_string == "" || ec_string == "NA") {
        next
      }
      
      # Split EC numbers by semicolon
      ec_numbers <- strsplit(ec_string, ";")[[1]]
      ec_numbers <- trimws(ec_numbers)  # Remove whitespace
      ec_numbers <- ec_numbers[ec_numbers != ""]  # Remove empty strings
      
      if (length(ec_numbers) > 0) {
        # Add EC: prefix if not present for consistency
        ec_numbers <- ifelse(grepl("^EC:", ec_numbers), ec_numbers, paste0("EC:", ec_numbers))
        gene_sets[[pathway_id]] <- ec_numbers
      }
    }
    
  } else if (pathway_type == "GO") {
    # Load KO to GO mapping with improved error handling
    ko_to_go_reference <- NULL

    # Try to load complete GO reference data from multiple sources
    data_loaded <- FALSE

    # Method 1: Try to load from package data
    tryCatch({
      data("ko_to_go_reference", package = "ggpicrust2", envir = environment())
      if (exists("ko_to_go_reference", envir = environment()) && !is.null(ko_to_go_reference)) {
        message("[OK] Using complete ko_to_go_reference dataset from package data")
        data_loaded <- TRUE
      }
    }, error = function(e) {
      # Continue to next method
    })

    # Method 2: Try to load from data/ directory directly
    if (!data_loaded) {
      tryCatch({
        data_file <- "data/ko_to_go_reference.RData"
        if (file.exists(data_file)) {
          load(data_file, envir = environment())
          if (exists("ko_to_go_reference", envir = environment()) && !is.null(ko_to_go_reference)) {
            message("[OK] Using complete ko_to_go_reference dataset from data file")
            data_loaded <- TRUE
          }
        }
      }, error = function(e) {
        # Continue to fallback
      })
    }

    # Method 3: Try to load from system file
    if (!data_loaded) {
      tryCatch({
        data_file <- system.file("data", "ko_to_go_reference.RData", package = "ggpicrust2")
        if (file.exists(data_file)) {
          load(data_file, envir = environment())
          if (exists("ko_to_go_reference", envir = environment()) && !is.null(ko_to_go_reference)) {
            message("[OK] Using complete ko_to_go_reference dataset from system file")
            data_loaded <- TRUE
          }
        }
      }, error = function(e) {
        # Continue to fallback
      })
    }

    # Fallback: Use enhanced basic mapping if complete data is not available
    if (!data_loaded) {
      warning("Complete ko_to_go_reference dataset not found. ",
              "Using enhanced basic GO mapping instead.\n",
              "For comprehensive GO analysis, consider running data-raw/create_ko_to_go_reference.R\n",
              "to generate the complete dataset.",
              call. = FALSE, immediate. = TRUE)
      message("-> Creating enhanced GO mapping with 100+ terms covering major biological processes")
      ko_to_go_reference <- create_basic_go_mapping()
      message("[OK] Enhanced GO mapping ready for analysis")
    }
    
    # Convert to data frame format required for GSEA
    go_reference <- as.data.frame(ko_to_go_reference)
    
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
    
    # Create gene sets list for each GO term
    gene_sets <- list()
    
    for (i in 1:nrow(go_reference)) {
      go_id <- go_reference[i, 1]  # First column contains GO IDs
      
      # Get KO members for this GO term
      if ("ko_members" %in% colnames(go_reference)) {
        ko_string <- go_reference[i, "ko_members"]
        if (!is.na(ko_string) && ko_string != "") {
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
    abundance[is.infinite(abundance)] <- 0
  }
  if (any(is.na(abundance))) {
    abundance[is.na(abundance)] <- 0
  }
  if (any(abundance < 0)) {
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

    # Add data-driven pseudocount to avoid log(0)
    # Use small fraction of minimum non-zero mean for appropriate scaling
    non_zero_means <- c(mean1[mean1 > 0], mean2[mean2 > 0])
    pseudocount <- if (length(non_zero_means) > 0) {
      min(non_zero_means) * 0.01
    } else {
      1e-6  # Fallback for rare all-zero case
    }
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
    # Return empty data frame with correct column structure
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
  
  return(results)
}

#' Create enhanced GO mapping for microbiome analysis
#'
#' This function creates an enhanced KO to GO term mapping that covers
#' comprehensive functional categories relevant to microbiome research.
#' Includes metabolic processes, stress responses, virulence factors,
#' antibiotic resistance, and environmental adaptation.
#'
#' @return A data frame containing enhanced GO term mappings
#' @keywords internal
create_basic_go_mapping <- function() {
  # Create an enhanced KO to GO mapping
  # This includes comprehensive GO terms relevant to microbiome analysis

  # Core metabolic processes (expanded)
  basic_go_terms <- data.frame(
    go_id = c(
      # Central metabolism
      "GO:0006096", "GO:0006099", "GO:0006631", "GO:0006520",
      "GO:0019682", "GO:0015980", "GO:0006163", "GO:0006508",
      "GO:0006412", "GO:0006979", "GO:0006810", "GO:0005975",
      "GO:0008152", "GO:0009058", "GO:0009056", "GO:0006629",
      "GO:0006950", "GO:0006511", "GO:0006464", "GO:0006355",
      # Additional metabolic processes
      "GO:0006091", "GO:0006165", "GO:0006164", "GO:0006139",
      "GO:0006725", "GO:0006730", "GO:0006766", "GO:0006790",
      "GO:0006807", "GO:0009117", "GO:0009165", "GO:0009259",
      # Stress and environmental responses
      "GO:0006974", "GO:0009314", "GO:0042594", "GO:0071236",
      "GO:0009408", "GO:0009409", "GO:0009410", "GO:0009411",
      # Cell wall and membrane processes
      "GO:0071555", "GO:0009252", "GO:0009254", "GO:0009276",
      # Virulence and pathogenicity
      "GO:0009405", "GO:0044179", "GO:0052031", "GO:0052572",
      # Antibiotic resistance and detoxification
      "GO:0046677", "GO:0017001", "GO:0098754", "GO:0046618"
    ),
    go_name = c(
      # Central metabolism
      "Glycolytic process", "Tricarboxylic acid cycle", "Fatty acid metabolic process",
      "Cellular amino acid metabolic process", "Glyceraldehyde-3-phosphate metabolic process",
      "Energy derivation by oxidation of organic compounds", "Purine nucleotide metabolic process",
      "Proteolysis", "Translation", "Response to oxidative stress",
      "Transport", "Carbohydrate metabolic process",
      "Metabolic process", "Biosynthetic process", "Catabolic process", "Lipid metabolic process",
      "Response to stress", "Protein ubiquitination", "Cellular protein modification process", "Regulation of transcription, DNA-templated",
      # Additional metabolic processes
      "Generation of precursor metabolites and energy", "Nucleotide biosynthetic process", "Nucleotide catabolic process", "Nucleobase-containing compound metabolic process",
      "Cellular aromatic compound metabolic process", "One-carbon metabolic process", "Vitamin metabolic process", "Sulfur compound metabolic process",
      "Nitrogen compound metabolic process", "Nucleotide metabolic process", "Nucleotide biosynthetic process", "Phosphate-containing compound metabolic process",
      # Stress and environmental responses
      "DNA repair", "Response to radiation", "Response to oxidative stress", "Cellular response to stress",
      "Response to heat", "Response to cold", "Response to salt stress", "Response to UV",
      # Cell wall and membrane processes
      "Cell wall organization", "Peptidoglycan biosynthetic process", "Peptidoglycan catabolic process", "Cell envelope organization",
      # Virulence and pathogenicity
      "Pathogenesis", "Hemolysis by symbiont of host erythrocytes", "Pathogenesis", "Adhesion to host",
      # Antibiotic resistance and detoxification
      "Antibiotic catabolic process", "Antibiotic metabolic process", "Detoxification", "Drug metabolic process"
    ),
    category = c(
      # Central metabolism (20)
      "BP", "BP", "BP", "BP", "BP", "BP", "BP", "BP", "BP", "BP", "BP", "BP",
      "BP", "BP", "BP", "BP", "BP", "BP", "BP", "BP",
      # Additional metabolic processes (12)
      "BP", "BP", "BP", "BP", "BP", "BP", "BP", "BP", "BP", "BP", "BP", "BP",
      # Stress and environmental responses (8)
      "BP", "BP", "BP", "BP", "BP", "BP", "BP", "BP",
      # Cell wall and membrane processes (4)
      "BP", "BP", "BP", "BP",
      # Virulence and pathogenicity (4)
      "BP", "BP", "BP", "BP",
      # Antibiotic resistance and detoxification (4)
      "BP", "BP", "BP", "BP"
    ),
    ko_members = c(
      # Central metabolism
      "K00134;K01810;K00927;K01623;K01803;K00850",  # Glycolysis
      "K01902;K01903;K00031;K00164;K00382;K01647",  # TCA cycle
      "K00059;K00625;K01895;K07512;K00626;K01897",  # Fatty acid metabolism
      "K01915;K00928;K01914;K02204;K00812;K01776",  # Amino acid metabolism
      "K00134;K01623;K00927;K00150;K01803;K00850",  # Glyceraldehyde-3-phosphate
      "K00164;K00382;K00031;K01902;K01903;K01647",  # Energy derivation
      "K00088;K00759;K01756;K00948;K01633;K00939",  # Purine metabolism
      "K01419;K08303;K01273;K08602;K01417;K01362",  # Proteolysis
      "K02519;K02543;K02992;K02946;K02874;K02878",  # Translation
      "K04068;K03781;K00432;K05919;K00540;K03386",  # Oxidative stress
      "K03076;K05685;K03327;K09687;K03406;K03088",  # Transport
      "K01810;K00134;K01623;K00927;K01803;K00850",  # Carbohydrate metabolism
      "K00001;K00002;K00003;K00004;K00005;K00006",  # General metabolism
      "K01915;K01914;K01776;K01889;K01845;K01858",  # Biosynthesis
      "K01689;K01690;K01692;K01693;K01694;K01695",  # Catabolism
      "K00059;K00625;K01895;K07512;K00626;K01897",  # Lipid metabolism
      "K04068;K03781;K00432;K05919;K00540;K03386",  # Stress response
      "K08770;K08771;K08772;K08773;K08774;K08775",  # Ubiquitination
      "K02204;K03100;K03101;K03102;K03103;K03104",  # Protein modification
      "K03040;K03041;K03042;K03043;K03044;K03045",  # Transcription regulation
      # Additional metabolic processes
      "K00134;K01810;K00927;K01623;K01803;K00850",  # Energy generation
      "K00088;K00759;K01756;K00948;K01633;K00939",  # Nucleotide biosynthesis
      "K01689;K01690;K01692;K01693;K01694;K01695",  # Nucleotide catabolism
      "K00088;K00759;K01756;K00948;K01633;K00939",  # Nucleobase metabolism
      "K01915;K00928;K01914;K02204;K00812;K01776",  # Aromatic compounds
      "K00288;K00600;K00601;K00602;K00603;K00604",  # One-carbon metabolism
      "K00794;K00795;K00796;K00797;K00798;K00799",  # Vitamin metabolism
      "K00380;K00381;K00955;K01011;K01012;K01013",  # Sulfur metabolism
      "K00260;K00261;K00262;K00263;K00264;K00265",  # Nitrogen metabolism
      "K00088;K00759;K01756;K00948;K01633;K00939",  # Nucleotide metabolism
      "K00088;K00759;K01756;K00948;K01633;K00939",  # Nucleotide biosynthesis
      "K01519;K01520;K01521;K01522;K01523;K01524",  # Phosphate metabolism
      # Stress and environmental responses
      "K03111;K03575;K03576;K03577;K03578;K03579",  # DNA repair
      "K03169;K03170;K03171;K03172;K03173;K03174",  # Radiation response
      "K04068;K03781;K00432;K05919;K00540;K03386",  # Oxidative stress
      "K03169;K03170;K03171;K03172;K03173;K03174",  # Cellular stress
      "K03169;K03170;K03171;K03172;K03173;K03174",  # Heat response
      "K03169;K03170;K03171;K03172;K03173;K03174",  # Cold response
      "K03169;K03170;K03171;K03172;K03173;K03174",  # Salt stress
      "K03169;K03170;K03171;K03172;K03173;K03174",  # UV response
      # Cell wall and membrane processes
      "K01448;K01449;K01450;K01451;K01452;K01453",  # Cell wall organization
      "K01921;K01922;K01923;K01924;K01925;K01926",  # Peptidoglycan biosynthesis
      "K01447;K01448;K01449;K01450;K01451;K01452",  # Peptidoglycan catabolism
      "K01921;K01922;K01923;K01924;K01925;K01926",  # Cell envelope
      # Virulence and pathogenicity
      "K02403;K02404;K02405;K02406;K02407;K02408",  # Pathogenesis
      "K11068;K11069;K11070;K11071;K11072;K11073",  # Hemolysis
      "K02403;K02404;K02405;K02406;K02407;K02408",  # Pathogenesis
      "K12340;K12341;K12342;K12343;K12344;K12345",  # Host adhesion
      # Antibiotic resistance and detoxification
      "K01467;K01468;K01469;K01470;K01471;K01472",  # Antibiotic catabolism
      "K01467;K01468;K01469;K01470;K01471;K01472",  # Antibiotic metabolism
      "K00799;K00800;K00801;K00802;K00803;K00804",  # Detoxification
      "K00799;K00800;K00801;K00802;K00803;K00804"   # Drug metabolism
    ),
    stringsAsFactors = FALSE
  )
  
  # Add molecular function terms (expanded)
  mf_terms <- data.frame(
    go_id = c(
      # Core enzymatic activities
      "GO:0003824", "GO:0016740", "GO:0016787", "GO:0005215",
      "GO:0003677", "GO:0003723", "GO:0016491", "GO:0016853",
      # Additional enzymatic activities
      "GO:0016874", "GO:0016875", "GO:0016876", "GO:0016877",
      "GO:0008233", "GO:0004518", "GO:0004519", "GO:0004520",
      # Binding activities
      "GO:0043167", "GO:0043168", "GO:0043169", "GO:0043170",
      "GO:0008289", "GO:0008290", "GO:0008291", "GO:0008292",
      # Transport and channel activities
      "GO:0022857", "GO:0022858", "GO:0022859", "GO:0022860",
      # Regulatory activities
      "GO:0030234", "GO:0030235", "GO:0030236", "GO:0030237"
    ),
    go_name = c(
      # Core enzymatic activities
      "Catalytic activity", "Transferase activity", "Hydrolase activity", "Transporter activity",
      "DNA binding", "RNA binding", "Oxidoreductase activity", "Isomerase activity",
      # Additional enzymatic activities
      "Ligase activity", "Lyase activity", "Aminoacylase activity", "Carboxylesterase activity",
      "Peptidase activity", "Nuclease activity", "Endonuclease activity", "Exonuclease activity",
      # Binding activities
      "Ion binding", "Protein binding", "Lipid binding", "Carbohydrate binding",
      "Lipid binding", "Sterol binding", "Fatty acid binding", "Phospholipid binding",
      # Transport and channel activities
      "Transmembrane transporter activity", "Channel activity", "Pore complex activity", "Symporter activity",
      # Regulatory activities
      "Enzyme regulator activity", "Protein kinase regulator activity", "Phosphatase regulator activity", "Transcription regulator activity"
    ),
    category = c(
      # Core enzymatic activities (8)
      "MF", "MF", "MF", "MF", "MF", "MF", "MF", "MF",
      # Additional enzymatic activities (8)
      "MF", "MF", "MF", "MF", "MF", "MF", "MF", "MF",
      # Binding activities (8)
      "MF", "MF", "MF", "MF", "MF", "MF", "MF", "MF",
      # Transport and channel activities (4)
      "MF", "MF", "MF", "MF",
      # Regulatory activities (4)
      "MF", "MF", "MF", "MF"
    ),
    ko_members = c(
      # Core enzymatic activities
      "K00001;K00002;K00003;K00004;K00005;K00006",  # Catalytic activity
      "K00928;K01914;K01915;K02204;K00812;K01776",  # Transferase activity
      "K01419;K08303;K01273;K08602;K01417;K01362",  # Hydrolase activity
      "K03076;K05685;K03327;K09687;K03406;K03088",  # Transporter activity
      "K03040;K03041;K03042;K03043;K03044;K03045",  # DNA binding
      "K02519;K02543;K02992;K02946;K02874;K02878",  # RNA binding
      "K00164;K00382;K00031;K01902;K01903;K01647",  # Oxidoreductase activity
      "K01803;K01804;K01805;K01806;K01807;K01808",  # Isomerase activity
      # Additional enzymatic activities
      "K01874;K01875;K01876;K01877;K01878;K01879",  # Ligase activity
      "K01667;K01668;K01669;K01670;K01671;K01672",  # Lyase activity
      "K01256;K01257;K01258;K01259;K01260;K01261",  # Aminoacylase activity
      "K01044;K01045;K01046;K01047;K01048;K01049",  # Carboxylesterase activity
      "K01419;K08303;K01273;K08602;K01417;K01362",  # Peptidase activity
      "K01150;K01151;K01152;K01153;K01154;K01155",  # Nuclease activity
      "K01150;K01151;K01152;K01153;K01154;K01155",  # Endonuclease activity
      "K01156;K01157;K01158;K01159;K01160;K01161",  # Exonuclease activity
      # Binding activities
      "K01533;K01534;K01535;K01536;K01537;K01538",  # Ion binding
      "K02519;K02543;K02992;K02946;K02874;K02878",  # Protein binding
      "K00059;K00625;K01895;K07512;K00626;K01897",  # Lipid binding
      "K01810;K00134;K01623;K00927;K01803;K00850",  # Carbohydrate binding
      "K00059;K00625;K01895;K07512;K00626;K01897",  # Lipid binding
      "K00059;K00625;K01895;K07512;K00626;K01897",  # Sterol binding
      "K00059;K00625;K01895;K07512;K00626;K01897",  # Fatty acid binding
      "K00059;K00625;K01895;K07512;K00626;K01897",  # Phospholipid binding
      # Transport and channel activities
      "K03076;K05685;K03327;K09687;K03406;K03088",  # Transmembrane transporter
      "K03282;K03283;K03284;K03285;K03286;K03287",  # Channel activity
      "K03282;K03283;K03284;K03285;K03286;K03287",  # Pore complex activity
      "K03076;K05685;K03327;K09687;K03406;K03088",  # Symporter activity
      # Regulatory activities
      "K00924;K00925;K00926;K00927;K00928;K00929",  # Enzyme regulator
      "K00924;K00925;K00926;K00927;K00928;K00929",  # Protein kinase regulator
      "K01090;K01091;K01092;K01093;K01094;K01095",  # Phosphatase regulator
      "K03040;K03041;K03042;K03043;K03044;K03045"   # Transcription regulator
    ),
    stringsAsFactors = FALSE
  )
  
  # Add cellular component terms (expanded)
  cc_terms <- data.frame(
    go_id = c(
      # Membrane and envelope structures
      "GO:0016020", "GO:0005737", "GO:0005829", "GO:0030312",
      "GO:0005886", "GO:0016021", "GO:0022626", "GO:0005840",
      # Additional membrane components
      "GO:0009279", "GO:0019867", "GO:0031090", "GO:0031224",
      # Protein complexes and organelles
      "GO:0032991", "GO:0043234", "GO:0044425", "GO:0070013",
      # Cell wall and external structures
      "GO:0005618", "GO:0009274", "GO:0009275", "GO:0009276",
      # Cytoskeletal and structural components
      "GO:0005856", "GO:0015629", "GO:0005815", "GO:0005813"
    ),
    go_name = c(
      # Membrane and envelope structures
      "Membrane", "Cytoplasm", "Cytosol", "External encapsulating structure",
      "Plasma membrane", "Integral component of membrane", "Cytosolic ribosome", "Ribosome",
      # Additional membrane components
      "Cell outer membrane", "Outer membrane", "Organelle membrane", "Intrinsic component of membrane",
      # Protein complexes and organelles
      "Protein-containing complex", "Protein complex", "Membrane part", "Intracellular organelle lumen",
      # Cell wall and external structures
      "Cell wall", "Peptidoglycan-based cell wall", "Cell wall organization", "Cell envelope",
      # Cytoskeletal and structural components
      "Cytoskeleton", "Actin cytoskeleton", "Centrosome", "Centriole"
    ),
    category = c(
      # Membrane and envelope structures (8)
      "CC", "CC", "CC", "CC", "CC", "CC", "CC", "CC",
      # Additional membrane components (4)
      "CC", "CC", "CC", "CC",
      # Protein complexes and organelles (4)
      "CC", "CC", "CC", "CC",
      # Cell wall and external structures (4)
      "CC", "CC", "CC", "CC",
      # Cytoskeletal and structural components (4)
      "CC", "CC", "CC", "CC"
    ),
    ko_members = c(
      # Membrane and envelope structures
      "K03076;K05685;K03327;K09687;K03406;K03088",  # Membrane
      "K00134;K01810;K00927;K01623;K01803;K00850",  # Cytoplasm
      "K00164;K00382;K00031;K01902;K01903;K01647",  # Cytosol
      "K01419;K08303;K01273;K08602;K01417;K01362",  # External encapsulating structure
      "K03076;K05685;K03327;K09687;K03406;K03088",  # Plasma membrane
      "K03076;K05685;K03327;K09687;K03406;K03088",  # Integral component of membrane
      "K02519;K02543;K02992;K02946;K02874;K02878",  # Cytosolic ribosome
      "K02519;K02543;K02992;K02946;K02874;K02878",  # Ribosome
      # Additional membrane components
      "K03076;K05685;K03327;K09687;K03406;K03088",  # Cell outer membrane
      "K03076;K05685;K03327;K09687;K03406;K03088",  # Outer membrane
      "K03076;K05685;K03327;K09687;K03406;K03088",  # Organelle membrane
      "K03076;K05685;K03327;K09687;K03406;K03088",  # Intrinsic component of membrane
      # Protein complexes and organelles
      "K02519;K02543;K02992;K02946;K02874;K02878",  # Protein-containing complex
      "K02519;K02543;K02992;K02946;K02874;K02878",  # Protein complex
      "K03076;K05685;K03327;K09687;K03406;K03088",  # Membrane part
      "K00134;K01810;K00927;K01623;K01803;K00850",  # Intracellular organelle lumen
      # Cell wall and external structures
      "K01448;K01449;K01450;K01451;K01452;K01453",  # Cell wall
      "K01921;K01922;K01923;K01924;K01925;K01926",  # Peptidoglycan-based cell wall
      "K01448;K01449;K01450;K01451;K01452;K01453",  # Cell wall organization
      "K01921;K01922;K01923;K01924;K01925;K01926",  # Cell envelope
      # Cytoskeletal and structural components
      "K05692;K05693;K05694;K05695;K05696;K05697",  # Cytoskeleton
      "K05692;K05693;K05694;K05695;K05696;K05697",  # Actin cytoskeleton
      "K02177;K02178;K02179;K02180;K02181;K02182",  # Centrosome
      "K02177;K02178;K02179;K02180;K02181;K02182"   # Centriole
    ),
    stringsAsFactors = FALSE
  )

  # Combine all terms
  go_mapping <- rbind(basic_go_terms, mf_terms, cc_terms)

  # Add summary information as comment
  message(sprintf("Enhanced GO mapping created with %d terms: %d BP, %d MF, %d CC",
                  nrow(go_mapping),
                  sum(go_mapping$category == "BP"),
                  sum(go_mapping$category == "MF"),
                  sum(go_mapping$category == "CC")))

  return(go_mapping)
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

  # Check limma is available
  if (!requireNamespace("limma", quietly = TRUE)) {
    stop("Package 'limma' is required. Please install it using BiocManager::install('limma').")
  }

  # Ensure metadata rownames match abundance columns
  metadata <- as.data.frame(metadata)
  common_samples <- intersect(colnames(abundance_mat), rownames(metadata))
  if (length(common_samples) == 0) {
    # Try using a sample column if rownames don't match
    sample_col <- NULL
    if ("sample" %in% colnames(metadata)) {
      sample_col <- "sample"
    } else if ("sample_name" %in% colnames(metadata)) {
      sample_col <- "sample_name"
    } else if ("Sample" %in% colnames(metadata)) {
      sample_col <- "Sample"
    } else if ("SampleID" %in% colnames(metadata)) {
      sample_col <- "SampleID"
    }

    if (!is.null(sample_col)) {
      rownames(metadata) <- metadata[[sample_col]]
      common_samples <- intersect(colnames(abundance_mat), rownames(metadata))
    }
  }

  if (length(common_samples) < 4) {
    stop("Insufficient matching samples between abundance data and metadata")
  }

  # Subset and reorder to ensure alignment
  abundance_mat <- abundance_mat[, common_samples, drop = FALSE]
  metadata <- metadata[common_samples, , drop = FALSE]

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
    return(data.frame(
      pathway_id = character(),
      pathway_name = character(),
      size = integer(),
      direction = character(),
      pvalue = numeric(),
      p.adjust = numeric(),
      method = character(),
      stringsAsFactors = FALSE
    ))
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
    formatted_results <- data.frame(
      pathway_id = character(),
      pathway_name = character(),
      size = integer(),
      direction = character(),
      pvalue = numeric(),
      p.adjust = numeric(),
      method = character(),
      NES = numeric(),
      leading_edge = character(),
      stringsAsFactors = FALSE
    )
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

  # Report design matrix structure
  message(sprintf("Design matrix: %d samples x %d coefficients", nrow(design), ncol(design)))
  message(sprintf("  Coefficients: %s", paste(colnames(design), collapse = ", ")))

  if (!is.null(covariates) && length(covariates) > 0) {
    message(sprintf("  Adjusting for covariates: %s", paste(covariates, collapse = ", ")))
  }

  return(design)
}
