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
  
  invisible(TRUE)
}

format_metadata_row_labels <- function(metadata, rows) {
  labels <- rownames(metadata)
  if (is.null(labels) || length(labels) != nrow(metadata)) {
    labels <- as.character(seq_len(nrow(metadata)))
  }
  labels <- labels[rows]
  bad_labels <- is.na(labels) | !nzchar(labels)
  labels[bad_labels] <- as.character(rows[bad_labels])
  labels
}

validate_gsea_group_labels_after_alignment <- function(metadata, group) {
  group_values <- metadata[[group]]
  group_chr <- as.character(group_values)
  invalid_group <- which(is.na(group_chr) | !nzchar(trimws(group_chr)))
  if (length(invalid_group) > 0) {
    stop(
      "Group column '", group, "' contains missing or empty values after sample alignment for sample(s): ",
      paste(head(format_metadata_row_labels(metadata, invalid_group), 5), collapse = ", "),
      ". Remove or impute these samples before GSEA.",
      call. = FALSE
    )
  }

  invisible(TRUE)
}

validate_gsea_group_after_alignment <- function(metadata, group) {
  validate_gsea_group_labels_after_alignment(metadata, group)

  group_values <- metadata[[group]]

  group_factor <- droplevels(factor(group_values))
  group_counts <- table(group_factor)
  if (length(group_counts) < 2) {
    found <- if (length(group_counts) == 0) {
      "none"
    } else {
      paste(names(group_counts), "=", as.integer(group_counts), collapse = ", ")
    }
    stop(
      "GSEA requires at least 2 groups with samples after sample alignment; found ",
      found,
      ".",
      call. = FALSE
    )
  }

  if (any(group_counts < 2)) {
    stop(
      "Each group must have at least 2 samples after sample alignment. Current group sizes: ",
      paste(names(group_counts), "=", as.integer(group_counts), collapse = ", "),
      ".",
      call. = FALSE
    )
  }

  invisible(TRUE)
}

resolve_preranked_comparison <- function(metadata, group, comparison = NULL) {
  if (!group %in% colnames(metadata)) {
    stop("Group column '", group, "' not found in metadata.", call. = FALSE)
  }

  group_values <- metadata[[group]]
  validate_gsea_group_labels_after_alignment(metadata, group)

  if (is.null(comparison)) {
    validate_gsea_group_after_alignment(metadata, group)
    group_factor <- droplevels(factor(group_values))
    group_levels <- levels(group_factor)

    if (length(group_levels) != 2) {
      group_counts <- table(group_factor)
      stop(
        "Preranked GSEA requires exactly two group levels when 'comparison' is NULL; found ",
        paste(names(group_counts), "=", as.integer(group_counts), collapse = ", "),
        ". Supply comparison = c(group1, group2) to define the ranking direction.",
        call. = FALSE
      )
    }

    return(group_levels)
  }

  if (!is.character(comparison) || length(comparison) != 2 ||
      any(is.na(comparison)) || any(!nzchar(comparison))) {
    stop(
      "'comparison' must be a length-2 character vector such as c('Treatment', 'Control').",
      call. = FALSE
    )
  }
  if (comparison[1] == comparison[2]) {
    stop("'comparison' must contain two distinct group levels.", call. = FALSE)
  }

  group_factor <- droplevels(factor(group_values))
  group_counts <- table(group_factor)
  available_groups <- names(group_counts)
  missing_groups <- setdiff(comparison, available_groups)
  if (length(missing_groups) > 0) {
    stop(
      "'comparison' group level(s) not found after sample alignment: ",
      paste(missing_groups, collapse = ", "),
      ". Available group levels: ",
      if (length(available_groups) > 0) paste(available_groups, collapse = ", ") else "none",
      ".",
      call. = FALSE
    )
  }

  selected_counts <- group_counts[comparison]
  if (any(selected_counts < 2)) {
    stop(
      "Each selected 'comparison' group must have at least 2 samples after sample alignment. Current group sizes: ",
      paste(comparison, "=", as.integer(selected_counts), collapse = ", "),
      ".",
      call. = FALSE
    )
  }

  comparison
}

validate_complete_design_variables <- function(metadata, group, covariates = NULL) {
  if (!is.character(group) || length(group) != 1 ||
      is.na(group) || !nzchar(group)) {
    stop("'group' must be a single non-empty character string.",
         call. = FALSE)
  }

  if (!is.null(covariates)) {
    if (!is.character(covariates) || anyNA(covariates) ||
        any(!nzchar(covariates))) {
      stop("'covariates' must be a character vector of non-empty metadata column names.",
           call. = FALSE)
    }
    if (anyDuplicated(covariates)) {
      duplicated_covariates <- unique(covariates[duplicated(covariates)])
      stop("'covariates' must contain unique column names. Duplicated: ",
           paste(duplicated_covariates, collapse = ", "),
           call. = FALSE)
    }
    if (group %in% covariates) {
      stop("'covariates' must not include the group column '", group, "'.",
           call. = FALSE)
    }
  }

  design_vars <- c(group, covariates)
  missing_cols <- setdiff(design_vars, colnames(metadata))
  if (length(missing_cols) > 0) {
    stop("Design variable(s) not found in metadata: ",
         paste(missing_cols, collapse = ", "),
         call. = FALSE)
  }

  model_data <- metadata[, design_vars, drop = FALSE]
  missing_matrix <- is.na(model_data)
  missing_rows <- which(rowSums(missing_matrix) > 0)
  if (length(missing_rows) > 0) {
    missing_vars <- colnames(model_data)[colSums(missing_matrix) > 0]
    stop(
      "Design variables contain missing values after sample alignment. ",
      "Column(s): ", paste(missing_vars, collapse = ", "),
      "; sample(s): ",
      paste(head(format_metadata_row_labels(metadata, missing_rows), 5), collapse = ", "),
      ". Remove or impute missing design variables before running camera/fry GSEA.",
      call. = FALSE
    )
  }

  invisible(TRUE)
}

#' Gene Set Enrichment Analysis for PICRUSt2 output
#'
#' This function performs Gene Set Enrichment Analysis (GSEA) on PICRUSt2 predicted functional data
#' to identify enriched pathways between different conditions.
#'
#' @param abundance A data frame containing gene/enzyme abundance data, with
#'   feature IDs in row names and samples as columns. Data frames may also
#'   provide a leading non-numeric feature ID column (for example \code{#NAME},
#'   \code{feature}, or \code{pathway}); it is converted to row names
#'   automatically. Values must be finite, non-missing, and
#'   non-negative count-like abundances; negative or non-finite values are
#'   rejected rather than being coerced to zero.
#'   For KEGG analysis: features should be KO IDs (e.g., K00001).
#'   For MetaCyc analysis: features should be EC numbers (e.g., EC:1.1.1.1 or 1.1.1.1), NOT pathway IDs.
#'   MetaCyc pathway-like identifiers are rejected because GSEA requires
#'   gene/enzyme-level input; use \code{\link{pathway_daa}} for pathway-level
#'   MetaCyc abundance tables.
#'   For GO analysis: features should be KO IDs that will be mapped to GO terms.
#'   NOTE: This function requires gene-level data, not pathway-level abundances.
#'   For pathway abundance analysis, use \code{\link{pathway_daa}} instead
#' @param metadata A data frame containing sample metadata. After sample
#'   alignment, all retained samples must have non-missing, non-empty values in
#'   the grouping column.
#' @param group A character string specifying the column name in metadata that
#'   contains the grouping variable
#' @param pathway_type A single character string specifying the pathway type:
#'   "KEGG", "MetaCyc", or "GO"
#' @param method A single character string specifying the GSEA method:
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
#'   for adjustment. Only supported when method is "camera" or "fry"; supplying
#'   covariates with preranked methods is an error because those methods use a
#'   precomputed rank vector rather than a design matrix. Default is NULL (no covariates).
#'   Covariate values must be complete for all aligned samples; rows with
#'   missing model variables are rejected rather than being silently dropped
#'   by \code{model.matrix()}. The resulting design matrix must also be finite
#'   and full rank; constant or fully confounded covariates are rejected.
#'   Example: \code{covariates = c("age", "sex", "BMI")}
#' @param contrast For "camera" or "fry" methods, specify the coefficient
#'   or contrast to test. Supplying \code{contrast} with preranked methods is
#'   an error; use \code{comparison} instead. Default \code{NULL} automatically tests the single
#'   non-reference group coefficient in two-group designs. Multi-group designs
#'   must specify this explicitly. A character value must exactly match a
#'   design column name or a non-reference group level; substring matching is
#'   not used. A numeric scalar is treated as a design-column index, and a
#'   numeric vector must have length equal to the number of design columns.
#'   Named numeric vectors are matched and reordered by design column names;
#'   unnamed numeric vectors are interpreted in design column order.
#' @param comparison For preranked methods only (\code{"fgsea"},
#'   \code{"GSEA"}, or \code{"clusterProfiler"}), an optional length-2
#'   character vector \code{c(group1, group2)} defining the ranking direction.
#'   Ranking statistics are calculated as group1 versus group2: positive
#'   \code{signal2noise}, \code{t_test}, \code{diff_abundance}, and
#'   \code{log2_ratio} values indicate higher abundance in \code{group1}
#'   (for \code{log2_ratio}, \code{log2(group1 / group2)}). If \code{NULL},
#'   exactly two aligned group levels must be present and their factor-level
#'   order is used. Multi-group preranked analyses must specify
#'   \code{comparison} explicitly.
#' @param inter.gene.cor Numeric value specifying the inter-gene correlation for camera method.
#'   Default is 0.01. Use NA to estimate correlation from data for each gene set.
#' @param rank_method A single character string specifying the ranking statistic for preranked methods
#'   (fgsea, GSEA, clusterProfiler): "signal2noise", "t_test", "log2_ratio", or "diff_abundance"
#' @param nperm An integer specifying the number of permutations (for clusterProfiler method only).
#'   The fgsea method uses adaptive multilevel splitting and does not require a fixed permutation count.
#' @param min_size An integer specifying the minimum gene set size
#' @param max_size An integer specifying the maximum gene set size
#' @param p_adjust_method A character string specifying the p-value adjustment method
#' @param p.adjust Deprecated. Use \code{p_adjust_method} instead.
#' @param seed An integer specifying the random seed for reproducibility
#' @param go_category A single character string specifying GO category to use.
#'   "all" (default) uses all categories present in the reference data.
#'   Valid categories are determined by the reference data (currently MF and CC).
#'   See \code{table(ko_to_go_reference$category)} for available categories.
#' @param organism Deprecated and has no effect. The KEGG and GO
#'   reference data bundled with ggpicrust2 are KO-based
#'   (organism-independent), so gene sets are returned in KO space
#'   regardless of this argument. Retained only for signature
#'   compatibility; passing any value other than the default
#'   \code{"ko"} emits a deprecation warning. Will be removed in a
#'   future release.
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
#'   For fgsea/clusterProfiler methods, additional columns include ES, NES,
#'   leading_edge, group1, and group2. Positive ES/NES values are in the
#'   \code{group1} versus \code{group2} direction. For camera/fry methods,
#'   limma does not return NES; ggpicrust2 retains a legacy \code{NES} column
#'   containing a signed \code{-log10(pvalue)} score for visualization
#'   compatibility and labels it with \code{score_type} and
#'   \code{score_label}.
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
#' For these methods, positive ES/NES values indicate gene-set enrichment near
#' the top of the ranked list. With \code{comparison = c(group1, group2)}, the
#' top of the list corresponds to features higher in \code{group1}; reverse
#' \code{comparison} to reverse the biological direction.
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
                        p_adjust_method = "BH",
                        seed = 42,
                        go_category = "all",
                        organism = "ko",
                        p.adjust = NULL,
                        comparison = NULL) {
  # Backward compatibility for deprecated parameter
  if (!is.null(p.adjust)) {
    warning("'p.adjust' parameter is deprecated. Use 'p_adjust_method' instead.", call. = FALSE)
    p_adjust_method <- p.adjust
  }
  
  # Normalize PICRUSt2-style input before validation so the feature ID column
  # is not counted as a sample or mixed into numeric matrix checks.
  abundance <- normalize_abundance_feature_ids(
    abundance,
    context = "pathway_gsea() abundance"
  )

  # Input validation using unified functions
  validate_abundance(abundance, min_samples = 4)
  validate_metadata(metadata)
  validate_group(metadata, group, min_groups = 2)
  
  validate_choice(pathway_type, c("KEGG", "MetaCyc", "GO"), "pathway_type")

  valid_methods <- c("camera", "fry", "fgsea", "GSEA", "clusterProfiler")
  validate_choice(method, valid_methods, "method")

  preranked_methods <- c("fgsea", "GSEA", "clusterProfiler")
  if (!is.null(comparison) && !method %in% preranked_methods) {
    stop(
      "'comparison' is only supported for preranked GSEA methods ('fgsea', 'GSEA', or 'clusterProfiler'). ",
      "Use 'contrast' to define camera/fry comparisons.",
      call. = FALSE
    )
  }
  if (!is.null(contrast) && method %in% preranked_methods) {
    stop(
      "'contrast' is only supported for limma-based GSEA methods ('camera' or 'fry'). ",
      "Use 'comparison' to define the ranking direction for preranked GSEA methods.",
      call. = FALSE
    )
  }

  validate_count_parameter(min_size, "min_size")
  validate_count_parameter(max_size, "max_size")
  if (min_size > max_size) {
    stop("'min_size' must be less than or equal to 'max_size'.", call. = FALSE)
  }

  validate_count_parameter(nperm, "nperm")

  if (!is.numeric(seed) || length(seed) != 1 ||
      is.na(seed) || !is.finite(seed) ||
      seed < 0 || seed != floor(seed)) {
    stop("'seed' must be a single non-negative integer.", call. = FALSE)
  }

  validate_p_adjust_method(p_adjust_method)

  # Validate rank_method only for preranked methods
  if (method %in% preranked_methods) {
    validate_choice(
      rank_method,
      c("signal2noise", "t_test", "log2_ratio", "diff_abundance"),
      "rank_method"
    )
  }

  # Validate covariates
  if (!is.null(covariates)) {
    if (method %in% preranked_methods) {
      stop("Covariates are only supported for 'camera' and 'fry' methods. ",
           "Preranked GSEA methods cannot use covariates in this wrapper; ",
           "remove 'covariates' or use method = 'camera'/'fry'.",
           call. = FALSE)
    } else {
      validate_complete_design_variables(metadata, group, covariates)
    }
  }

  # Validate inter.gene.cor for camera (correlation range: -1 to 1)
  if (method == "camera" && !is.null(inter.gene.cor)) {
    if (!is.numeric(inter.gene.cor) || length(inter.gene.cor) != 1 ||
        (!is.na(inter.gene.cor) &&
         (!is.finite(inter.gene.cor) || inter.gene.cor < -1 ||
          inter.gene.cor > 1))) {
      stop("inter.gene.cor must be a numeric value between -1 and 1, or NA")
    }
  }

  # GO category validation is deferred to prepare_gene_sets() where it is
  # checked against the actual reference data (data-driven, not hardcoded)

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
  if (method %in% preranked_methods) {
    message("Note: Preranked GSEA methods (fgsea, clusterProfiler) do not account for ",
            "inter-gene correlations, which may lead to unreliable p-values ",
            "(Wu et al., 2012). Consider using method='camera' or method='fry' for ",
            "more reliable statistical inference.")
  }
  
  # Set seed for reproducibility
  set.seed(seed)
  
  # Ensure abundance is a matrix with samples as columns
  abundance_mat <- as.matrix(abundance)
  validate_nonnegative_finite_matrix(abundance_mat, "abundance")
  validate_feature_rownames(abundance_mat, "abundance")
  
  # Input validation for MetaCyc pathway data
  if (pathway_type == "MetaCyc") {
    validate_metacyc_gsea_feature_ids(rownames(abundance_mat))
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

  if (method %in% preranked_methods && !is.null(comparison)) {
    comparison <- resolve_preranked_comparison(metadata, group, comparison)
  } else {
    validate_gsea_group_after_alignment(metadata, group)
  }
  if (method %in% c("camera", "fry")) {
    validate_complete_design_variables(metadata, group, covariates)
  }
  
  # Prepare gene sets
  gene_sets <- prepare_gene_sets(pathway_type, organism = organism, go_category = go_category)
  gene_sets <- validate_gene_sets(gene_sets, "pathway_gsea() gene_sets")

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
      p.adjust.method = p_adjust_method
    )

  } else if (method == "fgsea") {
    # =========================================================================
    # Preranked methods (fgsea)
    # =========================================================================

    # Calculate ranking metric for preranked methods
    ranked_list <- calculate_rank_metric(
      abundance_mat,
      metadata,
      group,
      method = rank_method,
      comparison = comparison
    )
    ranked_list <- validate_ranked_list(
      ranked_list,
      "pathway_gsea() ranked_list"
    )
    results <- run_fgsea(
      ranked_list,
      gene_sets,
      min_size,
      max_size,
      p_adjust_method = p_adjust_method
    )

  } else if (method == "GSEA" || method == "clusterProfiler") {
    # =========================================================================
    # clusterProfiler GSEA
    # =========================================================================

    # Calculate ranking metric for preranked methods
    ranked_list <- calculate_rank_metric(
      abundance_mat,
      metadata,
      group,
      method = rank_method,
      comparison = comparison
    )
    ranked_list <- validate_ranked_list(
      ranked_list,
      "pathway_gsea() ranked_list"
    )

    # Convert ranked list to format required by clusterProfiler
    gene_list <- sort(ranked_list, decreasing = TRUE)

    gene_sets_filtered <- filter_gene_sets_to_ranked_universe(
      gene_sets,
      names(gene_list),
      min_size,
      max_size
    )

    if (length(gene_sets_filtered) == 0) {
      warning(
        "No gene sets overlapped the ranked feature list after applying ",
        "min_size=", min_size, " and max_size=", max_size,
        ". Returning empty GSEA results.",
        call. = FALSE
      )
      results <- create_empty_gsea_result(method, full = TRUE)
    } else {
      # Run GSEA using clusterProfiler
      gsea_result <- clusterProfiler::GSEA(
        geneList = gene_list,
        TERM2GENE = data.frame(
          term = rep(names(gene_sets_filtered), lengths(gene_sets_filtered)),
          gene = unlist(gene_sets_filtered, use.names = FALSE),
          stringsAsFactors = FALSE
        ),
        minGSSize = min_size,
        maxGSSize = max_size,
        pvalueCutoff = 1,
        pAdjustMethod = p_adjust_method,
        nPermSimple = nperm,
        seed = seed
      )

      # Convert results to data frame. A valid gseaResult can contain zero
      # rows, in which case some optional columns (for example
      # core_enrichment) are absent. Preserve the standard ggpicrust2 schema
      # instead of returning a partial data frame.
      gsea_df <- if (is.null(gsea_result)) {
        data.frame()
      } else {
        as.data.frame(gsea_result)
      }

      if (nrow(gsea_df) > 0) {
        leading_edge <- if ("core_enrichment" %in% colnames(gsea_df)) {
          gsub("/", ";", gsea_df$core_enrichment, fixed = TRUE)
        } else {
          rep("", nrow(gsea_df))
        }

        # Rename columns to match fgsea output
        results <- data.frame(
          pathway_id = gsea_df$ID,
          pathway_name = gsea_df$Description,
          size = gsea_df$setSize,
          ES = gsea_df$enrichmentScore,
          NES = gsea_df$NES,
          pvalue = gsea_df$pvalue,
          p.adjust = gsea_df$p.adjust,
          leading_edge = leading_edge,
          stringsAsFactors = FALSE
        )
      } else {
        results <- create_empty_gsea_result(method, full = TRUE)
      }
    }
  }

  # Add method information to populated results
  # (empty results already have method from factory)
  if (nrow(results) > 0) {
    results$method <- method
    if (method %in% preranked_methods) {
      results$score_type <- "NES"
      results$score_label <- "Normalized Enrichment Score (NES)"
      rank_comparison <- attr(ranked_list, "comparison", exact = TRUE)
      if (length(rank_comparison) == 2) {
        results$group1 <- rank_comparison[1]
        results$group2 <- rank_comparison[2]
      }
    }
  }
  
  return(results)
}

validate_metacyc_gsea_feature_ids <- function(feature_names) {
  feature_names <- validate_nonempty_character_column(
    feature_names,
    "feature IDs",
    "pathway_gsea() MetaCyc abundance"
  )
  pathway_like <- is_pathway_id(feature_names)
  if (any(pathway_like)) {
    examples <- unique(feature_names[pathway_like])
    stop(
      "pathway_gsea(pathway_type = 'MetaCyc') requires gene/enzyme-level ",
      "EC abundance input, but the abundance row names contain pathway-like ",
      "identifier(s): ",
      paste(utils::head(examples, 5), collapse = ", "),
      ". Use EC_metagenome_out/pred_metagenome_unstrat.tsv for MetaCyc GSEA, ",
      "or use pathway_daa() for pathway-level MetaCyc abundance tables.",
      call. = FALSE
    )
  }
  invisible(TRUE)
}

#' Prepare gene sets for GSEA
#'
#' @param pathway_type A single character string specifying the pathway type:
#'   "KEGG", "MetaCyc", or "GO"
#' @param organism Deprecated and has no effect; gene sets are KO-based
#'   for both KEGG and GO. See \code{\link{pathway_gsea}} for details.
#'   Retained for signature compatibility only.
#' @param go_category A single character string specifying the GO category to use.
#'   "all" (default) uses all categories. Valid values are determined by
#'   the reference data; see \code{table(ko_to_go_reference$category)}.
#'
#' @return A list of pathway gene sets
#' @export
prepare_gene_sets <- function(pathway_type = "KEGG", organism = "ko", go_category = "all") {

  # Validate pathway_type
  valid_types <- c("KEGG", "MetaCyc", "GO")
  validate_choice(pathway_type, valid_types, "pathway_type")

  if (!is.character(go_category) || length(go_category) != 1 ||
      is.na(go_category) || !nzchar(go_category)) {
    stop("'go_category' must be a single non-empty character value.",
         call. = FALSE)
  }

  # Soft-deprecate `organism`. The KEGG and GO reference tables loaded
  # below -- `ko_to_kegg_reference` and `ko_to_go_reference` -- are
  # both keyed on KO identifiers and are organism-independent by
  # construction, so this argument never actually influenced either
  # branch. A caller that set `organism = "hsa"` expecting
  # human-specific pathways silently got the same KO gene sets as
  # `organism = "ko"`. Raise a visible deprecation warning at the
  # single choke point both `pathway_gsea()` and direct callers go
  # through, rather than letting the promise-implementation gap stay
  # silent. The default value is unchanged so existing callers that
  # relied on the (ignored) default keep working without a warning.
  if (!identical(organism, "ko")) {
    warning(sprintf(
      "'organism' argument is deprecated and has no effect (got '%s'). ",
      as.character(organism)[1]),
      "The KEGG and GO reference data bundled with ggpicrust2 are ",
      "KO-based (organism-independent), so all gene sets are returned ",
      "in KO space regardless of this argument. This parameter will ",
      "be removed in a future release.",
      call. = FALSE)
  }

  if (pathway_type == "KEGG") {
    # Load KEGG reference using unified loader
    ko_to_kegg_reference <- load_reference_data("ko_to_kegg")
    ko_to_kegg_reference <- filter_kegg_reference_to_pathways(ko_to_kegg_reference)

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

    # Derive valid categories from the actual reference data
    available_categories <- if ("category" %in% colnames(go_reference)) {
      unique(go_reference$category)
    } else {
      character(0)
    }

    # Validate go_category against what's actually in the data
    if (!is.null(go_category) && go_category != "all" &&
        !go_category %in% available_categories) {
      stop(sprintf(
        "go_category '%s' not found in reference data. Available: %s",
        go_category, paste(c(available_categories, "all"), collapse = ", ")),
        call. = FALSE)
    }

    # Filter by category
    if (!is.null(go_category) && go_category != "all" &&
        "category" %in% colnames(go_reference)) {
      go_reference <- go_reference[go_reference$category == go_category, ]
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
  
  validate_gene_sets(gene_sets, "prepare_gene_sets() gene_sets")
}

#' Validate a named gene-set list
#'
#' Gene-set names are pathway identifiers. Missing or duplicated names must be
#' rejected before passing the object to limma/fgsea because downstream code may
#' auto-repair duplicate row names (for example `set1` -> `set1.1`), which would
#' corrupt pathway IDs in the returned result table.
#'
#' @param gene_sets Named list of feature identifier vectors.
#' @param context Label used in error messages.
#' @return Sanitized gene-set list with character, unique member vectors.
#' @keywords internal
#' @noRd
validate_gene_sets <- function(gene_sets, context = "gene_sets") {
  if (!is.list(gene_sets)) {
    stop(context, " must be a named list of feature identifier vectors.",
         call. = FALSE)
  }

  set_names <- names(gene_sets)
  if (is.null(set_names) || length(set_names) != length(gene_sets)) {
    stop(context, " must have one name per gene set.", call. = FALSE)
  }
  set_names <- validate_nonempty_character_column(set_names, "names", context)
  if (anyDuplicated(set_names)) {
    duplicated_names <- unique(set_names[duplicated(set_names)])
    stop(
      context,
      " contains duplicated gene-set names: ",
      paste(utils::head(duplicated_names, 5), collapse = ", "),
      ". Gene-set names are used as pathway_id values and must be unique.",
      call. = FALSE
    )
  }

  sanitized <- vector("list", length(gene_sets))
  for (i in seq_along(gene_sets)) {
    members <- gene_sets[[i]]
    if (is.null(members) || !is.atomic(members)) {
      stop(
        context,
        " gene set '", set_names[i],
        "' must be an atomic vector of feature identifiers.",
        call. = FALSE
      )
    }
    members <- validate_nonempty_character_column(
      members,
      paste0("gene set '", set_names[i], "'"),
      context
    )
    sanitized[[i]] <- unique(members)
  }

  names(sanitized) <- set_names
  sanitized
}

#' Validate a named ranking-statistic vector for preranked GSEA
#'
#' @param ranked_list Named numeric vector of feature-level ranking statistics.
#' @param context Label used in error messages.
#' @param require_variation Logical. If TRUE, rejects ranking vectors where all
#'   finite statistics are tied.
#'
#' @return The validated ranked_list with normalized character names.
#' @keywords internal
#' @noRd
validate_ranked_list <- function(ranked_list,
                                 context = "ranked_list",
                                 require_variation = TRUE) {
  if (!is.atomic(ranked_list) || !is.numeric(ranked_list) ||
      !is.null(dim(ranked_list))) {
    stop(context, " must be a named numeric vector of ranking statistics.",
         call. = FALSE)
  }

  if (length(ranked_list) == 0) {
    stop(context, " must contain at least one ranking statistic.",
         call. = FALSE)
  }

  ranked_names <- names(ranked_list)
  if (is.null(ranked_names) || length(ranked_names) != length(ranked_list)) {
    stop(context, " must have one feature name per ranking statistic.",
         call. = FALSE)
  }
  ranked_names <- validate_nonempty_character_column(
    ranked_names,
    "names",
    context
  )
  if (anyDuplicated(ranked_names)) {
    duplicated_names <- unique(ranked_names[duplicated(ranked_names)])
    stop(
      context,
      " contains duplicated feature names: ",
      paste(utils::head(duplicated_names, 5), collapse = ", "),
      ". Feature names define the preranked GSEA universe and must be unique.",
      call. = FALSE
    )
  }

  if (anyNA(ranked_list) || any(!is.finite(ranked_list))) {
    stop(context, " must contain only finite, non-missing ranking statistics.",
         call. = FALSE)
  }

  if (isTRUE(require_variation) &&
      length(unique(as.numeric(ranked_list))) < 2) {
    stop(
      context,
      " must contain at least two distinct ranking statistics. ",
      "All features are tied, so preranked GSEA would be driven by arbitrary ",
      "input order rather than biological signal.",
      call. = FALSE
    )
  }

  names(ranked_list) <- ranked_names
  ranked_list
}

#' Filter gene sets to a preranked feature universe
#'
#' @param gene_sets Named list of gene sets.
#' @param universe Character vector of ranked feature identifiers.
#' @param min_size Minimum post-overlap gene set size.
#' @param max_size Maximum post-overlap gene set size.
#'
#' @return Named list of gene sets after universe intersection and size filtering.
#' @keywords internal
filter_gene_sets_to_ranked_universe <- function(gene_sets,
                                                universe,
                                                min_size,
                                                max_size) {
  universe <- as.character(universe)
  filtered <- lapply(gene_sets, function(genes) {
    unique(as.character(genes)[as.character(genes) %in% universe])
  })
  sizes <- lengths(filtered)
  filtered[sizes >= min_size & sizes <= max_size]
}

#' Calculate rank metric for GSEA
#'
#' @param abundance A matrix of abundance data
#' @param metadata A data frame of metadata
#' @param group A character string specifying the grouping variable
#' @param method A character string specifying the ranking method
#' @param comparison Optional length-2 character vector \code{c(group1, group2)}
#'   defining the ranking direction for preranked GSEA. Positive values indicate
#'   higher abundance in \code{group1}. If \code{NULL}, exactly two aligned
#'   group levels must be present and factor-level order is used.
#'
#' @return A named vector of ranking statistics
#' @details The abundance matrix must carry explicit, unique feature row
#'   names. The returned ranking vector is validated for unique feature names,
#'   finite numeric statistics, and at least two distinct values; a fully tied
#'   ranking is not meaningful for preranked GSEA because enrichment would then
#'   depend on arbitrary input order.
#' @keywords internal
calculate_rank_metric <- function(abundance, 
                                 metadata, 
                                 group, 
                                 method = "signal2noise",
                                 comparison = NULL) {
  
  # Ensure abundance is a matrix with samples as columns
  abundance <- normalize_abundance_feature_ids(
    abundance,
    context = "calculate_rank_metric() abundance"
  )
  abundance <- as.matrix(abundance)

  aligned <- align_samples(abundance, metadata, verbose = FALSE)
  abundance <- as.matrix(aligned$abundance)
  metadata <- aligned$metadata
  validate_feature_rownames(abundance, "calculate_rank_metric() abundance")

  # Extract group information after alignment so sample names come from the
  # authoritative abundance columns, not from metadata row names.
  if (is.null(comparison)) {
    validate_gsea_group_after_alignment(metadata, group)
  }
  Group <- factor(metadata[[group]])
  names(Group) <- colnames(abundance)

  validate_nonnegative_finite_matrix(abundance, "abundance")

  valid_rank_methods <- c("signal2noise", "t_test", "log2_ratio", "diff_abundance")
  validate_choice(method, valid_rank_methods, "method")

  group_levels <- resolve_preranked_comparison(metadata, group, comparison)

  keep_samples <- !is.na(Group) & Group %in% group_levels
  abundance <- abundance[, keep_samples, drop = FALSE]
  Group <- factor(Group[keep_samples], levels = group_levels)
  names(Group) <- colnames(abundance)

  # Split samples by group
  group1_samples <- names(Group)[Group == group_levels[1]]
  group2_samples <- names(Group)[Group == group_levels[2]]
  
  # Calculate ranking statistic based on method
  if (method == "signal2noise") {
    # Signal-to-noise ratio: (mean(group1) - mean(group2)) / (sd1 + sd2)
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
    
    for (i in seq_len(nrow(abundance))) {
      metric[i] <- tryCatch({
        t_test <- stats::t.test(abundance[i, group1_samples], abundance[i, group2_samples])
        as.numeric(t_test$statistic)
      }, error = function(e) {
        # Constant rows carry no ranking signal. Keep them in the ranked
        # vector as neutral rather than aborting the whole preranked GSEA run.
        if (grepl("constant|not enough", e$message, ignore.case = TRUE)) {
          0
        } else {
          stop(e)
        }
      })
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
  
  metric <- validate_ranked_list(
    metric,
    "calculate_rank_metric() ranked metric"
  )
  attr(metric, "comparison") <- group_levels
  return(metric)
}

#' Run fgsea using the recommended fgseaMultilevel method
#'
#' @param ranked_list A named vector of ranking statistics. Names define the
#'   feature universe and must be non-empty and unique; values must be finite
#'   and contain at least two distinct statistics.
#' @param gene_sets A list of pathway gene sets
#' @param min_size An integer specifying the minimum gene set size
#' @param max_size An integer specifying the maximum gene set size
#' @param p_adjust_method P-value adjustment method applied to fgsea raw p-values.
#'
#' @return A data frame of fgsea results
#' @keywords internal
run_fgsea <- function(ranked_list,
                     gene_sets,
                     min_size = 5,
                     max_size = 500,
                     p_adjust_method = "BH") {
  validate_p_adjust_method(p_adjust_method)
  gene_sets <- validate_gene_sets(gene_sets, "run_fgsea() gene_sets")
  ranked_list <- validate_ranked_list(ranked_list, "run_fgsea() ranked_list")

  # Use fgseaMultilevel (default when nperm is NOT passed).
  # fgseaMultilevel uses adaptive multilevel splitting Monte Carlo,
  # which provides accurate p-values without a fixed permutation count.
  # Passing nperm would force the deprecated fgseaSimple method.
  fgsea_result <- fgsea::fgsea(
    pathways = gene_sets,
    stats = ranked_list,
    minSize = min_size,
    maxSize = max_size
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
      p.adjust = stats::p.adjust(results$pval, method = p_adjust_method),
      leading_edge = sapply(results$leadingEdge, function(x) paste(x, collapse = ";")),
      stringsAsFactors = FALSE
    )
  } else {
    results <- create_empty_gsea_result("fgsea", full = TRUE)
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
#' @param contrast Contrast specification for camera/fry. See
#'   \code{\link{pathway_gsea}} for the public contract.
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
  # Samples are already aligned by the caller (pathway_gsea)
  abundance_mat <- as.matrix(abundance_mat)
  gene_sets <- validate_gene_sets(gene_sets, "run_limma_gsea() gene_sets")

  # Build design matrix
  design <- build_design_matrix(metadata, group, covariates)

  contrast_coef <- resolve_limma_contrast(design, metadata, group, contrast)

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
  # voom/limma do not accept NA, Inf, or negative values in count data.
  # Reject these inputs instead of coercing them to zero because such coercion
  # changes both the ranking statistics and the mean-variance relationship.
  validate_nonnegative_finite_matrix(abundance_mat, "abundance_mat")

  # Use limma-voom for count data transformation. This estimates the
  # mean-variance relationship and computes precision weights.
  # Pass raw counts to voom. limma's voom() applies its own 0.5 offset when
  # computing logCPM values; adding a pseudocount here would change library
  # sizes and distort the mean-variance trend.
  v <- tryCatch(
    limma::voom(abundance_mat, design, plot = FALSE),
    error = function(e) {
      stop(
        "limma::voom() failed, so ", method,
        " gene-set testing was not run. A log2 transform with unit weights ",
        "is not a statistically equivalent fallback because it discards ",
        "voom's observation-level mean-variance weights. Original error: ",
        conditionMessage(e),
        call. = FALSE
      )
    }
  )

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

    # Camera/fry do not produce a normalized enrichment score. Keep the legacy
    # NES column for visualization compatibility, but explicitly label it as a
    # signed p-value score so downstream plots do not present it as true NES.
    formatted_results$signed_log10_pvalue <- ifelse(
      formatted_results$direction == "Up",
      -log10(pmax(formatted_results$pvalue, 1e-300)),  # Positive for Up
      log10(pmax(formatted_results$pvalue, 1e-300))     # Negative for Down
    )
    formatted_results$NES <- formatted_results$signed_log10_pvalue
    formatted_results$score_type <- "signed_log10_pvalue"
    formatted_results$score_label <- "Signed -log10(p-value)"

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

  validate_complete_design_variables(metadata, group, covariates)

  # Ensure group is a factor
  metadata[[group]] <- factor(metadata[[group]])

  # Build formula from literal column names. This supports valid data-frame
  # columns such as "treatment group" without treating spaces or punctuation
  # as formula syntax.
  if (is.null(covariates) || length(covariates) == 0) {
    formula_obj <- build_one_sided_formula(group)
  } else {
    for (cov in covariates) {
      if (is.character(metadata[[cov]])) {
        metadata[[cov]] <- factor(metadata[[cov]])
      }
    }
    formula_obj <- build_one_sided_formula(c(group, covariates))
  }

  # Create design matrix
  design <- stats::model.matrix(formula_obj, data = metadata)
  # Clean up column names for readability before diagnostics are reported.
  colnames(design) <- gsub("\\(Intercept\\)", "Intercept", colnames(design))

  if (nrow(design) != nrow(metadata)) {
    stop(
      "Design matrix row count (", nrow(design),
      ") does not match metadata row count (", nrow(metadata),
      "). Missing design variables may have been dropped by model.matrix().",
      call. = FALSE
    )
  }

  if (any(!is.finite(design))) {
    stop(
      "Design matrix contains non-finite values. Check numeric covariates for ",
      "Inf, -Inf, or NaN before running camera/fry GSEA.",
      call. = FALSE
    )
  }

  qr_design <- qr(design)
  if (qr_design$rank < ncol(design)) {
    non_estimable <- colnames(design)[
      qr_design$pivot[seq.int(qr_design$rank + 1L, ncol(design))]
    ]
    stop(
      "Design matrix is rank deficient; camera/fry contrasts are not fully ",
      "estimable. This usually means a covariate is constant or perfectly ",
      "confounded with the group variable. Non-estimable column(s): ",
      paste(non_estimable, collapse = ", "),
      ".",
      call. = FALSE
    )
  }

  return(design)
}

#' Resolve a limma contrast for camera/fry
#'
#' @param design A design matrix from \code{model.matrix}
#' @param metadata A data frame containing sample metadata
#' @param group Character scalar naming the grouping variable
#' @param contrast NULL, a numeric contrast, or a character contrast. Named
#' numeric contrast vectors are aligned to \code{colnames(design)} before
#' being returned; unnamed vectors must already be in design column order.
#'
#' @return A column index, column name, or numeric contrast vector accepted by limma
#' @keywords internal
resolve_limma_contrast <- function(design, metadata, group, contrast = NULL) {
  design_cols <- colnames(design)
  group_levels <- levels(factor(metadata[[group]]))

  assign <- attr(design, "assign")
  if (is.null(assign)) {
    stop("Design matrix is missing term assignment metadata.", call. = FALSE)
  }

  terms_obj <- attr(design, "terms")
  if (!is.null(terms_obj)) {
    term_labels <- attr(terms_obj, "term.labels")
    group_term_index <- match(group, term_labels)
    if (is.na(group_term_index)) {
      group_term_index <- match(deparse(as.name(group)), term_labels)
    }
  } else {
    group_term_index <- 1L
  }
  if (is.na(group_term_index)) {
    group_cols <- grep(paste0("^", group), design_cols)
  } else {
    group_cols <- which(assign == group_term_index)
  }

  if (length(group_cols) == 0) {
    stop("Could not identify group coefficients in the design matrix. ",
         "Design columns: ", paste(design_cols, collapse = ", "),
         call. = FALSE)
  }

  if (is.null(contrast)) {
    if (length(group_cols) == 1) {
      return(group_cols)
    }
    stop("Multiple group coefficients are present (",
         paste(design_cols[group_cols], collapse = ", "),
         "). Specify 'contrast' explicitly for multi-group camera/fry analysis.",
         call. = FALSE)
  }

  if (is.numeric(contrast)) {
    if (length(contrast) == 1 && is.null(names(contrast))) {
      if (!is.finite(contrast) || contrast != as.integer(contrast) ||
          contrast < 1 || contrast > ncol(design)) {
        stop("Numeric contrast column index must be an integer between 1 and ",
             ncol(design), ".", call. = FALSE)
      }
      return(as.integer(contrast))
    }
    if (!is.null(names(contrast))) {
      contrast_names <- names(contrast)
      if (length(contrast_names) != length(contrast) ||
          any(is.na(contrast_names)) ||
          any(!nzchar(contrast_names))) {
        stop("Named numeric contrast vectors must have one non-empty name per value.",
             call. = FALSE)
      }
      if (anyDuplicated(contrast_names)) {
        duplicated_names <- unique(contrast_names[duplicated(contrast_names)])
        stop("Named numeric contrast vector contains duplicated coefficient name(s): ",
             paste(utils::head(duplicated_names, 5), collapse = ", "),
             call. = FALSE)
      }
      missing_names <- setdiff(design_cols, contrast_names)
      extra_names <- setdiff(contrast_names, design_cols)
      if (length(missing_names) > 0 || length(extra_names) > 0) {
        stop(
          "Named numeric contrast vector names must match the design columns exactly. ",
          if (length(missing_names) > 0) {
            paste0("Missing: ", paste(missing_names, collapse = ", "), ". ")
          } else "",
          if (length(extra_names) > 0) {
            paste0("Unexpected: ", paste(extra_names, collapse = ", "), ". ")
          } else "",
          "Design columns: ", paste(design_cols, collapse = ", "),
          call. = FALSE
        )
      }
      contrast <- contrast[design_cols]
    }
    if (length(contrast) != ncol(design)) {
      stop("Numeric contrast vector must have length equal to the number of ",
           "design columns (", ncol(design), ").", call. = FALSE)
    }
    if (any(!is.finite(contrast))) {
      stop("Numeric contrast vector must contain only finite values.",
           call. = FALSE)
    }
    return(contrast)
  }

  if (!is.character(contrast) || length(contrast) != 1) {
    stop("'contrast' must be NULL, a numeric column index, a numeric contrast ",
         "vector, or a single character string naming an exact design column ",
         "or non-reference group level.",
         call. = FALSE)
  }

  if (contrast %in% design_cols) {
    return(match(contrast, design_cols))
  }

  ref_level <- group_levels[1]
  non_ref_levels <- group_levels[-1]
  if (contrast == ref_level) {
    stop("Contrast '", contrast, "' is the reference level. Choose one of the ",
         "non-reference group levels (", paste(non_ref_levels, collapse = ", "),
         ") or supply a numeric contrast vector.",
         call. = FALSE)
  }
  if (contrast %in% non_ref_levels) {
    level_index <- match(contrast, non_ref_levels)
    if (level_index <= length(group_cols)) {
      return(group_cols[level_index])
    }
  }

  stop("Contrast '", contrast, "' was not found as an exact design column or ",
       "non-reference group level. Available design columns: ",
       paste(design_cols, collapse = ", "),
       ". Non-reference group levels: ", paste(non_ref_levels, collapse = ", "),
       ".",
       call. = FALSE)
}
