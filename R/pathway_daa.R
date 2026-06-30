#' @importFrom magrittr %>%
#' @importFrom stats wilcox.test setNames
NULL

#' Differential Abundance Analysis for Predicted Functional Pathways
#'
#' @name pathway_daa
#'
#' @description
#' Performs differential abundance analysis on predicted functional pathway data using various statistical methods.
#' This function supports multiple methods for analyzing differences in pathway abundance between groups,
#' including popular approaches like ALDEx2, DESeq2, edgeR, and others.
#'
#' @param abundance A data frame or matrix containing predicted functional pathway abundance,
#'        with pathways/features as rows and samples as columns.
#'        Data frames may also provide a leading non-numeric feature ID column
#'        (for example \code{#NAME}, \code{feature}, or \code{pathway}); it is
#'        converted to row names before sample alignment.
#'        Feature identifiers must be explicit, non-empty, and unique.
#'        The column names should match the sample names in metadata.
#'        Values should be finite, non-missing, non-negative counts or abundance
#'        measurements. Count-based backends that require or assume integer
#'        counts (ALDEx2, DESeq2, edgeR, and metagenomeSeq) round non-integer
#'        values with a warning before fitting.
#'
#' @param metadata A data frame or tibble containing sample information.
#'        Must include a 'sample' column with sample identifiers matching the column names in abundance data.
#'
#' @param group Character string specifying the column name in metadata that contains group information
#'        for differential abundance analysis. Group values must be non-missing
#'        and non-empty for all aligned/selected samples.
#'
#' @param daa_method Character string specifying the method for differential abundance analysis.
#'        Available choices are:
#'        \itemize{
#'          \item \code{"ALDEx2"}: ANOVA-Like Differential Expression tool
#'          \item \code{"DESeq2"}: Differential expression analysis based on negative binomial distribution
#'          \item \code{"edgeR"}: Exact test for differences between groups using negative binomial model
#'          \item \code{"limma voom"}: Limma-voom framework for RNA-seq analysis
#'          \item \code{"metagenomeSeq"}: Zero-inflated Gaussian mixture model
#'          \item \code{"LinDA"}: Linear models for differential abundance analysis
#'          \item \code{"Maaslin2"}: Multivariate Association with Linear Models
#'          \item \code{"Lefser"}: Linear discriminant analysis effect size
#'        }
#'        Default is "ALDEx2".
#'
#' @param select Character vector of unique sample names to include in the
#'        analysis. If NULL (default), all samples are included. The selected
#'        dataset must still contain at least four samples, at least two
#'        groups, and at least two samples per group.
#'
#' @param p_adjust_method Character string specifying the method for p-value adjustment.
#'        Choices are:
#'        \itemize{
#'          \item \code{"BH"}: Benjamini-Hochberg procedure (default)
#'          \item \code{"holm"}: Holm's step-down method
#'          \item \code{"bonferroni"}: Bonferroni correction
#'          \item \code{"hochberg"}: Hochberg's step-up method
#'          \item \code{"fdr"}: False Discovery Rate
#'          \item \code{"none"}: No adjustment
#'        }
#'
#' @param reference Character string specifying the reference level for the group comparison.
#'        If NULL (default), the first level is used as reference. When supplied,
#'        it must exactly match one observed group level after sample alignment
#'        and any \code{select} filtering.
#'
#' @param include_abundance_stats Logical value indicating whether to include
#'        abundance statistics (mean relative abundance and standard deviation
#'        per group) in the output. Default is FALSE. When the selected
#'        \code{daa_method} already provides a \code{log2_fold_change} column
#'        (ALDEx2 with effect size, DESeq2, edgeR, limma voom, LinDA, Maaslin2,
#'        metagenomeSeq), the method-native log2 fold change is preserved and
#'        the relative-abundance ratio is not recomputed.
#'
#' @param include_effect_size Logical value indicating whether to compute
#'        ALDEx2 effect size information via \code{ALDEx2::aldex.effect()}.
#'        When TRUE, adds \code{effect_size}, \code{diff_btw},
#'        \code{log2_fold_change}, \code{rab_all}, and \code{overlap} columns,
#'        aligning ALDEx2 output with the other DAA methods that return log2
#'        fold changes by default. Only applicable for two-group comparisons
#'        with the ALDEx2 method; ignored otherwise. Default is TRUE; set to
#'        FALSE to skip the extra \code{aldex.effect()} computation. For a
#'        two-group analysis, failure to compute or validate the requested
#'        effect-size output stops the analysis.
#'
#' @param p.adjust Deprecated. Use \code{p_adjust_method} instead.
#' @param .pre_aligned Internal logical. Set to TRUE only when the caller has
#'        already aligned abundance columns and metadata rows in identical
#'        sample order.
#' @param .sample_col Internal character. Sample identifier column used when
#'        \code{.pre_aligned = TRUE}.
#'
#' @param ... Reserved for future backend-specific parameters. Additional
#'        arguments are currently rejected rather than silently ignored, because
#'        ignored model/covariate arguments can make the fitted analysis differ
#'        from the analysis the user intended.
#'
#' @return A data frame containing the differential abundance analysis results. The structure of the results
#' depends on the chosen DAA method. For methods that support multi-group comparisons (like LinDA),
#' when there are more than two groups, the results will contain separate rows for each feature in each
#' pairwise comparison between the reference group and each non-reference group. The data frame includes
#' the following columns:
#' \itemize{
#'   \item \code{feature}: Feature/pathway identifier
#'   \item \code{method}: The DAA method used
#'   \item \code{group1}: Reference group
#'   \item \code{group2}: Comparison group
#'   \item \code{p_values}: P-values for the comparison
#'   \item \code{p_adjust}: Adjusted p-values
#'   \item \code{adj_method}: Method used for p-value adjustment
#' }
#' Method-native adjusted p-values are preserved when the backend provides
#' them directly (ALDEx2 \code{eBH}, DESeq2 \code{padj}, LinDA \code{padj},
#' and Maaslin2 \code{qval}). Other methods are adjusted by
#' \code{pathway_daa()} using \code{stats::p.adjust()} and
#' \code{p_adjust_method}. For wrapper-computed adjustments, p-values are
#' adjusted within each method and pairwise comparison when \code{method},
#' \code{group1}, and \code{group2} columns are available.
#'
#' Methods that fit a model on the abundance data (DESeq2, edgeR, limma voom,
#' LinDA, Maaslin2, metagenomeSeq) return a \code{log2_fold_change} column
#' computed in the method's own model space. ALDEx2 returns
#' \code{log2_fold_change} (plus \code{effect_size}, \code{diff_btw},
#' \code{rab_all}, \code{overlap}) when \code{include_effect_size = TRUE}
#' (the default), derived from \code{ALDEx2::aldex.effect()} in CLR space.
#' Lefser returns an \code{lda_score} column instead, which is its native
#' effect-size metric.
#'
#' When \code{include_abundance_stats = TRUE}, the following additional columns
#' are included:
#' \itemize{
#'   \item \code{mean_rel_abundance_group1}: Mean relative abundance for group1
#'   \item \code{sd_rel_abundance_group1}: Standard deviation of relative
#'         abundance for group1
#'   \item \code{mean_rel_abundance_group2}: Mean relative abundance for group2
#'   \item \code{sd_rel_abundance_group2}: Standard deviation of relative
#'         abundance for group2
#' }
#' A \code{log2_fold_change} column from relative abundance is only added when
#' the DAA method does not already provide one, to avoid conflating model-based
#' and ratio-based effect sizes. If the requested abundance statistics cannot
#' be calculated for every returned feature/group pair, the function fails
#' instead of returning a partially annotated result table.
#'
#' @examples
#' \donttest{
#' # Load example data
#' data(ko_abundance)
#' data(metadata)
#'
#' # Prepare abundance data
#' abundance_data <- as.data.frame(ko_abundance)
#' rownames(abundance_data) <- abundance_data[, "#NAME"]
#' abundance_data <- abundance_data[, -1]
#'
#' # Run differential abundance analysis using ALDEx2
#' results <- pathway_daa(
#'   abundance = abundance_data,
#'   metadata = metadata,
#'   group = "Environment"
#' )
#'
#' # Using a different method (DESeq2)
#' deseq_results <- pathway_daa(
#'   abundance = abundance_data,
#'   metadata = metadata,
#'   group = "Environment",
#'   daa_method = "DESeq2"
#' )
#'
#' # Create example data with more samples
#' abundance <- data.frame(
#'   sample1 = c(10, 20, 30),
#'   sample2 = c(20, 30, 40),
#'   sample3 = c(30, 40, 50),
#'   sample4 = c(40, 50, 60),
#'   sample5 = c(50, 60, 70),
#'   row.names = c("pathway1", "pathway2", "pathway3")
#' )
#'
#' metadata <- data.frame(
#'   sample = c("sample1", "sample2", "sample3", "sample4", "sample5"),
#'   group = c("control", "control", "treatment", "treatment", "treatment")
#' )
#'
#' # Run differential abundance analysis using ALDEx2
#' results <- pathway_daa(abundance, metadata, "group")
#'
#' # Using a different method (limma voom instead of DESeq2 for this small example)
#' limma_results <- pathway_daa(abundance, metadata, "group",
#'                             daa_method = "limma voom")
#'
#' # Analyze specific samples only
#' subset_results <- pathway_daa(abundance, metadata, "group",
#'                              select = c("sample1", "sample2", "sample3", "sample4"))
#'
#' # ALDEx2 returns effect size columns by default
#' # (effect_size, diff_btw, log2_fold_change, rab_all, overlap).
#' # Ranking by |log2_fold_change| is generally more biologically informative
#' # than ranking by p-value, especially for large datasets where small effects
#' # can reach statistical significance without being biologically meaningful.
#' aldex2_res <- pathway_daa(abundance, metadata, "group", daa_method = "ALDEx2")
#' head(aldex2_res)
#'
#' # Opt out of the extra aldex.effect() computation if only p-values are needed
#' aldex2_pvals_only <- pathway_daa(abundance, metadata, "group",
#'                                 daa_method = "ALDEx2",
#'                                 include_effect_size = FALSE)
#' }
#'
#' @references
#' \itemize{
#'   \item ALDEx2: Fernandes et al. (2014) Unifying the analysis of high-throughput sequencing datasets:
#'         characterizing RNA-seq, 16S rRNA gene sequencing and selective growth experiments by
#'         compositional data analysis. Microbiome.
#'   \item DESeq2: Love et al. (2014) Moderated estimation of fold change and dispersion for
#'         RNA-seq data with DESeq2. Genome Biology.
#'   \item edgeR: Robinson et al. (2010) edgeR: a Bioconductor package for differential expression
#'         analysis of digital gene expression data. Bioinformatics.
#'   \item limma-voom: Law et al. (2014) voom: precision weights unlock linear model analysis tools
#'         for RNA-seq read counts. Genome Biology.
#'   \item metagenomeSeq: Paulson et al. (2013) Differential abundance analysis for microbial
#'         marker-gene surveys. Nature Methods.
#'   \item Maaslin2: Mallick et al. (2021) Multivariable Association Discovery in Population-scale
#'         Meta-omics Studies.
#' }
NULL

#' Helper function to calculate abundance statistics for differential analysis
#'
#' This function calculates mean relative abundance, standard deviation, and log2 fold change
#' for each feature between two groups.
#'
#' @param abundance A matrix or data frame with features as rows and samples as columns
#' @param metadata A data frame containing sample information
#' @param group Character string specifying the group column name in metadata
#' @param features Character vector of feature names to calculate statistics for
#' @param group1 Character string specifying the first group name
#' @param group2 Character string specifying the second group name
#'
#' @return A data frame with columns:
#'   \item{feature}{Feature identifier}
#'   \item{mean_rel_abundance_group1}{Mean relative abundance for group1}
#'   \item{sd_rel_abundance_group1}{Standard deviation of relative abundance for group1}
#'   \item{mean_rel_abundance_group2}{Mean relative abundance for group2}
#'   \item{sd_rel_abundance_group2}{Standard deviation of relative abundance for group2}
#'   \item{log2_fold_change}{Log2 fold change (group2/group1)}
#'
#' @keywords internal
calculate_abundance_stats <- function(abundance, metadata, group, features, group1, group2) {
  # Validate inputs using unified functions
  abundance <- normalize_abundance_feature_ids(
    abundance,
    context = "calculate_abundance_stats() abundance"
  )
  validate_abundance(abundance)
  validate_metadata(metadata)
  validate_group(metadata, group)

  # Convert abundance to matrix if needed
  abundance_mat <- as.matrix(abundance)
  validate_feature_rownames(abundance_mat, "calculate_abundance_stats() abundance")

  # Align samples between abundance and metadata
  aligned <- align_samples(abundance_mat, metadata, verbose = FALSE)
  abundance_filtered <- aligned$abundance
  metadata_ordered <- aligned$metadata
  sample_col <- aligned$sample_col

  # Get group assignments
  group_assignments <- metadata_ordered[[group]]

  # Filter for the two groups of interest
  group1_samples <- colnames(abundance_filtered)[group_assignments == group1]
  group2_samples <- colnames(abundance_filtered)[group_assignments == group2]

  if (length(group1_samples) == 0) {
    stop("No samples found for group1: ", group1)
  }
  if (length(group2_samples) == 0) {
    stop("No samples found for group2: ", group2)
  }

  # Convert to relative abundance via the shared helper so zero-sum sample
  # columns fail fast here instead of producing NaN that later gets
  # dropped by `na.rm = TRUE` inside the mean()/sd() aggregations,
  # silently changing the effective sample size.
  relative_abundance <- compute_relative_abundance(
    abundance_filtered,
    context = "calculate_abundance_stats()"
  )

  # Filter for specified features
  if (!is.null(features)) {
    feature_ids <- validate_nonempty_character_column(
      features,
      "features",
      "calculate_abundance_stats()"
    )
    if (length(feature_ids) == 0) {
      stop("'features' must contain at least one feature identifier.",
           call. = FALSE)
    }
    # Multiple method rows can legitimately request the same feature for the
    # same group pair (for example ALDEx2 Welch and Wilcoxon outputs). The
    # abundance summary is feature-level, so calculate it once per feature
    # and let the caller merge it back to every DAA row that needs it.
    feature_ids <- unique(feature_ids)
    missing_features <- setdiff(feature_ids, rownames(relative_abundance))
    if (length(missing_features) > 0) {
      stop(
        "Feature(s) requested in calculate_abundance_stats() are missing ",
        "from abundance row names: ",
        paste(utils::head(missing_features, 5), collapse = ", "),
        if (length(missing_features) > 5) ", ..." else "",
        ". Use matching feature identifiers in abundance and DAA results.",
        call. = FALSE
      )
    }
    relative_abundance <- relative_abundance[feature_ids, , drop = FALSE]
  }

  # Per-feature mean and sd for the two groups of interest. Restrict the
  # matrix to group1 + group2 columns before aggregation so the helper does
  # not spend time on samples that wouldn't show up in the output anyway.
  # This is also the only entry point that computes the mean/sd inputs to
  # the log2 fold change: `pathway_errorbar()` reuses the same helper, so
  # any future change to the aggregation rule (e.g. na.rm strategy,
  # trimmed means, etc.) only needs to happen in one place.
  keep <- group_assignments %in% c(group1, group2)
  long_stats <- summarize_abundance_by_group(
    relative_abundance[, keep, drop = FALSE],
    group_assignments[keep]
  )

  feature_order <- rownames(relative_abundance)
  g1 <- long_stats[long_stats$group == group1, , drop = FALSE]
  g2 <- long_stats[long_stats$group == group2, , drop = FALSE]
  g1 <- g1[match(feature_order, g1$name), , drop = FALSE]
  g2 <- g2[match(feature_order, g2$name), , drop = FALSE]

  # Data-driven pseudocount must be computed from the full relative-
  # abundance matrix (same scope as before the refactor), not from the
  # per-group means -- otherwise the pseudocount would track fold-change
  # scale and inflate artificially for small-sample comparisons.
  pseudocount <- calculate_pseudocount(as.vector(relative_abundance))

  results <- data.frame(
    feature = feature_order,
    mean_rel_abundance_group1 = g1$mean,
    sd_rel_abundance_group1 = g1$sd,
    mean_rel_abundance_group2 = g2$mean,
    sd_rel_abundance_group2 = g2$sd,
    log2_fold_change = calculate_log2_fold_change(
      g1$mean, g2$mean, pseudocount = pseudocount
    ),
    stringsAsFactors = FALSE
  )

  return(results)
}

adjust_daa_p_values <- function(result, p_adjust_method) {
  if (is.null(result) || !"p_values" %in% colnames(result) || nrow(result) == 0) {
    return(result)
  }

  validate_probability_values(result$p_values, "p_values",
                              "DAA result", allow_na = TRUE)

  if ("p_adjust" %in% colnames(result)) {
    validate_probability_values(result$p_adjust, "p_adjust",
                                "DAA result", allow_na = TRUE)
    method_adj <- attr(result, "adj_method", exact = TRUE)
    if (is.null(method_adj) || length(method_adj) != 1 ||
        is.na(method_adj) || !nzchar(method_adj)) {
      method_adj <- "BH (method-specific)"
    }
    result$adj_method <- method_adj
    return(result)
  }

  result$p_adjust <- NA_real_
  adjust_by <- c("method", "group1", "group2")

  if (all(adjust_by %in% colnames(result))) {
    split_key <- do.call(
      interaction,
      c(result[adjust_by], list(drop = TRUE, lex.order = TRUE))
    )
    for (idx in split(seq_len(nrow(result)), split_key)) {
      result$p_adjust[idx] <- stats::p.adjust(
        result$p_values[idx],
        method = p_adjust_method
      )
    }
  } else {
    result$p_adjust <- stats::p.adjust(result$p_values, method = p_adjust_method)
  }

  validate_probability_values(result$p_adjust, "p_adjust",
                              "DAA result", allow_na = TRUE)
  result$adj_method <- p_adjust_method
  result
}

#' @rdname pathway_daa
#' @export
pathway_daa <- function(abundance, metadata, group, daa_method = "ALDEx2",
                       select = NULL, p_adjust_method = "BH", reference = NULL,
                       include_abundance_stats = FALSE, include_effect_size = TRUE,
                       p.adjust = NULL, .pre_aligned = FALSE,
                       .sample_col = NULL, ...) {
  # Backward compatibility for deprecated parameter
  if (!is.null(p.adjust)) {
    warning("'p.adjust' parameter is deprecated. Use 'p_adjust_method' instead.", call. = FALSE)
    p_adjust_method <- p.adjust
  }
  validate_p_adjust_method(p_adjust_method)

  # Single source of truth for supported DAA methods: this list drives
  # validation, the method->package mapping for optional-dependency
  # loading, and documents what the switch() below will accept.
  method_packages <- list(
    "ALDEx2" = "ALDEx2", "DESeq2" = "DESeq2", "edgeR" = "edgeR",
    "limma voom" = "limma", "metagenomeSeq" = "metagenomeSeq",
    "LinDA" = "MicrobiomeStat",
    "Maaslin2" = "Maaslin2", "Lefser" = "lefser"
  )

  # Validate daa_method up front. switch() below has no default, so an
  # unrecognized value silently yielded NULL; worse, both README and
  # some internal helpers historically referenced the misspellings
  # "linDA" and "Lefse", so users copy-pasted invalid values and got
  # a confusing NULL back. Fail fast with a message that shows the
  # canonical set.
  if (length(daa_method) != 1 || !is.character(daa_method) ||
      is.na(daa_method) || !nzchar(daa_method)) {
    stop("'daa_method' must be a single non-empty character string. ",
         "Supported methods: ",
         paste(names(method_packages), collapse = ", "), ".")
  }
  if (!daa_method %in% names(method_packages)) {
    # Common misspellings -> canonical form. Users frequently hit
    # "linDA" (from older docs) or "Lefse" (the LEfSe upstream name,
    # but our backend is the Bioconductor "Lefser" package); point
    # them at the right spelling instead of generic "not supported".
    suggestion <- switch(
      tolower(daa_method),
      "linda" = "LinDA",
      "lefse" = "Lefser",
      "aldex" = "ALDEx2",
      "aldex2" = "ALDEx2",
      "maaslin"  = "Maaslin2",
      "maaslin2" = "Maaslin2",
      "deseq" = "DESeq2",
      "deseq2" = "DESeq2",
      NULL
    )
    stop(
      sprintf("Unsupported daa_method: '%s'. ", daa_method),
      if (!is.null(suggestion)) sprintf("Did you mean '%s'? ", suggestion) else "",
      "Supported methods: ",
      paste(names(method_packages), collapse = ", "), "."
    )
  }

  extra_args <- list(...)
  if (length(extra_args) > 0) {
    extra_names <- names(extra_args)
    if (is.null(extra_names)) {
      extra_names <- rep("", length(extra_args))
    }
    unnamed <- !nzchar(extra_names)
    extra_names[unnamed] <- paste0("..", which(unnamed))
    stop(
      "pathway_daa() does not currently pass arguments in `...` to backend ",
      "DAA methods. Unsupported argument(s): ",
      paste(extra_names, collapse = ", "),
      ". Remove these arguments, or call the backend package directly for ",
      "covariate-adjusted or custom model analyses.",
      call. = FALSE
    )
  }

  require_package(method_packages[[daa_method]], purpose = daa_method)

  # Input validation using unified functions
  abundance <- normalize_abundance_feature_ids(
    abundance,
    context = "pathway_daa() abundance"
  )
  validate_abundance(abundance, min_samples = 4)
  validate_metadata(metadata)
  validate_group(metadata, group, min_groups = 2)

  # Alignment ownership is explicit. Standalone callers use the public
  # contract and get automatic sample matching here. ggpicrust2() aligns
  # once at the wrapper boundary because it also needs the same aligned
  # objects for plotting, then calls this function with `.pre_aligned = TRUE`
  # so we validate the invariant instead of repeating the alignment.
  .pre_aligned <- normalize_logical_flag(.pre_aligned, ".pre_aligned")
  if (.pre_aligned) {
    if (is.null(colnames(abundance))) {
      stop("Pre-aligned abundance data must have sample column names.", call. = FALSE)
    }
    metadata <- as.data.frame(metadata)
    sample_col <- .sample_col
    if (is.null(sample_col) || !sample_col %in% colnames(metadata)) {
      sample_col <- find_sample_column(metadata, colnames(abundance))
    }
    if (is.null(sample_col)) {
      stop("Pre-aligned metadata must contain sample identifiers matching abundance columns.",
           call. = FALSE)
    }
    if (sample_col == ".rownames") {
      metadata$.sample_id <- rownames(metadata)
      sample_col <- ".sample_id"
    }
    if (!identical(colnames(abundance), as.character(metadata[[sample_col]]))) {
      stop("Pre-aligned abundance and metadata are not in the same sample order.",
           call. = FALSE)
    }
    aligned <- list(
      abundance = abundance,
      metadata = metadata,
      sample_col = sample_col,
      n_samples = ncol(abundance)
    )
  } else {
    aligned <- align_samples(abundance, metadata, verbose = FALSE)
  }
  abundance <- aligned$abundance
  metadata <- tibble::as_tibble(aligned$metadata)
  sample_col <- aligned$sample_col

  if (aligned$n_samples < 4) {
    stop("At least 4 samples required for DAA (after alignment)")
  }

  # Handle sample selection
  if (!is.null(select)) {
    if (!is.character(select) || anyNA(select) || any(!nzchar(select))) {
      stop("'select' must be a character vector of non-empty sample names.",
           call. = FALSE)
    }
    if (anyDuplicated(select)) {
      stop("'select' must contain unique sample names.", call. = FALSE)
    }
    if (!all(select %in% colnames(abundance))) {
      stop("Some selected samples not in abundance data")
    }
    abundance <- abundance[, select, drop = FALSE]
    # Keep metadata rows in the exact order of `select` so row i of metadata
    # corresponds to column i of abundance. Simple %in% filtering preserves
    # the original order and can silently desynchronize Group vs samples.
    metadata <- metadata[match(select, metadata[[sample_col]]), , drop = FALSE]

    if (ncol(abundance) < 4) {
      stop("At least 4 samples required for DAA after `select` filtering; found ",
           ncol(abundance), ".", call. = FALSE)
    }
  }

  # Validate abundance matrix
  abundance_mat <- as.matrix(abundance)
  validate_daa_input(abundance_mat, method = daa_method)

  group_values <- as.character(metadata[[group]])
  invalid_group <- is.na(group_values) | !nzchar(trimws(group_values))
  if (any(invalid_group)) {
    bad_samples <- metadata[[sample_col]][invalid_group]
    stop(
      "Group column '", group,
      "' contains missing or empty values after sample alignment",
      if (!is.null(select)) " and `select` filtering" else "",
      ". Affected sample(s): ",
      paste(utils::head(bad_samples, 5), collapse = ", "),
      if (sum(invalid_group) > 5) ", ..." else "",
      ".",
      call. = FALSE
    )
  }

  Group <- factor(metadata[[group]])

  # Ensure factor levels only include groups present in the data
  # This is defensive programming to handle edge cases with unused factor levels
  Group <- droplevels(Group)
  Level <- levels(Group)
  length_Level <- length(Level)

  # Re-validate group counts after align_samples()/select have pruned samples.
  # The up-front validate_group() check sees the raw metadata, so sample
  # alignment or a narrow `select =` that removes every row of a level can
  # leave a single-group dataset that would crash inside the backends with a
  # less actionable error. Also require at least two samples per retained
  # group: the supported DAA backends estimate within-group variation or run
  # tests that are not meaningful for singleton groups.
  if (length_Level < 2) {
    stop(
      "DAA requires at least 2 groups with samples after alignment",
      if (!is.null(select)) " and `select` filtering" else "",
      "; found ", length_Level,
      if (length_Level == 1) paste0(" ('", Level, "')") else "",
      "."
    )
  }
  group_counts <- table(Group)
  small_groups <- names(group_counts)[group_counts < 2]
  if (length(small_groups) > 0) {
    stop(
      "DAA requires at least 2 samples per group after alignment",
      if (!is.null(select)) " and `select` filtering" else "",
      ". Group(s) with fewer than 2 samples: ",
      paste(small_groups, collapse = ", "),
      ".",
      call. = FALSE
    )
  }

  if (!is.null(reference)) {
    if (!is.character(reference) || length(reference) != 1 ||
        is.na(reference) || !nzchar(reference)) {
      stop("'reference' must be NULL or a single non-empty character string.",
           call. = FALSE)
    }
    if (!reference %in% Level) {
      stop(
        "Reference level '", reference, "' was not found in group '", group,
        "' after sample alignment",
        if (!is.null(select)) " and `select` filtering" else "",
        ". Available levels: ", paste(Level, collapse = ", "), ".",
        call. = FALSE
      )
    }
  }

  # Perform differential analysis
  result <- switch(
    daa_method,
    "ALDEx2" = perform_aldex2_analysis(abundance_mat, Group, Level, length_Level, reference, include_effect_size),
    "DESeq2" = perform_deseq2_analysis(abundance_mat, metadata, group, reference, Level, length_Level, p_adjust_method),
    "LinDA" = perform_linda_analysis(abundance, metadata, group, reference, Level, length_Level, p_adjust_method),
    "limma voom" = perform_limma_voom_analysis(abundance_mat, Group, reference, Level, length_Level),
    "edgeR" = perform_edger_analysis(abundance_mat, Group, reference, Level, length_Level),
    "metagenomeSeq" = perform_metagenomeseq_analysis(abundance_mat, metadata, group, reference, Level),
    "Maaslin2" = perform_maaslin2_analysis(abundance_mat, metadata, group, reference, Level, length_Level, p_adjust_method),
    "Lefser" = perform_lefser_analysis(abundance_mat, metadata, group, reference, Level)
  )


  result <- adjust_daa_p_values(result, p_adjust_method)

  # Add abundance statistics if requested
  if (include_abundance_stats && !is.null(result) && nrow(result) > 0) {
    # Check if we have the required columns for abundance stats
    if (all(c("feature", "group1", "group2") %in% colnames(result))) {
      tryCatch({
        # Get unique group pairs for calculation
        unique_comparisons <- unique(result[, c("group1", "group2")])

        # Calculate abundance stats for each unique comparison
        all_abundance_stats <- data.frame()

        for (i in seq_len(nrow(unique_comparisons))) {
          group1_name <- unique_comparisons$group1[i]
          group2_name <- unique_comparisons$group2[i]

          # Get features for this comparison
          comparison_features <- result$feature[
            result$group1 == group1_name & result$group2 == group2_name
          ]

          if (length(comparison_features) > 0) {
            abundance_stats <- calculate_abundance_stats(
              abundance = abundance_mat,
              metadata = metadata,
              group = group,
              features = comparison_features,
              group1 = group1_name,
              group2 = group2_name
            )

            # Add group information for merging
            abundance_stats$group1 <- group1_name
            abundance_stats$group2 <- group2_name

            all_abundance_stats <- rbind(all_abundance_stats, abundance_stats)
          }
        }

        # Merge abundance stats with results. If the DAA method already
        # provides a method-native log2_fold_change (e.g. ALDEx2 effect size
        # in CLR space, DESeq2 shrunk log2FC), keep it and drop the
        # relative-abundance-ratio version to avoid a .x/.y merge collision
        # and to avoid conflating two different effect-size definitions.
        if (nrow(all_abundance_stats) > 0) {
          if ("log2_fold_change" %in% colnames(result) &&
              "log2_fold_change" %in% colnames(all_abundance_stats)) {
            all_abundance_stats$log2_fold_change <- NULL
          }
          result <- merge(
            result,
            all_abundance_stats,
            by = c("feature", "group1", "group2"),
            all.x = TRUE
          )
        }
      }, error = function(e) {
        stop("Failed to calculate abundance statistics: ", e$message,
             call. = FALSE)
      })
    } else {
      stop(
        "Cannot calculate abundance statistics: required columns ",
        "(feature, group1, group2) not found in results.",
        call. = FALSE
      )
    }
  }

  return(result)
}

# Helper function: Perform ALDEx2 analysis
perform_aldex2_analysis <- function(abundance_mat, Group, Level, length_Level,
                                    reference = NULL,
                                    include_effect_size = TRUE) {
  # Filter zero-abundance features (ALDEx2 requirement)
  abundance_mat <- validate_daa_input(abundance_mat, method = "ALDEx2", filter_zero = TRUE)
  abundance_mat <- prepare_integer_counts(abundance_mat, "ALDEx2")

  # Resolve the reference level before converting the factor to ALDEx2's
  # numeric condition vector. aldex.effect() reports `diff.btw` in condition
  # order, so leaving the raw factor order in place made the user-supplied
  # reference argument silently ineffective and could flip the biological
  # interpretation of the reported effect size.
  ref_level <- if (!is.null(reference) && reference %in% Level) reference else Level[1]
  Group <- stats::relevel(factor(Group, levels = Level), ref = ref_level)
  Level <- levels(Group)
  length_Level <- length(Level)

  # Convert Group to numeric for ALDEx2
  Group <- as.numeric(Group)

  # Run analysis based on number of groups
  if (length_Level == 2) {
    ALDEx2_object <- tryCatch(
      ALDEx2::aldex.clr(abundance_mat, Group, mc.samples = 256, denom = "all", verbose = FALSE),
      error = function(e) stop("ALDEx2 CLR failed: ", e$message)
    )

    results <- tryCatch(
      ALDEx2::aldex.ttest(ALDEx2_object, paired.test = FALSE, verbose = FALSE),
      error = function(e) stop("ALDEx2 t-test failed: ", e$message)
    )

    # Get effect size results if requested
    effect_results <- NULL
    if (include_effect_size) {
      effect_results <- tryCatch(
        ALDEx2::aldex.effect(
          ALDEx2_object,
          verbose = FALSE
        ),
        error = function(e) {
          stop(
            "ALDEx2 effect-size calculation failed while ",
            "`include_effect_size = TRUE`: ", conditionMessage(e),
            ". Set `include_effect_size = FALSE` only if p-value-only output ",
            "is intended.",
            call. = FALSE
          )
        }
      )
      effect_results <- validate_and_align_aldex2_effect_results(
        effect_results,
        rownames(results)
      )
    }

    # Build result dataframe using original Level names
    # Use ALDEx2's pre-computed BH-corrected p-values (Monte Carlo-based)
    # These are more accurate than simple p.adjust() as they account for
    # the uncertainty in the CLR transformation
    base_df <- data.frame(
      feature = rep(rownames(results), 2),
      method = c(
        rep("ALDEx2_Welch's t test", nrow(results)),
        rep("ALDEx2_Wilcoxon rank test", nrow(results))
      ),
      group1 = rep(Level[1], 2 * nrow(results)),
      group2 = rep(Level[2], 2 * nrow(results)),
      p_values = c(results$we.ep, results$wi.ep),
      p_adjust = c(results$we.eBH, results$wi.eBH),
      stringsAsFactors = FALSE
    )

    # Replicate effect-size data for the Welch and Wilcoxon result blocks.
    if (!is.null(effect_results)) {
      base_df$effect_size <- rep(effect_results$effect, 2)
      base_df$diff_btw <- rep(effect_results$diff.btw, 2)
      base_df$log2_fold_change <- rep(effect_results$diff.btw, 2)
      base_df$rab_all <- rep(effect_results$rab.all, 2)
      base_df$overlap <- rep(effect_results$overlap, 2)
    }

    return(base_df)
    
  } else {
    # Multi-group analysis. ALDEx2::aldex.effect() only supports two-group
    # comparisons, so include_effect_size is silently ignored here regardless
    # of its value.
    ALDEx2_object <- tryCatch(
      ALDEx2::aldex.clr(abundance_mat, Group, mc.samples = 256, denom = "all", verbose = FALSE),
      error = function(e) stop("ALDEx2 CLR failed: ", e$message)
    )

    results <- ALDEx2::aldex.kw(ALDEx2_object)

    result_df <- data.frame(
      feature = rep(rownames(results), 2),
      method = c(rep(.daa_method_aldex2_kruskal_wallis, nrow(results)),
                 rep("ALDEx2_glm test", nrow(results))),
      p_values = c(results$kw.ep, results$glm.ep),
      p_adjust = c(results$kw.eBH, results$glm.eBH),
      stringsAsFactors = FALSE
    )

    # Add group columns
    group_cols <- matrix(rep(Level, each = nrow(result_df)), nrow = nrow(result_df))
    colnames(group_cols) <- paste0("group", seq_along(Level))
    result_df <- cbind(result_df[, c("feature", "method")], group_cols,
                       result_df[, c("p_values", "p_adjust"), drop = FALSE])

    return(result_df)
  }
}

validate_and_align_aldex2_effect_results <- function(effect_results,
                                                     feature_ids) {
  context <- "ALDEx2 effect-size output"
  if (!is.data.frame(effect_results)) {
    stop(context, " must be a data frame.", call. = FALSE)
  }

  feature_ids <- validate_nonempty_character_column(
    feature_ids,
    "expected feature identifiers",
    context
  )
  output_ids <- validate_nonempty_character_column(
    rownames(effect_results),
    "row names",
    context
  )
  if (anyDuplicated(output_ids)) {
    duplicated_ids <- unique(output_ids[duplicated(output_ids)])
    stop(
      context, " contains duplicated feature identifiers: ",
      paste(utils::head(duplicated_ids, 5), collapse = ", "),
      ".",
      call. = FALSE
    )
  }

  missing_ids <- setdiff(feature_ids, output_ids)
  unexpected_ids <- setdiff(output_ids, feature_ids)
  if (length(missing_ids) > 0 || length(unexpected_ids) > 0) {
    stop(
      context, " feature identifiers do not match the ALDEx2 test output",
      if (length(missing_ids) > 0) {
        paste0(
          "; missing: ",
          paste(utils::head(missing_ids, 5), collapse = ", ")
        )
      } else "",
      if (length(unexpected_ids) > 0) {
        paste0(
          "; unexpected: ",
          paste(utils::head(unexpected_ids, 5), collapse = ", ")
        )
      } else "",
      ".",
      call. = FALSE
    )
  }

  required_columns <- c("effect", "diff.btw", "rab.all", "overlap")
  missing_columns <- setdiff(required_columns, colnames(effect_results))
  if (length(missing_columns) > 0) {
    stop(
      context, " is missing required column(s): ",
      paste(missing_columns, collapse = ", "),
      ".",
      call. = FALSE
    )
  }

  effect_results <- effect_results[
    match(feature_ids, output_ids),
    required_columns,
    drop = FALSE
  ]
  for (column_name in required_columns) {
    validate_finite_numeric_values(
      effect_results[[column_name]],
      column_name,
      context,
      allow_na = FALSE
    )
  }
  validate_probability_values(
    effect_results$overlap,
    "overlap",
    context,
    allow_na = FALSE
  )
  rownames(effect_results) <- feature_ids
  effect_results
}

# Helper function: Perform DESeq2 analysis
perform_deseq2_analysis <- function(abundance_mat, metadata, group, reference,
                                    Level, length_Level, p_adjust_method) {
  counts <- prepare_integer_counts(as.matrix(abundance_mat), "DESeq2")
  feature_ids <- rownames(abundance_mat)

  # Resolve the reference level. When the user supplies one we honor it,
  # otherwise fall back to the first level so behavior stays identical to
  # the previous default.
  ref_level <- if (!is.null(reference) && reference %in% Level) reference else Level[1]
  non_ref_levels <- setdiff(Level, ref_level)
  factor_levels <- c(ref_level, non_ref_levels)

  # Ensure group column is a factor with reference in the first position so
  # DESeq2's default contrasts match the semantics expected downstream.
  metadata_group <- factor(metadata[[group]], levels = factor_levels)
  model_input <- prepare_model_group_column(metadata, group, metadata_group)
  metadata <- model_input$metadata
  model_group <- model_input$column

  # Create DESeqDataSet object
  se <- SummarizedExperiment::SummarizedExperiment(
    assays = list(counts = counts),
    colData = metadata
  )

  result <- tryCatch({
    dds <- DESeq2::DESeqDataSet(
      se,
      design = stats::as.formula(paste0("~", model_group))
    )

    dds <- DESeq2::estimateSizeFactors(dds)

    # Estimate dispersions with fallback chain: parametric -> local -> mean -> gene-wise
    # Note: fitType choice should NOT be based on sample size (DESeq2 documentation);
    # parametric is the default and recommended for most cases.
    dds <- tryCatch({
      DESeq2::estimateDispersions(dds, fitType = "parametric")
    }, error = function(e1) {
      tryCatch({
        DESeq2::estimateDispersions(dds, fitType = "local")
      }, error = function(e2) {
        tryCatch({
          DESeq2::estimateDispersions(dds, fitType = "mean")
        }, error = function(e3) {
          message("All dispersion fitting methods failed, using gene-wise estimates...")
          dds <- DESeq2::estimateDispersionsGeneEst(dds)
          DESeq2::dispersions(dds) <- SummarizedExperiment::mcols(dds)$dispGeneEst
          return(dds)
        })
      })
    })

    dds <- DESeq2::nbinomWaldTest(dds)

    # One row-block per non-reference level, identical shape to
    # edgeR's multi-group output.
    build_block <- function(lvl) {
      res <- DESeq2::results(
        dds,
        contrast = c(model_group, lvl, ref_level),
        pAdjustMethod = p_adjust_method
      )
      context <- paste0("DESeq2 results for '", lvl, "' vs '",
                        ref_level, "'")
      p_values <- validate_and_align_backend_result_values(
        values = res$pvalue,
        output_ids = rownames(res),
        feature_ids = feature_ids,
        value_name = "pvalue",
        context = context,
        value_type = "probability",
        allow_na = TRUE
      )
      p_adjust <- validate_and_align_backend_result_values(
        values = res$padj,
        output_ids = rownames(res),
        feature_ids = feature_ids,
        value_name = "padj",
        context = context,
        value_type = "probability",
        allow_na = TRUE
      )
      log2_fold_change <- validate_and_align_backend_result_values(
        values = res$log2FoldChange,
        output_ids = rownames(res),
        feature_ids = feature_ids,
        value_name = "log2FoldChange",
        context = context,
        value_type = "finite_numeric",
        allow_na = TRUE
      )
      data.frame(
        feature = feature_ids,
        method = "DESeq2",
        group1 = ref_level,
        group2 = lvl,
        p_values = p_values,
        p_adjust = p_adjust,
        log2_fold_change = log2_fold_change,
        stringsAsFactors = FALSE
      )
    }

    if (length_Level == 2) {
      build_block(non_ref_levels)
    } else {
      do.call(rbind, lapply(non_ref_levels, build_block))
    }
  }, error = function(e) {
    stop("DESeq2 analysis failed: ", e$message)
  })

  if (is.null(result)) {
    stop("DESeq2 analysis failed to produce results")
  }

  attr(result, "adj_method") <- paste0(p_adjust_method, " (method-specific)")
  return(result)
}

# Helper function: Perform limma voom analysis
perform_limma_voom_analysis <- function(abundance_mat, Group, reference, Level, length_Level) {
  
  # Ensure Group is a factor type
  Group <- factor(Group)
  if (!is.null(reference)) {
    Group <- relevel(Group, ref = reference)
  }

  # Create design matrix
  design <- model.matrix(~Group)

  # Perform voom transformation and analysis
  dge <- edgeR::DGEList(counts = abundance_mat, group = Group)
  dge <- edgeR::calcNormFactors(dge)
  v <- limma::voom(dge, design)
  fit <- limma::lmFit(v, design)
  fit <- limma::eBayes(fit)

  feature_ids <- rownames(abundance_mat)
  if (is.null(fit$p.value) || !is.matrix(fit$p.value) ||
      is.null(fit$coefficients) || !is.matrix(fit$coefficients)) {
    stop(
      "limma voom analysis failed: limma::eBayes() did not return ",
      "matrix p.value and coefficients components.",
      call. = FALSE
    )
  }
  if (ncol(fit$p.value) < length_Level ||
      ncol(fit$coefficients) < length_Level) {
    stop(
      "limma voom analysis failed: backend result has fewer coefficient ",
      "columns than the fitted group design.",
      call. = FALSE
    )
  }

  build_block <- function(coef_index, group1_label, group2_label) {
    context <- paste0("limma voom results for '", group2_label,
                      "' vs '", group1_label, "'")
    p_values <- validate_and_align_backend_result_values(
      values = fit$p.value[, coef_index],
      output_ids = rownames(fit$p.value),
      feature_ids = feature_ids,
      value_name = "p_values",
      context = context,
      value_type = "probability",
      allow_na = TRUE
    )
    log2_fold_change <- validate_and_align_backend_result_values(
      values = fit$coefficients[, coef_index],
      output_ids = rownames(fit$coefficients),
      feature_ids = feature_ids,
      value_name = "log2_fold_change",
      context = context,
      value_type = "finite_numeric",
      allow_na = TRUE
    )
    data.frame(
      feature = feature_ids,
      method = "limma voom",
      group1 = group1_label,
      group2 = group2_label,
      p_values = p_values,
      log2_fold_change = log2_fold_change,
      stringsAsFactors = FALSE
    )
  }

  # Extract results
  if (length_Level == 2) {
    # Two-group comparison - ensure consistent format with other methods
    group_levels <- levels(Group)
    results <- build_block(
      coef_index = 2,
      group1_label = group_levels[1],
      group2_label = group_levels[2]
    )
  } else {
    # Multi-group comparison handling.
    # fit$p.value columns [-1] correspond to levels(Group)[-1] in order
    # (since Group was releveled so the reference sits first). Build one
    # feature-aligned block per coefficient instead of flattening the matrix
    # by position, so backend row identifiers remain part of the contract.
    contrasts <- levels(Group)[-1]
    group1_label <- if (!is.null(reference)) reference else levels(Group)[1]
    results <- do.call(rbind, Map(
      function(coef_index, contrast_label) {
        build_block(
          coef_index = coef_index,
          group1_label = group1_label,
          group2_label = contrast_label
        )
      },
      seq_along(contrasts) + 1,
      contrasts
    ))
  }

  return(results)
}

# Helper function: Perform edgeR analysis
perform_edger_analysis <- function(abundance_mat, Group, reference, Level, length_Level) {

  # Resolve the reference level. Without this, edgeR's exactTest() used the
  # raw factor order (level 1 vs level 2), so the user-supplied `reference`
  # was silently ignored and the result labels were always `Level[1]` /
  # `Level[2]` regardless of intent. Relevel once so both the model and the
  # labels align with the requested contrast direction.
  ref_level <- if (!is.null(reference) && reference %in% Level) reference else Level[1]
  Group <- stats::relevel(factor(Group, levels = Level), ref = ref_level)
  Level <- levels(Group)

  # Create DGEList object
  counts <- prepare_integer_counts(abundance_mat, "edgeR")
  feature_ids <- rownames(abundance_mat)
  dge <- edgeR::DGEList(counts = counts, group = Group)
  dge <- edgeR::calcNormFactors(dge)
  dge <- edgeR::estimateCommonDisp(dge, verbose = TRUE)

  build_block <- function(et, group1_label, group2_label) {
    context <- paste0("edgeR exactTest results for '", group2_label,
                      "' vs '", group1_label, "'")
    if (is.null(et$table) || !is.data.frame(et$table) ||
        !all(c("PValue", "logFC") %in% colnames(et$table))) {
      stop(
        "edgeR analysis failed: exactTest() did not return a table with ",
        "PValue and logFC columns.",
        call. = FALSE
      )
    }
    p_values <- validate_and_align_backend_result_values(
      values = et$table$PValue,
      output_ids = rownames(et$table),
      feature_ids = feature_ids,
      value_name = "PValue",
      context = context,
      value_type = "probability",
      allow_na = FALSE
    )
    log2_fold_change <- validate_and_align_backend_result_values(
      values = et$table$logFC,
      output_ids = rownames(et$table),
      feature_ids = feature_ids,
      value_name = "logFC",
      context = context,
      value_type = "finite_numeric",
      allow_na = FALSE
    )
    data.frame(
      feature = feature_ids,
      method = "edgeR",
      group1 = group1_label,
      group2 = group2_label,
      p_values = p_values,
      log2_fold_change = log2_fold_change,
      stringsAsFactors = FALSE
    )
  }

  if (length_Level == 2) {
    # Two-group comparison: reference (level 1) vs non-reference (level 2).
    # exactTest(pair = c(1, 2)) is `log(level2 / level1)`, matching
    # group1 = ref, group2 = non-ref.
    et <- edgeR::exactTest(dge, pair = c(1, 2))
    results <- build_block(et, group1_label = Level[1], group2_label = Level[2])
  } else {
    # Multi-group: emit one block per (reference, non-reference) contrast so
    # the shape matches DESeq2 / limma voom / LinDA / Maaslin2. Previously
    # we enumerated every pair via combn(), which produced k*(k-1)/2 rows
    # per feature without any privileged reference -- inconsistent with the
    # other backends and with the documented `reference` semantics.
    non_ref_levels <- Level[-1]
    results_list <- lapply(non_ref_levels, function(lvl) {
      et <- edgeR::exactTest(dge, pair = c(Level[1], lvl))
      build_block(et, group1_label = Level[1], group2_label = lvl)
    })
    results <- do.call(rbind, results_list)
  }

  return(results)
}

# Helper function: Perform metagenomeSeq analysis
perform_metagenomeseq_analysis <- function(abundance_mat, metadata, group, reference, Level) {
  new_mrexperiment <- getExportedValue("metagenomeSeq", "newMRexperiment")
  cum_norm <- getExportedValue("metagenomeSeq", "cumNorm")
  cum_norm_stat_fast <- getExportedValue("metagenomeSeq", "cumNormStatFast")
  fit_feature_model <- getExportedValue("metagenomeSeq", "fitFeatureModel")

  # Convert metadata to data.frame. align_samples() has already ordered
  # metadata rows to match abundance_mat columns, so we set rownames from
  # the (authoritative) abundance column names rather than trusting a
  # hardcoded "sample" column that may not exist.
  metadata <- as.data.frame(metadata)
  rownames(metadata) <- colnames(abundance_mat)

  # Resolve the reference level and relevel the grouping factor so the
  # ~group model matrix uses the user-specified level as the intercept.
  # Without this, metagenomeSeq always contrasted against the natural
  # first level and our group1/group2 labels (Level[1]/Level[2]) were
  # fixed regardless of the `reference` argument, so flipping the
  # reference left both labels and coefficients unchanged.
  ref_level <- if (!is.null(reference) && reference %in% Level) reference else Level[1]
  metadata[[group]] <- stats::relevel(
    factor(metadata[[group]], levels = Level),
    ref = ref_level
  )
  Level <- levels(metadata[[group]])

  counts_all <- prepare_integer_counts(as.matrix(abundance_mat), "metagenomeSeq")

  # Single two-group fit. `fitFeatureModel()` is metagenomeSeq's documented
  # entry point for two-group comparisons (it tests a single coefficient and
  # returns one p-value per feature in `fit@pvalues`). For each pairwise
  # contrast we fit this two-group model on the subset of samples in the two
  # levels of interest and extract named vectors from the fitted model object.
  fit_pair <- function(counts_sub, meta_sub, group_col) {
    model_input <- prepare_model_group_column(
      meta_sub,
      group_col,
      meta_sub[[group_col]]
    )
    meta_model <- model_input$metadata
    model_group <- model_input$column
    context <- paste0(
      "metagenomeSeq comparison '",
      levels(meta_sub[[group_col]])[1],
      "' vs '",
      levels(meta_sub[[group_col]])[2],
      "'"
    )
    phenoData <- new("AnnotatedDataFrame",
                     data = meta_model,
                     varMetadata = data.frame(
                       labelDescription = colnames(meta_model),
                       row.names = colnames(meta_model)
                     ))
    obj <- try(
      new_mrexperiment(counts = counts_sub, phenoData = phenoData,
                       featureData = NULL, libSize = NULL, normFactors = NULL),
      silent = TRUE
    )
    if (inherits(obj, "try-error")) {
      stop("Failed to create MRexperiment object: ",
           attr(obj, "condition")$message)
    }
    p_norm <- resolve_metagenomeseq_css_quantile(
      obj = obj,
      cum_norm_stat_fast = cum_norm_stat_fast,
      context = context
    )
    obj <- cum_norm(obj, p = p_norm)

    mod <- stats::model.matrix(stats::as.formula(paste0("~", model_group)),
                               data = meta_model)
    fit <- fit_feature_model(obj, mod)

    extract_metagenomeseq_fit_values(
      fit = fit,
      feature_ids = rownames(counts_sub),
      context = context
    )
  }

  if (length(Level) == 2) {
    res <- fit_pair(counts_all, metadata, group)
    return(data.frame(
      feature = rownames(abundance_mat),
      method = "metagenomeSeq",
      group1 = Level[1],
      group2 = Level[2],
      p_values = res$p_values,
      log2_fold_change = res$log2_fold_change,
      stringsAsFactors = FALSE
    ))
  }

  # Multi-group: emit one (reference, non-reference) contrast block per
  # non-reference level. Previously the function built a full k-column
  # model matrix, called fit_feature_model() once, extracted coef = 2,
  # and hard-coded group1 = Level[1] / group2 = Level[2] -- which
  # silently dropped all contrasts beyond the first non-reference level
  # and shaped the output as if only two groups existed. fitFeatureModel
  # is also documented as a two-group entry point, so re-fitting per pair
  # on the subset of samples in the two levels of interest is the correct
  # multi-group handling and matches the shape returned by edgeR / LinDA
  # / Maaslin2 / DESeq2 / limma voom.
  non_ref_levels <- Level[-1]
  results_list <- lapply(non_ref_levels, function(lvl) {
    keep <- metadata[[group]] %in% c(Level[1], lvl)
    meta_sub <- metadata[keep, , drop = FALSE]
    # Relevel within the subset so the 2-column model matrix has the
    # reference level first and the current non-reference level second.
    meta_sub[[group]] <- factor(as.character(meta_sub[[group]]),
                                levels = c(Level[1], lvl))
    counts_sub <- counts_all[, keep, drop = FALSE]
    res <- fit_pair(counts_sub, meta_sub, group)
    data.frame(
      feature = rownames(abundance_mat),
      method = "metagenomeSeq",
      group1 = Level[1],
      group2 = lvl,
      p_values = res$p_values,
      log2_fold_change = res$log2_fold_change,
      stringsAsFactors = FALSE
    )
  })
  do.call(rbind, results_list)
}

resolve_metagenomeseq_css_quantile <- function(obj,
                                               cum_norm_stat_fast,
                                               context) {
  p_norm <- tryCatch(
    cum_norm_stat_fast(obj),
    error = function(e) {
      msg <- conditionMessage(e)
      if (grepl("one or zero features", msg, ignore.case = TRUE)) {
        stop(
          context,
          " failed to estimate metagenomeSeq CSS normalization quantile: ",
          "cumNormStatFast() requires at least two positive features in ",
          "every sample. Remove or reprocess samples with only one non-zero ",
          "feature before running metagenomeSeq.",
          call. = FALSE
        )
      }
      if (grepl("missing value where TRUE/FALSE needed", msg, fixed = TRUE)) {
        warning(
          context,
          " could not estimate a stable metagenomeSeq CSS normalization ",
          "quantile from cumNormStatFast(); using p = 0.5, matching ",
          "metagenomeSeq's internal default branch for low estimated ",
          "quantiles.",
          call. = FALSE
        )
        return(0.5)
      }
      stop(
        context,
        " failed to estimate metagenomeSeq CSS normalization quantile: ",
        msg,
        call. = FALSE
      )
    }
  )

  if (!is.numeric(p_norm) || length(p_norm) != 1 ||
      is.na(p_norm) || !is.finite(p_norm) ||
      p_norm <= 0 || p_norm > 1) {
    warning(
      context,
      " returned an invalid metagenomeSeq CSS normalization quantile; ",
      "using p = 0.5.",
      call. = FALSE
    )
    p_norm <- 0.5
  }

  unname(p_norm)
}

extract_metagenomeseq_fit_values <- function(fit, feature_ids, context) {
  p_values <- align_metagenomeseq_fit_vector(
    values = fit@pvalues,
    feature_ids = feature_ids,
    value_name = "p_values",
    context = context,
    fallback_names = fit@taxa
  )
  validate_probability_values(p_values, "p_values", context,
                              allow_na = FALSE)

  logfc <- NULL
  if (!is.null(fit@fitZeroLogNormal) &&
      !is.null(fit@fitZeroLogNormal$logFC)) {
    logfc <- fit@fitZeroLogNormal$logFC
  }
  logfc <- align_metagenomeseq_fit_vector(
    values = logfc,
    feature_ids = feature_ids,
    value_name = "logFC",
    context = context,
    fallback_names = fit@taxa
  )
  validate_finite_numeric_values(logfc, "logFC", context,
                                 allow_na = FALSE)

  list(p_values = p_values, log2_fold_change = logfc)
}

align_metagenomeseq_fit_vector <- function(values,
                                           feature_ids,
                                           value_name,
                                           context,
                                           fallback_names = NULL) {
  feature_ids <- validate_nonempty_character_column(
    feature_ids,
    "feature_ids",
    context
  )
  if (!is.numeric(values)) {
    stop(context, " did not return numeric metagenomeSeq ", value_name,
         " values.", call. = FALSE)
  }

  value_names <- names(values)
  invalid_names <- is.null(value_names) ||
    any(is.na(value_names)) ||
    any(!nzchar(value_names))
  if (invalid_names) {
    if (!is.null(fallback_names) && length(fallback_names) == length(values)) {
      value_names <- as.character(fallback_names)
    } else if (length(values) == length(feature_ids)) {
      value_names <- feature_ids
    } else {
      stop(
        context,
        " returned unnamed metagenomeSeq ", value_name,
        " values with length ", length(values),
        ", but ", length(feature_ids),
        " feature(s) were expected.",
        call. = FALSE
      )
    }
    names(values) <- value_names
  }

  if (anyDuplicated(value_names)) {
    duplicated_names <- unique(value_names[duplicated(value_names)])
    stop(
      context,
      " returned duplicated metagenomeSeq ", value_name,
      " feature identifiers: ",
      paste(utils::head(duplicated_names, 5), collapse = ", "),
      ".",
      call. = FALSE
    )
  }

  missing_features <- setdiff(feature_ids, value_names)
  if (length(missing_features) > 0) {
    stop(
      context,
      " did not return metagenomeSeq ", value_name,
      " values for feature(s): ",
      paste(utils::head(missing_features, 5), collapse = ", "),
      if (length(missing_features) > 5) ", ..." else "",
      ".",
      call. = FALSE
    )
  }

  values[feature_ids]
}

# Helper function: Perform Maaslin2 analysis
perform_maaslin2_analysis <- function(abundance_mat, metadata, group, reference,
                                      Level, length_Level, p_adjust_method) {

  # Prepare metadata first
  metadata <- as.data.frame(metadata)

  # Ensure metadata rownames match the column names of abundance_mat
  # This is crucial because the main function has already reordered metadata to match abundance
  rownames(metadata) <- colnames(abundance_mat)

  feature_map <- prepare_maaslin2_feature_mapping(rownames(abundance_mat))

  # Transpose abundance matrix (samples become rows, features become columns)
  abundance_mat_t <- t(abundance_mat)

  # Use a run-specific directory to avoid stale files from previous Maaslin2 runs.
  output_dir <- tempfile(pattern = "ggpicrust2_maaslin2_")
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  on.exit(unlink(output_dir, recursive = TRUE, force = TRUE), add = TRUE)

  # Maaslin2 emits an early logging::logwarn() before it sets up its own
  # log handler. If a prior Maaslin2 invocation in the same R session left
  # a writeToFile handler pointing at a now-deleted directory, that early
  # warning blows up with "cannot open the connection". Clear stale
  # handlers up front so repeat calls stay safe.
  if (requireNamespace("logging", quietly = TRUE)) {
    root_logger <- logging::getLogger()
    for (h_name in names(root_logger$handlers)) {
      try(logging::removeHandler(h_name), silent = TRUE)
    }
  }

  # Resolve the reference level up front so we can (a) pass it to Maaslin2
  # regardless of group count and (b) label group1 consistently. Without
  # forwarding reference in the 2-group case, Maaslin2 falls back to its
  # alphabetical-first-level default -- if the user asked for a different
  # reference, the resulting contrast direction silently disagrees with the
  # group1/group2 labels we emit.
  ref_level <- if (!is.null(reference) && reference %in% Level) reference else Level[1]
  maaslin_levels <- c(ref_level, setdiff(Level, ref_level))
  model_input <- prepare_model_group_column(
    metadata,
    group,
    factor(metadata[[group]], levels = maaslin_levels)
  )
  metadata <- model_input$metadata
  model_group <- model_input$column
  maaslin_reference <- paste0(model_group, ",", ref_level)

  # Run Maaslin2 analysis via dynamic lookup so the package remains optional
  maaslin2_fn <- getExportedValue("Maaslin2", "Maaslin2")
  fit_data <- maaslin2_fn(
    input_data = abundance_mat_t,
    input_metadata = metadata,
    output = output_dir,
    transform = "AST",
    fixed_effects = model_group,
    reference = maaslin_reference,
    correction = p_adjust_method,
    normalization = "TSS",
    standardize = TRUE,
    min_prevalence = 0.1,
    cores = 1,
    plot_heatmap = FALSE,
    plot_scatter = FALSE
  )

  # Read all results instead of significant results
  results_file <- file.path(output_dir, "all_results.tsv")

  if (!file.exists(results_file)) {
    stop("Maaslin2 analysis failed to produce results file")
  }

  maaslin2_results <- utils::read.table(results_file,
                                        header = TRUE,
                                        sep = "\t",
                                        stringsAsFactors = FALSE)

  maaslin2_results <- validate_maaslin2_results(
    maaslin2_results = maaslin2_results,
    group = model_group,
    feature_map = feature_map
  )

  results <- data.frame(
    feature = maaslin2_results$feature,
    method = "Maaslin2",
    group1 = ref_level,
    group2 = maaslin2_results$value,
    p_values = maaslin2_results$pval,
    p_adjust = maaslin2_results$qval,
    log2_fold_change = maaslin2_results$coef,
    stringsAsFactors = FALSE
  )

  attr(results, "adj_method") <- paste0(p_adjust_method, " (method-specific)")
  return(results)
}

prepare_maaslin2_feature_mapping <- function(feature_ids) {
  feature_ids <- validate_nonempty_character_column(
    feature_ids,
    "feature identifiers",
    "Maaslin2 input"
  )
  maaslin2_ids <- make.names(feature_ids)
  duplicated_ids <- unique(maaslin2_ids[duplicated(maaslin2_ids)])
  if (length(duplicated_ids) > 0) {
    examples <- vapply(
      utils::head(duplicated_ids, 5),
      function(id) {
        paste0(id, " <- ",
               paste(feature_ids[maaslin2_ids == id], collapse = ", "))
      },
      character(1)
    )
    stop(
      "Maaslin2 sanitizes feature identifiers with make.names(), which ",
      "would make some feature IDs ambiguous: ",
      paste(examples, collapse = "; "),
      ". Rename these features before running pathway_daa(daa_method = ",
      "'Maaslin2').",
      call. = FALSE
    )
  }
  stats::setNames(feature_ids, maaslin2_ids)
}

validate_maaslin2_results <- function(maaslin2_results, group, feature_map) {
  context <- "Maaslin2 all_results.tsv"
  if (!is.data.frame(maaslin2_results)) {
    stop(context, " must be a data frame.", call. = FALSE)
  }

  required_columns <- c("feature", "metadata", "value", "coef", "pval", "qval")
  missing_columns <- setdiff(required_columns, colnames(maaslin2_results))
  if (length(missing_columns) > 0) {
    stop(
      context,
      " is missing required column(s): ",
      paste(missing_columns, collapse = ", "),
      ". Expected columns include feature, metadata, value, coef, pval, and qval.",
      call. = FALSE
    )
  }

  metadata_ids <- validate_nonempty_character_column(
    maaslin2_results$metadata,
    "metadata",
    context
  )
  maaslin2_results <- maaslin2_results[metadata_ids == group, , drop = FALSE]
  if (nrow(maaslin2_results) == 0) {
    return(maaslin2_results[, required_columns, drop = FALSE])
  }

  output_features <- validate_nonempty_character_column(
    maaslin2_results$feature,
    "feature",
    context
  )
  unexpected_features <- setdiff(output_features, names(feature_map))
  if (length(unexpected_features) > 0) {
    stop(
      context,
      " contains feature identifier(s) that do not match the sanitized ",
      "Maaslin2 input feature IDs: ",
      paste(utils::head(unexpected_features, 5), collapse = ", "),
      if (length(unexpected_features) > 5) ", ..." else "",
      ".",
      call. = FALSE
    )
  }
  maaslin2_results$feature <- unname(feature_map[output_features])

  maaslin2_results$value <- validate_nonempty_character_column(
    maaslin2_results$value,
    "value",
    context
  )
  duplicate_pairs <- duplicated(paste(maaslin2_results$feature,
                                      maaslin2_results$value,
                                      sep = "\r"))
  if (any(duplicate_pairs)) {
    duplicated_rows <- paste0(
      maaslin2_results$feature[duplicate_pairs],
      " / ",
      maaslin2_results$value[duplicate_pairs]
    )
    stop(
      context,
      " contains duplicated feature/contrast rows: ",
      paste(utils::head(duplicated_rows, 5), collapse = ", "),
      if (sum(duplicate_pairs) > 5) ", ..." else "",
      ".",
      call. = FALSE
    )
  }

  validate_probability_values(maaslin2_results$pval, "pval", context,
                              allow_na = FALSE)
  validate_probability_values(maaslin2_results$qval, "qval", context,
                              allow_na = FALSE)
  validate_finite_numeric_values(maaslin2_results$coef, "coef", context,
                                 allow_na = FALSE)

  maaslin2_results
}

# Helper function: Perform Lefser analysis
perform_lefser_analysis <- function(abundance_mat, metadata, group, reference, Level) {

  # Check if we only have 2 groups (lefser requirement)
  if (length(Level) != 2) {
    stop("Lefser requires exactly 2 groups. Found ", length(Level), " groups.")
  }

  ref_level <- if (!is.null(reference) && reference %in% Level) reference else Level[1]
  metadata[[group]] <- stats::relevel(
    factor(metadata[[group]], levels = Level),
    ref = ref_level
  )
  Level <- levels(metadata[[group]])

  # Create SummarizedExperiment object
  se <- SummarizedExperiment::SummarizedExperiment(
    assays = list(counts = abundance_mat),
    colData = metadata
  )

  # lefser expects relative abundances scaled to one million per sample
  # (`lefser::relativeAb()` prepends this assay and `lefser()` uses assay 1
  # by default). Passing raw counts lets library-size differences drive
  # Kruskal-Wallis/Wilcoxon and LDA statistics, which is not the LEfSe data
  # scale and triggers a warning from lefser itself.
  se <- lefser::relativeAb(se)
  rel_abundance_mat <- SummarizedExperiment::assay(se, i = 1L)

  # Perform Lefser analysis to get LDA scores (effect sizes)
  lefser_results <- tryCatch({
    lefser::lefser(se, classCol = group)
  }, error = function(e) {
    stop("Lefser analysis failed: ", e$message, call. = FALSE)
  })

  # Get all features and group assignments
  all_features <- rownames(abundance_mat)
  group_vector <- metadata[[group]]

  # Calculate Kruskal-Wallis p-values for all features on the same
  # relative-abundance scale used by lefser. Without a subclass variable,
  # lefser filters features with Kruskal-Wallis before LDA; Wilcoxon is only
  # part of the subclass consistency step.
  p_values <- vapply(seq_len(nrow(abundance_mat)), function(i) {
    calculate_lefser_kruskal_p_value(
      feature_values = rel_abundance_mat[i, ],
      group_vector = group_vector,
      feature_id = all_features[i]
    )
  }, numeric(1))

  # Initialize results with Wilcoxon p-values
  results <- data.frame(
    feature = all_features,
    method = "Lefser",
    group1 = Level[1],
    group2 = Level[2],
    p_values = p_values,
    lda_score = rep(NA_real_, length(all_features)),  # LDA effect size from lefser
    stringsAsFactors = FALSE
  )

  # Add LDA scores for features identified by lefser
  if (!is.null(lefser_results) && nrow(lefser_results) > 0) {
    sig_indices <- match(lefser_results$features, all_features)
    valid_indices <- !is.na(sig_indices)
    results$lda_score[sig_indices[valid_indices]] <- lefser_results$scores[valid_indices]
    message(sprintf("Lefser identified %d features with LDA score > threshold", sum(valid_indices)))
  } else {
    message("No significant features found by Lefser LDA analysis.")
  }

  return(results)
}

calculate_lefser_kruskal_p_value <- function(feature_values,
                                             group_vector,
                                             feature_id = "feature") {
  feature_values <- as.numeric(feature_values)
  if (length(feature_values) != length(group_vector)) {
    stop(
      "Lefser Kruskal-Wallis test failed for feature '", feature_id,
      "': feature values length (", length(feature_values),
      ") does not match group vector length (", length(group_vector), ").",
      call. = FALSE
    )
  }

  if (anyNA(feature_values) || any(!is.finite(feature_values))) {
    stop(
      "Lefser Kruskal-Wallis test failed for feature '", feature_id,
      "': relative abundance values must be finite and non-missing.",
      call. = FALSE
    )
  }

  if (length(unique(feature_values)) <= 1) {
    return(1.0)
  }

  p_value <- tryCatch(
    stats::kruskal.test(feature_values ~ group_vector)$p.value,
    error = function(e) {
      stop(
        "Lefser Kruskal-Wallis test failed for feature '", feature_id,
        "': ", e$message,
        call. = FALSE
      )
    }
  )

  validate_probability_values(p_value, "p_value",
                              paste0("Lefser feature '", feature_id, "'"))
  p_value
}

# Helper function: Perform LinDA analysis
format_linda_output <- function(linda_output, group, reference, Level) {
  # Create an empty list to store results for each comparison
  results_list <- list()

  # Get all non-reference groups
  non_ref_groups <- Level[Level != reference]

  # Process each output element (each corresponds to a group comparison)
  for (i in seq_along(linda_output)) {
    output_name <- names(linda_output)[i]

    # Extract non-reference group name from output name
    group_prefix <- paste0(group, "")
    if (!is.null(output_name) && startsWith(output_name, group_prefix)) {
      comparison_group <- substring(output_name, nchar(group_prefix) + 1)
    } else {
      comparison_group <- if (i <= length(non_ref_groups)) non_ref_groups[i] else paste0("Group", i)
    }
    comparison_group <- as.character(comparison_group)[1]

    comparison_df <- linda_output[[i]]
    if (!is.data.frame(comparison_df) || nrow(comparison_df) == 0) {
      next
    }

    n_features <- nrow(comparison_df)
    feature_ids <- validate_linda_feature_ids(comparison_df,
                                              comparison_group)

    pvals <- extract_linda_output_column(
      comparison_df,
      column_name = "pvalue",
      comparison_group = comparison_group,
      n_features = n_features,
      probability = TRUE
    )
    lfc <- extract_linda_output_column(
      comparison_df,
      column_name = "log2FoldChange",
      comparison_group = comparison_group,
      n_features = n_features,
      probability = FALSE
    )
    padj <- extract_linda_output_column(
      comparison_df,
      column_name = "padj",
      comparison_group = comparison_group,
      n_features = n_features,
      probability = TRUE
    )

    # Create a results data frame for this comparison
    results_list[[length(results_list) + 1]] <- data.frame(
      feature = feature_ids,
      method = "LinDA",
      group1 = reference,
      group2 = comparison_group,
      p_values = as.numeric(pvals),
      p_adjust = as.numeric(padj),
      log2_fold_change = as.numeric(lfc),
      stringsAsFactors = FALSE
    )
  }

  # Combine all results
  if (length(results_list) > 0) {
    do.call(rbind, results_list)
  } else {
    create_empty_daa_result(include_log2fc = TRUE)
  }
}

validate_linda_feature_ids <- function(comparison_df, comparison_group) {
  feature_ids <- rownames(comparison_df)
  context <- paste0("LinDA output for comparison '", comparison_group, "'")
  feature_ids <- validate_nonempty_character_column(feature_ids,
                                                    "rownames",
                                                    context)
  if (anyDuplicated(feature_ids)) {
    duplicated_ids <- unique(feature_ids[duplicated(feature_ids)])
    stop(
      context, " contains duplicated feature identifiers: ",
      paste(utils::head(duplicated_ids, 5), collapse = ", "),
      call. = FALSE
    )
  }
  feature_ids
}

extract_linda_output_column <- function(comparison_df,
                                        column_name,
                                        comparison_group,
                                        n_features,
                                        probability = FALSE) {
  context <- paste0("LinDA output for comparison '", comparison_group, "'")
  if (!column_name %in% colnames(comparison_df)) {
    stop(
      context, " is missing required column '", column_name,
      "'. MicrobiomeStat::linda() output is expected to contain ",
      "'baseMean', 'log2FoldChange', 'lfcSE', 'stat', 'pvalue', ",
      "'padj', 'reject', and 'df'.",
      call. = FALSE
    )
  }

  values <- comparison_df[[column_name]]
  if (length(values) != n_features || is.matrix(values) || is.data.frame(values)) {
    stop(
      context, " column '", column_name, "' has invalid length or shape. ",
      "Expected one value per feature (", n_features, "), got length ",
      length(values), ".",
      call. = FALSE
    )
  }

  numeric_values <- suppressWarnings(as.numeric(values))
  if (anyNA(numeric_values) || any(!is.finite(numeric_values))) {
    stop(
      context, " column '", column_name,
      "' must contain finite numeric values without NA.",
      call. = FALSE
    )
  }

  if (probability) {
    validate_probability_values(numeric_values, column_name, context)
  }

  numeric_values
}

perform_linda_analysis <- function(abundance, metadata, group, reference, Level,
                                   length_Level, p_adjust_method) {
  # Filter zero-abundance features using unified validation
  feature.dat <- validate_daa_input(as.matrix(abundance), method = "LinDA", filter_zero = TRUE)
  meta.dat <- metadata

  # Handle reference level
  if (!is.null(reference) && !(reference %in% Level)) {
    warning(sprintf("Reference '%s' not in group levels; using '%s'", reference, Level[1]))
    reference <- Level[1]
  } else if (is.null(reference)) {
    reference <- Level[1]
  }

  # Relevel the grouping factor so that LinDA's `~ group` formula actually
  # treats the user-specified level as the reference. Without this, LinDA
  # silently keeps the natural first level as reference while our result
  # labeling uses `reference` -- producing rows with group1 == group2 and
  # statistics that don't reflect the requested contrast direction.
  meta.dat[[group]] <- stats::relevel(
    factor(meta.dat[[group]], levels = Level),
    ref = reference
  )
  # Keep Level in sync with the releveled factor so downstream lookups
  # (e.g. non-reference level enumeration) agree with the model fit.
  Level <- levels(meta.dat[[group]])

  model_input <- prepare_model_group_column(meta.dat, group, meta.dat[[group]])
  meta.dat <- model_input$metadata
  model_group <- model_input$column
  formula <- paste0("~ ", model_group)

  # Run LinDA
  linda_result <- tryCatch({
    # Run LinDA analysis
    linda_obj <- MicrobiomeStat::linda(
      feature.dat = feature.dat,
      meta.dat = meta.dat,
      formula = formula,
      feature.dat.type = "count",
      prev.filter = 0,
      mean.abund.filter = 0,
      adaptive = TRUE,
      zero.handling = "pseudo-count",
      pseudo.cnt = 0.5,
      p.adj.method = p_adjust_method,
      alpha = 0.05,
      n.cores = 1,
      verbose = TRUE    # Changed to TRUE to provide more detailed output
    )
    
    # Check the structure of linda_obj
    if (is.null(linda_obj) || is.null(linda_obj$output) || length(linda_obj$output) == 0) {
      message("LinDA analysis returned empty results. This might be due to insufficient variation in the data.")
      return(create_empty_daa_result(include_log2fc = TRUE))
    }
    
    format_linda_output(
      linda_output = linda_obj$output,
      group = model_group,
      reference = reference,
      Level = Level
    )
  }, error = function(e) {
    stop("LinDA analysis failed: ", e$message, call. = FALSE)
  })
  
  attr(linda_result, "adj_method") <- paste0(p_adjust_method, " (method-specific)")
  return(linda_result)
}
