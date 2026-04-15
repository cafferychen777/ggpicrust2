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
#'        The column names should match the sample names in metadata.
#'        Values should be counts or abundance measurements.
#'
#' @param metadata A data frame or tibble containing sample information.
#'        Must include a 'sample' column with sample identifiers matching the column names in abundance data.
#'
#' @param group Character string specifying the column name in metadata that contains group information
#'        for differential abundance analysis.
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
#' @param select Vector of sample names to include in the analysis.
#'        If NULL (default), all samples are included.
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
#'        If NULL (default), the first level is used as reference.
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
#'        FALSE to skip the extra \code{aldex.effect()} computation.
#'
#' @param p.adjust Deprecated. Use \code{p_adjust_method} instead.
#'
#' @param ... Additional arguments passed to the specific DAA method
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
#' and ratio-based effect sizes.
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
  validate_abundance(abundance)
  validate_metadata(metadata)
  validate_group(metadata, group)

  # Convert abundance to matrix if needed
  abundance_mat <- as.matrix(abundance)

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
    features_available <- intersect(features, rownames(relative_abundance))
    if (length(features_available) == 0) {
      stop("None of the specified features found in abundance data")
    }
    relative_abundance <- relative_abundance[features_available, , drop = FALSE]
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

#' @rdname pathway_daa
#' @export
pathway_daa <- function(abundance, metadata, group, daa_method = "ALDEx2",
                       select = NULL, p_adjust_method = "BH", reference = NULL,
                       include_abundance_stats = FALSE, include_effect_size = TRUE,
                       p.adjust = NULL, ...) {
  # Backward compatibility for deprecated parameter
  if (!is.null(p.adjust)) {
    warning("'p.adjust' parameter is deprecated. Use 'p_adjust_method' instead.", call. = FALSE)
    p_adjust_method <- p.adjust
  }

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

  require_package(method_packages[[daa_method]], purpose = daa_method)

  # Input validation using unified functions
  validate_abundance(abundance, min_samples = 4)
  validate_metadata(metadata)
  validate_group(metadata, group, min_groups = 2)

  # Align samples. This is the contract for STANDALONE callers of
  # pathway_daa() who may pass unaligned abundance/metadata. When the
  # wrapper ggpicrust2() calls us, inputs are already aligned; since
  # align_samples() is deterministic and idempotent, that pre-aligned
  # case is a cheap no-op here and not a source of drift. The two
  # call sites (this one and ggpicrust2.R) are invoked with identical
  # arguments so they cannot silently disagree on the aligned shape;
  # see the commentary at the ggpicrust2 call site for the full
  # rationale and the invariant that any future change to alignment
  # semantics must update both sites together.
  aligned <- align_samples(abundance, metadata, verbose = FALSE)
  abundance <- aligned$abundance
  metadata <- tibble::as_tibble(aligned$metadata)
  sample_col <- aligned$sample_col

  if (aligned$n_samples < 4) {
    stop("At least 4 samples required for DAA (after alignment)")
  }

  # Handle sample selection
  if (!is.null(select)) {
    if (!all(select %in% colnames(abundance))) {
      stop("Some selected samples not in abundance data")
    }
    abundance <- abundance[, select, drop = FALSE]
    # Keep metadata rows in the exact order of `select` so row i of metadata
    # corresponds to column i of abundance. Simple %in% filtering preserves
    # the original order and can silently desynchronize Group vs samples.
    metadata <- metadata[match(select, metadata[[sample_col]]), , drop = FALSE]
  }

  # Validate abundance matrix
  abundance_mat <- as.matrix(abundance)
  validate_daa_input(abundance_mat, method = daa_method)

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
  # less actionable error.
  if (length_Level < 2) {
    stop(
      "DAA requires at least 2 groups with samples after alignment",
      if (!is.null(select)) " and `select` filtering" else "",
      "; found ", length_Level,
      if (length_Level == 1) paste0(" ('", Level, "')") else "",
      "."
    )
  }

  # Perform differential analysis
  result <- switch(
    daa_method,
    "ALDEx2" = perform_aldex2_analysis(abundance_mat, Group, Level, length_Level, include_effect_size),
    "DESeq2" = perform_deseq2_analysis(abundance_mat, metadata, group, reference, Level, length_Level),
    "LinDA" = perform_linda_analysis(abundance, metadata, group, reference, Level, length_Level),
    "limma voom" = perform_limma_voom_analysis(abundance_mat, Group, reference, Level, length_Level),
    "edgeR" = perform_edger_analysis(abundance_mat, Group, reference, Level, length_Level),
    "metagenomeSeq" = perform_metagenomeseq_analysis(abundance_mat, metadata, group, reference, Level),
    "Maaslin2" = perform_maaslin2_analysis(abundance_mat, metadata, group, reference, Level, length_Level),
    "Lefser" = perform_lefser_analysis(abundance_mat, metadata, group, Level)
  )


  # Add multiple testing correction
  # Note: Some methods (e.g., ALDEx2) provide pre-computed BH-corrected p-values

  # which are more accurate as they account for Monte Carlo sampling uncertainty
  if (!is.null(result) && "p_values" %in% colnames(result) && nrow(result) > 0) {
    if ("p_adjust" %in% colnames(result)) {
      # p_adjust already computed by the method (e.g., ALDEx2)
      # Just add the adjustment method indicator
      result$adj_method <- "BH (method-specific)"
    } else {
      # Apply standard p-value adjustment for methods without pre-computed values
      result$p_adjust <- stats::p.adjust(result$p_values, method = p_adjust_method)
      result$adj_method <- p_adjust_method
    }
  }

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
        warning("Failed to calculate abundance statistics: ", e$message)
      })
    } else {
      warning("Cannot calculate abundance statistics: required columns (feature, group1, group2) not found in results")
    }
  }

  return(result)
}

# Helper function: Perform ALDEx2 analysis
perform_aldex2_analysis <- function(abundance_mat, Group, Level, length_Level, include_effect_size = FALSE) {
  # Filter zero-abundance features (ALDEx2 requirement)
  abundance_mat <- validate_daa_input(abundance_mat, method = "ALDEx2", filter_zero = TRUE)
  abundance_mat <- round(abundance_mat)

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
      effect_results <- tryCatch({
        # In newer versions of ALDEx2, aldex.effect() only needs the clr object
        # The conditions are already stored in the clr object
        ALDEx2::aldex.effect(
          ALDEx2_object,
          verbose = FALSE
        )
      }, error = function(e) {
        warning("Failed to calculate ALDEx2 effect sizes: ", e$message,
                ". Proceeding without effect size information.")
        return(NULL)
      })
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

    # Add effect size columns if available
    if (!is.null(effect_results)) {
      # Simple approach: add the most important effect size columns
      # Replicate effect size data for both test methods (Welch and Wilcoxon)
      n_features <- nrow(results)

      # Add basic effect size columns
      if ("effect" %in% colnames(effect_results)) {
        base_df$effect_size <- rep(effect_results$effect, 2)
      }
      if ("diff.btw" %in% colnames(effect_results)) {
        base_df$diff_btw <- rep(effect_results$diff.btw, 2)
        base_df$log2_fold_change <- rep(effect_results$diff.btw, 2)
      }
      if ("rab.all" %in% colnames(effect_results)) {
        base_df$rab_all <- rep(effect_results$rab.all, 2)
      }
      if ("overlap" %in% colnames(effect_results)) {
        base_df$overlap <- rep(effect_results$overlap, 2)
      }
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
      method = c(rep("ALDEx2_Kruskal-Wallace test", nrow(results)),
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

# Helper function: Perform DESeq2 analysis
perform_deseq2_analysis <- function(abundance_mat, metadata, group, reference, Level, length_Level) {
  counts <- round(as.matrix(abundance_mat))

  # Resolve the reference level. When the user supplies one we honor it,
  # otherwise fall back to the first level so behavior stays identical to
  # the previous default.
  ref_level <- if (!is.null(reference) && reference %in% Level) reference else Level[1]
  non_ref_levels <- setdiff(Level, ref_level)
  factor_levels <- c(ref_level, non_ref_levels)

  # Ensure group column is a factor with reference in the first position so
  # DESeq2's default contrasts match the semantics expected downstream.
  metadata[[group]] <- factor(metadata[[group]], levels = factor_levels)

  # Create DESeqDataSet object
  se <- SummarizedExperiment::SummarizedExperiment(
    assays = list(counts = counts),
    colData = metadata
  )

  result <- tryCatch({
    suppressWarnings({
      dds <- DESeq2::DESeqDataSet(se, design = as.formula(paste0("~", group)))

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
        res <- DESeq2::results(dds, contrast = c(group, lvl, ref_level))
        data.frame(
          feature = rownames(abundance_mat),
          method = "DESeq2",
          group1 = ref_level,
          group2 = lvl,
          p_values = res$pvalue,
          log2_fold_change = res$log2FoldChange,
          stringsAsFactors = FALSE
        )
      }

      if (length_Level == 2) {
        build_block(non_ref_levels)
      } else {
        do.call(rbind, lapply(non_ref_levels, build_block))
      }
    })
  }, error = function(e) {
    stop("DESeq2 analysis failed: ", e$message)
  })

  if (is.null(result)) {
    stop("DESeq2 analysis failed to produce results")
  }

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

  # Extract results
  if (length_Level == 2) {
    # Two-group comparison - ensure consistent format with other methods
    group_levels <- levels(Group)
    results <- data.frame(
      feature = rownames(abundance_mat),
      method = "limma voom",
      group1 = group_levels[1],
      group2 = group_levels[2],
      p_values = fit$p.value[,2],
      log2_fold_change = fit$coefficients[,2],
      stringsAsFactors = FALSE
    )
  } else {
    # Multi-group comparison handling.
    # fit$p.value columns [-1] correspond to levels(Group)[-1] in order
    # (since Group was releveled so the reference sits first). as.vector()
    # of that matrix stacks columns head-to-tail: all features for
    # contrast 1, then all features for contrast 2, etc. group2 must be
    # repeated *each = n_features* to stay aligned; plain recycling
    # produces interleaved B,C,B,C labels that scramble the contrast.
    contrasts <- levels(Group)[-1]
    n_features <- nrow(abundance_mat)
    group1_label <- if (!is.null(reference)) reference else levels(Group)[1]
    results <- data.frame(
      feature = rep(rownames(abundance_mat), length(contrasts)),
      method = "limma voom",
      group1 = group1_label,
      group2 = rep(contrasts, each = n_features),
      p_values = as.vector(fit$p.value[, -1]),
      log2_fold_change = as.vector(fit$coefficients[, -1]),
      stringsAsFactors = FALSE
    )
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
  dge <- edgeR::DGEList(counts = round(abundance_mat), group = Group)
  dge <- edgeR::calcNormFactors(dge)
  dge <- edgeR::estimateCommonDisp(dge, verbose = TRUE)

  if (length_Level == 2) {
    # Two-group comparison: reference (level 1) vs non-reference (level 2).
    # exactTest(pair = c(1, 2)) is `log(level2 / level1)`, matching
    # group1 = ref, group2 = non-ref.
    et <- edgeR::exactTest(dge, pair = c(1, 2))
    results <- data.frame(
      feature = rownames(abundance_mat),
      method = "edgeR",
      group1 = Level[1],
      group2 = Level[2],
      p_values = et$table$PValue,
      log2_fold_change = et$table$logFC,
      stringsAsFactors = FALSE
    )
  } else {
    # Multi-group: emit one block per (reference, non-reference) contrast so
    # the shape matches DESeq2 / limma voom / LinDA / Maaslin2. Previously
    # we enumerated every pair via combn(), which produced k*(k-1)/2 rows
    # per feature without any privileged reference -- inconsistent with the
    # other backends and with the documented `reference` semantics.
    non_ref_levels <- Level[-1]
    results_list <- lapply(non_ref_levels, function(lvl) {
      et <- edgeR::exactTest(dge, pair = c(Level[1], lvl))
      data.frame(
        feature = rownames(abundance_mat),
        method = "edgeR",
        group1 = Level[1],
        group2 = lvl,
        p_values = et$table$PValue,
        log2_fold_change = et$table$logFC,
        stringsAsFactors = FALSE
      )
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
  mrcoefs <- getExportedValue("metagenomeSeq", "MRcoefs")

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

  counts_all <- round(as.matrix(abundance_mat))

  # Single two-group fit. `fitFeatureModel()` is metagenomeSeq's documented
  # entry point for two-group comparisons (it tests a single coefficient and
  # returns one p-value per feature in `fit@pvalues`). For each pairwise
  # contrast we fit this two-group model on the subset of samples in the
  # two levels of interest and read `coef = 2` -- the non-reference dummy.
  # Returns a list of two length-nrow(counts_sub) vectors so the shape is
  # stable across contrasts regardless of log-fold-change extraction
  # success.
  fit_pair <- function(counts_sub, meta_sub, group_col) {
    phenoData <- new("AnnotatedDataFrame",
                     data = meta_sub,
                     varMetadata = data.frame(
                       labelDescription = colnames(meta_sub),
                       row.names = colnames(meta_sub)
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
    # cumNormStatFast() chooses a normalization quantile by comparing
    # per-sample CDFs; with very small or near-uniform matrices (e.g. the
    # minimum 4-sample input) the underlying quantile math returns NaN,
    # and metagenomeSeq then aborts inside an `if (x <= 0.5)` check with
    # "missing value where TRUE/FALSE needed". Precompute the factor
    # ourselves and fall back to metagenomeSeq's documented default
    # (p = 0.5) whenever it is degenerate.
    p_norm <- tryCatch(cum_norm_stat_fast(obj), error = function(e) NA_real_)
    if (!is.finite(p_norm)) p_norm <- 0.5
    obj <- cum_norm(obj, p = p_norm)

    mod <- stats::model.matrix(as.formula(paste0("~", group_col)),
                               data = meta_sub)
    fit <- fit_feature_model(obj, mod)

    coef_table <- tryCatch(
      mrcoefs(fit, coef = 2),
      error = function(e) {
        warning("Failed to extract coefficients from metagenomeSeq: ",
                e$message)
        NULL
      }
    )
    log2fc <- if (!is.null(coef_table) && "logFC" %in% colnames(coef_table)) {
      coef_table$logFC
    } else {
      rep(NA_real_, nrow(counts_sub))
    }
    list(p_values = fit@pvalues, log2_fold_change = log2fc)
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

# Helper function: Perform Maaslin2 analysis
perform_maaslin2_analysis <- function(abundance_mat, metadata, group, reference, Level, length_Level) {

  # Prepare metadata first
  metadata <- as.data.frame(metadata)

  # Ensure metadata rownames match the column names of abundance_mat
  # This is crucial because the main function has already reordered metadata to match abundance
  rownames(metadata) <- colnames(abundance_mat)

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
  maaslin_reference <- paste0(group, ",", ref_level)

  # Run Maaslin2 analysis via dynamic lookup so the package remains optional
  maaslin2_fn <- getExportedValue("Maaslin2", "Maaslin2")
  fit_data <- maaslin2_fn(
    input_data = abundance_mat_t,
    input_metadata = metadata,
    output = output_dir,
    transform = "AST",
    fixed_effects = group,
    reference = maaslin_reference,
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

  # Keep only rows that correspond to contrasts on `group`. Maaslin2 emits
  # one row per (feature, non-reference level) for a categorical fixed effect,
  # which gives us the multi-group shape natively: N_features rows in the
  # 2-group case, N_features * (k - 1) in a k-group case.
  maaslin2_results <- maaslin2_results[maaslin2_results$metadata == group, , drop = FALSE]

  # Maaslin2 replaces hyphens (-) with dots (.) in feature names, so we
  # translate Maaslin2's feature column back to the original IDs rather
  # than trying to match original IDs into Maaslin2 output (which loses
  # multi-group rows via match()'s first-hit semantics).
  original_features <- rownames(abundance_mat)
  translated_features <- original_features[match(maaslin2_results$feature,
                                                 gsub("-", ".", original_features))]
  # Fall back to the raw Maaslin2 name if the reverse map missed.
  translated_features[is.na(translated_features)] <-
    maaslin2_results$feature[is.na(translated_features)]

  results <- data.frame(
    feature = translated_features,
    method = "Maaslin2",
    group1 = ref_level,
    group2 = maaslin2_results$value,
    p_values = maaslin2_results$pval,
    log2_fold_change = maaslin2_results$coef,
    stringsAsFactors = FALSE
  )

  return(results)
}

# Helper function: Perform Lefser analysis
perform_lefser_analysis <- function(abundance_mat, metadata, group, Level) {

  # Check if we only have 2 groups (lefser requirement)
  if (length(Level) != 2) {
    stop("Lefser requires exactly 2 groups. Found ", length(Level), " groups.")
  }

  # Create SummarizedExperiment object
  se <- SummarizedExperiment::SummarizedExperiment(
    assays = list(counts = abundance_mat),
    colData = metadata
  )

  # Perform Lefser analysis to get LDA scores (effect sizes)
  lefser_results <- tryCatch({
    lefser::lefser(se, classCol = group)
  }, error = function(e) {
    message("Lefser analysis failed: ", e$message)
    NULL
  })

  # Get all features and group assignments
  all_features <- rownames(abundance_mat)
  group_vector <- metadata[[group]]

  # Calculate Wilcoxon p-values for all features
  # This is consistent with LEfSe methodology which uses Kruskal-Wallis/Wilcoxon internally
  p_values <- sapply(seq_len(nrow(abundance_mat)), function(i) {
    feature_values <- abundance_mat[i, ]
    group1_values <- feature_values[group_vector == Level[1]]
    group2_values <- feature_values[group_vector == Level[2]]

    # Handle edge cases
    if (length(unique(c(group1_values, group2_values))) <= 1) {
      return(1.0)  # No variation, not significant
    }

    tryCatch({
      wilcox.test(group1_values, group2_values, exact = FALSE)$p.value
    }, error = function(e) {
      1.0  # Return 1.0 if test fails
    })
  })

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
    feature_ids <- rownames(comparison_df)
    if (is.null(feature_ids)) {
      feature_ids <- paste0("feature_", seq_len(n_features))
    }

    pvals <- comparison_df$pvalue
    if (is.null(pvals)) {
      warning("LinDA output missing 'pvalue' column; filling with NA for this comparison.")
      pvals <- rep(NA_real_, n_features)
    } else if (length(pvals) != n_features) {
      warning(
        "LinDA output has inconsistent 'pvalue' length (",
        length(pvals), " vs ", n_features, "); filling with NA for this comparison."
      )
      pvals <- rep(NA_real_, n_features)
    }

    lfc <- comparison_df$log2FoldChange
    if (is.null(lfc)) {
      warning("LinDA output missing 'log2FoldChange' column; filling with NA for this comparison.")
      lfc <- rep(NA_real_, n_features)
    } else if (length(lfc) != n_features) {
      warning(
        "LinDA output has inconsistent 'log2FoldChange' length (",
        length(lfc), " vs ", n_features, "); filling with NA for this comparison."
      )
      lfc <- rep(NA_real_, n_features)
    }

    # Create a results data frame for this comparison
    results_list[[length(results_list) + 1]] <- data.frame(
      feature = feature_ids,
      method = "LinDA",
      group1 = reference,
      group2 = comparison_group,
      p_values = as.numeric(pvals),
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

perform_linda_analysis <- function(abundance, metadata, group, reference, Level, length_Level) {
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

  formula <- paste0("~ ", group)

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
      p.adj.method = "BH",
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
      group = group,
      reference = reference,
      Level = Level
    )
  }, error = function(e) {
    # Catch and report errors, but return an empty data frame instead of NULL
    message(sprintf("Error in LinDA analysis: %s", e$message))
    create_empty_daa_result(include_log2fc = TRUE)
  })
  
  return(linda_result)
}
