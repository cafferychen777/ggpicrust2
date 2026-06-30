#' Compare Metagenome Results
#'
#' @param metagenomes A list of metagenome matrices with rows as KOs and columns
#' as samples. Each matrix must have unique, non-empty feature row names and
#' sample column names, and finite non-negative numeric abundance values.
#' Each matrix in the list should correspond to a different metagenome.
#' @param names A unique, non-empty character vector of names for the
#' metagenomes in the same order as in the `metagenomes` list.
#' @param daa_method Character. Paired differential abundance method.
#'   Choices are \code{"ALDEx2"} (paired ALDEx2 t and Wilcoxon tests) or
#'   \code{"paired Wilcoxon"} (paired Wilcoxon signed-rank tests on
#'   sample-wise relative abundances). Methods that model the metagenomes as
#'   independent groups are not supported because all matrices represent the
#'   same aligned biological samples.
#' @param p_adjust_method A character specifying the method for p-value adjustment.
#' Possible choices are: "BH" (Benjamini-Hochberg), "holm", "bonferroni", "hochberg", "fdr", and "none".
#' The default is "BH".
#' @param reference Optional metagenome name used as the first group in
#'   pairwise DAA comparisons. Other metagenome pairs are still compared.
#' @param p.adjust Deprecated. Use \code{p_adjust_method} instead.
#' @param correlation_permutations Non-negative integer. Number of joint
#'   sample-label permutations used to test the median per-feature Spearman
#'   correlation. Use 0 to skip correlation p-values. Default 999.
#' @param correlation_seed Non-negative integer used to generate correlation
#'   permutations reproducibly. The caller's random-number state is restored.
#' @param correlation_p_adjust_method P-value adjustment method for the unique
#'   off-diagonal metagenome correlation tests. Default \code{"BH"}.
#' @name compare_metagenome_results
#' @return A list containing three elements:
#' \itemize{
#' \item "daa": paired differential abundance results for every metagenome pair.
#' \item "correlation": a list containing \code{cor_matrix},
#' \code{p_matrix}, \code{p_adjust_matrix}, and \code{n_features_matrix}.
#' Correlations are medians of finite per-feature Spearman correlations.
#' P-values use joint sample-label permutation and are adjusted across unique
#' off-diagonal metagenome pairs. Diagonal p-values are \code{NA}.
#' \item "heatmap": a ComplexHeatmap object visualizing the correlation matrix.
#' Use \code{print()} or \code{draw()} to display it.
#' }
#'
#' @details
#' Metagenome matrices are aligned to the same feature and sample identifiers.
#' Because each matrix measures the same biological samples, DAA must retain
#' this pairing. Independent-group DAA methods would treat repeated
#' measurements as independent replicates and are therefore rejected.
#'
#' Correlation inference permutes the sample columns of one metagenome jointly
#' across all features. This preserves within-metagenome feature dependence
#' while breaking only the cross-metagenome sample correspondence. Monte Carlo
#' p-values use the plus-one correction \code{(b + 1) / (B + 1)}.
#'
#' @examples
#' \donttest{
#' library(dplyr)
#' library(ComplexHeatmap)
#' # Generate example data
#' set.seed(123)
#' # First metagenome
#' metagenome1 <- abs(matrix(rnorm(1000), nrow = 100, ncol = 10))
#' rownames(metagenome1) <- paste0("KO", 1:100)
#' colnames(metagenome1) <- paste0("sample", 1:10)
#' # Second metagenome
#' metagenome2 <- abs(matrix(rnorm(1000), nrow = 100, ncol = 10))
#' rownames(metagenome2) <- paste0("KO", 1:100)
#' colnames(metagenome2) <- paste0("sample", 1:10)
#' # Put the metagenomes into a list
#' metagenomes <- list(metagenome1, metagenome2)
#' # Define names
#' names <- c("metagenome1", "metagenome2")
#' # Call the function
#' results <- compare_metagenome_results(
#'   metagenomes,
#'   names,
#'   daa_method = "paired Wilcoxon",
#'   correlation_permutations = 99
#' )
#' # Print the correlation matrix
#' print(results$correlation$cor_matrix)
#' # Display the heatmap
#' print(results$heatmap)
#' }
#' @references
#' Fernandes AD, Macklaim JM, Linn TG, Reid G, Gloor GB. Unifying the
#' analysis of high-throughput sequencing datasets: characterizing RNA-seq,
#' 16S rRNA gene sequencing and selective growth experiments by compositional
#' data analysis. Microbiome. 2014;2:15.
#'
#' Phipson B, Smyth GK. Permutation P-values Should Never Be Zero:
#' Calculating Exact P-values When Permutations Are Randomly Drawn.
#' Statistical Applications in Genetics and Molecular Biology. 2010;9(1).
#' @importFrom stats median wilcox.test
#' @export
compare_metagenome_results <- function(metagenomes, names, daa_method = "ALDEx2",
                                       p_adjust_method = "BH", reference = NULL,
                                       p.adjust = NULL,
                                       correlation_permutations = 999,
                                       correlation_seed = 123,
                                       correlation_p_adjust_method = "BH") {
  # Backward compatibility for deprecated parameter
  if (!is.null(p.adjust)) {
    warning("'p.adjust' parameter is deprecated. Use 'p_adjust_method' instead.", call. = FALSE)
    p_adjust_method <- p.adjust
  }
  validate_p_adjust_method(p_adjust_method)
  validate_p_adjust_method(correlation_p_adjust_method,
                           "correlation_p_adjust_method")
  validate_choice(daa_method, c("ALDEx2", "paired Wilcoxon"), "daa_method")
  validate_count_parameter(correlation_permutations,
                           "correlation_permutations",
                           allow_zero = TRUE)
  validate_count_parameter(correlation_seed, "correlation_seed",
                           allow_zero = TRUE)
  correlation_permutations <- as.integer(correlation_permutations)
  correlation_seed <- as.integer(correlation_seed)

  if (!is.list(metagenomes) || length(metagenomes) < 2) {
    stop("'metagenomes' must be a list containing at least two metagenome matrices.")
  }
  if (length(metagenomes) != length(names)) {
    stop("The length of 'metagenomes' must match the length of 'names'")
  }
  if (!is.character(names)) {
    stop("'names' must be a character vector.", call. = FALSE)
  }
  names <- validate_nonempty_character_column(names, "names", "names")
  if (anyDuplicated(names)) {
    duplicated_names <- unique(names[duplicated(names)])
    stop("'names' contains duplicated labels: ",
         paste(utils::head(duplicated_names, 5), collapse = ", "),
         call. = FALSE)
  }
  if (!is.null(reference)) {
    if (!is.character(reference) || length(reference) != 1 ||
        is.na(reference) || !nzchar(trimws(reference))) {
      stop("'reference' must be NULL or a single non-empty metagenome name.",
           call. = FALSE)
    }
    if (!reference %in% names) {
      stop("'reference' must match one of 'names': ",
           paste(names, collapse = ", "), ".", call. = FALSE)
    }
  }

  # Align every metagenome on the shared feature set AND the shared sample
  # set before anything else. Both downstream steps index by position:
  #
  #   * `cbind()` concatenates feature-by-position for the DAA matrix.
  #   * the per-feature Spearman correlation indexes
  #     `metagenomes[[i]][k, ]` vs `metagenomes[[j]][k, ]` -- i.e. compares
  #     a feature's abundance pattern *across samples* between two
  #     metagenomes, which is only biologically meaningful if column k
  #     refers to the same biological sample in both matrices.
  #
  # So two alignments are needed in parallel:
  #   1. Intersect row names (features). Two metagenomes with identical
  #      feature names in different orders used to produce correlations
  #      that could flip sign vs the correct by-name alignment.
  #   2. Intersect column names (samples). Two metagenomes with identical
  #      sample names in different orders used to correlate "sample 1 of
  #      metagenome A" against "sample 5 of metagenome B" and produce
  #      meaningless negative correlations. Mismatched column counts also
  #      used to fall through to `stats::cor()` and abort mid-loop with
  #      "incompatible dimensions" instead of failing at the boundary.
  metagenomes <- lapply(seq_along(metagenomes), function(i) {
    m <- metagenomes[[i]]
    context <- sprintf("metagenomes[[%s]]", names[i])
    if (!is.matrix(m) && !is.data.frame(m)) {
      stop(context, " must be a matrix or data frame.", call. = FALSE)
    }
    m <- as.matrix(m)
    validate_feature_rownames(m, context)
    sample_ids <- colnames(m)
    if (is.null(sample_ids)) {
      stop(context, " must have column names (sample identifiers) ",
           "so samples can be aligned across metagenomes.",
           call. = FALSE)
    }
    sample_ids <- validate_nonempty_character_column(sample_ids,
                                                     "colnames",
                                                     context)
    if (anyDuplicated(sample_ids)) {
      duplicated_samples <- unique(sample_ids[duplicated(sample_ids)])
      stop(context, " column names contain duplicated sample identifiers: ",
           paste(utils::head(duplicated_samples, 5), collapse = ", "),
           ".",
           call. = FALSE)
    }
    validate_nonnegative_finite_matrix(m, context, check_duplicates = FALSE)
    m
  })

  require_package("circlize", "comparison heatmap")
  require_package("ComplexHeatmap", "comparison heatmap")

  if (any(vapply(metagenomes, function(m) is.null(rownames(m)), logical(1)))) {
    stop("Every element of 'metagenomes' must have row names (feature identifiers) ",
         "so features can be aligned across metagenomes.")
  }
  if (any(vapply(metagenomes, function(m) is.null(colnames(m)), logical(1)))) {
    stop("Every element of 'metagenomes' must have column names (sample identifiers) ",
         "so samples can be aligned across metagenomes.")
  }
  shared_features <- Reduce(intersect, lapply(metagenomes, rownames))
  if (length(shared_features) == 0) {
    stop("No shared feature identifiers (row names) across the provided metagenomes.")
  }
  shared_samples <- Reduce(intersect, lapply(metagenomes, colnames))
  if (length(shared_samples) == 0) {
    stop("No shared sample identifiers (column names) across the provided metagenomes. ",
         "Per-sample correlation requires the same biological samples to appear ",
         "in every metagenome.")
  }
  if (length(shared_samples) < 3) {
    stop(
      "At least 3 shared samples are required for per-feature Spearman ",
      "correlation and paired metagenome comparison; found ",
      length(shared_samples), ".",
      call. = FALSE
    )
  }
  # Warn if alignment drops samples. A partial overlap is almost always a
  # user mistake (the function is designed for parallel quantifications on
  # the same samples); proceeding silently would hide it.
  for (nm_i in seq_along(metagenomes)) {
    dropped <- setdiff(colnames(metagenomes[[nm_i]]), shared_samples)
    if (length(dropped) > 0) {
      warning(sprintf(
        "Metagenome '%s': %d sample(s) dropped by cross-metagenome intersection (%s).",
        names[nm_i], length(dropped),
        paste(utils::head(dropped, 5),
              collapse = ", ")),
        call. = FALSE)
    }
  }
  metagenomes <- lapply(metagenomes,
                        function(m) m[shared_features, shared_samples, drop = FALSE])

  # Perform pairwise DAA while retaining the shared-sample pairing.
  daa_results <- run_paired_metagenome_daa(
    metagenomes = metagenomes,
    names = names,
    daa_method = daa_method,
    p_adjust_method = p_adjust_method,
    reference = reference
  )

  # Compute per-feature Spearman correlations between metagenomes. Safe
  # to index by both row position (shared feature order) and column
  # position (shared sample order) now that every metagenome has been
  # aligned in both dimensions above.
  n_metagenomes <- length(names)
  cor_matrix <- diag(1, nrow = n_metagenomes, ncol = n_metagenomes)
  p_matrix <- matrix(NA_real_, nrow = n_metagenomes, ncol = n_metagenomes)
  p_adjust_matrix <- matrix(NA_real_, nrow = n_metagenomes,
                            ncol = n_metagenomes)
  n_features_matrix <- matrix(NA_integer_, nrow = n_metagenomes,
                              ncol = n_metagenomes)
  spearman_inputs <- lapply(metagenomes, prepare_spearman_rank_matrix)
  diag(n_features_matrix) <- vapply(
    spearman_inputs,
    function(x) sum(x$row_norm > 0),
    integer(1)
  )
  permutations <- generate_sample_permutations(
    n_samples = length(shared_samples),
    n_permutations = correlation_permutations,
    seed = correlation_seed
  )

  for (i in seq_len(n_metagenomes - 1L)) {
    for (j in seq.int(i + 1L, n_metagenomes)) {
      pair_result <- summarize_spearman_pair(
        spearman_inputs[[i]],
        spearman_inputs[[j]]
      )
      if (pair_result$n_features == 0) {
        stop(sprintf(
          "No finite per-feature Spearman correlations could be computed for metagenomes '%s' and '%s'. This usually means every shared feature is constant across the aligned samples in at least one metagenome.",
          names[i], names[j]
        ), call. = FALSE)
      }
      cor_matrix[i, j] <- pair_result$median_correlation
      cor_matrix[j, i] <- pair_result$median_correlation
      n_features_matrix[i, j] <- pair_result$n_features
      n_features_matrix[j, i] <- pair_result$n_features

      if (correlation_permutations > 0) {
        p_value <- permutation_spearman_p_value(
          spearman_inputs[[i]],
          spearman_inputs[[j]],
          pair_result = pair_result,
          permutations = permutations
        )
        p_matrix[i, j] <- p_value
        p_matrix[j, i] <- p_value
      }
    }
  }

  pair_indices <- which(upper.tri(p_matrix), arr.ind = TRUE)
  raw_pair_p <- p_matrix[pair_indices]
  finite_pair_p <- is.finite(raw_pair_p)
  adjusted_pair_p <- rep(NA_real_, length(raw_pair_p))
  adjusted_pair_p[finite_pair_p] <- stats::p.adjust(
    raw_pair_p[finite_pair_p],
    method = correlation_p_adjust_method
  )
  for (k in seq_len(nrow(pair_indices))) {
    i <- pair_indices[k, 1]
    j <- pair_indices[k, 2]
    p_adjust_matrix[i, j] <- adjusted_pair_p[k]
    p_adjust_matrix[j, i] <- adjusted_pair_p[k]
  }

  rownames(cor_matrix) <- names
  colnames(cor_matrix) <- names
  rownames(p_matrix) <- names
  colnames(p_matrix) <- names
  rownames(p_adjust_matrix) <- names
  colnames(p_adjust_matrix) <- names
  rownames(n_features_matrix) <- names
  colnames(n_features_matrix) <- names

  cor_results <- list(
    cor_matrix = cor_matrix,
    p_matrix = p_matrix,
    p_adjust_matrix = p_adjust_matrix,
    n_features_matrix = n_features_matrix,
    test = "Joint sample-label permutation of median per-feature Spearman correlation"
  )

  # Build heatmap object (return it, don't print — let the user decide)
  color_mapping <- circlize::colorRamp2(c(-1, 0, 1), c("#0571b0", "white", "#ca0020"))
  heatmap_obj <- ComplexHeatmap::Heatmap(cor_matrix, name = "correlation",
    show_row_names = TRUE, show_column_names = TRUE, col = color_mapping)

  list(daa = daa_results, correlation = cor_results, heatmap = heatmap_obj)
}

#' Run paired DAA across metagenome pairs
#'
#' @noRd
run_paired_metagenome_daa <- function(metagenomes, names, daa_method,
                                      p_adjust_method, reference = NULL) {
  comparison_order <- seq_along(names)
  if (!is.null(reference)) {
    reference_index <- match(reference, names)
    comparison_order <- c(reference_index,
                          setdiff(comparison_order, reference_index))
  }
  comparisons <- utils::combn(comparison_order, 2, simplify = FALSE)

  results <- lapply(comparisons, function(comparison) {
    group1_index <- comparison[1]
    group2_index <- comparison[2]
    if (identical(daa_method, "ALDEx2")) {
      run_paired_aldex2_comparison(
        metagenomes[[group1_index]],
        metagenomes[[group2_index]],
        group1 = names[group1_index],
        group2 = names[group2_index],
        p_adjust_method = p_adjust_method
      )
    } else {
      run_paired_wilcoxon_comparison(
        metagenomes[[group1_index]],
        metagenomes[[group2_index]],
        group1 = names[group1_index],
        group2 = names[group2_index],
        p_adjust_method = p_adjust_method
      )
    }
  })

  do.call(rbind, results)
}

#' Run paired Wilcoxon tests on relative abundance
#'
#' @noRd
run_paired_wilcoxon_comparison <- function(group1_abundance,
                                           group2_abundance,
                                           group1,
                                           group2,
                                           p_adjust_method) {
  group1_relative <- compute_relative_abundance(
    group1_abundance,
    context = paste0("metagenome '", group1, "'")
  )
  group2_relative <- compute_relative_abundance(
    group2_abundance,
    context = paste0("metagenome '", group2, "'")
  )

  test_results <- lapply(seq_len(nrow(group1_relative)), function(i) {
    values1 <- group1_relative[i, ]
    values2 <- group2_relative[i, ]
    paired_difference <- values2 - values1

    if (all(paired_difference == 0)) {
      p_value <- 1
      statistic <- 0
    } else {
      test <- suppressWarnings(stats::wilcox.test(
        values2,
        values1,
        paired = TRUE,
        exact = FALSE
      ))
      p_value <- test$p.value
      statistic <- unname(test$statistic)
      if (!is.finite(p_value)) {
        stop(
          "Paired Wilcoxon test returned a non-finite p-value for feature '",
          rownames(group1_relative)[i], "' in comparison '", group1,
          "' vs '", group2, "'.",
          call. = FALSE
        )
      }
    }

    data.frame(
      p_values = p_value,
      statistic = statistic,
      median_difference = stats::median(paired_difference),
      log2_fold_change = calculate_log2_fold_change(
        mean(values1),
        mean(values2),
        reference_values = c(values1, values2)
      ),
      stringsAsFactors = FALSE
    )
  })
  test_results <- do.call(rbind, test_results)

  data.frame(
    feature = rownames(group1_relative),
    method = "Paired Wilcoxon signed-rank test",
    group1 = group1,
    group2 = group2,
    p_values = test_results$p_values,
    p_adjust = stats::p.adjust(test_results$p_values,
                               method = p_adjust_method),
    adj_method = p_adjust_method,
    statistic = test_results$statistic,
    median_difference = test_results$median_difference,
    log2_fold_change = test_results$log2_fold_change,
    stringsAsFactors = FALSE
  )
}

#' Run paired ALDEx2 tests for one metagenome pair
#'
#' @noRd
run_paired_aldex2_comparison <- function(group1_abundance,
                                         group2_abundance,
                                         group1,
                                         group2,
                                         p_adjust_method) {
  require_package("ALDEx2", "paired metagenome comparison")

  paired_counts <- cbind(group1_abundance, group2_abundance)
  paired_counts <- prepare_integer_counts(paired_counts, "paired ALDEx2")
  n_samples <- ncol(group1_abundance)
  colnames(paired_counts) <- c(
    paste0("group1_", colnames(group1_abundance)),
    paste0("group2_", colnames(group2_abundance))
  )
  conditions <- c(rep("group1", n_samples), rep("group2", n_samples))

  clr <- tryCatch(
    ALDEx2::aldex.clr(
      paired_counts,
      conditions,
      mc.samples = 256,
      denom = "all",
      verbose = FALSE
    ),
    error = function(e) {
      stop("Paired ALDEx2 CLR failed for '", group1, "' vs '", group2,
           "': ", conditionMessage(e), call. = FALSE)
    }
  )
  tests <- tryCatch(
    ALDEx2::aldex.ttest(clr, paired.test = TRUE, verbose = FALSE),
    error = function(e) {
      stop("Paired ALDEx2 tests failed for '", group1, "' vs '", group2,
           "': ", conditionMessage(e), call. = FALSE)
    }
  )
  effects <- tryCatch(
    ALDEx2::aldex.effect(clr, paired.test = TRUE, verbose = FALSE),
    error = function(e) {
      stop("Paired ALDEx2 effect estimation failed for '", group1,
           "' vs '", group2, "': ", conditionMessage(e), call. = FALSE)
    }
  )
  effects <- validate_and_align_aldex2_effect_results(
    effects,
    rownames(tests)
  )

  raw_p_values <- c(tests$we.ep, tests$wi.ep)
  adjusted_p_values <- if (p_adjust_method %in% c("BH", "fdr")) {
    c(tests$we.eBH, tests$wi.eBH)
  } else {
    c(
      stats::p.adjust(tests$we.ep, method = p_adjust_method),
      stats::p.adjust(tests$wi.ep, method = p_adjust_method)
    )
  }
  n_features <- nrow(tests)

  data.frame(
    feature = rep(rownames(tests), 2),
    method = c(
      rep("ALDEx2 paired t test", n_features),
      rep("ALDEx2 paired Wilcoxon signed-rank test", n_features)
    ),
    group1 = group1,
    group2 = group2,
    p_values = raw_p_values,
    p_adjust = adjusted_p_values,
    adj_method = if (p_adjust_method %in% c("BH", "fdr")) {
      "BH (ALDEx2 expected)"
    } else {
      p_adjust_method
    },
    effect_size = rep(effects$effect, 2),
    diff_btw = rep(effects$diff.btw, 2),
    log2_fold_change = rep(effects$diff.btw, 2),
    rab_all = rep(effects$rab.all, 2),
    overlap = rep(effects$overlap, 2),
    stringsAsFactors = FALSE
  )
}

#' Prepare row-wise rank matrices for Spearman correlation
#'
#' @noRd
prepare_spearman_rank_matrix <- function(abundance) {
  ranked <- t(apply(
    abundance,
    1,
    rank,
    ties.method = "average"
  ))
  if (nrow(abundance) == 1) {
    ranked <- matrix(ranked, nrow = 1,
                     dimnames = list(rownames(abundance), colnames(abundance)))
  }
  centered <- sweep(ranked, 1, rowMeans(ranked), "-")
  list(
    centered = centered,
    row_norm = sqrt(rowSums(centered^2))
  )
}

#' Summarize finite feature-wise Spearman correlations
#'
#' @noRd
summarize_spearman_pair <- function(group1_ranks, group2_ranks) {
  denominator <- group1_ranks$row_norm * group2_ranks$row_norm
  valid <- is.finite(denominator) & denominator > 0
  if (!any(valid)) {
    return(list(
      median_correlation = NA_real_,
      n_features = 0L,
      valid = valid,
      denominator = numeric(0)
    ))
  }

  feature_correlations <- rowSums(
    group1_ranks$centered[valid, , drop = FALSE] *
      group2_ranks$centered[valid, , drop = FALSE]
  ) / denominator[valid]
  feature_correlations <- pmax(-1, pmin(1, feature_correlations))

  list(
    median_correlation = stats::median(feature_correlations),
    n_features = sum(valid),
    valid = valid,
    denominator = denominator[valid]
  )
}

#' Generate reproducible sample permutations without changing caller RNG state
#'
#' @noRd
generate_sample_permutations <- function(n_samples, n_permutations, seed) {
  if (n_permutations == 0) {
    return(list())
  }

  had_random_seed <- exists(".Random.seed", envir = .GlobalEnv,
                            inherits = FALSE)
  if (had_random_seed) {
    old_random_seed <- get(".Random.seed", envir = .GlobalEnv,
                           inherits = FALSE)
  }
  on.exit({
    if (had_random_seed) {
      assign(".Random.seed", old_random_seed, envir = .GlobalEnv)
    } else if (exists(".Random.seed", envir = .GlobalEnv,
                      inherits = FALSE)) {
      rm(".Random.seed", envir = .GlobalEnv)
    }
  }, add = TRUE)

  set.seed(seed)
  lapply(seq_len(n_permutations), function(i) sample.int(n_samples))
}

#' Compute a joint sample-label permutation p-value
#'
#' @noRd
permutation_spearman_p_value <- function(group1_ranks, group2_ranks,
                                         pair_result, permutations) {
  valid <- pair_result$valid
  group1_centered <- group1_ranks$centered[valid, , drop = FALSE]
  group2_centered <- group2_ranks$centered[valid, , drop = FALSE]
  denominator <- pair_result$denominator

  permutation_statistics <- vapply(permutations, function(permutation) {
    feature_correlations <- rowSums(
      group1_centered *
        group2_centered[, permutation, drop = FALSE]
    ) / denominator
    feature_correlations <- pmax(-1, pmin(1, feature_correlations))
    stats::median(feature_correlations)
  }, numeric(1))

  tolerance <- sqrt(.Machine$double.eps)
  exceedances <- sum(
    abs(permutation_statistics) >=
      abs(pair_result$median_correlation) - tolerance
  )
  (exceedances + 1) / (length(permutations) + 1)
}
