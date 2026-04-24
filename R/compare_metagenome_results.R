#' Compare Metagenome Results
#'
#' @param metagenomes A list of metagenomes matrices with rows as KOs and columns as samples.
#' Each matrix in the list should correspond to a different metagenome.
#' @param names A vector of names for the metagenomes in the same order as in the `metagenomes` list.
#' @param daa_method A character specifying the method for differential abundance analysis (DAA).
#' Possible choices are: "ALDEx2", "DESeq2", "edgeR", "limma voom", "metagenomeSeq", "LinDA",
#' "Maaslin2", and "Lefser". The default is "ALDEx2".
#' @param p_adjust_method A character specifying the method for p-value adjustment.
#' Possible choices are: "BH" (Benjamini-Hochberg), "holm", "bonferroni", "hochberg", "fdr", and "none".
#' The default is "BH".
#' @param reference A character specifying the reference group level for DAA.
#' This parameter is used when there are more than two groups. The default is NULL.
#' @param p.adjust Deprecated. Use \code{p_adjust_method} instead.
#' @name compare_metagenome_results
#' @return A list containing three elements:
#' \itemize{
#' \item "daa": a data frame of results from the `pathway_daa` function
#' containing the differential abundance analysis results.
#' \item "correlation": a list with two elements: "cor_matrix" and "p_matrix", which are
#' matrices of per-feature median Spearman correlation coefficients and their
#' corresponding p-values, respectively, between every pair of metagenomes.
#' \item "heatmap": a ComplexHeatmap object visualizing the correlation matrix.
#' Use \code{print()} or \code{draw()} to display it.
#' }
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
#' results <- compare_metagenome_results(metagenomes, names, daa_method = "LinDA")
#' # Print the correlation matrix
#' print(results$correlation$cor_matrix)
#' # Display the heatmap
#' print(results$heatmap)
#' }
#' @importFrom stats cor median wilcox.test
#' @export
compare_metagenome_results <- function(metagenomes, names, daa_method = "ALDEx2",
                                       p_adjust_method = "BH", reference = NULL,
                                       p.adjust = NULL) {
  # Backward compatibility for deprecated parameter
  if (!is.null(p.adjust)) {
    warning("'p.adjust' parameter is deprecated. Use 'p_adjust_method' instead.", call. = FALSE)
    p_adjust_method <- p.adjust
  }

  if (length(metagenomes) != length(names)) {
    stop("The length of 'metagenomes' must match the length of 'names'")
  }
  require_package("circlize", "comparison heatmap")
  require_package("ComplexHeatmap", "comparison heatmap")

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

  # Concatenate metagenomes with pseudonym column names. After the
  # alignment above, every matrix has identical rownames AND identical
  # colnames in the same order, so cbind() is guaranteed to stack the
  # same feature across columns, and the downstream per-feature
  # correlation compares the same biological samples across matrices.
  pseudonym_metagenomes <- lapply(seq_along(metagenomes), function(i) {
    new_metagenome <- metagenomes[[i]]
    colnames(new_metagenome) <- paste0(names[i], "_", colnames(new_metagenome))
    return(new_metagenome)
  })
  concatenated_metagenomes <- do.call(cbind, pseudonym_metagenomes)

  # Create a new metadata with original sample names
  new_metadata <- data.frame(sample = colnames(concatenated_metagenomes),
                             group = rep(names, sapply(metagenomes, ncol)))

  # Perform DAA
  daa_results <- pathway_daa(abundance = concatenated_metagenomes,
                             metadata = new_metadata,
                             group = "group",
                             daa_method = daa_method,
                             p_adjust_method = p_adjust_method,
                             reference = reference)

  # Compute per-feature Spearman correlations between metagenomes. Safe
  # to index by both row position (shared feature order) and column
  # position (shared sample order) now that every metagenome has been
  # aligned in both dimensions above.
  n_metagenomes <- length(names)
  cor_matrix <- matrix(NA, nrow = n_metagenomes, ncol = n_metagenomes)
  p_matrix <- matrix(NA, nrow = n_metagenomes, ncol = n_metagenomes)

  for (i in seq_along(names)) {
    for (j in seq(from = i, to = n_metagenomes)) {
      feature_cors <- vapply(seq_along(shared_features), function(k) {
        stats::cor(metagenomes[[i]][k, ], metagenomes[[j]][k, ], method = "spearman")
      }, numeric(1))
      cor_matrix[i, j] <- median(feature_cors, na.rm = TRUE)
      cor_matrix[j, i] <- cor_matrix[i, j]

      # Test if the per-feature correlations differ from 0 as a group.
      # `exact = FALSE` forces the normal approximation instead of letting
      # wilcox.test() try (and fail on ties) to compute an exact p-value.
      # Ties are endemic to Spearman correlations -- any pair of features
      # with the same rank pattern yields identical coefficients -- so the
      # default `exact = NULL` would spam "cannot compute exact p-value
      # with ties" warnings that are purely an internal implementation
      # detail and would drown out user-facing warnings (e.g. dropped
      # samples from the cross-metagenome intersection just above). The
      # numerical p-value is unchanged whenever ties exist, and this
      # function's typical n (feature count, usually thousands) puts the
      # normal approximation well within its accurate regime anyway.
      p_matrix[i, j] <- stats::wilcox.test(feature_cors, mu = 0,
                                           exact = FALSE)$p.value
      p_matrix[j, i] <- p_matrix[i, j]
    }
  }

  rownames(cor_matrix) <- names
  colnames(cor_matrix) <- names
  rownames(p_matrix) <- names
  colnames(p_matrix) <- names

  cor_results <- list(cor_matrix = cor_matrix, p_matrix = p_matrix)

  # Build heatmap object (return it, don't print — let the user decide)
  color_mapping <- circlize::colorRamp2(c(-1, 0, 1), c("#0571b0", "white", "#ca0020"))
  heatmap_obj <- ComplexHeatmap::Heatmap(cor_matrix, name = "correlation",
    show_row_names = TRUE, show_column_names = TRUE, col = color_mapping)

  list(daa = daa_results, correlation = cor_results, heatmap = heatmap_obj)
}
