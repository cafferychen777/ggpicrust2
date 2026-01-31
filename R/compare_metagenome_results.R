#' Compare Metagenome Results
#'
#' @param metagenomes A list of metagenomes matrices with rows as KOs and columns as samples.
#' Each matrix in the list should correspond to a different metagenome.
#' @param names A vector of names for the metagenomes in the same order as in the `metagenomes` list.
#' @param daa_method A character specifying the method for differential abundance analysis (DAA).
#' Possible choices are: "ALDEx2", "DESeq2", "edgeR", "limma voom", "metagenomeSeq", "LinDA",
#' "Maaslin2", and "Lefse". The default is "ALDEx2".
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

  # Concatenate metagenomes with pseudonym column names
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

  # Compute per-feature Spearman correlations between metagenomes
  # For each pair, compute correlation per feature row, then summarize with median
  n_metagenomes <- length(names)
  cor_matrix <- matrix(NA, nrow = n_metagenomes, ncol = n_metagenomes)
  p_matrix <- matrix(NA, nrow = n_metagenomes, ncol = n_metagenomes)

  for (i in seq_along(names)) {
    for (j in seq(from = i, to = n_metagenomes)) {
      # Per-feature correlation: correlate sample profiles for each feature
      feature_cors <- sapply(seq_len(nrow(metagenomes[[i]])), function(k) {
        stats::cor(metagenomes[[i]][k, ], metagenomes[[j]][k, ], method = "spearman")
      })
      cor_matrix[i, j] <- median(feature_cors, na.rm = TRUE)
      cor_matrix[j, i] <- cor_matrix[i, j]

      # Test if median correlation differs from 0
      p_matrix[i, j] <- stats::wilcox.test(feature_cors, mu = 0)$p.value
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
