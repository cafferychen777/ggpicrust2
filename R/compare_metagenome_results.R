#' Compare Metagenome Results
#'
#' @param metagenomes A list of metagenomes matrices with rows as KOs and columns as samples.
#' Each matrix in the list should correspond to a different metagenome.
#' @param names A vector of names for the metagenomes in the same order as in the `metagenomes` list.
#' @param daa_method A character specifying the method for differential abundance analysis (DAA).
#' Possible choices are: "ALDEx2", "DESeq2", "edgeR", "limma voom", "metagenomeSeq", "LinDA",
#' "Maaslin2", and "Lefse". The default is "ALDEx2".
#' @param p.adjust A character specifying the method for p-value adjustment.
#' Possible choices are: "BH" (Benjamini-Hochberg), "holm", "bonferroni", "hochberg", "fdr", and "none".
#' The default is "BH".
#' @param reference A character specifying the reference group level for DAA.
#' This parameter is used when there are more than two groups. The default is NULL.
#' @name compare_metagenome_results
#' @return A list containing two elements:
#' \itemize{
#' \item "daa": a list of results from the `pathway_daa` function. Each result is a data frame
#' containing the differential abundance analysis results with columns for the feature ID,
#' the test statistic, the raw p-value, and the adjusted p-value.
#' \item "correlation": a list with two elements: "cor_matrix" and "p_matrix", which are
#' matrices of Spearman correlation coefficients and their corresponding p-values, respectively,
#' between every pair of metagenomes.
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
#' # Print the p-value matrix
#' print(results$correlation$p_matrix)
#' }
#' @export
utils::globalVariables(c("cor.test","Heatmap"))
compare_metagenome_results <- function(metagenomes, names, daa_method = "ALDEx2", p.adjust = "BH", reference = NULL) {
  if(length(metagenomes) != length(names)){
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
                             p.adjust = p.adjust,
                             reference = reference)

  # Compute Spearman correlations and p-values
  cor_matrix <- matrix(NA, nrow = length(names), ncol = length(names))
  p_matrix <- matrix(NA, nrow = length(names), ncol = length(names))
  for (i in 1:length(names)) {
    for (j in i:length(names)) {
      cor_test <- cor.test(t(metagenomes[[i]]), t(metagenomes[[j]]), method = "spearman")
      cor_matrix[i, j] <- cor_test$estimate
      cor_matrix[j, i] <- cor_test$estimate
      p_matrix[i, j] <- cor_test$p.value
      p_matrix[j, i] <- cor_test$p.value
    }
  }

  # Add rownames and colnames
  rownames(cor_matrix) <- names
  colnames(cor_matrix) <- names
  rownames(p_matrix) <- names
  colnames(p_matrix) <- names

  # Create a list containing correlation matrix and p-value matrix
  cor_results <- list(cor_matrix = cor_matrix, p_matrix = p_matrix)

  color_mapping <- circlize::colorRamp2(c(-1, 0, 1), c("#0571b0", "white", "#ca0020"))
  print(Heatmap(cor_matrix, name = "correlation", show_row_names = TRUE, show_column_names = TRUE, col = color_mapping))

  # Return a list containing DAA results and correlation results
  return(list(daa = daa_results, correlation = cor_results))
}
