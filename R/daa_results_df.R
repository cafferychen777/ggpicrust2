#' DAA Results Dataset
#'
#' This dataset is the result of processing 'kegg_abundance' through the 'LinDA' method in the 'pathway_daa' function.
#' It includes information about the feature, groups compared, p values, and method used.
#'
#' @format A data frame with columns:
#' \describe{
#'   \item{adj_method}{Method used for p-value adjustment.}
#'   \item{feature}{The feature (pathway) being compared.}
#'   \item{group1}{The first group in the comparison.}
#'   \item{group2}{The second group in the comparison.}
#'   \item{method}{The method used for the comparison.}
#'   \item{p_adjust}{The adjusted p-value from the comparison.}
#'   \item{p_values}{The raw p-value from the comparison.}
#' }
#' @source From ggpicrust2 package demonstration.
#' @references Douglas GM, Maffei VJ, Zaneveld J, Yurgel SN, Brown JR, Taylor CM, Huttenhower C, Langille MGI. PICRUSt2 for prediction of metagenome functions. Nat Biotechnol. 2020.
"daa_results_df"
