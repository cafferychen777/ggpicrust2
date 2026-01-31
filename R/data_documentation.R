#' KO to KEGG Pathway Reference Data
#'
#' A comprehensive mapping between KEGG Orthology (KO) identifiers and KEGG pathways.
#' This dataset contains mappings covering 532 pathways and 23,466 unique KO IDs,
#' filtered to include only real KEGG pathway maps (5-digit IDs).
#'
#' @format A data frame with 9 variables:
#' \describe{
#'   \item{pathway_id}{KEGG pathway identifier (e.g., "ko00010")}
#'   \item{pathway_number}{KEGG pathway number}
#'   \item{pathway_name}{Full name of the pathway}
#'   \item{ko_id}{KEGG Orthology identifier (e.g., "K00001")}
#'   \item{ko_description}{Description of the KO}
#'   \item{ec_number}{EC number associated with the KO (if applicable)}
#'   \item{level1}{KEGG pathway hierarchy Level 1 classification}
#'   \item{level2}{KEGG pathway hierarchy Level 2 classification}
#'   \item{level3}{KEGG pathway hierarchy Level 3 classification}
#' }
#'
#' @details
#' This reference data is used by the \code{\link{ko2kegg_abundance}} function to convert
#' KO abundance data to KEGG pathway abundance. The data is stored internally and does not
#' require internet connectivity to use.
#'
#' The dataset covers major KEGG pathway categories including:
#' \itemize{
#'   \item Metabolism
#'   \item Genetic Information Processing
#'   \item Environmental Information Processing
#'   \item Cellular Processes
#'   \item Organismal Systems
#'   \item Human Diseases
#' }
#'
#' @source KEGG database (\url{https://www.kegg.jp/})
#'
#' @examples
#' # Load the reference data
#' data(ko_to_kegg_reference)
#'
#' # View structure
#' str(ko_to_kegg_reference)
#'
#' # Get unique pathways
#' unique_pathways <- unique(ko_to_kegg_reference$pathway_id)
#' length(unique_pathways)
#'
#' # Find KOs for a specific pathway
#' glycolysis_kos <- ko_to_kegg_reference[ko_to_kegg_reference$pathway_id == "ko00010", ]
#' head(glycolysis_kos)
#'
#' @seealso \code{\link{ko2kegg_abundance}} for converting KO abundance to pathway abundance
#'
"ko_to_kegg_reference"
