#' Convert KO abundance in picrust2 export files to KEGG pathway abundance
#'
#' This function takes a file containing KO abundance in picrust2 export format and converts it to KEGG pathway abundance.
#' The file should be in .tsv, .txt or .csv format.
#'
#' @param file A character, the address of the file containing KO abundance in picrust2 export format
#'
#' @return A data frame, kegg_abundance, with KEGG pathway abundance values
#'
#' @examples
#' ko2kegg_abundance(file = "path/to/pred_metagenome_unstrat.tsv")
#'
#' @export
ko2kegg_abundance <- function(file) {
  
  # Check if file has a .tsv, .txt or .csv extension
  if (!grepl("\\.(tsv|txt|csv)$", file)) {
    stop("Error: Please input file as .tsv, .txt or .csv\nThe best input file is what you get from picrust2 output file 'pred_metagenome_unstrat.tsv'")
  }
  
  # Read the file using the readr package
  abundance <- readr::read_delim(file, delim = "\t", trim_ws = TRUE)

#Load the KEGG reference data using the system.file function
kegg_reference_file <- system.file("extdata", "kegg_reference.RData", package = "ggpicrust2")
load(kegg_reference_file)

#Get the sample names and KEGG pathway names
sample_names <- colnames(abundance)[-1]
kegg_names <- ko_to_kegg_reference[, 1]

#Initialize a matrix to store the KEGG pathway abundance
kegg_abundance <- matrix(NA, nrow = nrow(kegg_names), ncol = length(sample_names))
colnames(kegg_abundance) <- sample_names
rownames(kegg_abundance) <- as.matrix(kegg_names)

#Use the dplyr package to calculate the KEGG pathway abundance
library(dplyr)
kegg_abundance <- as.data.frame(kegg_abundance)
kegg_abundance <- kegg_names %>%
mutate(abundance = map2_dbl(., sample_names, ~ {
kegg_name <- .x
sample_name <- .y
# Get the KO codes associated with the current KEGG pathway
ko_to_kegg <- ko_to_kegg_reference[ko_to_kegg_reference[, 1] == kegg_name, -1]
ko_to_kegg <- ko_to_kegg[!is.na(ko_to_kegg)]
# Calculate the abundance of the current KEGG pathway in the current sample
sum(abundance[as.matrix(abundance[, 1]) %in% ko_to_kegg, sample_name])
})) %>%
filter(abundance != 0)

#Return the updated kegg abundance data frame
return(kegg_abundance)
}
