#' Convert KO abundance in picrust2 export files to KEGG pathway abundance
#'
#' This function takes a file containing KO (KEGG Orthology) abundance data in picrust2 export format and converts it to KEGG pathway abundance data.
#' The input file should be in .tsv, .txt, or .csv format.
#'
#' @param file A character string representing the file path of the input file containing KO abundance data in picrust2 export format. The input file should have KO identifiers in the first column and sample identifiers in the first row. The remaining cells should contain the abundance values for each KO-sample pair.
#'
#' @return
#' A data frame with KEGG pathway abundance values. Rows represent KEGG pathways, identified by their KEGG pathway IDs. Columns represent samples, identified by their sample IDs from the input file. Each cell contains the abundance of a specific KEGG pathway in a given sample, calculated by summing the abundances of the corresponding KOs in the input file.
#'
#' @export
#' @examples
#' \dontrun{
#' library(ggpicrust2)
#' # Prepare an example input file
#' input_file <- system.file("extdata", "example_ko_abundance.tsv", package = "ggpicrust2")
#'
#' # Run ko2kegg_abundance function
#' kegg_abundance <- ko2kegg_abundance(file = input_file)
#' }
ko2kegg_abundance <- function (file)
{
    file_format <- substr(file, nchar(file) - 3, nchar(file))
    switch(file_format, .txt = abundance <- readr::read_delim(file,
        delim = "\t", escape_double = FALSE, trim_ws = TRUE),
        .tsv = abundance <- readr::read_delim(file, delim = "\t", escape_double = FALSE,
            trim_ws = TRUE), .csv = abundance <- readr::read_delim(file,
            delim = "\t", escape_double = FALSE, trim_ws = TRUE),
        stop("Error: Please input file as .tsv, .txt or .csv\nThe best input file is what you get from picrust2 output file 'pred_metagenome_unstrat.tsv'"))
    message("Calculation may take a long time, please be patient.")
    load(system.file("extdata", "kegg_reference.RData", package = "ggpicrust2"))
    sample_names <- colnames(abundance)[-1]
    kegg_names <- ko_to_kegg_reference[, 1]
    kegg_abundance <- matrix(NA, nrow = nrow(kegg_names), ncol = length(sample_names))
    colnames(kegg_abundance) <- sample_names
    rownames(kegg_abundance) <- as.matrix(kegg_names)
    for (i in seq_len(nrow(kegg_abundance))) {
        for (j in seq_len(ncol(kegg_abundance))) {
            kegg_name <- rownames(kegg_abundance)[i]
            sample_name <- colnames(kegg_abundance)[j]
            ko_to_kegg <- ko_to_kegg_reference[ko_to_kegg_reference[,
                1] == kegg_name, -1]
            ko_to_kegg <- ko_to_kegg[!is.na(ko_to_kegg)]
            kegg_abundance[i, j] <- sum(abundance[as.matrix(abundance[,
                1]) %in% ko_to_kegg, sample_name])
        }
    }
    kegg_abundance <- kegg_abundance[rowSums(kegg_abundance) !=
        0, ]
    message("The kegg pathway with zero abundance in all the different samples has been removed.")
    kegg_abundance <- as.data.frame(kegg_abundance)
    return(kegg_abundance)
}
