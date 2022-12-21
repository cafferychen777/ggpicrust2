#file="/Users/apple/Microbiome/C9orf72/Code And Data/picrust2_out/KO_metagenome_out/pred_metagenome_unstrat.tsv/pred_metagenome_unstrat.tsv"
ko2kegg_abundance <- function(file){
  file_format <- substr(file, nchar(file)-3, nchar(file))
  switch(file_format,
         ".txt" = abundance <- read_delim(file, delim = "\t", escape_double = FALSE, trim_ws = TRUE),
         ".tsv" = abundance <- read_delim(file, delim = "\t", escape_double = FALSE, trim_ws = TRUE),
         ".csv" = abundance <- read_delim(file, delim = "\t", escape_double = FALSE, trim_ws = TRUE),
         stop("Error: Please input file as .tsv, .txt or .csv\nThe best input file is what you get from picrust2 output file 'pred_metagenome_unstrat.tsv'"))
  message("Calculation may take a long time, please be patient.")
  library(readr)
  library(dplyr)
  library(tidyr)
  library(tibble)
  load("/Users/apple/Microbiome/ggpicrust2/ggpicrust2/data/kegg_reference.RData")
  sample_names <- colnames(abundance)[-1]
  kegg_names <- ko_to_kegg_reference[,1]
  kegg_abundance <- matrix(NA,nrow = nrow(kegg_names),ncol= length(sample_names))
  colnames(kegg_abundance) <- sample_names
  rownames(kegg_abundance) <- as.matrix(kegg_names)
  for (i in 1:nrow(kegg_abundance)){
    for (j in 1:ncol(kegg_abundance)){
      kegg_name <- rownames(kegg_abundance)[i]
      sample_name <- colnames(kegg_abundance)[j]
      ko_to_kegg <- ko_to_kegg_reference[ko_to_kegg_reference[,1]==kegg_name,-1]
      ko_to_kegg <- ko_to_kegg[!is.na(ko_to_kegg)]
      kegg_abundance[i,j] <- sum(abundance[as.matrix(abundance[,1]) %in% ko_to_kegg,sample_name])
    }
  }
  kegg_abundance <- kegg_abundance[rowSums(kegg_abundance) != 0, ]
  message("The kegg pathway with zero abundance in all the different samples has been removed.")
  kegg_abundance <- as.data.frame(kegg_abundance)
  return(kegg_abundance)
}