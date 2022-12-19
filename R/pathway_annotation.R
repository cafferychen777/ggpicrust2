#file=.tsv or .txt. or .csv
#pathway="KO", "EC" and "MetaCyc"
#file="/Users/apple/Microbiome/C9orf72/Code And Data/picrust2_out/EC_metagenome_out/pred_metagenome_unstrat.tsv/pred_metagenome_unstrat.tsv"
#file="/Users/apple/Microbiome/C9orf72/Code And Data/picrust2_out/KO_metagenome_out/pred_metagenome_unstrat.tsv/pred_metagenome_unstrat.tsv"
#file="/Users/apple/Microbiome/C9orf72/Code And Data/picrust2_out/pathways_out/path_abun_unstrat.tsv/path_abun_unstrat.tsv"
pathway_annotation <- function(file,pathway){
  library(dplyr)
  library(tibble)
  library(readr)
  file_format <- substr(file, nchar(file)-3, nchar(file))
  switch(file_format,
         ".txt" = abundance <- read_delim(file, delim = "\t", escape_double = FALSE, trim_ws = TRUE),
         ".tsv" = abundance <- read_delim(file, delim = "\t", escape_double = FALSE, trim_ws = TRUE),
         ".csv" = abundance <- read_delim(file, delim = "\t", escape_double = FALSE, trim_ws = TRUE),
          stop("Error: Please input file as .tsv, .txt or .csv\nThe best input file is what you get from picrust2 output file 'pred_metagenome_unstrat.tsv'"))
  abundance <- abundance %>% add_column(description = rep(NA,length = nrow(abundance)),.after = 1)
  switch(pathway,
         "KO" = {
           load("/Users/apple/Microbiome/ggpicrust2/ggpicrust2/data/KO_reference.RData")
           for (i in 1:nrow(abundance)) {
             abundance[i, 2] <-
               KO_reference[KO_reference[, 1] %in% abundance[i, 1], 5][1]
           }
         }
         ,
         "EC" = {
           load("/Users/apple/Microbiome/ggpicrust2/ggpicrust2/data/EC_reference.RData")
           for (i in 1:nrow(abundance)) {
             abundance[i, 2] <-
               EC_reference[EC_reference[, 1] %in% abundance[i, 1], 2]
           }
           message("EC description may appear to be duplicated")
         },
         "MetaCyc" = {
           load("/Users/apple/Microbiome/ggpicrust2/ggpicrust2/data/MetaCyc_reference.RData")
           for (i in 1:nrow(abundance)) {
             abundance[i, 2] <-
               MetaCyc_reference[MetaCyc_reference[, 1] %in% abundance[i, 1], 2]
           }
         },
         stop("Only provide 'KO', 'EC' and 'MetaCyc' pathway"))
  return(abundance)
}



