#' Title
#'
#' @param file
#' @param pathway
#' @param daa_results_df
#'
#' @return
#' @export
#'
#' @examples
pathway_annotation <-
  function(file = NULL,
           pathway = NULL,
           daa_results_df = NULL,
           ko_to_kegg = FALSE) {
    if (is.null(file) & is.null(daa_results_df)) {
      stop("Please input the picrust2 output or results of pathway_daa daa_results_df")
    }
    if (!is.null(file)) {
      file_format <- substr(file, nchar(file) - 3, nchar(file))
      switch(
        file_format,
        ".txt" = abundance <-
          read_delim(
            file,
            delim = "\t",
            escape_double = FALSE,
            trim_ws = TRUE
          ),
        ".tsv" = abundance <-
          read_delim(
            file,
            delim = "\t",
            escape_double = FALSE,
            trim_ws = TRUE
          ),
        ".csv" = abundance <-
          read_delim(
            file,
            delim = "\t",
            escape_double = FALSE,
            trim_ws = TRUE
          ),
        stop(
          "Error: Please input file as .tsv, .txt or .csv\nThe best input file is what you get from picrust2 output file 'pred_metagenome_unstrat.tsv'"
        )
      )
      abundance <-
        abundance %>% add_column(description = rep(NA, length = nrow(abundance)),
                                 .after = 1)
      switch(
        pathway,
        "KO" = {
          load(system.file("extdata", "KO_reference.RData", package = "ggpicrust2"))
          for (i in seq_len(nrow(abundance))) {
            abundance[i, 2] <-
              KO_reference[KO_reference[, 1] %in% abundance[i, 1], 5][1]
          }
        },
        "EC" = {
          load(system.file("extdata", "EC_reference.RData", package = "ggpicrust2"))
          for (i in seq_len(nrow(abundance))) {
            abundance[i, 2] <-
              EC_reference[EC_reference[, 1] %in% abundance[i, 1], 2]
          }
          message("EC description may appear to be duplicated")
        },
        "MetaCyc" = {
          load(system.file("extdata", "MetaCyc_reference.RData", package = "ggpicrust2"))
          for (i in seq_len(nrow(abundance))) {
            abundance[i, 2] <-
              MetaCyc_reference[MetaCyc_reference[, 1] %in% abundance[i, 1], 2]
          }
        },
        stop("Only provide 'KO', 'EC' and 'MetaCyc' pathway")
      )
      return(abundance)
    }
    if (!is.null(daa_results_df)) {
      if (ko_to_kegg == FALSE) {
        daa_results_df$description <- "nonsense"
        switch(
          pathway,
          "KO" = {
            load(system.file("extdata", "KO_reference.RData", package = "ggpicrust2"))
            for (i in seq_len(nrow(daa_results_df))) {
              daa_results_df[i,]$description <-
                KO_reference[KO_reference[, 1] %in% daa_results_df[i,]$feature, 5][1]
            }
          },
          "EC" = {
            load(system.file("extdata", "EC_reference.RData", package = "ggpicrust2"))
            for (i in seq_len(nrow(daa_results_df))) {
              daa_results_df[i,]$description <-
                EC_reference[EC_reference[, 1] %in% daa_results_df[i,]$feature, 2]
            }
            message("EC description may appear to be duplicated")
          },
          "MetaCyc" = {
            load(system.file("extdata", "MetaCyc_reference.RData", package = "ggpicrust2"))
            for (i in seq_len(nrow(daa_results_df))) {
              daa_results_df[i,]$description <-
                EC_reference[EC_reference[, 1] %in% daa_results_df[i,]$feature, 2]
            }
          },
          stop("Only provide 'KO', 'EC' and 'MetaCyc' pathway")
        )
        return(daa_results_df)
      }else {daa_results_filtered_df <-
        daa_results_df[daa_results_df$p_adjust < 0.05,]
      daa_results_filtered_df$pathway_name <- "nonsense"
      daa_results_filtered_df$pathway_description <- "nonsense"
      daa_results_filtered_df$pathway_class <- "nonsense"
      daa_results_filtered_df$pathway_map <- "nonsense"
      keggGet_results <- list()
      message(
        "We are connecting to the KEGG database to get the latest results, please wait patiently."
      )
      for (i in seq_len(nrow(daa_results_filtered_df))) {
        keggGet_results[[i]] <- keggGet(daa_results_filtered_df$feature[i])
        daa_results_filtered_df[i,]$pathway_name <-
          keggGet_results[[i]][[1]]$NAME
        daa_results_filtered_df[i,]$pathway_description <-
          keggGet_results[[i]][[1]]$DESCRIPTION
        daa_results_filtered_df[i,]$pathway_class <-
          keggGet_results[[i]][[1]]$CLASS
        daa_results_filtered_df[i,]$pathway_map <-
          keggGet_results[[i]][[1]]$PATHWAY_MAP
      }
      daa_results_filtered_annotation_df <- daa_results_filtered_df
      return(daa_results_filtered_annotation_df)
          }
    }
  }
