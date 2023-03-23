#' Pathway information annotation of "EC", "KO", "MetaCyc" pathway
#'
#' @param file A character, address to store picrust2 export files
#' @param pathway A character, consisting of "EC", "KO", "MetaCyc"
#' @param daa_results_df A data frame, output of pathway_daa
#' @param ko_to_kegg A character, decide if convert ko abundance to kegg pathway abundance
#'
#' @return A data frame containing pathway annotation information. The data frame has the following columns:
#' \itemize{
#'   \item \code{feature}: The feature ID of the pathway (e.g., KO, EC, or MetaCyc ID).
#'   \item \code{description}: The description or name of the pathway.
#'   \item Other columns depending on the input parameters and type of pathway.
#' }
#' If \code{ko_to_kegg} is set to TRUE, the output data frame will also include the following columns:
#' \itemize{
#'   \item \code{pathway_name}: The name of the KEGG pathway.
#'   \item \code{pathway_description}: The description of the KEGG pathway.
#'   \item \code{pathway_class}: The class of the KEGG pathway.
#'   \item \code{pathway_map}: The KEGG pathway map ID.
#' }
#'
#' @export
#'
#' @examples
#' \dontrun{
#' # Prepare the required input files and data frames
#' # Then, you can use the pathway_annotation function as follows:
#'
#' result <- pathway_annotation(file = "path/to/picrust2/export/file.txt",
#'                              pathway = "KO",
#'                              daa_results_df = NULL,
#'                              ko_to_kegg = FALSE)
#'}
pathway_annotation <-
  function(file = NULL,
           pathway = NULL,
           daa_results_df = NULL,
           ko_to_kegg = FALSE) {
    if (is.null(file) && is.null(daa_results_df)) {
      stop("Please input the picrust2 output or results of pathway_daa daa_results_df")
    }
    if (!is.null(file)) {
      file_format <- substr(file, nchar(file) - 3, nchar(file))
      switch(file_format,
        ".txt" = abundance <-
          readr::read_delim(
            file,
            delim = "\t",
            escape_double = FALSE,
            trim_ws = TRUE
          ),
        ".tsv" = abundance <-
          readr::read_delim(
            file,
            delim = "\t",
            escape_double = FALSE,
            trim_ws = TRUE
          ),
        ".csv" = abundance <-
          readr::read_delim(
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
        abundance %>% tibble::add_column(
          description = rep(NA, length = nrow(abundance)),
          .after = 1
        )
      switch(pathway,
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
        daa_results_df$description <- NA
        switch(pathway,
          "KO" = {
            load(system.file("extdata", "KO_reference.RData", package = "ggpicrust2"))
            for (i in seq_len(nrow(daa_results_df))) {
              daa_results_df[i, ]$description <-
                KO_reference[KO_reference[, 1] %in% daa_results_df[i, ]$feature, 5][1]
            }
          },
          "EC" = {
            load(system.file("extdata", "EC_reference.RData", package = "ggpicrust2"))
            for (i in seq_len(nrow(daa_results_df))) {
              daa_results_df[i, ]$description <-
                EC_reference[EC_reference[, 1] %in% daa_results_df[i, ]$feature, 2]
            }
            message("EC description may appear to be duplicated")
          },
          "MetaCyc" = {
            load(system.file("extdata", "MetaCyc_reference.RData", package = "ggpicrust2"))
            for (i in seq_len(nrow(daa_results_df))) {
              daa_results_df[i, ]$description <-
                EC_reference[EC_reference[, 1] %in% daa_results_df[i, ]$feature, 2]
            }
          },
          stop("Only provide 'KO', 'EC' and 'MetaCyc' pathway")
        )
        return(daa_results_df)
      } else {
        daa_results_filtered_df <-
          daa_results_df[daa_results_df$p_adjust < 0.05, ]
        if (nrow(daa_results_filtered_df) == 0 ){
          stop("There are no statistically significant biomarker")
        }
        daa_results_filtered_df$pathway_name <- NA
        daa_results_filtered_df$pathway_description <- NA
        daa_results_filtered_df$pathway_class <- NA
        daa_results_filtered_df$pathway_map <- NA
        keggGet_results <- list()
        message(
          "We are connecting to the KEGG database to get the latest results, please wait patiently."
        )
        if (nrow(daa_results_filtered_df) > 100) {
          message(
            "The pathways with statistically significance are too many. The database cannot execute the query request in a single time. Please seperate it."
          )
        }
        if (nrow(daa_results_filtered_df) <= 10) {
          for (i in seq_len(nrow(daa_results_filtered_df))) {
            a <- 0
            repeat {
              tryCatch(
                {
                  keggGet_results[[i]] <-
                    KEGGREST::keggGet(daa_results_filtered_df$feature[i])
                  a <- 1
                },
                error = function(e) {
                }
              )
              if (a == 1) {
                break
              }
            }
            daa_results_filtered_df[i, ]$pathway_name <-
              keggGet_results[[i]][[1]]$NAME
            daa_results_filtered_df[i, ]$pathway_description <-
              keggGet_results[[i]][[1]]$DESCRIPTION
            daa_results_filtered_df[i, ]$pathway_class <-
              keggGet_results[[i]][[1]]$CLASS
            daa_results_filtered_df[i, ]$pathway_map <-
              keggGet_results[[i]][[1]]$PATHWAY_MAP
          }
        }
        if (nrow(daa_results_filtered_df) > 10 &&
          nrow(daa_results_filtered_df) < 99) {
          n <-
            length(c(seq(
              10, nrow(daa_results_filtered_df), 10
            ), nrow(daa_results_filtered_df)))
          j <- 1
          seq <-
            c(seq(10, nrow(daa_results_filtered_df), 10), nrow(daa_results_filtered_df))
          for (i in seq) {
            if (i %% 10 == 0) {
              keggGet_results[[j]] <-
                KEGGREST::keggGet(daa_results_filtered_df$feature[seq(i -
                  9, i, 1)])
              # for (k in seq(i - 9, i, 1)) {
              #   if (k %% 10 == 0) {
              #     daa_results_filtered_df[k, ]$pathway_name <-
              #       keggGet_results[[j]][[10]]$NAME
              #     daa_results_filtered_df[i, ]$pathway_description <-
              #       keggGet_results[[j]][[10]]$DESCRIPTION
              #     daa_results_filtered_df[i, ]$pathway_class <-
              #       keggGet_results[[j]][[10]]$CLASS
              #     daa_results_filtered_df[i, ]$pathway_map <-
              #       keggGet_results[[j]][[10]]$PATHWAY_MAP
              #   } else{
              #     daa_results_filtered_df[k, ]$pathway_name <-
              #       keggGet_results[[j]][[k %% 10]]$NAME
              #     daa_results_filtered_df[i, ]$pathway_description <-
              #       keggGet_results[[j]][[k %% 10]]$DESCRIPTION
              #     daa_results_filtered_df[i, ]$pathway_class <-
              #       keggGet_results[[j]][[k %% 10]]$CLASS
              #     daa_results_filtered_df[i, ]$pathway_map <-
              #       keggGet_results[[j]][[k %% 10]]$PATHWAY_MAP
              #   }
              # }
            } else {
              keggGet_results[[j]] <-
                KEGGREST::keggGet(daa_results_filtered_df$feature[seq(nrow(daa_results_filtered_df) %/% 10 *
                  10 + 1, i, 1)])
              # for (k in seq(nrow(daa_results_filtered_df) %/% 10 * 10 +
              #               1, i, 1)) {
              #   daa_results_filtered_df[k, ]$pathway_name <-
              #     keggGet_results[[j]][[k %% 10]]$NAME
              #   daa_results_filtered_df[i, ]$pathway_description <-
              #     keggGet_results[[j]][[k %% 10]]$DESCRIPTION
              #   daa_results_filtered_df[i, ]$pathway_class <-
              #     keggGet_results[[j]][[k %% 10]]$CLASS
              #   daa_results_filtered_df[i, ]$pathway_map <-
              #     keggGet_results[[j]][[k %% 10]]$PATHWAY_MAP
              # }
            }
            j <- j + 1
          }
          for (k in 1:n) {
            w <- length(keggGet_results[[k]])
            for (j in 1:w) {
              daa_results_filtered_df[daa_results_filtered_df$feature == keggGet_results[[k]][[j]]$ENTRY, ]$pathway_name <-
                keggGet_results[[k]][[j]]$NAME[1]
              daa_results_filtered_df[daa_results_filtered_df$feature ==
                keggGet_results[[k]][[j]]$ENTRY, ]$pathway_description <-
                keggGet_results[[k]][[j]]$DESCRIPTION[1]
              daa_results_filtered_df[daa_results_filtered_df$feature == keggGet_results[[k]][[j]]$ENTRY, ]$pathway_class <-
                keggGet_results[[k]][[j]]$CLASS[1]
              daa_results_filtered_df[daa_results_filtered_df$feature == keggGet_results[[k]][[j]]$ENTRY, ]$pathway_map <-
                keggGet_results[[k]][[j]]$PATHWAY_MAP[1]
            }
          }
        }
        daa_results_filtered_annotation_df <-
          daa_results_filtered_df
        return(daa_results_filtered_annotation_df)
      }
    }
  }
