#' Pathway information annotation of "EC", "KO", "MetaCyc" pathway
#'
#' This function has two primary use cases:
#' 1. Annotating pathway information using the output file from PICRUSt2.
#' 2. Annotating pathway information from the output of `pathway_daa` function, and converting KO abundance to KEGG pathway abundance when `ko_to_kegg` is set to TRUE.
#'
#' @param file A character, address to store PICRUSt2 export files. Provide this parameter when using the function for the first use case.
#' @param pathway A character, consisting of "EC", "KO", "MetaCyc"
#' @param daa_results_df A data frame, output of pathway_daa. Provide this parameter when using the function for the second use case.
#' @param ko_to_kegg A logical, decide if convert KO abundance to KEGG pathway abundance. Default is FALSE. Set to TRUE when using the function for the second use case.
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
#' # Use case 1: Annotating pathway information using the output file from PICRUSt2
#' result1 <- pathway_annotation(file = "path/to/picrust2/export/file.txt",
#'                               pathway = "KO",
#'                               daa_results_df = NULL,
#'                               ko_to_kegg = FALSE)
#'
#' # Use case 2: Annotating pathway information from the output of pathway_daa function
#' # and converting KO abundance to KEGG pathway abundance
#' # This use case will be demonstrated using both a hypothetical example, and a real dataset.
#'
#' ## Hypothetical example
#' hypothetical_daa_results_df <- data.frame() # Replace this with your actual data frame
#' result2 <- pathway_annotation(file = NULL,
#'                               pathway = "KO",
#'                               daa_results_df = hypothetical_daa_results_df,
#'                               ko_to_kegg = TRUE)
#'
#' ## Real dataset example
#' # Load the real dataset
#' data(daa_results_df)
#' result3 <- pathway_annotation(file = NULL,
#'                               pathway = "KO",
#'                               daa_results_df = daa_results_df,
#'                               ko_to_kegg = TRUE)
#' }
pathway_annotation <-
  function(file = NULL,
           pathway = NULL,
           daa_results_df = NULL,
           ko_to_kegg = FALSE) {

    message("Starting pathway annotation...")

    if (is.null(file) && is.null(daa_results_df)) {
      stop("Please input the picrust2 output or results of pathway_daa daa_results_df")
    }
    if (!is.null(file)) {
      message("Reading the input file...")
      file_format <- substr(file, nchar(file) - 3, nchar(file))
      switch(file_format,
             ".txt" = {
               message("Loading .txt file...")
               abundance <- readr::read_delim(
                 file,
                 delim = "\t",
                 escape_double = FALSE,
                 trim_ws = TRUE
               )
               message(".txt file successfully loaded.")
             },
             ".tsv" = {
               message("Loading .tsv file...")
               abundance <- readr::read_delim(
                 file,
                 delim = "\t",
                 escape_double = FALSE,
                 trim_ws = TRUE
               )
               message(".tsv file successfully loaded.")
             },
             ".csv" = {
               message("Loading .csv file...")
               abundance <- readr::read_delim(
                 file,
                 delim = "\t",
                 escape_double = FALSE,
                 trim_ws = TRUE
               )
               message(".csv file successfully loaded.")
             },
             stop(
               "Invalid file format. Please input file in .tsv, .txt or .csv format. The best input file format is the output file from PICRUSt2, specifically 'pred_metagenome_unstrat.tsv'."
             )
      )
      abundance <-
        abundance %>% tibble::add_column(
          description = rep(NA, length = nrow(abundance)),
          .after = 1
        )
      switch(pathway,
             "KO" = {
               message("Loading KO reference data...")
               load(system.file("extdata", "KO_reference.RData", package = "ggpicrust2"))
               message("Annotating abundance data with KO reference...")
               for (i in seq_len(nrow(abundance))) {
                 abundance[i, 2] <- KO_reference[KO_reference[, 1] %in% abundance[i, 1], 5][1]
               }
               message("Abundance data annotation with KO reference completed.")
             },
             "EC" = {
               message("Loading EC reference data...")
               load(system.file("extdata", "EC_reference.RData", package = "ggpicrust2"))
               message("Annotating abundance data with EC reference...")
               for (i in seq_len(nrow(abundance))) {
                 abundance[i, 2] <- EC_reference[EC_reference[, 1] %in% abundance[i, 1], 2]
               }
               message("Abundance data annotation with EC reference completed.")
               message("Note: EC description may appear to be duplicated due to shared EC numbers across different reactions.")
             },
             "MetaCyc" = {
               message("Loading MetaCyc reference data...")
               load(system.file("extdata", "MetaCyc_reference.RData", package = "ggpicrust2"))
               message("Annotating abundance data with MetaCyc reference...")
               for (i in seq_len(nrow(abundance))) {
                 abundance[i, 2] <- MetaCyc_reference[MetaCyc_reference[, 1] %in% abundance[i, 1], 2]
               }
               message("Abundance data annotation with MetaCyc reference completed.")
             },
             stop("Invalid pathway option. Please provide one of the following options: 'KO', 'EC', 'MetaCyc'.")
      )
      return(abundance)
    }
    if (!is.null(daa_results_df)) {
      message("DAA results data frame is not null. Proceeding...")
      if (ko_to_kegg == FALSE) {
        message("KO to KEGG is set to FALSE. Proceeding with standard workflow...")
        daa_results_df$description <- NA
        switch(pathway,
               "KO" = {
                 message("Loading KO reference data...")
                 load(system.file("extdata", "KO_reference.RData", package = "ggpicrust2"))
                 for (i in seq_len(nrow(daa_results_df))) {
                   daa_results_df[i, ]$description <-
                     KO_reference[KO_reference[, 1] %in% daa_results_df[i, ]$feature, 5][1]
                 }
               },
               "EC" = {
                 message("Loading EC reference data...")
                 load(system.file("extdata", "EC_reference.RData", package = "ggpicrust2"))
                 for (i in seq_len(nrow(daa_results_df))) {
                   daa_results_df[i, ]$description <-
                     EC_reference[EC_reference[, 1] %in% daa_results_df[i, ]$feature, 2]
                 }
                 message("EC description may appear to be duplicated")
               },
               "MetaCyc" = {
                 message("Loading MetaCyc reference data...")
                 load(system.file("extdata", "MetaCyc_reference.RData", package = "ggpicrust2"))
                 for (i in seq_len(nrow(daa_results_df))) {
                   daa_results_df[i, ]$description <-
                     MetaCyc_reference[MetaCyc_reference[, 1] %in% daa_results_df[i, ]$feature, 2]
                 }
               },
               stop("Only provide 'KO', 'EC' and 'MetaCyc' pathway")
        )
        message("Returning DAA results data frame...")
        return(daa_results_df)
      } else {
        message("KO to KEGG is set to TRUE. Proceeding with KEGG pathway annotations...")
        daa_results_filtered_df <- daa_results_df[daa_results_df$p_adjust < 0.05, ]
        if (nrow(daa_results_filtered_df) == 0) {
          stop(
            "No statistically significant biomarkers found. 'Statistically significant biomarkers' refer to those biomarkers that demonstrate a significant difference in expression between different groups, as determined by a statistical test (p_adjust < 0.05 in this case).\n",
            "You might consider re-evaluating your experiment design or trying alternative statistical analysis methods. Consult with a biostatistician or a data scientist if you are unsure about the next steps."
          )
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
            "Too many statistically significant pathways. The database cannot handle the query all at once. Please break it down into smaller queries."
          )
        }
        if (nrow(daa_results_filtered_df) <= 10) {
          message("Processing pathways individually...")
          for (i in seq_len(nrow(daa_results_filtered_df))) {
            message("Beginning annotation for pathway ", i, " of ", nrow(daa_results_filtered_df), "...")
            a <- 0
            start_time <- Sys.time() # start timer
            repeat {
              tryCatch(
                {
                  keggGet_results[[i]] <- KEGGREST::keggGet(daa_results_filtered_df$feature[i])
                  a <- 1
                },
                error = function(e) {
                  message("An error occurred. Retrying...")
                }
              )
              if (a == 1) {
                break
              }
            }
            end_time <- Sys.time() # end timer
            time_taken <- end_time - start_time # calculate time taken
            message("Annotated pathway ", i, " of ", nrow(daa_results_filtered_df), ". Time taken: ", round(time_taken, 2), " seconds.")
            daa_results_filtered_df[i, ]$pathway_name <- keggGet_results[[i]][[1]]$NAME
            daa_results_filtered_df[i, ]$pathway_description <- keggGet_results[[i]][[1]]$DESCRIPTION[1]
            daa_results_filtered_df[i, ]$pathway_class <- keggGet_results[[i]][[1]]$CLASS
            daa_results_filtered_df[i, ]$pathway_map <- keggGet_results[[i]][[1]]$PATHWAY_MAP
          }
        }
        if (nrow(daa_results_filtered_df) > 10 && nrow(daa_results_filtered_df) < 99) {
          message("Processing pathways in chunks...")
          start_time <- Sys.time() # start timer
          n <- length(c(seq(10, nrow(daa_results_filtered_df), 10), nrow(daa_results_filtered_df)))
          j <- 1
          seq <- c(seq(10, nrow(daa_results_filtered_df), 10), nrow(daa_results_filtered_df))
          for (i in seq) {
            if (i %% 10 == 0) {
              keggGet_results[[j]] <- KEGGREST::keggGet(daa_results_filtered_df$feature[seq(i - 9, i, 1)])
            } else {
              keggGet_results[[j]] <- KEGGREST::keggGet(daa_results_filtered_df$feature[seq(nrow(daa_results_filtered_df) %/% 10 * 10 + 1, i, 1)])
            }
            j <- j + 1
          }
          end_time <- Sys.time() # end timer
          time_taken <- end_time - start_time # calculate time taken
          message("Finished processing chunks. Time taken: ", round(time_taken, 2), " seconds.")

          message("Finalizing pathway annotations...")
          start_time <- Sys.time() # start timer
          for (k in 1:n) {
            w <- length(keggGet_results[[k]])
            for (j in 1:w) {
              daa_results_filtered_df[daa_results_filtered_df$feature == keggGet_results[[k]][[j]]$ENTRY, ]$pathway_name <- keggGet_results[[k]][[j]]$NAME[1]
              daa_results_filtered_df[daa_results_filtered_df$feature == keggGet_results[[k]][[j]]$ENTRY, ]$pathway_description <- keggGet_results[[k]][[j]]$DESCRIPTION[1]
              daa_results_filtered_df[daa_results_filtered_df$feature == keggGet_results[[k]][[j]]$ENTRY, ]$pathway_class <- keggGet_results[[k]][[j]]$CLASS[1]
              daa_results_filtered_df[daa_results_filtered_df$feature == keggGet_results[[k]][[j]]$ENTRY, ]$pathway_map <- keggGet_results[[k]][[j]]$PATHWAY_MAP[1]
            }
          }
          end_time <- Sys.time() # end timer
          time_taken <- end_time - start_time # calculate time taken
          message("Finished finalizing pathway annotations. Time taken: ", round(time_taken, 2), " seconds.")
        }
        daa_results_filtered_annotation_df <-
          daa_results_filtered_df
        message("Returning DAA results filtered annotation data frame...")
        return(daa_results_filtered_annotation_df)
      }
    }
  }
