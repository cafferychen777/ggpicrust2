#' function which integrates pathway name/description annotations, ten of the most advanced differential abundance (DA) methods, and visualization of DA results.
#'
#' @param file A character, address to store picrust2 export files
#' @param metadata A tibble, consisting of samples information
#' @param group A character, group name
#' @param pathway A character, consisting of "EC", "KO", "MetaCyc"
#' @param daa_method A character, choosing the da method
#' @param ko_to_kegg A charachter to control if converting ko abundance to kegg abundance
#' @param p.adjust A character, the method of adjust p
#' @param order A character to control the main plot rows' order
#' @param p_values_bar A character to control if the main plot has the p_values_bar
#' @param x_lab A character to control x_lab name
#' @param select A vector consisting of pathway names
#' @param reference A character, several of da methods need a reference group level
#' @param colors A vector consisting of colors number
#'
#' @return daa.results.df
#' @export
#'
#' @examples
ggpicrust2 <-
  function(file,
           metadata,
           group,
           pathway,
           daa_method = "ALDEx2",
           ko_to_kegg = FALSE,
           p.adjust = "BH",
           order = "group",
           p_values_bar = TRUE,
           x_lab = "pathway_name",
           select = NULL,
           reference = NULL,
           colors = NULL)
  {
    switch(ko_to_kegg,
           "TRUE" = {
             if (ko_to_kegg == TRUE) {
               abundance <- ko2kegg_abundance(file)
             }
             daa_results_df <-
               pathway_daa(
                 abundance = abundance,
                 metadata = metadata,
                 group = group,
                 daa_method = daa_method,
                 select = select,
                 p.adjust = p.adjust,
                 reference = reference
               )
             if (x_lab == "pathway_name") {
               daa_results_df  <-
                 pathway_annotation(daa_results_df = daa_results_df,ko_to_kegg = TRUE)
             }
             j <- 1
             for (i in unique(daa_results_df$method)) {
               daa_sub_method_results_df <-
                 daa_results_df[daa_results_df[, "method"] == i, ]
               combination_bar_plot <-
                 pathway_errorbar(
                   abundance,
                   daa_sub_method_results_df,
                   metadata[, group],
                   ko_to_kegg = ko_to_kegg,
                   order = "pathway_class",
                   colors = colors,
                   x_lab = x_lab
                 )
               print(combination_bar_plot)
               message(paste0("No.", j, " plot is method ", i))
             }
             return(daa_results_df)
           },
           "FALSE" = {
             abundance <-
               read_delim(
                 file,
                 delim = "\t",
                 escape_double = FALSE,
                 trim_ws = TRUE
               )
             abundance <- column_to_rownames(abundance, var = "function")
             daa_results_df <-
               pathway_daa(
                 abundance = abundance,
                 metadata = metadata,
                 group = group,
                 daa_method = daa_method,
                 select = select,
                 p.adjust = p.adjust,
                 reference = reference
               )
             daa_results_df <-
               pathway_annotation(
                 pathway = pathway,
                 ko_to_kegg = FALSE,
                 daa_results_df = daa_results_df
               )
             for (i in unique(daa_results_df$method)) {
               daa_sub_method_results_df <-
                 daa_results_df[daa_results_df[, "method"] == i, ]
               combination_bar_plot <-
                 pathway_errorbar(
                   abundance,
                   daa_sub_method_results_df,
                   metadata[, group],
                   ko_to_kegg = ko_to_kegg,
                   order = order,
                   colors = colors,
                   x_lab = x_lab
                 )
               print(combination_bar_plot)
               message(paste0("No.", j, " plot is method ", i))
             }
             return(daa_results_df)
           })
  }
