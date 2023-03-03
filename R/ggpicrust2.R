# This function integrates pathway name/description annotations, ten of the most advanced differential abundance (DA) methods, and visualization of DA results.
#
# @param file A character, the file path to store picrust2 export files
# @param metadata A tibble, consisting of sample information
# @param group A character, name of the group
# @param pathway A character, consisting of "EC", "KO", "MetaCyc"
# @param daa_method A character, the chosen differential abundance analysis (DA) method
# @param ko_to_kegg A character to control the conversion of KO abundance to KEGG abundance
# @param p.adjust A character, the method to adjust p-values
# @param order A character to control the order of the main plot rows
# @param p_values_bar A character to control if the main plot has the p_values bar
# @param x_lab A character to control the x-axis label name
# @param select A vector consisting of pathway names to be selected
# @param reference A character, a reference group level for several DA methods
# @param colors A vector consisting of colors number
#
# @return daa.results.df, a dataframe of DA results
# @export
#
# @examples
ggpicrust2 <- function(file,
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
                       colors = NULL) {
  # switch function to choose between converting KO abundance to KEGG abundance or not
  switch(ko_to_kegg,
         "TRUE" = {
           # if the ko_to_kegg argument is set to TRUE, call the ko2kegg_abundance function to convert abundance
           if (ko_to_kegg == TRUE) {
             abundance <- ko2kegg_abundance(file)
           }
           # call pathway_daa function to perform DA analysis
           daa_results_df <- pathway_daa(abundance = abundance,
                                         metadata = metadata,
                                         group = group,
                                         daa_method = daa_method,
                                         select = select,
                                         p.adjust = p.adjust,
                                         reference = reference)
           # if the x_lab argument is set to "pathway_name", call the pathway_annotation function to add pathway name/description annotations
           if (x_lab == "pathway_name") {
             daa_results_df  <- pathway_annotation(daa_results_df = daa_results_df, ko_to_kegg = TRUE)
           }
           # loop over unique values of the "method" column in the daa_results_df dataframe
           j <- 1
           for (i in unique(daa_results_df$method)) {
             # select rows from the daa_results_df dataframe that have the current method value
             daa_sub_method_results_df <- daa_results_df[daa_results_df[, "method"] == i, ]
                     combination_bar_plot <-
                 pathway_errorbar(
                   abundance, # Input abundance data
                   daa_sub_method_results_df, # Results data frame from differential abundance analysis
                   metadata[, group], # Metadata with group information
                   ko_to_kegg = ko_to_kegg, # Option to convert KO numbers to KEGG pathway names
                   order = "pathway_class", # Order to arrange the pathways in the plot
                   colors = colors, # Colors to use in the plot
                   x_lab = x_lab # Label for the x-axis
                 )
               print(combination_bar_plot) # Print the combination bar plot
               message(paste0("No.", j, " plot is method ", i)) # Message indicating the plot and the method used
             }
             return(daa_results_df) # Return the results data frame
           },
           "FALSE" = {
             abundance <- read_delim(
                 file, # Input file
                 delim = "\t", # Delimiter in the file
                 escape_double = FALSE, # Option to escape double quotes
                 trim_ws = TRUE # Option to trim white spaces
               )
             abundance <- column_to_rownames(abundance, var = "function") # Convert a column to row names
             daa_results_df <- pathway_daa(
                 abundance = abundance, # Input abundance data
                 metadata = metadata, # Metadata
                 group = group, # Group information
                 daa_method = daa_method, # Differential abundance analysis method
                 select = select, # Option to select pathways with significant differences
                 p.adjust = p.adjust, # Method to adjust p-values
                 reference = reference # Option to use a reference group
               )
             daa_results_df <-
               pathway_annotation(
                 pathway = pathway, # Option to use custom pathway names
                 ko_to_kegg = FALSE, # Option to convert KO numbers to KEGG pathway names
                 daa_results_df = daa_results_df # Results data frame from differential abundance analysis
               )
             for (i in unique(daa_results_df$method)) { # Loop through unique methods
               daa_sub_method_results_df <-
                 daa_results_df[daa_results_df[, "method"] == i, ] # Subset results for each method
               combination_bar_plot <-
                 pathway_errorbar(
                   abundance, # Input abundance data
                   daa_sub_method_results_df, # Results data frame from differential abundance analysis
                   metadata[, group], # Metadata with group information
                   ko_to_kegg = ko_to_kegg, # Option to convert KO numbers to KEGG pathway names
                   order = order, # Order to arrange the pathways in the plot
                   colors = colors, # Colors to use in the plot
                   x_lab = x_lab # Label for the x-axis
                 )
               print(combination_bar_plot) # Print the combination bar plot
               message(paste0("No.", j, " plot is method ", i)) # Message indicating the plot and the method used
             }
             return(daa_results_df) # Return the results data frame
           })
                       }