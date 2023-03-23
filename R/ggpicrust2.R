#' This function integrates pathway name/description annotations, ten of the most advanced differential abundance (DA) methods, and visualization of DA results.
#'
#' @param file A character, the file path to store picrust2 export files
#' @param metadata A tibble, consisting of sample information
#' @param group A character, name of the group
#' @param pathway A character, consisting of "EC", "KO", "MetaCyc"
#' @param daa_method A character, the chosen differential abundance analysis (DA) method
#' @param ko_to_kegg A character to control the conversion of KO abundance to KEGG abundance
#' @param p.adjust A character, the method to adjust p-values
#' @param order A character to control the order of the main plot rows
#' @param p_values_bar A character to control if the main plot has the p_values bar
#' @param x_lab A character to control the x-axis label name
#' @param select A vector consisting of pathway names to be selected
#' @param reference A character, a reference group level for several DA methods
#' @param colors A vector consisting of colors number
#'
#' @return daa.results.df, a dataframe of DA results
#' @value
#' A list of sub-lists, each containing a ggplot2 plot (`plot`) and a dataframe of differential abundance results (`results`) for a specific DA method.
#' Each plot visualizes the differential abundance results of a specific DA method, and the corresponding dataframe contains the results used to create the plot.
#' @export
#' @examples
#' \dontrun{
#' # Load necessary data: abundance data and metadata
#' abundance_file <- "path/to/your/abundance_file.tsv"
#' metadata <- read.csv("path/to/your/metadata.csv")
#'
#' # Run ggpicrust2 with desired parameters
#' results <- ggpicrust2(file = abundance_file,
#'                       metadata = metadata,
#'                       group = "your_group_column",
#'                       pathway = "KO",
#'                       daa_method = "ALDEx2")
#'
#' # Access the plot and results dataframe for the first DA method
#' example_plot <- results[[1]]$plot
#' example_results <- results[[1]]$results
#'}
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
  # 创建一个空list
  plot_result_list <- list()

  switch(ko_to_kegg,
         "TRUE" = {
           plot_result_list <- list()
           abundance <- ko2kegg_abundance(file)
           daa_results_df <- pathway_daa(
             abundance = abundance,
             metadata = metadata,
             group = group,
             daa_method = daa_method,
             select = select,
             p.adjust = p.adjust,
             reference = reference
           )
           if (sum(as.numeric(daa_results_df$p_adjust <= 0.05)) == 0){
             stop("There are no statistically significant biomarkers")
           }
           if (x_lab == "pathway_name") {
             daa_results_df  <-
               pathway_annotation(daa_results_df = daa_results_df, ko_to_kegg = TRUE)
           }
           j <- 1
           for (i in unique(daa_results_df$method)) {
             daa_sub_method_results_df <-
               daa_results_df[daa_results_df[, "method"] == i,]
             daa_sub_method_results_df_sorted <- data.frame()
             if (select == "NULL"){
               daa_sub_method_results_df_sorted <- daa_sub_method_results_df
             }else if (select == "Top 10"){
               # 对 daa_sub_method_results_df 按照 p_adjust 进行排序，自小向大
               daa_sub_method_results_df_sorted <- daa_sub_method_results_df[order(daa_sub_method_results_df$p_adjust),]
               # 保留 p_adjust 最小的十条记录
               daa_sub_method_results_df_sorted <- daa_sub_method_results_df_sorted[1:10,]
             }else if (select == "Top 20"){
               # 对 daa_sub_method_results_df 按照 p_adjust 进行排序，自小向大
               daa_sub_method_results_df_sorted <- daa_sub_method_results_df[order(daa_sub_method_results_df$p_adjust),]
               # 保留 p_adjust 最小的二十条记录
               daa_sub_method_results_df_sorted <- daa_sub_method_results_df_sorted[1:20,]
             }else if (select == "Top 30"){
               # 对 daa_sub_method_results_df 按照 p_adjust 进行排序，自小向大
               daa_sub_method_results_df_sorted <- daa_sub_method_results_df[order(daa_sub_method_results_df$p_adjust),]
               # 保留 p_adjust 最小的三十条记录
               daa_sub_method_results_df_sorted <- daa_sub_method_results_df_sorted[1:30,]
             }
             combination_bar_plot <-
               pathway_errorbar(
                 abundance = abundance,
                 daa_results_df = daa_sub_method_results_df_sorted,
                 Group = metadata[, group],
                 ko_to_kegg = ko_to_kegg,
                 order = "pathway_class",
                 colors = colors,
                 x_lab = x_lab
               )
             # 创建一个子list，包含一个combination_bar_plot和对应的daa_results_df子集
             sub_list <- list(plot = combination_bar_plot, results = daa_sub_method_results_df)
             # 将子list添加到主list中
             plot_result_list[[j]] <- sub_list
             j <- j + 1
           }
           return(plot_result_list)
         },
         "FALSE" = {
           plot_result_list <- list()
           abundance <- readr::read_delim(
             file,
             delim = "\t",
             escape_double = FALSE,
             trim_ws = TRUE
           )
           abundance <-
             tibble::column_to_rownames(abundance, var = "function")
           daa_results_df <- pathway_daa(
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
           j <- 1
           for (i in unique(daa_results_df$method)) {
             daa_sub_method_results_df <-
               daa_results_df[daa_results_df[, "method"] == i,]
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
             # 创建一个子list
             sub_list <- list(plot = combination_bar_plot, results = daa_sub_method_results_df)
             # 将子list添加到主list中
             plot_result_list[[j]] <- sub_list
             j <- j + 1
           }
           return(plot_result_list)
         })
}
