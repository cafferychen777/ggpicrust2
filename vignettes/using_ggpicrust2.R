## ----include = FALSE----------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----ggpicrust2(), eval = FALSE-----------------------------------------------
#  # If you want to analyze the abundance of KEGG pathways instead of KO within the pathway, please set `ko_to_kegg` to TRUE.
#  # KEGG pathways typically have more descriptive explanations.
#  
#  library(readr)
#  library(ggpicrust2)
#  library(tibble)
#  library(tidyverse)
#  library(ggprism)
#  library(patchwork)
#  
#  # Load necessary data: abundance data and metadata
#  abundance_file <- "path/to/your/abundance_file.tsv"
#  metadata <- read_delim(
#      "path/to/your/metadata.txt",
#      delim = "\t",
#      escape_double = FALSE,
#      trim_ws = TRUE
#  )
#  
#  # Run ggpicrust2 with input file path
#  results_file_input <- ggpicrust2(file = abundance_file,
#                                   metadata = metadata,
#                                   group = "your_group_column", # For example dataset, group = "Environment"
#                                   pathway = "KO",
#                                   daa_method = "LinDA",
#                                   ko_to_kegg = TRUE,
#                                   order = "pathway_class",
#                                   p_values_bar = TRUE,
#                                   x_lab = "pathway_name")
#  
#  # Run ggpicrust2 with imported data.frame
#  abundance_data <- read_delim(abundance_file, delim = "\t", col_names = TRUE, trim_ws = TRUE)
#  
#  # Run ggpicrust2 with input data
#  results_data_input <- ggpicrust2(data = abundance_data,
#                                   metadata = metadata,
#                                   group = "your_group_column", # For example dataset, group = "Environment"
#                                   pathway = "KO",
#                                   daa_method = "LinDA",
#                                   ko_to_kegg = TRUE,
#                                   order = "pathway_class",
#                                   p_values_bar = TRUE,
#                                   x_lab = "pathway_name")
#  
#  # Access the plot and results dataframe for the first DA method
#  example_plot <- results_file_input[[1]]$plot
#  example_results <- results_file_input[[1]]$results
#  
#  # Use the example data in ggpicrust2 package
#  data(ko_abundance)
#  data(metadata)
#  results_file_input <- ggpicrust2(data = ko_abundance,
#                                   metadata = metadata,
#                                   group = "Environment",
#                                   pathway = "KO",
#                                   daa_method = "LinDA",
#                                   ko_to_kegg = TRUE,
#                                   order = "pathway_class",
#                                   p_values_bar = TRUE,
#                                   x_lab = "pathway_name")
#  
#  # Analyze the EC or MetaCyc pathway
#  data(metacyc_abundance)
#  results_file_input <- ggpicrust2(data = metacyc_abundance,
#                                   metadata = metadata,
#                                   group = "Environment",
#                                   pathway = "MetaCyc",
#                                   daa_method = "LinDA",
#                                   ko_to_kegg = FALSE,
#                                   order = "group",
#                                   p_values_bar = TRUE,
#                                   x_lab = "description")
#  results_file_input[[1]]$plot
#  results_file_input[[1]]$results

## ----alternative, eval = FALSE------------------------------------------------
#  library(readr)
#  library(ggpicrust2)
#  library(tibble)
#  library(tidyverse)
#  library(ggprism)
#  library(patchwork)
#  
#  # If you want to analyze KEGG pathway abundance instead of KO within the pathway, turn ko_to_kegg to TRUE.
#  # KEGG pathways typically have more explainable descriptions.
#  
#  # Load metadata as a tibble
#  # data(metadata)
#  metadata <- read_delim("path/to/your/metadata.txt", delim = "\t", escape_double = FALSE, trim_ws = TRUE)
#  
#  # Load KEGG pathway abundance
#  # data(kegg_abundance)
#  kegg_abundance <- ko2kegg_abundance("path/to/your/pred_metagenome_unstrat.tsv")
#  
#  # Perform pathway differential abundance analysis (DAA) using ALDEx2 method
#  # Please change group to "your_group_column" if you are not using example dataset
#  daa_results_df <- pathway_daa(abundance = kegg_abundance, metadata = metadata, group = "Environment", daa_method = "ALDEx2", select = NULL, reference = NULL)
#  
#  # Filter results for ALDEx2_Welch's t test method
#  # Please check the unique(daa_results_df$method) and choose one
#  daa_sub_method_results_df <- daa_results_df[daa_results_df$method == "ALDEx2_Wilcoxon rank test", ]
#  
#  # Annotate pathway results using KO to KEGG conversion
#  daa_annotated_sub_method_results_df <- pathway_annotation(pathway = "KO", daa_results_df = daa_sub_method_results_df, ko_to_kegg = TRUE)
#  
#  # Generate pathway error bar plot
#  # Please change Group to metadata$your_group_column if you are not using example dataset
#  p <- pathway_errorbar(abundance = kegg_abundance, daa_results_df = daa_annotated_sub_method_results_df, Group = metadata$Environment, p_values_threshold = 0.05, order = "pathway_class", select = NULL, ko_to_kegg = TRUE, p_value_bar = TRUE, colors = NULL, x_lab = "pathway_name")
#  
#  # If you want to analyze EC, MetaCyc, and KO without conversions, turn ko_to_kegg to FALSE.
#  
#  # Load metadata as a tibble
#  # data(metadata)
#  metadata <- read_delim("path/to/your/metadata.txt", delim = "\t", escape_double = FALSE, trim_ws = TRUE)
#  
#  # Load KO abundance as a data.frame
#  # data(ko_abundance)
#  ko_abundance <- read.delim("path/to/your/pred_metagenome_unstrat.tsv")
#  
#  # Perform pathway DAA using ALDEx2 method
#  # Please change column_to_rownames() to the feature column if you are not using example dataset
#  # Please change group to "your_group_column" if you are not using example dataset
#  daa_results_df <- pathway_daa(abundance = ko_abundance %>% column_to_rownames("#NAME"), metadata = metadata, group = "Environment", daa_method = "ALDEx2", select = NULL, reference = NULL)
#  
#  # Filter results for ALDEx2_Kruskal-Wallace test method
#  daa_sub_method_results_df <- daa_results_df[daa_results_df$method == "ALDEx2_Wilcoxon rank test", ]
#  
#  # Annotate pathway results without KO to KEGG conversion
#  daa_annotated_sub_method_results_df <- pathway_annotation(pathway = "KO", daa_results_df = daa_sub_method_results_df, ko_to_kegg = FALSE)
#  
#  # Generate pathway error bar plot
#  # Please change column_to_rownames() to the feature column
#  # Please change Group to metadata$your_group_column if you are not using example dataset
#  p <- pathway_errorbar(abundance = ko_abundance %>% column_to_rownames("#NAME"), daa_results_df = daa_annotated_sub_method_results_df, Group = metadata$Environment, p_values_threshold = 0.05, order = "group",
#  select = daa_annotated_sub_method_results_df %>% arrange(p_adjust) %>% slice(1:20) %>% dplyr::select(feature) %>% pull(),
#  ko_to_kegg = FALSE,
#  p_value_bar = TRUE,
#  colors = NULL,
#  x_lab = "description")
#  
#  # Workflow for MetaCyc Pathway and EC
#  
#  # Load MetaCyc pathway abundance and metadata
#  data("metacyc_abundance")
#  data("metadata")
#  
#  # Perform pathway DAA using LinDA method
#  # Please change column_to_rownames() to the feature column if you are not using example dataset
#  # Please change group to "your_group_column" if you are not using example dataset
#  metacyc_daa_results_df <- pathway_daa(abundance = metacyc_abundance %>% column_to_rownames("pathway"), metadata = metadata, group = "Environment", daa_method = "LinDA")
#  
#  # Annotate MetaCyc pathway results without KO to KEGG conversion
#  metacyc_daa_annotated_results_df <- pathway_annotation(pathway = "MetaCyc", daa_results_df = daa_results_df, ko_to_kegg = FALSE)
#  
#  # Generate pathway error bar plot
#  # Please change column_to_rownames() to the feature column
#  # Please change Group to metadata$your_group_column if you are not using example dataset
#  pathway_errorbar(abundance = metacyc_abundance %>% column_to_rownames("pathway"), daa_results_df = metacyc_daa_annotated_results_df, Group = metadata$Environment, ko_to_kegg = FALSE, p_values_threshold = 0.05, order = "group", select = NULL, p_value_bar = TRUE, colors = NULL, x_lab = "description")
#  
#  # Generate pathway heatmap
#  # Please change column_to_rownames() to the feature column if you are not using example dataset
#  # Please change group to "your_group_column" if you are not using example dataset
#  feature_with_p_0.05 <- metacyc_daa_results_df %>% filter(p_adjust < 0.05)
#  pathway_heatmap(abundance = metacyc_abundance %>% filter(pathway %in% feature_with_p_0.05$feature) %>% column_to_rownames("pathway"), metadata = metadata, group = "Environment")
#  
#  # Generate pathway PCA plot
#  # Please change column_to_rownames() to the feature column if you are not using example dataset
#  # Please change group to "your_group_column" if you are not using example dataset
#  pathway_pca(abundance = metacyc_abundance %>% column_to_rownames("pathway"), metadata = metadata, group = "Environment")
#  
#  # Run pathway DAA for multiple methods
#  # Please change column_to_rownames() to the feature column if you are not using example dataset
#  # Please change group to "your_group_column" if you are not using example dataset
#  methods <- c("ALDEx2", "DESeq2", "edgeR")
#  daa_results_list <- lapply(methods, function(method) {
#    pathway_daa(abundance = metacyc_abundance %>% column_to_rownames("pathway"), metadata = metadata, group = "Environment", daa_method = method)
#  })
#  
#  # Compare results across different methods
#  comparison_results <- compare_daa_results(daa_results_list = daa_results_list, method_names = c("ALDEx2_Welch's t test", "ALDEx2_Wilcoxon rank test", "DESeq2", "edgeR"))
#  

## ----ko2kegg_abundance sample,echo = TRUE,eval=FALSE--------------------------
#  # Sample usage of the ko2kegg_abundance function
#  
#  library(ggpicrust2)
#  
#  # Assume that the KO abundance table is stored in a file named "ko_abundance.tsv"
#  ko_abundance_file <- "ko_abundance.tsv"
#  
#  # Convert KO abundance to KEGG pathway abundance
#  kegg_abundance <- ko2kegg_abundance(file = ko_abundance_file)
#  
#  # Alternatively, if the KO abundance data is already loaded as a data frame named "ko_abundance"
#  data("ko_abundance")
#  kegg_abundance <- ko2kegg_abundance(data = ko_abundance)
#  
#  # The resulting kegg_abundance data frame can now be used for further analysis and visualization.
#  

## ----pathway_daa sample,echo = TRUE,eval=FALSE--------------------------------
#  # The abundance table is recommended to be a data.frame rather than a tibble.
#  # The abundance table should have feature names or pathway names as row names, and sample names as column names.
#  # You can use the output of ko2kegg_abundance
#  ko_abundance_file <- "path/to/your/pred_metagenome_unstrat.tsv"
#  kegg_abundance <- ko2kegg_abundance(ko_abundance_file) # Or use data(kegg_abundance)
#  
#  metadata <- read_delim("path/to/your/metadata.txt", delim = "\t", escape_double = FALSE, trim_ws = TRUE)
#  
#  # The default DAA method is "ALDEx2"
#  # Please change group to "your_group_column" if you are not using example dataset
#  daa_results_df <- pathway_daa(abundance = kegg_abundance, metadata = metadata, group = "Environment", daa_method = "linDA", select = NULL, p.adjust = "BH", reference = NULL)
#  
#  # If you have more than 3 group levels and want to use the LinDA, limma voom, or Maaslin2 methods, you should provide a reference.
#  metadata <- read_delim("path/to/your/metadata.txt", delim = "\t", escape_double = FALSE, trim_ws = TRUE)
#  
#  # Please change group to "your_group_column" if you are not using example dataset
#  daa_results_df <- pathway_daa(abundance = kegg_abundance, metadata = metadata, group = "Group", daa_method = "LinDA", select = NULL, p.adjust = "BH", reference = "Harvard BRI")
#  
#  # Other example
#  data("metacyc_abundance")
#  data("metadata")
#  metacyc_daa_results_df <- pathway_daa(abundance = metacyc_abundance %>% column_to_rownames("pathway"), metadata = metadata, group = "Environment", daa_method = "LinDA", select = NULL, p.adjust = "BH", reference = NULL)

## ----compare_daa_results sample,echo = TRUE,eval=FALSE------------------------
#  library(ggpicrust2)
#  library(tidyverse)
#  data("metacyc_abundance")
#  data("metadata")
#  
#  # Run pathway_daa function for multiple methods
#  # Please change column_to_rownames() to the feature column if you are not using example dataset
#  # Please change group to "your_group_column" if you are not using example dataset
#  methods <- c("ALDEx2", "DESeq2", "edgeR")
#  daa_results_list <- lapply(methods, function(method) {
#    pathway_daa(abundance = metacyc_abundance %>% column_to_rownames("pathway"), metadata = metadata, group = "Environment", daa_method = method)
#  })
#  
#  method_names <- c("ALDEx2_Welch's t test","ALDEx2_Wilcoxon rank test","DESeq2", "edgeR")
#  # Compare results across different methods
#  comparison_results <- compare_daa_results(daa_results_list = daa_results_list, method_names = method_names)

## ----pathway_annotation sample,echo = TRUE,eval=FALSE-------------------------
#  
#  # Make sure to check if the features in `daa_results_df` correspond to the selected pathway
#  
#  # Annotate KEGG Pathway
#  data("kegg_abundance")
#  data("metadata")
#  # Please change group to "your_group_column" if you are not using example dataset
#  daa_results_df <- pathway_daa(abundance = kegg_abundance, metadata = metadata, group = "Environment", daa_method = "LinDA")
#  daa_annotated_results_df <- pathway_annotation(pathway = "KO", daa_results_df = daa_results_df, ko_to_kegg = TRUE)
#  
#  # Annotate KO
#  data("ko_abundance")
#  data("metadata")
#  # Please change column_to_rownames() to the feature column if you are not using example dataset
#  # Please change group to "your_group_column" if you are not using example dataset
#  daa_results_df <- pathway_daa(abundance = ko_abundance %>% column_to_rownames("#NAME"), metadata = metadata, group = "Environment", daa_method = "LinDA")
#  daa_annotated_results_df <- pathway_annotation(pathway = "KO", daa_results_df = daa_results_df, ko_to_kegg = FALSE)
#  
#  # Annotate KEGG
#  # daa_annotated_results_df <- pathway_annotation(pathway = "EC", daa_results_df = daa_results_df, ko_to_kegg = FALSE)
#  
#  # Annotate MetaCyc Pathway
#  data("metacyc_abundance")
#  data("metadata")
#  # Please change column_to_rownames() to the feature column if you are not using example dataset
#  # Please change group to "your_group_column" if you are not using example dataset
#  metacyc_daa_results_df <- pathway_daa(abundance = metacyc_abundance %>% column_to_rownames("pathway"), metadata = metadata, group = "Environment", daa_method = "LinDA")
#  metacyc_daa_annotated_results_df <- pathway_annotation(pathway = "MetaCyc", daa_results_df = metacyc_daa_results_df, ko_to_kegg = FALSE)

## ----pathway_errorbar sample,echo = TRUE,eval=FALSE---------------------------
#  data("ko_abundance")
#  data("metadata")
#  kegg_abundance <- ko2kegg_abundance(data = ko_abundance) # Or use data(kegg_abundance)
#  # Please change group to "your_group_column" if you are not using example dataset
#  daa_results_df <- pathway_daa(kegg_abundance, metadata = metadata, group = "Environment", daa_method = "LinDA")
#  daa_annotated_results_df <- pathway_annotation(pathway = "KO", daa_results_df = daa_results_df, ko_to_kegg = TRUE)
#  # Please change Group to metadata$your_group_column if you are not using example dataset
#  p <- pathway_errorbar(abundance = kegg_abundance,
#             daa_results_df = daa_annotated_results_df,
#             Group = metadata$Environment,
#             ko_to_kegg = TRUE,
#             p_values_threshold = 0.05,
#             order = "pathway_class",
#             select = NULL,
#             p_value_bar = TRUE,
#             colors = NULL,
#             x_lab = "pathway_name")
#  
#  # If you want to analysis the EC. MetaCyc. KO without conversions.
#  data("metacyc_abundance")
#  data("metadata")
#  metacyc_daa_results_df <- pathway_daa(abundance = metacyc_abundance %>% column_to_rownames("pathway"), metadata = metadata, group = "Environment", daa_method = "LinDA")
#  metacyc_daa_annotated_results_df <- pathway_annotation(pathway = "MetaCyc", daa_results_df = metacyc_daa_results_df, ko_to_kegg = FALSE)
#  p <- pathway_errorbar(abundance = metacyc_abundance %>% column_to_rownames("pathway"),
#             daa_results_df = metacyc_daa_annotated_results_df,
#             Group = metadata$Environment,
#             ko_to_kegg = FALSE,
#             p_values_threshold = 0.05,
#             order = "group",
#             select = NULL,
#             p_value_bar = TRUE,
#             colors = NULL,
#             x_lab = "description")

## ----echo = TRUE,eval=FALSE---------------------------------------------------
#  # Create example functional pathway abundance data
#  abundance_example <- matrix(rnorm(30), nrow = 3, ncol = 10)
#  colnames(abundance_example) <- paste0("Sample", 1:10)
#  rownames(abundance_example) <- c("PathwayA", "PathwayB", "PathwayC")
#  
#  # Create example metadata
#  # Please change your sample id's column name to sample_name
#  metadata_example <- data.frame(sample_name = colnames(abundance_example),
#                                 group = factor(rep(c("Control", "Treatment"), each = 5)))
#  
#  # Create a heatmap
#  pathway_heatmap(abundance_example, metadata_example, "group")

## ----echo = TRUE,eval=FALSE---------------------------------------------------
#  # Load the data
#  data("metacyc_abundance")
#  
#  # Load the metadata
#  data("metadata")
#  
#  # Perform differential abundance analysis
#  metacyc_daa_results_df <- pathway_daa(
#    abundance = metacyc_abundance %>% column_to_rownames("pathway"),
#    metadata = metadata,
#    group = "Environment",
#    daa_method = "LinDA"
#  )
#  
#  # Annotate the results
#  annotated_metacyc_daa_results_df <- pathway_annotation(
#    pathway = "MetaCyc",
#    daa_results_df = metacyc_daa_results_df,
#    ko_to_kegg = FALSE
#  )
#  
#  # Filter features with p < 0.05
#  feature_with_p_0.05 <- metacyc_daa_results_df %>%
#    filter(p_adjust < 0.05)
#  
#  # Create the heatmap
#  pathway_heatmap(
#    abundance = metacyc_abundance %>%
#      right_join(
#        annotated_metacyc_daa_results_df %>% select(all_of(c("feature","description"))),
#        by = c("pathway" = "feature")
#      ) %>%
#      filter(pathway %in% feature_with_p_0.05$feature) %>%
#      select(-"pathway") %>%
#      column_to_rownames("description"),
#    metadata = metadata,
#    group = "Environment"
#  )

## ----echo = TRUE,eval=FALSE---------------------------------------------------
#  # Create example functional pathway abundance data
#  abundance_example <- matrix(rnorm(30), nrow = 3, ncol = 10)
#  colnames(abundance_example) <- paste0("Sample", 1:10)
#  rownames(kegg_abundance_example) <- c("PathwayA", "PathwayB", "PathwayC")
#  
#  # Create example metadata
#  metadata_example <- data.frame(sample_name = colnames(kegg_abundance_example),
#                                  group = factor(rep(c("Control", "Treatment"), each = 5)))
#  # Perform PCA and create visualizations
#  pathway_pca(abundance = abundance_example, metadata = metadata_example, "group")

## ----echo = TRUE,eval=FALSE---------------------------------------------------
#  # Create example functional pathway abundance data
#  data("metacyc_abundance")
#  data("metadata")
#  
#  pathway_pca(abundance = metacyc_abundance %>% column_to_rownames("pathway"), metadata = metadata, group = "Environment")

## ----compare_metagenome_results sample,echo = TRUE,eval=FALSE-----------------
#  library(ComplexHeatmap)
#  set.seed(123)
#  # First metagenome
#  metagenome1 <- abs(matrix(rnorm(1000), nrow = 100, ncol = 10))
#  rownames(metagenome1) <- paste0("KO", 1:100)
#  colnames(metagenome1) <- paste0("sample", 1:10)
#  # Second metagenome
#  metagenome2 <- abs(matrix(rnorm(1000), nrow = 100, ncol = 10))
#  rownames(metagenome2) <- paste0("KO", 1:100)
#  colnames(metagenome2) <- paste0("sample", 1:10)
#  # Put the metagenomes into a list
#  metagenomes <- list(metagenome1, metagenome2)
#  # Define names
#  names <- c("metagenome1", "metagenome2")
#  # Call the function
#  results <- compare_metagenome_results(metagenomes, names)
#  # Print the correlation matrix
#  print(results$correlation$cor_matrix)
#  # Print the p-value matrix
#  print(results$correlation$p_matrix)

