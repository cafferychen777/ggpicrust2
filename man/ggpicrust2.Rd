% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ggpicrust2.R
\name{ggpicrust2}
\alias{ggpicrust2}
\title{This function integrates pathway name/description annotations, ten of the most advanced differential abundance (DA) methods, and visualization of DA results.}
\usage{
ggpicrust2(
  file = NULL,
  data = NULL,
  metadata,
  group,
  pathway,
  daa_method = "ALDEx2",
  ko_to_kegg = FALSE,
  p.adjust = "BH",
  order = "group",
  p_values_bar = TRUE,
  x_lab = NULL,
  select = NULL,
  reference = NULL,
  colors = NULL
)
}
\arguments{
\item{file}{A character string representing the file path of the input file containing KO abundance data in picrust2 export format. The input file should have KO identifiers in the first column and sample identifiers in the first row. The remaining cells should contain the abundance values for each KO-sample pair.}

\item{data}{An optional data.frame containing KO abundance data in the same format as the input file. If provided, the function will use this data instead of reading from the file. By default, this parameter is set to NULL.}

\item{metadata}{A tibble, consisting of sample information}

\item{group}{A character, name of the group}

\item{pathway}{A character, consisting of "EC", "KO", "MetaCyc"}

\item{daa_method}{a character specifying the method for differential abundance analysis, default is "ALDEx2", choices are:
- "ALDEx2": ANOVA-Like Differential Expression tool for high throughput sequencing data
- "DESeq2": Differential expression analysis based on the negative binomial distribution using DESeq2
- "edgeR": Exact test for differences between two groups of negative-binomially distributed counts using edgeR
- "limma voom": Limma-voom framework for the analysis of RNA-seq data
- "metagenomeSeq": Fit logistic regression models to test for differential abundance between groups using metagenomeSeq
- "LinDA": Linear models for differential abundance analysis of microbiome compositional data
- "Maaslin2": Multivariate Association with Linear Models (MaAsLin2) for differential abundance analysis}

\item{ko_to_kegg}{A character to control the conversion of KO abundance to KEGG abundance}

\item{p.adjust}{a character specifying the method for p-value adjustment, default is "BH", choices are:
- "BH": Benjamini-Hochberg correction
- "holm": Holm's correction
- "bonferroni": Bonferroni correction
- "hochberg": Hochberg's correction
- "fdr": False discovery rate correction
- "none": No p-value adjustment.}

\item{order}{A character to control the order of the main plot rows}

\item{p_values_bar}{A character to control if the main plot has the p_values bar}

\item{x_lab}{A character to control the x-axis label name, you can choose from "feature","pathway_name" and "description"}

\item{select}{A vector consisting of pathway names to be selected}

\item{reference}{A character, a reference group level for several DA methods}

\item{colors}{A vector consisting of colors number}
}
\value{
daa.results.df, a dataframe of DA results
A list of sub-lists, each containing a ggplot2 plot (`plot`) and a dataframe of differential abundance results (`results`) for a specific DA method.
Each plot visualizes the differential abundance results of a specific DA method, and the corresponding dataframe contains the results used to create the plot.
}
\description{
This function integrates pathway name/description annotations, ten of the most advanced differential abundance (DA) methods, and visualization of DA results.
}
\examples{
\dontrun{
# Load necessary data: abundance data and metadata
abundance_file <- "path/to/your/abundance_file.tsv"
metadata <- read.csv("path/to/your/metadata.csv")

# Run ggpicrust2 with input file path
results_file_input <- ggpicrust2(file = abundance_file,
                                 metadata = metadata,
                                 group = "your_group_column",
                                 pathway = "KO",
                                 daa_method = "LinDA",
                                 ko_to_kegg = "TRUE",
                                 order = "pathway_class",
                                 p_values_bar = TRUE,
                                 x_lab = "pathway_name")

# Run ggpicrust2 with imported data.frame
abundance_data <- read_delim(abundance_file, delim="\t", col_names=TRUE, trim_ws=TRUE)

# Run ggpicrust2 with input data
results_data_input <- ggpicrust2(data = abundance_data,
                                 metadata = metadata,
                                 group = "your_group_column",
                                 pathway = "KO",
                                 daa_method = "LinDA",
                                 ko_to_kegg = "TRUE",
                                 order = "pathway_class",
                                 p_values_bar = TRUE,
                                 x_lab = "pathway_name")

# Access the plot and results dataframe for the first DA method
example_plot <- results_file_input[[1]]$plot
example_results <- results_file_input[[1]]$results

# Use the example data in ggpicrust2 package
data(ko_abundance)
data(metadata)
results_file_input <- ggpicrust2(data = ko_abundance,
                                 metadata = metadata,
                                 group = "Environment",
                                 pathway = "KO",
                                 daa_method = "LinDA",
                                 ko_to_kegg = TRUE,
                                 order = "pathway_class",
                                 p_values_bar = TRUE,
                                 x_lab = "pathway_name")
# Analyze the EC or MetaCyc pathway
data(metacyc_abundance)
results_file_input <- ggpicrust2(data = metacyc_abundance,
                                 metadata = metadata,
                                 group = "Environment",
                                 pathway = "MetaCyc",
                                 daa_method = "LinDA",
                                 ko_to_kegg = FALSE,
                                 order = "group",
                                 p_values_bar = TRUE,
                                 x_lab = "description")
}
}
