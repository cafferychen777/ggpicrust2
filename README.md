
<!-- README.md is generated from README.Rmd. Please edit that file -->

# ggpicrust2 Documentation

*ggpicrust2* is a comprehensive package that integrates pathway
name/description annotations, ten of the most advanced differential
abundance (DA) methods, and visualization of DA results. It offers a
comprehensive solution for analyzing and interpreting the results of
PICRUSt2 functional prediction in a seamless and intuitive way. Whether
you are a researcher, data scientist, or bioinformatician, *ggpicrust2*
can help you better understand the underlying biological processes and
mechanisms at play in your PICRUSt2 output data. So if you are
interested in exploring the output data of PICRUSt2, *ggpicrust2* is the
tool you need.

[![CRAN
version](https://www.r-pkg.org/badges/version/ggpicrust2)](https://CRAN.R-project.org/package=ggpicrust2)
[![Downloads](https://cranlogs.r-pkg.org/badges/grand-total/ggpicrust2)](https://CRAN.R-project.org/package=ggpicrust2)
[![License:
MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/license/mit/)

## News

Update alert! 🚨 ggpicrust2 Ver1.6.3 is now available. To get the latest
version, simply use

``` r
devtools::install_github('cafferychen777/ggpicrust2')
```

Don’t miss out on the major bug fixes and improvements. 😀

## Table of Contents

- [Citation](#citation)
- [Installation](#installation)
- [Stay Updated](#stay-updated)
- [Workflow](#workflow)
- [Output](#output)
- [Function Details](#function-details)
  - [ko2kegg_abundance()](#ko2kegg_abundance)
  - [pathway_daa()](#pathway_daa)
  - [compare_daa_results()](#compare_daa_results)
  - [pathway_annotation()](#pathway_annotation)
  - [pathway_errorbar()](#pathway_errorbar)
  - [pathway_heatmap()](#pathway_heatmap)
  - [pathway_pca()](#pathway_pca)
  - [compare_metagenome_results()](#compare_metagenome_results)
- [FAQ](#faq)
- [Author’s Other Projects](#authors-other-projects)

## Citation

If you use *ggpicrust2* in your research, please cite the following
paper:

Chen Yang, Aaron Burberry, Xuan Cao, Jiahao Mai, Liangliang Zhang.
(2023). [ggpicrust2: an R package for PICRUSt2 predicted functional
profile analysis and visualization.](https://arxiv.org/abs/2303.10388)
arXiv preprint arXiv:2303.10388.

BibTeX entry: @misc{yang2023ggpicrust2, title={ggpicrust2: an R package
for PICRUSt2 predicted functional profile analysis and visualization},
author={Chen Yang and Aaron Burberry and Jiahao Mai and Xuan Cao and
Liangliang Zhang}, year={2023}, eprint={2303.10388},
archivePrefix={arXiv}, primaryClass={stat.AP} }

ResearchGate preprint link: [Click
here](https://www.researchgate.net/publication/369380349_ggpicrust2_an_R_package_for_PICRUSt2_predicted_functional_profile_analysis_and_visualization)

## Installation

You can install the stable version of *ggpicrust2* from CRAN with:

``` r
install.packages("ggpicrust2")
```

To install the latest development version of *ggpicrust2* from GitHub,
you can use:

``` r
# Install the devtools package if not already installed
# install.packages("devtools")

# Install ggpicrust2 from GitHub
devtools::install_github("cafferychen777/ggpicrust2")
```

## Dependent CRAN Packages

| Package        | Description                                                                                                                       |
|----------------|-----------------------------------------------------------------------------------------------------------------------------------|
| ALDEx2         | Differential abundance analysis of taxonomic and functional features                                                              |
| aplot          | Create interactive plots                                                                                                          |
| dplyr          | A fast consistent tool for working with data frame like objects both in memory and out of memory both in memory and out of memory |
| edgeR          | Empirical Analysis of Digital Gene Expression Data in R                                                                           |
| ggplot2        | An implementation of the Grammar of Graphics in R                                                                                 |
| grid           | A rewrite of the graphics layout capabilities of R                                                                                |
| limma          | Linear Models for Microarray and RNA-Seq Data                                                                                     |
| Maaslin2       | Tools for microbiome analysis                                                                                                     |
| metagenomeSeq  | Statistical analysis for sparse high-throughput sequencing                                                                        |
| MicrobiomeStat | Statistical analysis of microbiome data                                                                                           |
| readr          | Read rectangular data (csv tsv fwf) into R                                                                                        |
| stats          | The R Stats Package                                                                                                               |
| tibble         | Simple Data Frames                                                                                                                |
| tidyr          | Easily tidy data with spread() and gather() functions                                                                             |
| ggprism        | Interactive 3D plots with ‘prism’ graphics                                                                                        |
| cowplot        | Streamlined Plot Theme and Plot Annotations for ‘ggplot2’                                                                         |
| KEGGREST       | R Interface to KEGG REST API                                                                                                      |
| ggforce        | Easily add secondary axes, zooms, and image overlays to ‘ggplot2’                                                                 |
| ggplotify      | Convert complex plots into ‘grob’ or ‘ggplot’ objects                                                                             |
| magrittr       | A Forward-Pipe Operator for R                                                                                                     |
| utils          | The R Utils Package                                                                                                               |

## Dependent Bioconductor Packages

| Package              | Description                                                           |
|----------------------|-----------------------------------------------------------------------|
| phyloseq             | Handling and analysis of high-throughput microbiome census data       |
| SummarizedExperiment | SummarizedExperiment container for storing data and metadata together |
| Biobase              | Base functions for Bioconductor                                       |
| devtools             | Tools to make developing R packages easier                            |
| ComplexHeatmap       | Making Complex Heatmaps in R                                          |
| BiocGenerics         | S4 generic functions for Bioconductor                                 |
| BiocManager          | Access the Bioconductor Project Package Repositories                  |

``` r
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

pkgs <- c("phyloseq", "SummarizedExperiment", "Biobase", "devtools", 
          "ComplexHeatmap", "BiocGenerics", "BiocManager")

BiocManager::install(pkgs)
```

## Stay Updated

Follow me on Twitter for the latest *ggpicrust2* updates:
[![](https://img.shields.io/twitter/follow/CafferyYang?style=social)](https://twitter.com/CafferyYang)

## Workflow

The easiest way to analyze the PICRUSt2 output is using ggpicrust2()
function. The entire pipeline can be run with ggpicrust2() function.

ggpicrust2() integrates ko abundance to kegg pathway abundance
conversion, annotation of pathway, differential abundance (DA) analysis,
DA results visualization. When you have trouble running ggpicrust2(),
you can debug it by running a separate function, which will greatly
increase the speed of your analysis and visualization.

![](https://github.moeyy.xyz/https://raw.githubusercontent.com/cafferychen777/ggpicrust2_paper/main/paper_figure/workflow_1.png)

![](https://github.moeyy.xyz/https://raw.githubusercontent.com/cafferychen777/ggpicrust2_paper/main/paper_figure/workflow_2.png)

### ggpicrust2()

You can download the example dataset from the provided [Google Drive
link](https://drive.google.com/drive/folders/1on4RKgm9NkaBCykMCCRvVJuEJeNVVqAF?usp=share_link)
or use the dataset included in the package.

``` r
# If you want to analyze the abundance of KEGG pathways instead of KO within the pathway, please set `ko_to_kegg` to TRUE.
# KEGG pathways typically have more descriptive explanations.

library(readr)
library(ggpicrust2)
library(tibble)
library(tidyverse)
library(ggprism)
library(patchwork)

# Load necessary data: abundance data and metadata
abundance_file <- "path/to/your/abundance_file.tsv"
metadata <- read_delim(
    "path/to/your/metadata.txt",
    delim = "\t",
    escape_double = FALSE,
    trim_ws = TRUE
)

# Run ggpicrust2 with input file path
results_file_input <- ggpicrust2(file = abundance_file,
                                 metadata = metadata,
                                 group = "your_group_column",
                                 pathway = "KO",
                                 daa_method = "LinDA",
                                 ko_to_kegg = TRUE,
                                 order = "pathway_class",
                                 p_values_bar = TRUE,
                                 x_lab = "pathway_name")

# Run ggpicrust2 with imported data.frame
abundance_data <- read_delim(abundance_file, delim = "\t", col_names = TRUE, trim_ws = TRUE)

# Run ggpicrust2 with input data
results_data_input <- ggpicrust2(data = abundance_data,
                                 metadata = metadata,
                                 group = "your_group_column",
                                 pathway = "KO",
                                 daa_method = "LinDA",
                                 ko_to_kegg = TRUE,
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
```

### If an error occurs with ggpicrust2, please use the following workflow.

``` r
#If you want to analyze kegg pathway abundance instead of ko within the pathway. You should turn ko_to_kegg to TRUE.
#The kegg pathway typically have the more explainable description.
#metadata should be tibble.
library(readr)
library(ggpicrust2)
library(tibble)
library(tidyverse)
library(ggprism)
library(patchwork)
metadata <-
  read_delim(
    "path/to/your/metadata.txt",
    delim = "\t",
    escape_double = FALSE,
    trim_ws = TRUE
  )

kegg_abundance <-
  ko2kegg_abundance(
    "path/to/your/pred_metagenome_unstrat.tsv"
  )

group <- "Environment"

daa_results_df <-
  pathway_daa(
    abundance = kegg_abundance,
    metadata = metadata,
    group = group,
    daa_method = "ALDEx2",
    select = NULL,
    reference = NULL
  )

#if you are using LinDA, limme voom and Maaslin2, please specify a reference just like followings code.
# daa_results_df <- pathway_daa(abundance = kegg_abundance, metadata = metadata, group = group, daa_method = "LinDA", select = NULL, reference = "Harvard BRI")

daa_sub_method_results_df <-
  daa_results_df[daa_results_df$method == "ALDEx2_Kruskal-Wallace test", ]

daa_annotated_sub_method_results_df <-
  pathway_annotation(pathway = "KO",
                     daa_results_df = daa_sub_method_results_df,
                     ko_to_kegg = TRUE)

Group <-
  metadata$Environment # column which you are interested in metadata

# select parameter format in pathway_error() is c("ko00562", "ko00440", "ko04111", "ko05412", "ko00310", "ko04146", "ko00600", "ko04142", "ko00604", "ko04260", "ko04110", "ko04976", "ko05222", "ko05416", "ko00380", "ko05322", "ko00625", "ko00624", "ko00626", "ko00621")

daa_results_list <-
  pathway_errorbar(
    abundance = kegg_abundance,
    daa_results_df = daa_annotated_sub_method_results_df,
    Group = Group,
    p_values_threshold = 0.05,
    order = "pathway_class",
    select = NULL,
    ko_to_kegg = TRUE,
    p_value_bar = TRUE,
    colors = NULL,
    x_lab = "pathway_name"
  )
# x_lab can be set "feature" which represents kegg pathway id.

# If you want to analysis the EC. MetaCyc. KO without conversions. You should turn ko_to_kegg to FALSE.
library(readr)
library(ggpicrust2)
library(tibble)
library(tidyverse)
library(ggprism)
library(patchwork)
metadata <-
  read_delim(
    "path/to/your/metadata.txt",
    delim = "\t",
    escape_double = FALSE,
    trim_ws = TRUE
  )

#ko_abundance should be data.frame instead of tibble
#ID should be in the first column
ko_abundance <-
  read.delim(
    "path/to/your/pred_metagenome_unstrat.tsv"
  )
#sometimes there are function, pathway or something else. Change function. when needed
rownames(ko_abundance) <- ko_abundance[,1]
ko_abundance <- ko_abundance[, -1]

group <- "Environment"

daa_results_df <-
  pathway_daa(
    abundance = ko_abundance,
    metadata = metadata,
    group = group,
    daa_method = "ALDEx2",
    select = NULL,
    reference = NULL
  )

#if you are using LinDA, limme voom and Maaslin2, please specify a reference just like followings code.
# daa_results_df <- pathway_daa(abundance = kegg_abundance, metadata = metadata, group = group, daa_method = "LinDA", select = NULL, reference = "Harvard BRI")

daa_sub_method_results_df <-
  daa_results_df[daa_results_df$method == "ALDEx2_Kruskal-Wallace test", ]

daa_annotated_sub_method_results_df <-
  pathway_annotation(pathway = "KO",
                     daa_results_df = daa_sub_method_results_df,
                     ko_to_kegg = FALSE)

Group <-
  metadata$Environment # column which you are interested in metadata

# select parameter format in pathway_error() is c("K00001", "K00002", "K00003", "K00004")

daa_results_list <-
  pathway_errorbar(
    abundance = ko_abundance,
    daa_results_df = daa_annotated_sub_method_results_df,
    Group = Group,
    p_values_threshold = 0.05,
    order = "pathway_class",
    select = NULL,
    ko_to_kegg = FALSE,
    p_value_bar = TRUE,
    colors = NULL,
    x_lab = "description"
  )
```

## Output

The typical output of the ggpicrust2 is like this.

![](https://raw.iqiq.io/cafferychen777/ggpicrust2_paper/main/paper_figure/figure1.jpg)

## function details

### ko2kegg_abundance()

KEGG Orthology(KO) is a classification system developed by the Kyoto
Encyclopedia of Genes and Genomes (KEGG) data-base(Kanehisa et al.,
2022). It uses a hierarchical structure to classify enzymes based on the
reactions they catalyze. To better understand pathways’ role in
different groups and classify the pathways, the KO abundance table needs
to be converted to KEGG pathway abundance. But PICRUSt2 removes the
function from PICRUSt. ko2kegg_abundance() can help convert the table.

``` r
# Sample usage of the ko2kegg_abundance function

# Assume that the KO abundance table is stored in a file named "ko_abundance.tsv"

ko_abundance_file <- "ko_abundance.tsv"

kegg_abundance <- ko2kegg_abundance(ko_abundance_file)

# The resulting kegg_abundance data frame can now be used for further analysis and visualization.
```

### pathway_daa()

Differential abundance(DA) analysis plays a major role in PICRUSt2
downstream analysis. pathway_daa() integrates almost all DA methods
applicable to the predicted functional profile which there excludes
ANCOM and ANCOMBC. It includes
[ALDEx2](https://bioconductor.riken.jp/packages/3.9/bioc/html/ALDEx2.html)(Fernandes
et al., 2013),
[DESeq2](https://bioconductor.org/packages/release/bioc/html/DESeq2.html)(Love
et al., 2014),
[Maaslin2](https://www.bioconductor.org/packages/release/bioc/html/Maaslin2.html)(Mallick
et al., 2021),
[LinDA](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-022-02655-5)(Zhou
et al., 2022),
[edgeR](https://bioconductor.org/packages/release/bioc/html/edgeR.html)(Robinson
et al., 2010) , [limma
voom](https://ucdavis-bioinformatics-training.github.io/2018-June-RNA-Seq-Workshop/thursday/DE.html)(Ritchie
et al., 2015),
[metagenomeSeq](https://www.bioconductor.org/packages/release/bioc/html/metagenomeSeq.html#:~:text=metagenomeSeq%20is%20designed%20to%20address,the%20testing%20of%20feature%20correlations.)(Paulson
et al., 2013),
[Lefser](http://bioconductor.org/packages/release/bioc/html/lefser.html)(Segata
et al., 2011).

``` r
#the abundance table is better to be data.frame rather than tibble
#you can use ko2_kegg_abundance output
abundance <- ko2kegg_abundance(ko_abundance_file_path)
metadata <-
  read_delim(
    "~/Microbiome/C9orf72/Code And Data/new_metadata.txt",
    delim = "\t",
    escape_double = FALSE,
    trim_ws = TRUE
  )
group <- "Environment"
#daa_method default is "ALDex2"
daa_results_df <- pathway_daa(abundance = abundance,
           metadata = metadata,
           group = group ,
           daa_method = "limma voom",
           select = NULL,
           p.adjust = "BH",
           reference = NULL)
#If you group levels >3 and want to use the LinDA, limme voom and Maaslin2, you should give a reference.
metadata <-
  read_delim(
    "~/Microbiome/C9orf72/Code And Data/new_metadata.txt",
    delim = "\t",
    escape_double = FALSE,
    trim_ws = TRUE
  )
group <- "Group"
daa_results_df <- pathway_daa(abundance = abundance,
           metadata = metadata,
           group = group ,
           daa_method = "limma voom",
           select = NULL,
           p.adjust = "BH",
           reference = "Harvard BRI")
```

### compare_daa_results() &

``` r
library(ggpicrust2)
library(tidyverse)
data("metacyc_abundance")
data("metadata")

# Run pathway_daa function for multiple methods
methods <- c("ALDEx2", "DESeq2", "edgeR")
daa_results_list <- lapply(methods, function(method) {
  pathway_daa(abundance = metacyc_abundance %>% column_to_rownames("pathway"), metadata = metadata, group = "Environment", daa_method = method)
})

# Compare results across different methods
comparison_results <- compare_daa_results(daa_results_list = daa_results_list, method_names = c("ALDEx2_Welch's t test","ALDEx2_Wilcoxon rank test","DESeq2", "edgeR"))
```

### pathway_annotation()

**If you are in China and you are using kegg pathway annotation, Please
make sure your internet can break through the firewall.**

``` r
daa_annotated_results_df <- pathway_annotation(pathway = "KO", daa_results_df = daa_results_df, ko_to_kegg = TRUE)

daa_annotated_results_df <- pathway_annotation(pathway = "KO", daa_results_df = daa_results_df, ko_to_kegg = FALSE)

daa_annotated_results_df <- pathway_annotation(pathway = "EC", daa_results_df = daa_results_df, ko_to_kegg = FALSE)

daa_annotated_results_df <- pathway_annotation(pathway = "MetaCyc", daa_results_df = daa_results_df, ko_to_kegg = FALSE)
```

### pathway_errorbar()

``` r
pathway_errorbar(abundance = abundance,
           daa_results_df = daa_annotated_results_df,
           Group = metadata$Environment,
           ko_to_kegg = TRUE,
           p_values_threshold = 0.05,
           order = "pathway_class",
           select = NULL,
           p_value_bar = TRUE,
           colors = NULL,
           x_lab = "pathway_name")

# If you want to analysis the EC. MetaCyc. KO without conversions. 
pathway_errorbar(abundance = abundance,
           daa_results_df = daa_annotated_results_df,
           Group = metadata$Environment,
           ko_to_kegg = FALSE,
           p_values_threshold = 0.05,
           order = "group",
           select = NULL,
           p_value_bar = TRUE,
           colors = NULL,
           x_lab = "description")
```

### pathway_heatmap()

In this section, we will demonstrate how to create a pathway heatmap
using the `pathway_heatmap` function in the ggpicrust2 package. This
function visualizes the relative abundance of pathways in different
samples.

First, we will create some example functional pathway abundance data and
metadata:

``` r
# Create example functional pathway abundance data
abundance_example <- matrix(rnorm(30), nrow = 10, ncol = 3)
rownames(abundance_example) <- paste0("Sample", 1:10)
colnames(abundance_example) <- c("PathwayA", "PathwayB", "PathwayC")

# Create example metadata
# Please change your sample id's column name to sample_name
metadata_example <- data.frame(sample_name = rownames(abundance_example),
                               group = factor(rep(c("Control", "Treatment"), each = 5)))
```

Next, we will use the pathway_heatmap function to create a heatmap:

``` r
# Create a heatmap
heatmap_plot <- ggpicrust2::pathway_heatmap(t(abundance_example), metadata_example, "group")
```

Finally, we will display the resulting heatmap:

``` r
print(heatmap_plot)
```

### pathway_pca()

In this section, we will demonstrate how to perform Principal Component
Analysis (PCA) on functional pathway abundance data and create
visualizations of the PCA results using the `pathway_pca` function in
the ggpicrust2 package.

First, we will create some example functional pathway abundance data and
metadata:

``` r
# Create example functional pathway abundance data
abundance_example <- data.frame(A = rnorm(10), B = rnorm(10), C = rnorm(10))

# Create example metadata
metadata_example <- tibble::tibble(sample_id = 1:10,
                                   group = factor(rep(c("Control", "Treatment"), each = 5)))
```

Next, we will use the pathway_pca function to perform PCA and create
visualizations:

``` r
# Perform PCA and create visualizations
pca_plot <- ggpicrust2::pathway_pca(t(abundance_example), metadata_example, "group")
```

Finally, we will display the resulting PCA plot:

``` r
print(pca_plot)
```

### compare_metagenome_results

``` r
set.seed(123)
# First metagenome
metagenome1 <- abs(matrix(rnorm(1000), nrow = 100, ncol = 10))
rownames(metagenome1) <- paste0("KO", 1:100)
colnames(metagenome1) <- paste0("sample", 1:10)
# Second metagenome
metagenome2 <- abs(matrix(rnorm(1000), nrow = 100, ncol = 10))
rownames(metagenome2) <- paste0("KO", 1:100)
colnames(metagenome2) <- paste0("sample", 1:10)
# Put the metagenomes into a list
metagenomes <- list(metagenome1, metagenome2)
# Define names
names <- c("metagenome1", "metagenome2")
# Call the function
results <- compare_metagenome_results(metagenomes, names)
# Print the correlation matrix
print(results$correlation$cor_matrix)
# Print the p-value matrix
print(results$correlation$p_matrix)
```

## Share

\[![Twitter](https://twitter.com/intent/tweet?url=https%253A%252F%252Fgithub.com%252Fcafferychen777%252Fggpicrust2&text=Check%20out%20this%20awesome%20package%20on%20GitHub!)

[![Facebook](https://img.shields.io/badge/Share_on-Facebook-1877F2?logo=facebook&style=social)](https://www.facebook.com/sharer/sharer.php?u=https%3A%2F%2Fgithub.com%2Fcafferychen777%2Fggpicrust2&quote=Check%20out%20this%20awesome%20package%20on%20GitHub%21)

[![LinkedIn](https://img.shields.io/badge/Share_on-LinkedIn-0077B5?logo=linkedin&style=social)](https://www.linkedin.com/shareArticle?mini=true&url=https%3A%2F%2Fgithub.com%2Fcafferychen777%2Fggpicrust2&title=Check%20out%20this%20awesome%20package%20on%20GitHub%21)

[![Reddit](https://img.shields.io/badge/Share_on-Reddit-FF4500?logo=reddit&style=social)](https://www.reddit.com/submit?url=https%3A%2F%2Fgithub.com%2Fcafferychen777%2Fggpicrust2&title=Check%20out%20this%20awesome%20package%20on%20GitHub%21)

## FAQ

### Issue 1: pathway_errorbar error

When using `pathway_errorbar` with the following parameters:

``` r
pathway_errorbar(abundance = abundance,
                 daa_results_df = daa_results_df,
                 Group = metadata$Environment,
                 ko_to_kegg = TRUE,
                 p_values_threshold = 0.05,
                 order = "pathway_class",
                 select = NULL,
                 p_value_bar = TRUE,
                 colors = NULL,
                 x_lab = "pathway_name")
```

You may encounter an error:

    Error in `ggplot_add()`:
    ! Can't add `e2` to a <ggplot> object.
    Run `rlang::last_trace()` to see where the error occurred.

Make sure you have the `patchwork` package loaded:

``` r
library(patchwork)
```

### Issue 2: guide_train.prism_offset_minor error

You may encounter an error with `guide_train.prism_offset_minor`:

    Error in guide_train.prism_offset_minor(guide, panel_params[[aesthetic]]) : 
      No minor breaks exist, guide_prism_offset_minor needs minor breaks to work

Ensure that the `ggprism` package is loaded:

``` r
library(ggprism)
```

### Issue 3: SSL certificate problem

When encountering the following error:

    SSL peer certificate or SSH remote key was not OK: [rest.kegg.jp] SSL certificate problem: certificate has expired

If you are in China, make sure your computer network can bypass the
firewall.

### Issue 4: Bad Request (HTTP 400)

When encountering the following error:

    Error in .getUrl(url, .flatFileParser) : Bad Request (HTTP 400).

Please restart R session.

### Issue 5: Error in grid.Call(C_textBounds, as.graphicsAnnot(xlabel),x$x, x$y, :

When encountering the following error:

    Error in grid.Call(C_textBounds, as.graphicsAnnot(xlabel),x$x, x$y, :

Please having some required fonts installed. You can refer to this
[thread](https://stackoverflow.com/questions/71362738/r-error-in-grid-callc-textbounds-as-graphicsannotxlabel-xx-xy-polygo).

### Issue 6: Visualization becomes cluttered when there are more than 30 features of statistical significance.

When faced with this issue, consider the following solutions:

**Solution 1: Utilize the ‘select’ parameter**

The ‘select’ parameter allows you to specify which features you wish to
visualize. Here’s an example of how you can apply this in your code:

    ggpicrust2::pathway_errorbar(
      abundance = kegg_abundance,
      daa_results_df = daa_results_df_annotated,
      Group = metadata$Day,
      p_values_threshold = 0.05,
      order = "pathway_class",
      select = c("ko05340", "ko00564", "ko00680", "ko00562", "ko03030", "ko00561", "ko00440", "ko00250", "ko00740", "ko04940", "ko00010", "ko00195", "ko00760", "ko00920", "ko00311", "ko00310", "ko04146", "ko00600", "ko04141", "ko04142", "ko00604", "ko04260", "ko00909", "ko04973", "ko00510", "ko04974"),
      ko_to_kegg = TRUE,
      p_value_bar = FALSE,
      colors = NULL,
      x_lab = "pathway_name"
    )

**Solution 2: Limit to the Top 20 features**

If there are too many significant features to visualize effectively, you
might consider limiting your visualization to the top 20 features with
the smallest adjusted p-values:

    daa_results_df_annotated <- daa_results_df_annotated[!is.na(daa_results_df_annotated$pathway_name),]

    daa_results_df_annotated$p_adjust <- round(daa_results_df_annotated$p_adjust,5)

    low_p_feature <- daa_results_df_annotated[order(daa_results_df_annotated$p_adjust), ]$feature[1:20]


    p <- ggpicrust2::pathway_errorbar(
      abundance = kegg_abundance,
      daa_results_df = daa_results_df_annotated,
      Group = metadata$Day,
      p_values_threshold = 0.05,
      order = "pathway_class",
      select = low_p_feature,
      ko_to_kegg = TRUE,
      p_value_bar = FALSE,
      colors = NULL,
      x_lab = "pathway_name"
    )

### Issue 7: There are no statistically significant biomarkers

If you are not finding any statistically significant biomarkers in your
analysis, there could be several reasons for this:

1.  **The true difference between your groups is small or
    non-existent.** If the microbial communities or pathways you’re
    comparing are truly similar, then it’s correct and expected that you
    won’t find significant differences.

2.  **Your sample size might be too small to detect the differences.**
    Statistical power, the ability to detect differences if they exist,
    increases with sample size.

3.  **The variation within your groups might be too large.** If there’s
    a lot of variation in microbial communities within a single group,
    it can be hard to detect differences between groups.

Here are a few suggestions:

1.  **Increase your sample size**: If possible, adding more samples to
    your analysis can increase your statistical power, making it easier
    to detect significant differences.

2.  **Decrease intra-group variation**: If there’s a lot of variation
    within your groups, consider whether there are outliers or subgroups
    that are driving this variation. You might need to clean your data,
    or to stratify your analysis to account for these subgroups.

3.  **Change your statistical method or adjust parameters**: Depending
    on the nature of your data and your specific question, different
    statistical methods might be more or less powerful. If you’re
    currently using a parametric test, consider using a non-parametric
    test, or vice versa. Also, consider whether adjusting the parameters
    of your current test might help.

Remember, not finding significant results is also a result and can be
informative, as it might indicate that there are no substantial
differences between the groups you’re studying. It’s important to
interpret your results in the context of your specific study and not to
force statistical significance where there isn’t any.

With these strategies, you should be able to create a more readable and
informative visualization, even when dealing with a large number of
significant features.

## Author’s other projects

[MicrobiomeGallery](https://a95dps-caffery-chen.shinyapps.io/MicrobiomeGallery/),
a web-based platform for sharing microbiome data visualization code and
datasets. （It’s still under construction）

![](https://raw.githubusercontent.com/cafferychen777/ggpicrust2_paper/main/paper_figure/MicrobiomeGallery_preview.jpg)
