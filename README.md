
<!-- README.md is generated from README.Rmd. Please edit that file -->

# ggpicrust2 Community User Tutorial Competition

As an open source developer, I am excited to announce a competition for
the best tutorials on using *ggpicrust2* - an R package designed for
analyzing and interpreting the results of PICRUSt2 functional
prediction.

We welcome submissions in any format, and the first place winner will
receive a full set of *ggpicrust2* community swag, including a logo
t-shirt, tote bag, and hat. The second, third, and fourth place winners
will each receive a cultural t-shirt with the *ggpicrust2* logo.

The competition will run from April 3 to June 3, and to enter, simply
use Rstudio’s rpub function to publish your tutorial online, and then
send an email to [Caffery](cafferychen7850@gmail.com) with the link.

Whether you are a researcher, data scientist, or bioinformatician,
*ggpicrust2* can help you better understand the underlying biological
processes and mechanisms at play in your PICRUSt2 output data. So, join
us in building an even better *ggpicrust2* community!

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

## Table of Contents

- [Citation](#citation)
- [Installation](#installation)
- [Stay Updated](#stay-updated)
- [Workflow](#workflow)
- [Output](#output)
- [Function Details](#function-details)
  - [ko2kegg_abundance()](#ko2kegg_abundance)
  - [pathway_daa()](#pathway_daa)
  - [pathway_annotation()](#pathway_annotation)
  - [pathway_errorbar()](#pathway_errorbar)
  - [pathway_heatmap()](#pathway_heatmap)
  - [pathway_pca()](#pathway_pca)
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

You can also find the paper on
[ResearchGate](https://www.researchgate.net/publication/369380349_ggpicrust2_an_R_package_for_PICRUSt2_predicted_functional_profile_analysis_and_visualization).

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
increase the speed of your analysis and visualization

![](https://github.moeyy.xyz/https://raw.githubusercontent.com/cafferychen777/ggpicrust2_paper/main/paper_figure/workflow_1.png)

![](https://github.moeyy.xyz/https://raw.githubusercontent.com/cafferychen777/ggpicrust2_paper/main/paper_figure/workflow_2.png)

### ggpicrust2()

``` r
#If you want to analysis kegg pathway abundance instead of ko within the pathway. You should turn ko_to_kegg to TRUE.
#The kegg pathway typically have the more explainable description.
metadata <-
  read_delim(
    "~/Microbiome/C9orf72/Code And Data/new_metadata.txt",
    delim = "\t",
    escape_double = FALSE,
    trim_ws = TRUE
  )
group <- "Enviroment"
daa_results_list <-
  ggpicrust2(
    file = "/Users/apple/Microbiome/C9orf72/Code And Data/picrust2_out/KO_metagenome_out/pred_metagenome_unstrat.tsv/pred_metagenome_unstrat.tsv",
    metadata = metadata,
    group = "Enviroment",
    pathway = "KO",
    daa_method = "LinDA",
    p_values_bar = TRUE,
    p.adjust = "BH",
    ko_to_kegg = TRUE,
    order = "pathway_class"
    select = NULL,
    reference = NULL # If your metadata[,group] has more than two levels, please specify a reference.
  )



#If you want to analysis the EC. MetaCyc. KO without conversions. You should turn ko_to_kegg to FALSE.
metadata <-
  read_delim(
    "~/Microbiome/C9orf72/Code And Data/new_metadata.txt",
    delim = "\t",
    escape_double = FALSE,
    trim_ws = TRUE
  )
group <- "Enviroment"
daa_results_list <-
  ggpicrust2(
    file = "//Users/apple/Microbiome/C9orf72/Code And Data/picrust2_out/EC_metagenome_out/pred_metagenome_unstrat.tsv/pred_metagenome_unstrat.tsv",
    metadata = metadata,
    group = "Enviroment",
    pathway = "EC",
    daa_method = "LinDA",
    p_values_bar = TRUE,
    order = "group",
    ko_to_kegg = FALSE,
    x_lab = "description",
    p.adjust = "BH",
    select = NULL,
    reference = NULL
  )
```

### If an error occurs with ggpicrust2, please use the following workflow.

``` r
#If you want to analysis kegg pathway abundance instead of ko within the pathway. You should turn ko_to_kegg to TRUE.
#The kegg pathway typically have the more explainable description.
#metadata should be tibble.
metadata <-
  read_delim(
    "~/Microbiome/C9orf72/Code And Data/new_metadata.txt",
    delim = "\t",
    escape_double = FALSE,
    trim_ws = TRUE
  )

kegg_abundance <-
  ko2kegg_abundance(
    "/Users/apple/Downloads/pred_metagenome_unstrat.tsv/pred_metagenome_unstrat.tsv"
  )

group <- "Enviroment"

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
# daa_results_df <- pathway_daa(abundance = kegg_abundance, metadata = metadata, group = group, daa_method = "ALDEx2", select = NULL, reference = "Harvard BRI")

daa_sub_method_results_df <-
  daa_results_df[daa_results_df$method == "ALDEx2_Kruskal-Wallace test", ]

daa_annotated_sub_method_results_df <-
  pathway_annotation(pathway = "KO",
                     daa_results_df = daa_sub_method_results_df,
                     ko_to_kegg = TRUE)

Group <-
  metadata$Enviroment # column which you are interested in metadata

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
metadata <-
  read_delim(
    "~/Microbiome/C9orf72/Code And Data/new_metadata.txt",
    delim = "\t",
    escape_double = FALSE,
    trim_ws = TRUE
  )

#ko_abundance should be data.frame instead of tibble
#ID should be in the first column
ko_abundance <-
  read.delim(
    "/Users/apple/Downloads/pred_metagenome_unstrat.tsv/pred_metagenome_unstrat.tsv"
  )
#sometimes there are function, pathway or something else. Change function. when needed
rownames(ko_abundance) <- ko_abundance$function.
ko_abundance <- ko_abundance[, -1]

group <- "Enviroment"

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
# daa_results_df <- pathway_daa(abundance = kegg_abundance, metadata = metadata, group = group, daa_method = "ALDEx2", select = NULL, reference = "Harvard BRI")

daa_sub_method_results_df <-
  daa_results_df[daa_results_df$method == "ALDEx2_Kruskal-Wallace test", ]

daa_annotated_sub_method_results_df <-
  pathway_annotation(pathway = "KO",
                     daa_results_df = daa_sub_method_results_df,
                     ko_to_kegg = FALSE)

Group <-
  metadata$Enviroment # column which you are interested in metadata

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

![](https://cdn.jsdelivr.ren/gh/cafferychen777/ggpicrust2_paper@main/paper_figure/figure1.jpg)

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
[DEseq2](https://bioconductor.org/packages/release/bioc/html/DESeq2.html)(Love
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
[lefser](http://bioconductor.org/packages/release/bioc/html/lefser.html)(Segata
et al., 2011).

``` r
#the abundance table is better to be data.frame rather than tibble
#you can use ko2_kegg_abundance output
abundance <- ko2kegg_abundance(ko_abundance_file)
metadata <-
  read_delim(
    "~/Microbiome/C9orf72/Code And Data/new_metadata.txt",
    delim = "\t",
    escape_double = FALSE,
    trim_ws = TRUE
  )
group <- "Enviroment"
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

### pathway_annotation()

**If you are in China and you are using kegg pathway annotation, Please
make sure your internet can break through the firewall.**

``` r
daa_results_df <- pathway_annotation(pathway = "KO", daa_results_df = daa_results_df, ko_to_kegg = TRUE)
```

### pathway_errorbar()

``` r
pathway_errorbar(abundance = abundance,
           daa_results_df = daa_results_df,
           Group = metadata$Enviroment,
           ko_to_kegg = TRUE,
           p_values_threshold = 0.05,
           order = "pathway_class",
           select = NULL,
           p_value_bar = TRUE,
           colors = NULL,
           x_lab = NULL)
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

## Share

[![Twitter](https://img.shields.io/twitter/url?url=https%253A%252F%252Fgithub.com%252Fcafferychen777%252Fggpicrust2&style=social)](https://twitter.com/intent/tweet?url=https%3A%2F%2Fgithub.com%2Fcafferychen777%2Fggpicrust2&text=Check%20out%20this%20awesome%20package%20on%20GitHub%21)

[![Facebook](https://img.shields.io/badge/Share_on-Facebook-1877F2?logo=facebook&style=social)](https://www.facebook.com/sharer/sharer.php?u=https%3A%2F%2Fgithub.com%2Fcafferychen777%2Fggpicrust2&quote=Check%20out%20this%20awesome%20package%20on%20GitHub%21)

[![LinkedIn](https://img.shields.io/badge/Share_on-LinkedIn-0077B5?logo=linkedin&style=social)](https://www.linkedin.com/shareArticle?mini=true&url=https%3A%2F%2Fgithub.com%2Fcafferychen777%2Fggpicrust2&title=Check%20out%20this%20awesome%20package%20on%20GitHub%21)

[![Reddit](https://img.shields.io/badge/Share_on-Reddit-FF4500?logo=reddit&style=social)](https://www.reddit.com/submit?url=https%3A%2F%2Fgithub.com%2Fcafferychen777%2Fggpicrust2&title=Check%20out%20this%20awesome%20package%20on%20GitHub%21)

## Author’s other projects

[MicrobiomeGallery](https://a95dps-caffery-chen.shinyapps.io/MicrobiomeGallery/),
a web-based platform for sharing microbiome data visualization code and
datasets. （It’s still under construction）

![](https://raw.githubusercontent.com/cafferychen777/ggpicrust2_paper/main/paper_figure/MicrobiomeGallery_preview.jpg)
