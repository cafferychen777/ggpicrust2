
<!-- README.md is generated from README.Rmd. Please edit that file -->

# ggpicrust2

ggpicrust2 is a comprehensive package that integrates pathway
name/description annotations, ten of the most advanced differential
abundance (DA) methods, and visualization of DA results. It offers a
comprehensive solution for analyzing and interpreting the results of
picrust2 functional prediction in a seamless and intuitive way. Whether
you are a researcher, data scientist, or bioinformatician, ggpicrust2
can help you better understand the underlying biological processes and
mechanisms at play in your picrust2 output data. So if you are
interested in exploring the output data of picrust2, ggpicrust2 is the
tool you need.

## Installation

You can install the development version of ggpicrust2 from
[GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("cafferychen777/ggpicrust2")
```

We are actively preparing to upload the package to bioconductor.

## Workflow

The easiest way to analyze the picrust2 output is using ggpicrust2()
function. The entire pipeline can be run with ggpicrust2() function.

ggpicrust2() integrates ko abundance to kegg pathway abundance
conversion, annotation of pathway, differential abundance (DA) analysis,
DA results visualization.

![](https://cafferychen777.github.io/ggpicrust2/reference/figures/ggpicrust2%20fig1.jpeg)

![](https://cafferychen777.github.io/ggpicrust2/reference/figures/ggpicrust2%20fig2.jpeg)

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
daa_results_df <-
  ggpicrust2(
    file = "/Users/apple/Microbiome/C9orf72/Code And Data/picrust2_out/KO_metagenome_out/pred_metagenome_unstrat.tsv/pred_metagenome_unstrat.tsv",
    metadata = metadata,
    group = "Enviroment",
    pathway = "KO",
    daa_method = "ALDEx2",
    order = "pathway_class",
    ko_to_kegg = TRUE,
    x_lab = "pathway_name",
    p.adjust = "BH",
    select = NULL,
    reference = NULL
  )
#The visualization will be published in viewer.

#If you want to analysis the EC. MetaCyc. KO without conversions. You should turn ko_to_kegg to FALSE.
metadata <-
  read_delim(
    "~/Microbiome/C9orf72/Code And Data/new_metadata.txt",
    delim = "\t",
    escape_double = FALSE,
    trim_ws = TRUE
  )
group <- "Enviroment"
daa_results_df <-
  ggpicrust2(
    file = "//Users/apple/Microbiome/C9orf72/Code And Data/picrust2_out/EC_metagenome_out/pred_metagenome_unstrat.tsv/pred_metagenome_unstrat.tsv",
    metadata = metadata,
    group = "Enviroment",
    pathway = "EC",
    daa_method = "ALDEx2",
    order = "pathway_class",
    ko_to_kegg = FALSE,
    x_lab = "description",
    p.adjust = "BH",
    select = NULL,
    reference = NULL
  )
#The visualization will be published in viewer.
```

## Output

The typical output of the ggpicrust2 is like this.

![](https://cafferychen777.github.io/ggpicrust2/reference/figures/pathway_errorbar.jpg)

You’ll still need to render `README.Rmd` regularly, to keep `README.md`
up-to-date. `devtools::build_readme()` is handy for this. You could also
use GitHub Actions to re-render `README.Rmd` every time you push. An
example workflow can be found here:
<https://github.com/r-lib/actions/tree/v1/examples>.

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
ANCOM and ANCOMBC. It includes ALDEx2(Fernandes et al., 2013),
DEseq2(Love et al., 2014), Maaslin2(Mallick et al., 2021), Lin-DA(Zhou
et al., 2022), edgeR(Robinson et al., 2010) , limma voom(Ritchie et al.,
2015), metagenomeSeq(Paulson et al., 2013), lefser(Segata et al., 2011).

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
#If you group levels >3 and want to use the LinDA or limma voom, you should give a reference.
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

### pathway_annotation

``` r
daa_results_df <- pathway_annotation(pathway = "KO", daa_results_df = daa_results_df, ko_to_kegg = TRUE)
```

### pathway_errorbar

``` r
pathway_errorbar(abundance = abundance,
           daa_results_df = daa_results_df,
           Group = metadata$Enviroment,
           ko_to_kegg = TRUE,
           p_values_threshold = 0.05,
           order = "group",
           select = NULL,
           p_value_bar = TRUE,
           colors = NULL,
           x_lab = NULL)
```

### pathway_heatmap

pathway_heatmap() can visualize the patterns in PICRUSt2 output data,
which can be useful for identifying trends or highlighting areas of
in-terest.

``` r
abundance <- ko2kegg_abundance(ko_abundance_file)
metadata <-
  read_delim(
    "~/Microbiome/C9orf72/Code And Data/new_metadata.txt",
    delim = "\t",
    escape_double = FALSE,
    trim_ws = TRUE
  )
group <- "Enviroment"
pathway_heatmap(abundance = abundance, metadata = metadata, group = group)
```

### pathway_pca()

pathway_pca() can show the difference after dimensional reduction by
PCA.

``` r
# generate example abundance matrix and metadata dataframe
abundance <- matrix(rnorm(200), ncol = 20)
metadata <- data.frame(Group = rep(c("A", "B"), each = 10),
                       Treatment = rep(c("X", "Y"), 10))

# perform PCA and create plot
pathway_pca_plot <- pathway_pca(abundance = abundance, metadata = metadata, group = "Group")
print(pathway_pca_plot)
```
