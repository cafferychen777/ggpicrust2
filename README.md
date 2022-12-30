
<!-- README.md is generated from README.Rmd. Please edit that file -->

# ggpicrust2

<!-- badges: start -->

[![Codecov test
coverage](https://codecov.io/gh/cafferychen777/ggpicrust2/branch/main/graph/badge.svg)](https://app.codecov.io/gh/cafferychen777/ggpicrust2?branch=main)
[![R-CMD-check](https://github.com/cafferychen777/ggpicrust2/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/cafferychen777/ggpicrust2/actions/workflows/R-CMD-check.yaml)

<!-- badges: end -->

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

![](/Users/apple/Microbiome/ggpicrust2/ggpicrust2/docs/reference/figures/pathway_errorbar.pdf)

You’ll still need to render `README.Rmd` regularly, to keep `README.md`
up-to-date. `devtools::build_readme()` is handy for this. You could also
use GitHub Actions to re-render `README.Rmd` every time you push. An
example workflow can be found here:
<https://github.com/r-lib/actions/tree/v1/examples>.

You can also embed plots, for example:

In that case, don’t forget to commit and push the resulting figure
files, so they display on GitHub and CRAN.
