# ggpicrust2 2.1.1

## Bug Fixes

* Improved error handling for KEGG database connections (Issue #138):
  - Added robust error handling for HTTP 404 errors when connecting to the KEGG database
  - Function now continues processing other KO IDs even when some IDs return HTTP 404 errors
  - Added detailed logging about which KO IDs were not found
  - Only throws a fatal error when all KO IDs fail to be processed
  - Includes summary statistics about successful, not found, and error counts

* Fixed the LinDA analysis for multi-group comparisons (Issue #144):
  - Modified the `perform_linda_analysis` function to handle multi-group comparisons correctly
  - The function now creates separate result entries for each pairwise comparison
  - Added the `log2FoldChange` column to the results for effect size information

# ggpicrust2 2.1.0

## Major Changes

* Added Gene Set Enrichment Analysis (GSEA) functionality with the following new functions:
  - `pathway_gsea()`: Performs GSEA analysis, supporting KEGG, MetaCyc, and GO pathways
  - `visualize_gsea()`: Creates visualizations of GSEA results, including enrichment plots, dotplots, barplots, network plots, and heatmaps
  - `compare_gsea_daa()`: Compares GSEA and Differential Abundance Analysis (DAA) results
  - `gsea_pathway_annotation()`: Adds pathway annotations to GSEA results
  - `ggpicrust2_extended()`: Provides integrated analysis combining ggpicrust2 and GSEA functionality

* Improved network and heatmap visualization capabilities with richer parameter options and better error handling

* Added preliminary support for MetaCyc and GO pathways

* Fixed various bugs and optimized code structure

# ggpicrust2 2.0.1

## Major Changes

* Fixed a bug in the pathway_annotation function, resolving issues caused by changes in the KEGG API response structure

# ggpicrust2 2.0.0

## 主要变更

* Refactored the package dependencies, moving most Bioconductor packages from Imports to Suggests, reducing mandatory dependencies.

* Added conditional checks to ensure that packages are only used when they are available.

* Fixed the function export issue, ensuring that all public API functions are correctly exported.

* Updated the example code, improving its stability and compatibility.

* Updated the documentation format to comply with the latest CRAN standards.

* Fixed invalid URL links in the README.

* Optimized code quality, removing unused variables.

* Added missing import declarations to ensure package integrity.

# ggpicrust2 1.7.5

# ggpicrust2 1.7.2

# ggpicrust2 1.7.1

# ggpicrust2 1.7.0

# ggpicrust2 1.6.6

# ggpicrust2 1.6.5

# ggpicrust2 1.6.4

# ggpicrust2 1.6.3

# ggpicrust2 1.6.2

# ggpicrust2 1.6.1

# ggpicrust2 1.6.0

# ggpicrust2 1.5.1

# ggpicrust2 1.5.0

# ggpicrust2 1.4.12

# ggpicrust2 1.4.11

# ggpicrust2 1.4.10

# ggpicrust2 1.4.9

# ggpicrust2 1.4.8

# ggpicrust2 1.4.7

# ggpicrust2 1.4.6

# ggpicrust2 1.4.5

# ggpicrust2 1.4.4

# ggpicrust2 1.4.3

# ggpicrust2 1.4.2

# ggpicrust2 1.4.1

# ggpicrust2 1.4.0

# ggpicrust2 1.3.1

# ggpicrust2 1.3.0

# ggpicrust2 1.2.3

# ggpicrust2 1.2.2

# ggpicrust2 1.2.1

# ggpicrust2 1.2.0

# ggpicrust2 1.1.1

# ggpicrust2 1.1.0

# ggpicrust2 1.0.1

# ggpicrust2 1.0.0

# ggpicrust2 0.2.0

* Added a `NEWS.md` file to track changes to the package.
