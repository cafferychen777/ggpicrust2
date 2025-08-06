# ggpicrust2 2.4.0

## New Features

### Enhanced Multi-Grouping Support for pathway_heatmap()

* **Added secondary_groups parameter (#171)**:
  - Enables multi-level grouping in pathway heatmaps
  - Supports nested faceting with multiple grouping variables
  - Maintains full backward compatibility with existing code
  - Automatic color generation for complex grouping structures

* **Deprecated facet_by parameter**:
  - Replaced with more flexible secondary_groups parameter
  - Shows deprecation warning with migration guidance
  - Will be removed in future major version

* **Enhanced Examples and Documentation**:
  - Added comprehensive examples for multi-grouping usage
  - Migration guide from facet_by to secondary_groups
  - Updated function documentation with new parameter descriptions

* **Comprehensive Testing**:
  - Added unit tests for all new functionality
  - Integration tests with real data
  - Backward compatibility validation
  - Error handling and edge case testing

## Bug Fixes

* **Improved color allocation for multi-level grouping**:
  - Smart color generation based on grouping complexity
  - Automatic color repetition when insufficient colors provided
  - Better handling of color assignment in nested structures

---

# ggpicrust2 2.3.3

## Major Bug Fixes and Improvements

### DAA Methods Comprehensive Fixes

* **Fixed limma voom compatibility with compare_daa_results (#163)**:
  - Added missing group1 and group2 columns for two-group comparisons
  - Ensures consistent output format with other DAA methods
  - Resolves "Unknown comparison type" error when using limma voom with compare_daa_results

* **Fixed Maaslin2 multiple critical issues (#164)**:
  - Fixed sample name matching between abundance matrix and metadata
  - Implemented intelligent feature name matching to handle Maaslin2's hyphen-to-dot conversion
  - Resolves "Unable to find samples in data and metadata files" error
  - Eliminates NA p-values caused by feature name mismatches
  - Robust handling of various feature naming conventions (hyphens, dots, underscores)

* **Enhanced DESeq2 robustness**:
  - Added fallback mechanism for dispersion estimation failures
  - Uses gene-wise estimates when standard dispersion estimation fails
  - Improved compatibility with small sample sizes and low-variance data

* **Standardized Lefser implementation**:
  - Added lefser package to method detection list
  - Standardized output format to match other DAA methods
  - Returns results for all features, not just significant ones
  - Converts effect scores to p-values for consistency

### Comprehensive Testing and Validation

* **100% success rate achieved for all 8 DAA methods**:
  - ALDEx2, DESeq2, edgeR, limma voom, metagenomeSeq, Maaslin2, LinDA, Lefser
  - All methods now fully compatible with compare_daa_results function
  - Consistent output format across all methods
  - Extensive testing with various feature naming patterns and edge cases

### Technical Improvements

* Enhanced error handling and robustness across all DAA methods
* Improved feature name matching algorithms
* Better sample-metadata alignment validation
* Defensive programming against edge cases with factor levels

# ggpicrust2 2.3.2

## Bug Fixes

* Fixed MetaCyc pathway annotation NA description issue (#154):
  - Standardized MetaCyc reference data column names from 'X1'/'X2' to 'id'/'description'
  - Applied fix in all three loading paths of load_reference_data function
  - Resolves column name mismatch that caused all MetaCyc annotations to return NA
  - Tested with 100% success rate on sample data
  - Maintains backward compatibility with KO and EC pathway types

# ggpicrust2 2.3.1

## Bug Fixes

* Fixed MetaCyc reference data loading issue:
  - Enhanced the file search mechanism in the `load_reference_data` function
  - Added multiple search paths for reference data files
  - Improved error messages with more diagnostic information
  - Fixed "Reference data file not found" error that some users encountered

# ggpicrust2 2.3.0

## Bug Fixes

* Fixed NA handling in pathway_errorbar function:
  - Improved handling of NA values in the feature column
  - Added robust error checking for group ordering option
  - Prevents errors when processing MetaCyc pathway data with missing annotations

* Enhanced pathway_pca function to handle zero variance data:
  - Automatically detects and filters out columns (samples) with zero variance
  - Automatically detects and filters out rows (pathways) with zero variance
  - Updates metadata to match remaining samples after filtering
  - Provides informative warnings about removed samples/pathways
  - Improves error handling with clear diagnostic messages

* Fixed file extension handling in ko2kegg_abundance function:
  - Resolved "the condition has length > 1" error when processing files
  - Improved extension detection for different file types
  - Enhanced robustness for handling various input formats

# ggpicrust2 2.2.2

## Bug Fixes

* Fixed column name compatibility issues in reference data files:
  - Standardized column names in KO_reference.RData and EC_reference.RData
  - Resolved warnings when using pathway_annotation() with ko_to_kegg=FALSE
  - Improved backward compatibility with existing code

# ggpicrust2 2.2.1

## Bug Fixes

* Added aggregate_by_group parameter to pathway_heatmap function:
  - Allows displaying representative samples (e.g., mean) for each group instead of all individual samples
  - Improved handling of NA values in heatmap visualization
  - Added aggregate_fun parameter for customizing the aggregation function

# ggpicrust2 2.1.4

## Reference Data Updates

* Updated reference databases for improved pathway annotation:
  - EC reference data updated from 3,180 to 8,371 entries (163% increase)
  - KO reference data updated from 23,917 to 27,531 unique KO IDs (15.4% increase)
  - These updates provide more comprehensive and accurate pathway annotations

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
