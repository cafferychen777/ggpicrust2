# ggpicrust2 2.5.8

## Bug Fixes

* Fixed CRAN pretest failure in `tests/testthat/test-pathway_daa.R` by splitting
  default/core method coverage from optional extended methods that require
  non-mainstream dependencies (e.g., `Maaslin2`).
* Kept extended DAA method coverage available behind
  `GGPICRUST2_RUN_EXTENDED_DAA_TESTS=true`.
* Corrected `metacyc_reference` documentation to match actual data columns
  (`id`, `description`), resolving `codoc` mismatch warnings.

# ggpicrust2 2.5.7

## Bug Fixes

* Added a regression test for `pathway_errorbar()` to explicitly cover the
  `ko_to_kegg = TRUE` + `order = "pathway_class"` path.
* This guards against reintroducing the historical
  `tibble::column_to_rownames()` / `Can't find column '.'` failure mode.

# ggpicrust2 2.5.6

## Breaking Changes

* **Removed `ggpicrust2_extended()` function**:
  - This wrapper function provided minimal value over calling `ggpicrust2()` and `pathway_gsea()` separately
  - Users can achieve the same functionality by calling the individual functions directly
  - This change reduces maintenance burden and improves code clarity

## Major Features

### Covariate Adjustment & Improved Statistical Methods for pathway_gsea() (#193)

* **Added limma camera and fry methods as new GSEA options**:
  - `method = "camera"` (now default): Competitive gene set test using limma's camera function
  - `method = "fry"`: Fast rotation gene set test (self-contained)
  - Both methods account for inter-gene correlations, providing more reliable p-values than preranked GSEA

* **Added covariate adjustment support**:
  - New `covariates` parameter to adjust for confounding factors (age, sex, BMI, etc.)
  - Covariates are incorporated into the design matrix for proper statistical adjustment
  - Essential for microbiome studies where host factors can confound results

* **New parameters**:
  - `covariates`: Character vector of covariate column names from metadata
  - `contrast`: For multi-group comparisons, specify the contrast to test
  - `inter.gene.cor`: Inter-gene correlation for camera method (default: 0.01)

* **Scientific background**:
  - Wu et al. (2012) demonstrated that preranked GSEA methods can produce "spectacularly wrong p-values" due to not accounting for inter-gene correlations
  - The camera and fry methods from limma address this limitation
  - Reference: Wu, D., & Smyth, G. K. (2012). Nucleic Acids Research, 40(17), e133.

* **Backward compatibility**:
  - Existing `fgsea` and `clusterProfiler` methods remain available
  - Added informational message when using preranked methods about p-value reliability

* **New internal functions**:
  - `run_limma_gsea()`: Core implementation for camera/fry methods
  - `build_design_matrix()`: Constructs design matrix with covariates

* **Updated documentation and vignettes**:
  - Method selection guide with comparison table
  - Covariate adjustment examples
  - Updated gsea_analysis.Rmd vignette

### New Visualization Functions

* **Added `pathway_volcano()` function**:
  - Creates publication-quality volcano plots for differential abundance analysis
  - Visualizes both statistical significance (-log10 p-value) and effect size (log2 fold change)
  - Smart label placement using ggrepel to avoid overlapping labels
  - Color-coded significance categories (Up/Down/Not Significant)
  - Automatic handling of NA pathway names (won't display "NA" labels)
  - Handles infinite p-values (when p = 0) gracefully
  - Customizable thresholds, colors, and appearance

* **Added `pathway_ridgeplot()` function**:
  - Creates ridge plots (joy plots) for GSEA results interpretation
  - Shows distribution of gene abundances/fold changes within enriched pathways
  - Color-coded by enrichment direction (Up/Down)
  - Automatic pathway-KO mapping using built-in ko_to_kegg_reference data
  - Supports KEGG and GO pathway types
  - Helps identify whether pathways are predominantly up- or down-regulated
  - Requires ggridges package (added to Suggests)

* **New dependencies added to Suggests**:
  - `ggridges`: Required for ridge plot visualization
  - `ggrepel`: Used for smart label placement in volcano plots

---

# ggpicrust2 2.5.5

## Major Features

### Prokaryote-Specific Pathway Filtering (#191)

* **New `filter_for_prokaryotes` parameter in `ko2kegg_abundance()`**:
  - Defaults to TRUE, automatically filtering out eukaryote-specific pathways
  - Removes biologically irrelevant pathways from bacterial/archaeal analysis:
    - Cancer pathways (overview and specific types)
    - Neurodegenerative diseases (Alzheimer's, Parkinson's, etc.)
    - Substance dependence (addiction pathways)
    - Cardiovascular diseases
    - Endocrine and metabolic diseases (human-specific)
    - Immune diseases (human-specific)
    - Organismal systems (immune, nervous, endocrine, digestive, etc.)
  - Retains prokaryote-relevant pathways:
    - All Metabolism pathways
    - Infectious disease: bacterial (Salmonella, E. coli, Tuberculosis, etc.)
    - Drug resistance: antimicrobial (antibiotic resistance)
    - Genetic/Environmental Information Processing
    - Cellular Processes
  - Set `filter_for_prokaryotes = FALSE` for eukaryotic analysis or to include all pathways
  - Reduces pathway count from ~370 to ~290 for typical bacterial analyses

### Local KEGG Database Implementation (#113)

* **Replaced KEGG API dependency with local database**:
  - Implemented comprehensive local KO-to-KEGG pathway mapping (61,655 mappings)
  - Covers 557 pathways and 27,127 KO IDs
  - Eliminates dependency on external KEGG API
  - Provides 100% data coverage with 0% missing values (vs 84% NA in previous format)
  - Significantly improved performance and reliability

* **Enhanced data structure**:
  - Migrated from wide format (306×326 matrix) to long format (61,655×9 table)
  - Added rich pathway metadata including hierarchical classification (Level1-3)
  - Includes pathway names, KO descriptions, and EC numbers
  - Optimized with fast lookup index for O(N) performance

* **Updated functions**:
  - `ko2kegg_abundance()`: Now uses internal database instead of KEGG API
  - `pathway_gsea()`: Updated to use long-format data for gene set enrichment
  - Both functions maintain backward compatibility
  - Added comprehensive input validation and error handling

* **PICRUSt 2.6.2 compatibility**:
  - Automatic detection and cleaning of "ko:" prefix in KO IDs
  - Handles both old (K##### format) and new (ko:K##### format) PICRUSt2 outputs
  - Seamless migration path for existing users

* **Data quality improvements**:
  - Validates KO ID format (K##### pattern)
  - Detects negative abundance values
  - Reports missing values with detailed statistics
  - Identifies all-zero KOs across samples
  - Checks for duplicate column names

* **Performance**:
  - Processing speed: 2,145-6,452 KOs/second
  - Handles large datasets (1000+ KOs, 50+ samples) efficiently
  - Progress bar for long-running operations

* **Testing**:
  - Added comprehensive test suite with 64 tests
  - Covers data structure, functionality, performance, and edge cases
  - All tests passing with 100% coverage of new features

This resolves Discussion #113 and provides a robust, API-independent solution for KEGG pathway analysis.

# ggpicrust2 2.5.4

## Improvements

### Enhanced Error Messages for pathway_annotation() (#142)

* **Significantly improved diagnostic messages when no significant pathways are found**:
  - Provides clear explanation when p_adjust < 0.05 filter removes all pathways
  - Shows detailed statistics (total features, significant count, minimum p-value)
  - Lists possible reasons (sample size, effect size, variability, method choice)
  - Offers concrete recommendations (data quality checks, alternative methods)
  - Suggests using ko_to_kegg = FALSE for local annotation without KEGG API

* **Enhanced error messages for KEGG API failures**:
  - Detailed diagnostic information when all queries fail
  - Distinguishes between "not found (HTTP 404)" and "network errors"
  - Lists affected KO IDs for troubleshooting
  - Provides specific recommendations based on error type
  - Suggests local annotation as reliable alternative

* **Updated documentation**:
  - Clarified p_adjust < 0.05 filtering behavior in function description
  - Added note about NA columns when no significant pathways exist
  - Improved parameter descriptions for ko_to_kegg parameter
  - Enhanced return value documentation with filtering behavior

* **User experience improvements**:
  - All messages follow R package standards (warning(), stop() instead of cat())
  - No emoji characters (professional text output)
  - Clear formatting for readability
  - Actionable information to help users resolve issues quickly

This addresses Discussion #142 where users received NA annotations without understanding why.

## Bug Fixes

### pathway_gsea() MetaCyc Support Fix (#174)

* **Fixed undefined organism parameter bug**:
  - Added `organism = "ko"` parameter with default value
  - Prevents "object 'organism' not found" error
  - Properly passes organism to prepare_gene_sets() function

* **Added input validation for MetaCyc data**:
  - Detects when MetaCyc pathway IDs are provided instead of EC numbers
  - Issues clear warning directing users to use pathway_daa() for pathway-level data
  - Helps users understand the difference between gene-level and pathway-level analysis
  
* **Enhanced documentation**:
  - Clarified that GSEA requires gene-level data (EC numbers for MetaCyc)
  - Added cross-reference to pathway_daa() for pathway abundance analysis
  - Improved parameter descriptions to prevent data type confusion

* **Scientific integrity maintained**:
  - Ensures correct analysis method for each data type
  - Prevents misleading results from incorrect data usage
  - Guides users to appropriate functions based on their data

### Critical annotation_custom() Fix (#184)

* **Fixed annotation_custom() parameter type issue in pathway_errorbar()**:
  - Removed unit object wrappers from `annotation_custom()` position parameters
  - Now uses numeric values directly for xmin, xmax, ymin, ymax parameters
  - Resolves "no applicable method for 'rescale' applied to an object of class 'c('simpleUnit', 'unit', 'unit_v2')'" error
  - Fixes pathway class background color rendering failures

* **Root cause identified and resolved**:
  - `annotation_custom()` expects numeric values, not unit objects for position parameters
  - Previous fix attempt (changing ggplot2::unit to grid::unit) was incorrect
  - Both ggplot2::unit() and grid::unit() return identical objects - the issue was using unit objects at all
  - Theme-related unit usage (legend.key.size, plot.margin) remains unchanged and correct

This fix resolves the critical rendering issue where users encountered errors when generating
pathway error bar plots with pathway class backgrounds, particularly with R 4.4+ and ggplot2 4.0.0.

# ggpicrust2 2.5.3

## Previous Release
* Initial attempt to fix issue #184 (superseded by v2.5.4)

# ggpicrust2 2.5.2

## Bug Fixes

### Critical Edge Case Resolution

* **Fixed pathway_errorbar empty data handling**:
  - `pathway_errorbar()` now returns NULL instead of crashing when no annotation data is available
  - Main `ggpicrust2()` function gracefully handles NULL plot objects
  - Users still receive complete results data even when plots cannot be generated
  - Improved warning messages for better user experience

* **Enhanced robustness for datasets with no significant pathways**:
  - Prevents crashes when all pathways have p_adjust > 0.05
  - Returns meaningful results with empty annotation columns for visualization
  - Maintains data integrity throughout the analysis pipeline
  - Works seamlessly with all PICRUSt2 versions including 2.6.2

* **Function signature consistency**:
  - Fixed `pathway_annotation()` function calls in main function
  - Ensures proper parameter passing throughout the workflow
  - Resolves compatibility issues between internal functions

These fixes ensure the package works reliably with all types of microbiome data, 
including edge cases where no statistically significant pathways are found.

# ggpicrust2 2.5.1

## Bug Fixes

### PICRUSt 2.6.2 Compatibility (#174) - Complete Resolution

* **Added automatic KO ID format detection and conversion**:
  - Automatically detects PICRUSt 2.6.2 format with "ko:" prefixes
  - Transparently removes "ko:" prefixes during data loading
  - Maintains full backward compatibility with PICRUSt 2.5.2 format
  - Eliminates "subscript out of bounds" errors caused by format mismatches

* **Enhanced core functions for seamless compatibility**:
  - Updated `ko2kegg_abundance()` with automatic format detection
  - Updated `ggpicrust2()` main function with compatibility layer
  - Added clear informational messages about format conversion
  - Zero manual preprocessing required for users

* **Comprehensive testing and validation**:
  - Tested with real PICRUSt 2.6.2 output files (5,000+ KO features)
  - Verified compatibility with all major DAA methods
  - Confirmed backward compatibility with existing workflows
  - Performance optimized for large datasets

# ggpicrust2 2.5.0

## Bug Fixes

### PICRUSt 2.6.2 Compatibility (#174)

* **Fixed compatibility issues with PICRUSt 2.6.2 output**:
  - Added comprehensive data validation for PICRUSt compatibility
  - Improved zero-abundance data filtering with fallback strategies
  - Enhanced error handling with specific PICRUSt version guidance
  - Added graceful handling of sparse data scenarios

* **Enhanced ALDEx2 error detection**:
  - Better error messages for insufficient or invalid data
  - Improved handling of edge cases (single features, all-zero data)
  - Added PICRUSt version-specific troubleshooting guidance

* **Improved LinDA analysis robustness**:
  - Better filtering of zero-abundance features before analysis
  - Enhanced data validation and error reporting
  - Added warnings for very sparse data scenarios

* **Enhanced ko2kegg_abundance function**:
  - Added fallback strategies for zero-abundance KO data
  - Improved compatibility warnings and error messages
  - Better handling of PICRUSt format variations

* **Added comprehensive compatibility guide**:
  - Created detailed troubleshooting documentation
  - Provided alternative analysis strategies
  - Added version-specific recommendations

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
