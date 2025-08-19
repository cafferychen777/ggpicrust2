# PICRUSt Version Compatibility Guide

## Issue Description

Some users have reported compatibility issues when using ggpicrust2 with PICRUSt 2.6.2 output files, while the same analysis works fine with PICRUSt 2.5.2 output. The typical error pattern includes:

```
The kegg pathway with zero abundance in all the different samples has been removed.
0 pathways filtered (prevalence < 1)
0 features are filtered!
The filtered data has 0 samples and 0 features that will be tested!
Fit linear models ...
Error in x[, ii] : subscript out of bounds
```

## Root Cause

The issue stems from differences in output format between PICRUSt versions, which can lead to:
1. Overly aggressive filtering of zero-abundance pathways
2. Data format inconsistencies that cause all features to be filtered out
3. Empty data matrices that cause downstream analysis methods (especially LinDA) to fail

## Solutions

### 1. Automatic Handling (Recommended)

The latest version of ggpicrust2 includes improved compatibility handling:

```r
# The package now automatically detects and handles PICRUSt version issues
library(ggpicrust2)

# Your analysis should work with both PICRUSt 2.5.2 and 2.6.2 outputs
results <- ggpicrust2(
  file = "pred_metagenome_unstrat.tsv",
  metadata = metadata,
  group = "your_group_column",
  pathway = "KO",
  daa_method = "ALDEx2",  # Try ALDEx2 first if LinDA fails
  ko_to_kegg = TRUE
)
```

### 2. Alternative DAA Methods

If you encounter issues with LinDA, try other methods:

```r
# ALDEx2 is generally more robust to data format issues
daa_results <- pathway_daa(
  abundance = abundance,
  metadata = metadata,
  group = "your_group",
  daa_method = "ALDEx2"  # Instead of "LinDA"
)

# Or try DESeq2
daa_results <- pathway_daa(
  abundance = abundance,
  metadata = metadata,
  group = "your_group",
  daa_method = "DESeq2"
)
```

### 3. Data Preprocessing

If issues persist, try preprocessing your data:

```r
# Load your PICRUSt output
abundance_data <- read.delim("pred_metagenome_unstrat.tsv", 
                            check.names = FALSE, 
                            row.names = 1)

# Remove features with zero abundance across all samples
abundance_filtered <- abundance_data[rowSums(abundance_data) > 0, ]

# Check data quality
cat("Original features:", nrow(abundance_data), "\n")
cat("Filtered features:", nrow(abundance_filtered), "\n")
cat("Proportion retained:", nrow(abundance_filtered)/nrow(abundance_data), "\n")

# Proceed with analysis
daa_results <- pathway_daa(
  abundance = abundance_filtered,
  metadata = metadata,
  group = "your_group",
  daa_method = "ALDEx2"
)
```

### 4. Fallback to PICRUSt 2.5.2

If compatibility issues persist, consider using PICRUSt 2.5.2 output:

1. Re-run your analysis with PICRUSt 2.5.2
2. Use the 2.5.2 output with ggpicrust2
3. Report the issue to the ggpicrust2 developers with your specific data

## Troubleshooting

### Error: "All abundance values are zero"
- Check your data loading: ensure the file is read correctly
- Verify column names match between abundance data and metadata
- Try different file reading parameters (e.g., `sep="\t"`, `check.names=FALSE`)

### Error: "subscript out of bounds"
- This typically occurs with LinDA method on filtered data
- Switch to ALDEx2 or DESeq2 methods
- Check that your metadata sample names match abundance data columns

### Warning: "X% of features have zero abundance"
- This is normal for PICRUSt data, but >90% may indicate format issues
- Verify your PICRUSt analysis parameters
- Consider using less stringent filtering

## Reporting Issues

If you continue to experience problems, please report them with:
1. PICRUSt version used
2. ggpicrust2 version
3. Sample of your data (first few rows/columns)
4. Complete error message
5. Your analysis code

## Version Information

This compatibility guide applies to:
- ggpicrust2 version 2.4.1+
- PICRUSt2 versions 2.5.2 and 2.6.2
- R version 4.0+
