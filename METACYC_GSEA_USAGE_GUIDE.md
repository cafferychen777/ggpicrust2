# MetaCyc GSEA Implementation - Usage Guide

## Overview

The ggpicrust2 package now provides full support for MetaCyc pathway gene set enrichment analysis (GSEA). This implementation allows you to analyze EC enzyme abundance data against MetaCyc metabolic pathways.

## Key Features

✅ **Complete MetaCyc Integration**: 50+ MetaCyc pathways with EC number mappings  
✅ **Full GSEA Workflow**: Compatible with existing ggpicrust2 GSEA functions  
✅ **Pathway Annotations**: Automatic pathway name/description lookup  
✅ **Validation & Error Handling**: Robust input validation and error messages  
✅ **Multiple GSEA Methods**: Support for fgsea and clusterProfiler  

## Installation & Requirements

```r
# Required packages
library(ggpicrust2)
library(fgsea)  # For GSEA analysis

# Verify MetaCyc support
gene_sets <- prepare_gene_sets(pathway_type = "MetaCyc")
length(gene_sets)  # Should show ~50 pathways
```

## Data Format Requirements

### EC Abundance Data
Your abundance data must have:
- **Row names**: EC numbers in format `"EC:X.X.X.X"` (e.g., "EC:1.1.1.1")
- **Columns**: Sample identifiers
- **Values**: Numerical abundance/expression values

```r
# Example EC abundance data structure
head(rownames(ec_abundance))
# [1] "EC:1.1.1.1"  "EC:1.2.1.12" "EC:2.7.1.1"  "EC:2.7.1.11"
# [2] "EC:4.1.2.13" "EC:1.2.1.59"
```

### Metadata
Standard ggpicrust2 metadata format:
```r
metadata <- data.frame(
  sample_id = colnames(ec_abundance),
  condition = rep(c("Control", "Treatment"), each = 5),
  stringsAsFactors = FALSE
)
rownames(metadata) <- metadata$sample_id
```

## Basic Usage

### 1. Prepare Gene Sets
```r
# Load MetaCyc pathways
metacyc_gene_sets <- prepare_gene_sets(pathway_type = "MetaCyc")

# Inspect available pathways
cat("Available MetaCyc pathways:", length(metacyc_gene_sets), "\n")
head(names(metacyc_gene_sets))
```

### 2. Run GSEA Analysis
```r
# Run MetaCyc GSEA
gsea_results <- pathway_gsea(
  abundance = ec_abundance,
  metadata = metadata,
  group = "condition",
  pathway_type = "MetaCyc",
  method = "fgsea",
  nperm = 1000,
  min_size = 5,
  max_size = 100
)
```

### 3. Annotate Results
```r
# Add pathway descriptions
annotated_results <- gsea_pathway_annotation(
  gsea_results = gsea_results,
  pathway_type = "MetaCyc"
)

# View top results
head(annotated_results[, c("pathway_id", "pathway_name", "NES", "pvalue")])
```

## Available MetaCyc Pathways

The implementation includes major metabolic pathways:

### Central Carbon Metabolism
- **1CMET2-PWY**: N10-formyl-tetrahydrofolate biosynthesis
- **ANAGLYCOLYSIS-PWY**: Glycolysis III (from glucose)
- **CALVIN-PWY**: Calvin-Benson cycle
- **TCA**: TCA cycle (incomplete pathway data)

### Amino Acid Metabolism
- **ARGSYN-PWY**: Arginine biosynthesis I (via L-ornithine)
- **ARG+POLYAMINE-SYN**: Superpathway of arginine and polyamine biosynthesis
- **BRANCHED-CHAIN-AA-SYN-PWY**: Branched-chain amino acid biosynthesis

### Cofactor Biosynthesis
- **BIOTIN-BIOSYNTHESIS-PWY**: Biotin biosynthesis I
- **COA-PWY**: CoA biosynthesis I
- **COBALSYN-PWY**: Adenosylcobalamin biosynthesis

### And Many More...
Total: 50+ pathways with 5-20 EC numbers each

## Advanced Usage

### Custom Parameters
```r
# Fine-tune GSEA parameters
gsea_results <- pathway_gsea(
  abundance = ec_abundance,
  metadata = metadata,
  group = "condition",
  pathway_type = "MetaCyc",
  method = "fgsea",
  rank_method = "signal2noise",    # or "t_test", "log2_ratio", "diff_abundance"
  nperm = 1000,
  min_size = 3,                    # Minimum pathway size
  max_size = 200,                  # Maximum pathway size
  p.adjust = "BH",                 # Multiple testing correction
  seed = 42                        # For reproducibility
)
```

### Using clusterProfiler
```r
# Alternative GSEA method
gsea_results <- pathway_gsea(
  abundance = ec_abundance,
  metadata = metadata,
  group = "condition",
  pathway_type = "MetaCyc",
  method = "clusterProfiler",
  min_size = 5,
  max_size = 100
)
```

### Quality Control
```r
# Check pathway coverage
gene_sets <- prepare_gene_sets(pathway_type = "MetaCyc")
pathway_sizes <- sapply(gene_sets, length)
summary(pathway_sizes)

# Check EC number overlap with your data
your_ecs <- rownames(ec_abundance)
pathway_ecs <- unique(unlist(gene_sets))
overlap <- intersect(your_ecs, pathway_ecs)
cat("EC overlap:", length(overlap), "/", length(your_ecs), "ECs\n")
```

## Troubleshooting

### Common Issues

1. **"No pathways tested" error**
   - Check EC number format: must be "EC:X.X.X.X"
   - Ensure sufficient overlap between your ECs and pathway ECs

2. **"MetaCyc reference not found" error**
   - Reinstall/update ggpicrust2 package
   - Check that extdata files are properly installed

3. **Warning about pathway sizes**
   - Adjust min_size and max_size parameters
   - Some pathways may be too small/large for your analysis

### Data Validation
```r
# Validate your data before GSEA
cat("EC format check:\n")
ec_pattern <- "^EC:[0-9]+\\.[0-9]+\\.[0-9]+\\.[0-9]+$"
valid_ecs <- grepl(ec_pattern, rownames(ec_abundance))
cat("Valid ECs:", sum(valid_ecs), "/", length(valid_ecs), "\n")

# Check metadata alignment
cat("Sample alignment check:\n")
overlap_samples <- intersect(colnames(ec_abundance), rownames(metadata))
cat("Aligned samples:", length(overlap_samples), "\n")
```

## Comparison with KEGG

| Feature | MetaCyc | KEGG |
|---------|---------|------|
| **Input** | EC numbers | KO numbers |
| **Pathways** | 50+ curated | 100+ comprehensive |
| **Focus** | Enzymatic reactions | Gene functions |
| **Detail Level** | High biochemical detail | Broader functional categories |

## Example Complete Workflow

```r
library(ggpicrust2)
library(fgsea)

# 1. Load your data
# ec_abundance <- read.csv("your_ec_abundance.csv", row.names = 1)
# metadata <- read.csv("your_metadata.csv", row.names = 1)

# 2. Validate data format
cat("Validating data...\n")
stopifnot(all(grepl("^EC:", rownames(ec_abundance))))
stopifnot(all(colnames(ec_abundance) %in% rownames(metadata)))

# 3. Run MetaCyc GSEA
cat("Running MetaCyc GSEA...\n")
results <- pathway_gsea(
  abundance = ec_abundance,
  metadata = metadata,
  group = "condition",
  pathway_type = "MetaCyc",
  method = "fgsea"
)

# 4. Annotate and summarize
cat("Annotating results...\n")
annotated <- gsea_pathway_annotation(results, pathway_type = "MetaCyc")

# 5. Export results
write.csv(annotated, "metacyc_gsea_results.csv")

cat("Analysis complete!\n")
cat("Significant pathways (p < 0.05):", sum(annotated$pvalue < 0.05), "\n")
```

## Citation

When using MetaCyc GSEA functionality, please cite:

```
Yang, C. et al. ggpicrust2: an R package for PICRUSt2 predicted functional profile analysis and visualization. Bioinformatics 39, btad470 (2023).
```

And consider citing the MetaCyc database:
```
Caspi, R. et al. The MetaCyc database of metabolic pathways and enzymes - a 2019 update. Nucleic Acids Res. 48, D445-D453 (2020).
```