# Group Size Validation Implementation

## Summary

Added comprehensive group size validation to improve statistical reliability awareness in GSEA functions, following Linus's principle: "fail fast with clear reasons."

## Changes Made

### 1. New Function: `validate_group_sizes()`

**Location**: `/Users/apple/Research/ggpicrust2/code/ggpicrust2/R/pathway_gsea.R` (lines 1-46)

**Purpose**: Validates group sizes and balance for statistical reliability

**Functionality**:
- **Error for insufficient samples**: Stops execution if any group has <2 samples
- **Warning for small groups**: Warns if any group has <3 samples (power issues)
- **Warning for severe imbalance**: Warns if two-group ratio >3:1 (bias concern)
- **Warning for multiple groups**: Warns about pairwise comparison behavior

**Example Messages**:
```
ERROR: "Groups with <2 samples detected in 'treatment' (A : 1; B : 1). Statistical comparison impossible."
WARNING: "Small group sizes in 'treatment' (B : 2). Statistical power severely limited. Recommend n>=3 per group."
WARNING: "Severe group imbalance in 'treatment' (ratio 4.5:1). Results may be biased."
```

### 2. Integration Points

**pathway_gsea() function** (line 225):
```r
# Group size validation - critical for statistical reliability
validate_group_sizes(Group, group)
```

**calculate_rank_metric() function** (line 304):
```r
# Group size validation - critical for statistical reliability
validate_group_sizes(Group, group)
```

### 3. Sample Matching Fix

Improved sample matching logic in both functions:
```r
common_samples <- intersect(colnames(abundance), names(Group))
abundance <- abundance[, common_samples]
Group <- Group[common_samples]
```

## Testing Verification

Comprehensive testing confirms:
- ✅ Valid balanced groups work silently
- ✅ Small groups generate appropriate warnings
- ✅ Severely imbalanced groups generate warnings
- ✅ Insufficient groups stop execution with clear errors
- ✅ Multiple groups generate guidance warnings
- ✅ Integration with both main functions works correctly

## Statistical Guidance Provided

1. **Minimum Requirements**: At least 2 samples per group for any statistical comparison
2. **Recommended Size**: At least 3 samples per group for reliable statistical power
3. **Balance Considerations**: Avoid >3:1 ratios between groups to prevent bias
4. **Multi-group Awareness**: Alerts users about pairwise comparison behavior

## Implementation Philosophy

Following Linus's core principles:
- **Good Taste**: One unified validation function eliminates special cases
- **Simplicity**: Clear, direct messages with specific guidance
- **Practical**: Addresses real statistical reliability concerns
- **Non-breaking**: Users get warnings but can proceed if they understand the risks

## Files Modified

1. `/Users/apple/Research/ggpicrust2/code/ggpicrust2/R/pathway_gsea.R`
   - Added `validate_group_sizes()` function
   - Modified `pathway_gsea()` function
   - Modified `calculate_rank_metric()` function

## Usage Impact

- **Existing workflows continue to work**
- **Users get informed about statistical reliability issues**
- **Clear guidance provided for improving study design**
- **No performance impact (simple table() operations)**

This implementation ensures users are aware of statistical limitations while maintaining workflow continuity.