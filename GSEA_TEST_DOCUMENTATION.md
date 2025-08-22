# Comprehensive GSEA Mathematical Test Suite Documentation

## Overview

This document describes the comprehensive test suite created for validating the mathematical correctness and statistical accuracy of GSEA calculation functions in `pathway_gsea.R`.

## Test Files Created

### 1. `test-pathway_gsea_comprehensive_mathematical.R`
**Focus**: Mathematical correctness and algorithmic implementation

**Key Test Categories**:
- **Ranking Metric Validation**: Tests mathematical formulas for all ranking methods
  - Signal-to-noise ratio: `(mean1-mean2)/(sd1+sd2)`
  - T-test statistic computation
  - Log2 fold change calculations: `log2(mean1/mean2)`
  - Difference in abundance: `mean1-mean2`

- **Statistical Properties**: Validates statistical behavior
  - Features with known differences rank appropriately
  - Ranking consistency across methods
  - Statistical power analysis with different sample sizes

- **Edge Case Handling**: Tests numerical stability
  - Zero variance handling (adds small constant 0.00001)
  - Extremely small/large values
  - Scale invariance properties

- **Parameter Validation**: Comprehensive input validation
  - Invalid data types and formats
  - Missing required parameters
  - Boundary conditions

### 2. `test-pathway_gsea_integration_validation.R` 
**Focus**: Method comparison and result standardization

**Key Test Categories**:
- **Method Comparison**: fgsea vs clusterProfiler integration
  - Result format consistency across methods
  - Column name mapping validation
  - Statistical value comparisons

- **Result Standardization**: Ensures consistent output structure
  - Standard column names: `pathway_id`, `pathway_name`, `size`, `ES`, `NES`, `pvalue`, `p.adjust`, `leading_edge`, `method`
  - Data type consistency across methods
  - Leading edge gene extraction accuracy

- **Parameter Handling**: Tests parameter consistency
  - `nperm`, `min_size`, `max_size`, `p.adjust` method handling
  - Large gene set collections performance
  - Empty gene set robustness

### 3. `test-pathway_gsea_statistical_validation.R`
**Focus**: Statistical accuracy and boundary conditions

**Key Test Categories**:
- **Statistical Validity**: Tests statistical properties
  - Ranking metrics produce valid statistical orderings
  - Known enrichment patterns are detected correctly
  - P-value uniformity under null hypothesis
  - Multiple testing correction accuracy

- **Mathematical Edge Cases**: Boundary condition testing
  - Perfect group separation
  - Identical means with different variances
  - Extreme variance ratios
  - Constant data handling

- **Power Analysis**: Statistical power validation
  - Power increases with sample size
  - Power increases with effect size
  - Reasonable power for large effects

- **Enrichment Score Validation**: Tests ES calculation properties
  - Correct signs for enriched/depleted pathways
  - NES normalization properties
  - Leading edge consistency

### 4. `test-pathway_gsea_performance_stress.R`
**Focus**: Performance, scalability, and robustness

**Key Test Categories**:
- **Scalability Testing**: Performance with large datasets
  - Execution time scaling with data size
  - Memory usage efficiency
  - Large gene set collection handling

- **Stress Testing**: Robustness under adverse conditions
  - High sparsity data (90% zeros)
  - Problematic data patterns (identical values, outliers)
  - Malformed input handling

- **Consistency Testing**: Reproducibility validation
  - Repeated calculations produce identical results
  - Memory cleanup after large calculations
  - Thread safety considerations

## Mathematical Validation Approach

### 1. Deterministic Test Data Generation
```r
create_deterministic_test_data <- function(n_features, n_samples, effect_size) {
  # Creates abundance matrix with known differential patterns
  # First half of features have true differential abundance
  # Second half are null (no difference)
  # Uses log-normal distribution to mimic microbiome data
}
```

### 2. Known Enrichment Validation
- Creates pathways with known enriched/depleted gene sets
- Validates that enrichment scores have correct signs
- Tests that true differential features rank higher than null features

### 3. Statistical Property Verification
- Tests correlation between effect size and ranking statistics
- Validates power analysis results
- Checks p-value distribution under null hypothesis

## Key Mathematical Validations

### Signal-to-Noise Ratio
```r
# Formula: (mean1 - mean2) / (sd1 + sd2)
# Handles zero standard deviation by adding small constant
# Scale-invariant property validated
```

### T-Test Statistic
```r
# Uses built-in t.test() function
# Validates against manual calculations
# Handles constant data appropriately
```

### Log2 Fold Change
```r
# Formula: log2(mean1 / mean2)
# Adds small constant to avoid division by zero
# Scale-invariant property validated
```

### Difference in Abundance
```r
# Formula: mean1 - mean2
# Scales proportionally with data magnitude
# Simple but effective for abundance differences
```

## Statistical Validation Framework

### 1. Known Effect Size Testing
- Creates data with predetermined fold changes (2x, 4x, etc.)
- Validates ranking metrics detect these differences
- Tests statistical significance calculations

### 2. Null Distribution Testing
- Creates data with no true differences
- Validates p-values are approximately uniform
- Tests false positive rate control

### 3. Power Analysis
- Tests with varying sample sizes (6, 10, 20, 40)
- Tests with varying effect sizes (0.5, 1.0, 2.0)
- Validates power increases appropriately

## Performance Benchmarks

### Scalability Targets
- **Small datasets**: 50 features × 20 samples → < 1 second
- **Medium datasets**: 500 features × 50 samples → < 5 seconds  
- **Large datasets**: 2000 features × 100 samples → < 10 seconds

### Memory Efficiency
- Memory increase should not exceed 5× data size
- Proper cleanup after large calculations
- No memory leaks in repeated calculations

## Integration with Existing Tests

### Complementary to Existing Tests
The comprehensive test suite complements the existing `test-pathway_gsea.R` by:
- Adding mathematical validation (missing in original)
- Testing statistical properties extensively
- Adding performance and stress testing
- Providing edge case coverage

### Mock Strategy
Uses sophisticated mocking to:
- Test internal functions in isolation
- Validate result format standardization
- Enable testing without external dependencies
- Provide deterministic test results

## Usage Instructions

### Running Individual Test Suites
```r
# Mathematical correctness tests
testthat::test_file("tests/testthat/test-pathway_gsea_comprehensive_mathematical.R")

# Integration and method comparison tests  
testthat::test_file("tests/testthat/test-pathway_gsea_integration_validation.R")

# Statistical validation tests
testthat::test_file("tests/testthat/test-pathway_gsea_statistical_validation.R")

# Performance and stress tests
testthat::test_file("tests/testthat/test-pathway_gsea_performance_stress.R")
```

### Running All GSEA Tests
```r
devtools::test(filter = "pathway_gsea")
```

## Expected Outcomes

### Mathematical Accuracy
- All ranking metrics produce mathematically correct results
- Statistical properties are validated
- Edge cases are handled gracefully

### Method Consistency
- fgsea and clusterProfiler produce comparable results
- Result formats are standardized across methods
- Parameter handling is consistent

### Performance
- Scales reasonably with data size (sub-quadratic)
- Memory usage is efficient
- Robust to problematic data patterns

### Statistical Validity
- Proper p-value distributions
- Correct multiple testing adjustments
- Appropriate statistical power characteristics

## Notes for Future Maintenance

### Test Data Generation
- Use fixed seeds for reproducibility
- Create realistic microbiome-like distributions
- Include both positive and negative controls

### Mock Strategy
- Mock external dependencies (fgsea, clusterProfiler)
- Provide deterministic results for validation
- Test internal function behavior independently

### Performance Monitoring
- Monitor execution times as codebase evolves
- Watch for memory leaks in repeated calculations
- Validate scaling behavior with larger datasets

This comprehensive test suite ensures the mathematical correctness, statistical validity, and robust performance of the GSEA implementation in ggpicrust2.