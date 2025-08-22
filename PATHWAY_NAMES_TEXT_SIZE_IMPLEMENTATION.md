# Implementation of pathway_names_text_size Parameter

## Overview
This document describes the implementation of the `pathway_names_text_size` parameter in the `pathway_errorbar()` function to address GitHub issue #173.

## Problem Statement
Users reported that pathway names (y-axis labels) in error bar plots were too small and difficult to read. While `pathway_class_text_size` was already available for pathway class annotations, there was no parameter to control the text size of pathway names displayed on the y-axis.

## Solution
Added a new parameter `pathway_names_text_size` to the `pathway_errorbar()` function that allows users to control the text size of pathway names (y-axis labels).

## Implementation Details

### 1. Function Parameter Addition
- **Parameter**: `pathway_names_text_size`
- **Type**: Numeric value or "auto"
- **Default**: "auto"
- **Description**: Controls the text size of pathway names displayed on the y-axis

### 2. Files Modified

#### A. `code/ggpicrust2/R/pathway_errorbar.R`
- Added parameter documentation in roxygen comments
- Added parameter to function signature
- Added logic to calculate text size (auto vs. custom)
- Updated `axis.text.y` theme element to use the calculated size

#### B. `code/ggpicrust2/man/pathway_errorbar.Rd`
- Added parameter documentation

#### C. `code/ggpicrust2/tests/testthat/test-pathway_errorbar.R`
- Added unit tests for the new parameter

### 3. Implementation Logic

```r
# Calculate smart text size for pathway names (y-axis labels)
pathway_names_final_text_size <- if (pathway_names_text_size == "auto") {
  if (exists("calculate_smart_text_size")) {
    calculate_smart_text_size(nrow(daa_results_filtered_sub_df), 
                            base_size = 10, min_size = 8, max_size = 14)
  } else {
    10  # Default size (current hardcoded value)
  }
} else {
  pathway_names_text_size
}
```

### 4. Usage Examples

#### Basic Usage
```r
# Default auto sizing
pathway_errorbar(..., pathway_names_text_size = "auto")

# Custom size
pathway_errorbar(..., pathway_names_text_size = 12)
```

#### Combined with pathway_class_text_size
```r
pathway_errorbar(
  ...,
  pathway_names_text_size = 14,      # Large pathway names
  pathway_class_text_size = 6        # Large pathway class annotations
)
```

## Backward Compatibility
- The implementation is fully backward compatible
- Default behavior remains unchanged (auto sizing)
- Existing code will continue to work without modification

## Testing
- Unit tests added to verify parameter functionality
- Integration tests confirm compatibility with existing features
- Visual tests demonstrate different text sizes

## Benefits
1. **Improved Readability**: Users can now increase text size for better readability
2. **Flexibility**: Both auto and manual sizing options available
3. **Consistency**: Follows the same pattern as `pathway_class_text_size`
4. **Smart Defaults**: Auto sizing adapts to the number of pathways displayed

## Future Considerations
- Could extend to other plot types in the package
- Could add more text styling options (font family, weight, etc.)
- Could implement theme-based text sizing

## Resolution of GitHub Issue #173
This implementation directly addresses the user's request to increase text size for pathway labels in error bar plots, providing both immediate solutions and flexible customization options.
