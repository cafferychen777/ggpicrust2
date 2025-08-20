#!/usr/bin/env Rscript
# Test the reference data loading issue

cat("=== Testing reference data loading scenarios ===\n")

# Test 1: Simulate the package environment issue
cat("\n1. Testing load_reference_data function behavior\n")

# Source the function
source('R/pathway_annotation.R')

# Test what happens when ggpicrust2 package namespace doesn't exist
tryCatch({
  result <- load_reference_data("MetaCyc")
  cat("SUCCESS: Reference data loaded\n")
  cat("Structure:", str(result))
}, error = function(e) {
  cat("ERROR loading reference data:\n")
  cat(e$message, "\n")
})

# Test 2: Check if file exists where expected
cat("\n2. Testing file paths\n")

ref_file <- "MetaCyc_reference.RData"
ref_data_name <- "MetaCyc_reference"

# Test inst/extdata path (which should work in our dev environment)
ref_path1 <- file.path("inst/extdata", ref_file)
cat(sprintf("Path 1 (%s): %s\n", ref_path1, ifelse(file.exists(ref_path1), "EXISTS", "NOT FOUND")))

# Test system.file path (which would work in installed package)
ref_path2 <- system.file("extdata", ref_file, package = "ggpicrust2", mustWork = FALSE)
cat(sprintf("Path 2 (system.file): %s\n", ifelse(ref_path2 != "" && file.exists(ref_path2), "EXISTS", "NOT FOUND")))

# Test inst/extdata with system.file (alternative path in code)
ref_path3 <- system.file("inst/extdata", ref_file, package = "ggpicrust2", mustWork = FALSE)
cat(sprintf("Path 3 (system.file inst/extdata): %s\n", ifelse(ref_path3 != "" && file.exists(ref_path3), "EXISTS", "NOT FOUND")))

# Test 3: Manual loading to verify the actual content
cat("\n3. Manual reference data loading\n")

if (file.exists(ref_path1)) {
  cat("Loading reference data manually...\n")
  load(ref_path1)
  if (exists("MetaCyc_reference")) {
    cat("SUCCESS: MetaCyc_reference loaded\n")
    cat(sprintf("Dimensions: %s\n", paste(dim(MetaCyc_reference), collapse=" x ")))
    cat(sprintf("Column names: %s\n", paste(colnames(MetaCyc_reference), collapse=", ")))
    
    # Test the annotation process with some known pathways
    test_pathways <- c("1CMET2-PWY", "ALL-CHORISMATE-PWY", "NONEXISTENT-PWY")
    
    # Apply the standardization logic
    ref_data <- MetaCyc_reference
    if (all(c("X1", "X2") %in% colnames(ref_data))) {
      colnames(ref_data) <- c("id", "description")
    }
    
    matches <- match(test_pathways, ref_data$id)
    cat("\nTest annotation results:\n")
    for(i in 1:length(test_pathways)) {
      desc <- if(is.na(matches[i])) "NOT FOUND" else ref_data$description[matches[i]]
      cat(sprintf("  %s -> %s\n", test_pathways[i], desc))
    }
    
  } else {
    cat("ERROR: MetaCyc_reference object not found in file\n")
  }
} else {
  cat("ERROR: Reference file not found\n")
}

cat("\n=== Test completed ===\n")