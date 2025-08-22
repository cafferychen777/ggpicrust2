#!/usr/bin/env Rscript
# Debug the function directly

library(ggpicrust2)

# Create test data
test_abundance <- matrix(rnorm(100), nrow = 10, ncol = 10)
rownames(test_abundance) <- paste0("Gene", 1:10)
colnames(test_abundance) <- paste0("Sample", 1:10)

test_metadata <- data.frame(
  Environment = rep(c("Group1", "Group2"), 5),
  row.names = colnames(test_abundance)
)

cat("Test abundance matrix:\n")
print(dim(test_abundance))
print(head(rownames(test_abundance)))

cat("\nTest metadata:\n")
print(dim(test_metadata))
print(head(rownames(test_metadata)))

# Test the function directly
cat("\nTesting calculate_rank_metric function...\n")
result <- ggpicrust2:::calculate_rank_metric(
  abundance = test_abundance,
  metadata = test_metadata,
  group = "Environment",
  method = "signal2noise"
)

cat("Result class:", class(result), "\n")
cat("Result length:", length(result), "\n")
cat("Has names:", !is.null(names(result)), "\n")
cat("Names:", head(names(result)), "\n")
cat("Values:", head(result), "\n")

# Test if names exist and are correct
if (!is.null(names(result))) {
  cat("Names match rownames:", identical(names(result), rownames(test_abundance)), "\n")
} else {
  cat("ERROR: Names are NULL!\n")
}