#' Skip tests that are too slow for regular testing
#' @param threshold Time threshold in seconds (default: 1)
skip_if_testing_is_slow <- function(threshold = 1) {
  if (!identical(Sys.getenv("TEST_SLOW"), "true")) {
    testthat::skip("Skipping slow tests")
  }
} 