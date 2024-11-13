.onLoad <- function(libname, pkgname) {
  op <- options()
  op.ggpicrust2 <- list(
    ggpicrust2.cache_enabled = TRUE,
    ggpicrust2.verbose = TRUE,
    ggpicrust2.max_retries = 3,
    ggpicrust2.timeout = 300
  )
  toset <- !(names(op.ggpicrust2) %in% names(op))
  if (any(toset)) options(op.ggpicrust2[toset])
  
  invisible()
} 