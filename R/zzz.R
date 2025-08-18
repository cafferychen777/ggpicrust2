#' @importFrom grDevices col2rgb
#' @importFrom methods new
#' @importFrom stats as.formula model.matrix relevel p.adjust reorder var prcomp sd
#' @importFrom utils head data

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

.onAttach <- function(libname, pkgname) {
  packageStartupMessage("Loading required package: ggpicrust2")
  packageStartupMessage("To cite ggpicrust2 in publications use:")
  packageStartupMessage("Chen Yang, Jiahao Mai, Xuan Cao, Aaron Burberry, Fabio Cominelli, Liangliang Zhang, ggpicrust2: an R package for PICRUSt2 predicted functional profile analysis and visualization, Bioinformatics, Volume 39, Issue 8, August 2023, btad470, https://doi.org/10.1093/bioinformatics/btad470")
} 