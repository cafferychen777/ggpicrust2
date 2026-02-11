\name{metacyc_reference}
\alias{metacyc_reference}
\docType{data}
\title{
MetaCyc Pathway Reference Dataset
}
\description{
A reference dataset mapping MetaCyc pathway identifiers to their descriptions.
Used internally by ggpicrust2 for annotating MetaCyc pathway analysis results.
}
\usage{data("metacyc_reference")}
\format{
  A data frame with 2714 observations and the following columns:
  \describe{
    \item{\code{id}}{Character. MetaCyc pathway identifier (e.g., "GLYCOLYSIS", "TCA")}
    \item{\code{description}}{Character. Human-readable pathway description}
  }
}
\source{
MetaCyc database (\url{https://metacyc.org})
}
\examples{
data("metacyc_reference")
head(metacyc_reference)
}
\keyword{datasets}
