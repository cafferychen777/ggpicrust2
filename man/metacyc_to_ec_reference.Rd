\name{metacyc_to_ec_reference}
\alias{metacyc_to_ec_reference}
\docType{data}
\title{
MetaCyc Pathway to EC Number Mapping Dataset
}
\description{
A reference dataset mapping MetaCyc pathway identifiers to their associated
Enzyme Commission (EC) numbers. Used internally by ggpicrust2 for MetaCyc
pathway analysis, enabling the mapping between EC-level functional predictions
and MetaCyc pathways.
}
\usage{data("metacyc_to_ec_reference")}
\format{
  A data frame with 575 observations and the following columns:
  \describe{
    \item{\code{pathway}}{Character. MetaCyc pathway identifier (e.g., "1CMET2-PWY")}
    \item{\code{ec_numbers}}{Character. Semicolon-separated list of EC numbers associated with the pathway}
  }
}
\source{
MetaCyc database (\url{https://metacyc.org})
}
\examples{
data("metacyc_to_ec_reference")
head(metacyc_to_ec_reference)

# Count EC numbers per pathway
ec_counts <- sapply(strsplit(metacyc_to_ec_reference$ec_numbers, ";"), length)
summary(ec_counts)
}
\keyword{datasets}
