\name{ec_reference}
\alias{ec_reference}
\docType{data}
\title{
EC Number Reference Dataset
}
\description{
A reference dataset mapping Enzyme Commission (EC) numbers to their descriptions.
Used internally by ggpicrust2 for annotating enzyme-level functional predictions.
}
\usage{data("ec_reference")}
\format{
  A data frame with 8405 observations and the following columns:
  \describe{
    \item{\code{id}}{Character. EC number in the format "EC:X.X.X.X"}
    \item{\code{description}}{Character. Human-readable enzyme name/description}
  }
}
\source{
KEGG REST API (\url{https://rest.kegg.jp})
}
\examples{
data("ec_reference")
head(ec_reference)
}
\keyword{datasets}
