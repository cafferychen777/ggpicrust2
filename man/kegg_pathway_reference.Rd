\name{kegg_pathway_reference}
\alias{kegg_pathway_reference}
\docType{data}
\title{
KEGG Pathway Name Reference Dataset
}
\description{
A reference dataset mapping KEGG pathway IDs to their human-readable names.
Used internally by ggpicrust2 for pathway annotation in DAA and GSEA results.
}
\usage{data("kegg_pathway_reference")}
\format{
  A data frame with 505 observations and the following columns:
  \describe{
    \item{\code{pathway}}{Character. KEGG pathway ID in the format "koXXXXX" (e.g., "ko00010")}
    \item{\code{pathway_name}}{Character. Human-readable pathway name (e.g., "Glycolysis / Gluconeogenesis")}
  }
}
\source{
KEGG REST API (\url{https://rest.kegg.jp})
}
\examples{
data("kegg_pathway_reference")
head(kegg_pathway_reference)
}
\keyword{datasets}
