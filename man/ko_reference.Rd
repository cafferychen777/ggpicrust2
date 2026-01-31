\name{ko_reference}
\alias{ko_reference}
\docType{data}
\title{
KEGG Orthology (KO) Reference Dataset
}
\description{
A comprehensive reference dataset mapping KEGG Orthology (KO) identifiers to their
pathway classifications and descriptions. Each KO entry can appear in multiple rows
if it belongs to multiple pathways.
}
\usage{data("ko_reference")}
\format{
  A data frame with 58693 observations and the following columns:
  \describe{
    \item{\code{id}}{Character. KO identifier (e.g., "K00001")}
    \item{\code{PathwayL1}}{Character. Top-level KEGG pathway category (e.g., "Metabolism")}
    \item{\code{PathwayL2}}{Character. Second-level pathway category (e.g., "Carbohydrate metabolism")}
    \item{\code{Pathway}}{Character. Specific pathway name with ID (e.g., "Glycolysis / Gluconeogenesis [PATH:ko00010]")}
    \item{\code{description}}{Character. KO entry description with gene name and EC number}
  }
}
\source{
KEGG REST API (\url{https://rest.kegg.jp})
}
\examples{
data("ko_reference")
head(ko_reference)

# Check pathway hierarchy
table(ko_reference$PathwayL1)
}
\keyword{datasets}
