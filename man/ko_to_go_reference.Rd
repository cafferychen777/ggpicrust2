\name{ko_to_go_reference}
\alias{ko_to_go_reference}
\docType{data}
\title{
KO to GO Reference Mapping Dataset
}
\description{
A comprehensive reference dataset that maps KEGG Orthology (KO) identifiers to Gene Ontology (GO) terms. This dataset enables GO pathway analysis in ggpicrust2 by providing the necessary mappings between functional predictions and GO biological processes, molecular functions, and cellular components.
}
\usage{data("ko_to_go_reference")}
\format{
  A data frame with the following columns:
  \describe{
    \item{\code{go_id}}{Character. GO term identifier in the format "GO:XXXXXXX"}
    \item{\code{go_name}}{Character. Human-readable name of the GO term}
    \item{\code{category}}{Character. GO category: "BP" (Biological Process), "MF" (Molecular Function), or "CC" (Cellular Component)}
    \item{\code{ko_members}}{Character. Semicolon-separated list of KO identifiers associated with this GO term}
  }
}
\details{
This dataset contains mappings between KO identifiers and GO terms, enabling comprehensive GO pathway analysis in microbiome functional profiling studies. The dataset includes:

\itemize{
  \item \strong{Biological Processes (BP):} Metabolic pathways, stress responses, cell wall processes, virulence factors, and environmental adaptations
  \item \strong{Molecular Functions (MF):} Enzymatic activities, binding functions, transport activities, and regulatory functions
  \item \strong{Cellular Components (CC):} Membrane structures, protein complexes, cell wall components, and cytoskeletal elements
}

The dataset is created by combining information from multiple sources:
\itemize{
  \item KEGG database API for KO to GO mappings
  \item Enhanced basic mapping covering major microbiome-relevant functions
  \item Manual curation for microbiome-specific processes
}

Each GO term includes a list of associated KO identifiers that are known to be involved in that biological process, molecular function, or cellular component. This enables pathway-level analysis of PICRUSt2 functional predictions.
}
\source{
The dataset is compiled from:
\itemize{
  \item KEGG database (https://www.kegg.jp/)
  \item Gene Ontology database (http://geneontology.org/)
  \item Enhanced basic mapping for microbiome analysis
}
}
\references{
\itemize{
  \item Kanehisa, M., & Goto, S. (2000). KEGG: kyoto encyclopedia of genes and genomes. Nucleic acids research, 28(1), 27-30.
  \item Ashburner, M., et al. (2000). Gene ontology: tool for the unification of biology. Nature genetics, 25(1), 25-29.
  \item Chen Yang, et al. (2023). ggpicrust2: an R package for PICRUSt2 predicted functional profile analysis and visualization. Bioinformatics, 39(8), btad470.
}
}
\examples{
# Load the dataset
data("ko_to_go_reference")

# Explore the dataset structure
head(ko_to_go_reference)
str(ko_to_go_reference)

# Check the distribution of GO categories
table(ko_to_go_reference$category)

# Find GO terms related to metabolism
metabolism_terms <- ko_to_go_reference[grepl("metabolic", ko_to_go_reference$go_name, ignore.case = TRUE), ]
head(metabolism_terms)

# Get KO members for a specific GO term
glycolysis_term <- ko_to_go_reference[ko_to_go_reference$go_id == "GO:0006096", ]
if (nrow(glycolysis_term) > 0) {
  ko_list <- strsplit(glycolysis_term$ko_members, ";")[[1]]
  cat("KO identifiers in glycolytic process:", paste(ko_list, collapse = ", "))
}

# Use in pathway analysis
\dontrun{
library(ggpicrust2)
library(tibble)

# Load example data
data("ko_abundance")
data("metadata")

# Perform GO pathway GSEA analysis
gsea_results <- pathway_gsea(
  abundance = ko_abundance \%>\% column_to_rownames("#NAME"),
  metadata = metadata \%>\% column_to_rownames("sample_name"),
  group = "Environment",
  method = "fgsea",
  pathway_type = "GO",
  go_category = "BP",
  rank_method = "signal2noise"
)

# View results
head(gsea_results)
}
}
\seealso{
\code{\link{pathway_gsea}}, \code{\link{ko_abundance}}, \code{\link{metadata}}
}
\keyword{datasets}
\keyword{microbiome}
\keyword{functional-analysis}
\keyword{gene-ontology}
