\name{exportCoMorTrax}
\docType{methods}
\alias{exportCoMorTrax}
\alias{exportCoMorTrax,CoMorTraxObject-method}
\title{Export CoMorTrax Results to Disk}
\description{
  Writes all analysis results stored in a \code{\link{CoMorTraxObject}}
  to a structured output directory. Results are exported as CSV files
  (for interoperability with downstream tools) and optionally as a
  compressed RDS file (for reloading into R).
}
\usage{
\S4method{exportCoMorTrax}{CoMorTraxObject}(object, outDir,
    format = c("csv","rds"), compress = TRUE, verbose = TRUE)
}
\arguments{
  \item{object}{A \code{\link{CoMorTraxObject}}.}
  \item{outDir}{Character. Path to the output directory. Will be created
    if it does not exist.}
  \item{format}{Character vector. Which formats to export. Any combination
    of \code{"csv"} and \code{"rds"}.}
  \item{compress}{Logical. Use \code{xz} compression for the RDS file.
    Default \code{TRUE}.}
  \item{verbose}{Logical. Default \code{TRUE}.}
}
\value{
  Invisibly returns a character vector of output file paths. Called for
  its side-effect of writing files.

  The output directory structure is:
  \preformatted{
  <outDir>/
  ├── interaction_effects/
  │   ├── interaction_table.csv
  │   ├── synergistic_genes.csv
  │   └── antagonistic_genes.csv
  ├── differential_expression/
  │   └── <comparison>_<cellType>.csv
  ├── pathway_enrichment/
  │   └── <geneSet>_<database>.csv
  ├── cell_vulnerability/
  │   └── vulnerability_scores.csv
  ├── networks/
  │   └── network_rewiring.csv
  ├── classifier/
  │   └── classifier_performance.csv
  └── CoMorTrax_object.rds   (if format includes "rds")
  }
}
\author{Ahmed Mokhtar Ramzy Salem}
\seealso{\code{\link{summaryCoMorTrax}}, \code{\link{validateCoMorTrax}}}
