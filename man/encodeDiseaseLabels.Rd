\name{encodeDiseaseLabels}
\docType{methods}
\alias{encodeDiseaseLabels}
\alias{encodeDiseaseLabels,CoMorTraxObject-method}
\alias{listComorbidityComparisons}
\title{Encode Multi-Label Disease States}
\description{
  Constructs the binary disease matrix and disease-group membership lists
  that drive all downstream comorbidity analyses. Supports three input
  modes: (1) multiple binary indicator columns, (2) a single categorical
  column with multi-label values separated by a delimiter, or (3) an
  explicit per-cell disease vector supplied by the caller.
}
\usage{
\S4method{encodeDiseaseLabels}{CoMorTraxObject}(object, diseaseCol,
    separator = ";", controlLabel = "control",
    minCellsPerGroup = 5, verbose = TRUE)

listComorbidityComparisons(object)
}
\arguments{
  \item{object}{A \code{\link{CoMorTraxObject}} (must contain a
    \code{counts} assay).}
  \item{diseaseCol}{Character vector. Either (a) the names of multiple
    binary (0/1) columns in \code{colData}, one per disease, or (b) the
    name of a single categorical column whose values are multi-label
    strings such as \code{"AD;T2D"}.}
  \item{separator}{Character. Delimiter used to split multi-label values
    when \code{diseaseCol} is a single column. Default \code{";"}.}
  \item{controlLabel}{Character. The value (or absence of all disease flags)
    that identifies control cells. Default \code{"control"}.}
  \item{minCellsPerGroup}{Integer. Minimum number of cells required for a
    disease group to be retained. Default 5.}
  \item{verbose}{Logical. Print progress and group-size messages.
    Default \code{TRUE}.}
}
\value{
  The input \code{\link{CoMorTraxObject}} updated with:
  \itemize{
    \item \code{diseaseMatrix(object)} — a \code{DataFrame} (cells × diseases)
      of binary indicators.
    \item \code{diseaseGroups(object)} — a named \code{list} mapping each
      disease state (including all comorbid combinations) to the integer
      indices of member cells.
    \item An additional \code{colData} column \code{disease_state} encoding
      each cell's combined disease label.
  }

  \code{listComorbidityComparisons()} returns a \code{character} vector of
  all valid comparison strings (e.g. \code{"AD_T2D_vs_Control"}) that can
  be passed to \code{\link{runComorbidityDE}}.
}
\author{Ahmed Mokhtar Ramzy Salem}
\seealso{
  \code{\link{runComorbidityDE}},
  \code{\link{quantifyInteractionEffects}}
}
\examples{
library(CoMorTrax)
set.seed(1)
sce <- simulateCoMorTrax(n_cells = 200, n_genes = 100,
                          diseases = c("AD","T2D"))
cmt <- createCoMorTraxObject(sce)
cmt <- encodeDiseaseLabels(cmt,
    diseaseCol   = c("disease_AD", "disease_T2D"),
    controlLabel = "control",
    verbose      = FALSE)
head(diseaseMatrix(cmt))
listComorbidityComparisons(cmt)
}
