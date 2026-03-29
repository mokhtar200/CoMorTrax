\name{CoMorTrax-accessors}
\docType{methods}
\alias{diseaseMatrix}
\alias{diseaseMatrix<-}
\alias{diseaseGroups}
\alias{diseaseGroups<-}
\alias{deResults}
\alias{deResults<-}
\alias{interactionScores}
\alias{interactionScores<-}
\alias{pathwayResults}
\alias{pathwayResults<-}
\alias{cellVulnerability}
\alias{cellVulnerability<-}
\alias{networks}
\alias{networks<-}
\alias{classifier}
\alias{classifier<-}
\title{Accessor and Replacement Methods for CoMorTraxObject}
\description{
  Getter and setter methods for all analysis-result slots of a
  \code{\link{CoMorTraxObject}}. Each getter returns the content of the
  corresponding slot; each setter (the \code{<-} form) replaces it.
  These slots are normally populated by the dedicated analysis functions
  and should rarely need direct assignment.
}
\usage{
diseaseMatrix(x)
diseaseMatrix(x) <- value

diseaseGroups(x)
diseaseGroups(x) <- value

deResults(x)
deResults(x, comparison = NULL)
deResults(x) <- value

interactionScores(x)
interactionScores(x) <- value

pathwayResults(x)
pathwayResults(x) <- value

cellVulnerability(x)
cellVulnerability(x) <- value

networks(x)
networks(x) <- value

classifier(x)
classifier(x) <- value
}
\arguments{
  \item{x}{A \code{\link{CoMorTraxObject}}.}
  \item{comparison}{Character. Optional comparison name (e.g.
    \code{"AD_T2D_vs_Control"}) to retrieve a specific DE result list.
    When \code{NULL} (default), the full nested list is returned.}
  \item{value}{Replacement value. Must be of the correct type for the
    slot (see \code{\link{CoMorTraxObject-class}} for slot types).}
}
\value{
  Getters return the slot contents:
  \describe{
    \item{\code{diseaseMatrix}}{A \code{DataFrame} (cells × diseases) of
      binary indicators.}
    \item{\code{diseaseGroups}}{A named \code{list} of cell-index vectors
      per disease state.}
    \item{\code{deResults}}{A nested \code{list} of DE result
      \code{data.frame}s.}
    \item{\code{interactionScores}}{A \code{DataFrame} with interaction
      delta statistics per gene.}
    \item{\code{pathwayResults}}{A nested \code{list} of enrichment
      results.}
    \item{\code{cellVulnerability}}{A \code{DataFrame} of vulnerability
      scores per cell type.}
    \item{\code{networks}}{A named \code{list} of \pkg{igraph} objects,
      one per disease state.}
    \item{\code{classifier}}{A \code{list} with the trained model,
      CV metrics, and feature importance.}
  }

  Setters return the modified \code{\link{CoMorTraxObject}} invisibly.
}
\author{Ahmed Mokhtar Ramzy Salem}
\seealso{
  \code{\link{CoMorTraxObject-class}},
  \code{\link{createCoMorTraxObject}}
}
\examples{
library(CoMorTrax)
set.seed(1)
sce <- simulateCoMorTrax(n_cells = 100, n_genes = 80,
                          diseases = c("AD","T2D"))
cmt <- createCoMorTraxObject(sce)

# Disease matrix is empty before encoding
dim(diseaseMatrix(cmt))

# Encode and inspect
cmt <- encodeDiseaseLabels(cmt,
    diseaseCol = c("disease_AD","disease_T2D"),
    verbose = FALSE)
head(diseaseMatrix(cmt))
names(diseaseGroups(cmt))
}
