\name{CoMorTraxObject-class}
\docType{class}
\alias{CoMorTraxObject-class}
\alias{CoMorTraxObject}
\title{CoMorTraxObject: Comorbidity-Aware Single-Cell Experiment Container}
\description{
  \code{CoMorTraxObject} is the central S4 data structure of the
  \pkg{CoMorTrax} package. It extends \code{\link[SingleCellExperiment]{SingleCellExperiment}}
  and stores both raw/processed single-cell data and the results of every
  analytical module (disease encoding, differential expression, interaction
  effects, pathway enrichment, cell-type vulnerability, network modeling,
  and machine learning classification).
}
\section{Slots}{
  \describe{
    \item{\code{diseaseMatrix}:}{\code{DataFrame} — binary cells × diseases matrix.
      Each column is a disease indicator (0/1). Populated by
      \code{\link{encodeDiseaseLabels}}.}
    \item{\code{diseaseGroups}:}{\code{list} — named list mapping each
      disease state to its corresponding cell indices. Populated by
      \code{\link{encodeDiseaseLabels}}.}
    \item{\code{deResults}:}{\code{list} — nested list of differential
      expression results. Top level is cell type; second level is comparison
      name. Each element is a \code{data.frame} with columns
      \code{gene}, \code{logFC}, \code{pvalue}, \code{padj},
      \code{significant}, \code{direction}.}
    \item{\code{interactionScores}:}{\code{DataFrame} — interaction effect
      table with columns \code{gene}, \code{cellType}, \code{observed},
      \code{expected}, \code{interactionDelta}, \code{permPvalue},
      \code{permPadj}, \code{interactionClass}.}
    \item{\code{pathwayResults}:}{\code{list} — pathway enrichment results
      per cell type and gene set (synergistic / antagonistic / additive).}
    \item{\code{cellVulnerability}:}{\code{DataFrame} — vulnerability scores
      with columns \code{cellType}, \code{nSynGenes}, \code{nAntGenes},
      \code{pathwayScore}, \code{meanDelta}, \code{vulnerabilityScore},
      \code{rank}.}
    \item{\code{networks}:}{\code{list} — per-state igraph network objects
      plus community and centrality annotations.}
    \item{\code{networkDiff}:}{\code{list} — network rewiring statistics
      between pairs of disease states.}
    \item{\code{classifier}:}{\code{list} — trained classifier model,
      cross-validation metrics, feature importance, and per-cell predictions.}
    \item{\code{spatialResults}:}{\code{list} — spatial transcriptomics
      comorbidity scores per spot.}
    \item{\code{params}:}{\code{list} — all parameters used during object
      creation and analysis.}
    \item{\code{coMorTraxVersion}:}{\code{character} — package version at
      time of object creation.}
  }
}
\section{Extends}{
  \code{\link[SingleCellExperiment]{SingleCellExperiment}}
}
\author{Ahmed Mokhtar Ramzy Salem}
\seealso{
  \code{\link{createCoMorTraxObject}},
  \code{\link{preprocessCoMorTrax}},
  \code{\link{encodeDiseaseLabels}},
  \code{\link{runComorbidityDE}},
  \code{\link{quantifyInteractionEffects}}
}
\examples{
library(CoMorTrax)
set.seed(1)
sce <- simulateCoMorTrax(n_cells = 100, n_genes = 80,
                          diseases = c("AD","T2D"))
cmt <- createCoMorTraxObject(sce)
cmt
}
\keyword{classes}
