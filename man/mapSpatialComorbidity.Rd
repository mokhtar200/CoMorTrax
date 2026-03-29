\name{mapSpatialComorbidity}
\docType{methods}
\alias{mapSpatialComorbidity}
\alias{mapSpatialComorbidity,CoMorTraxObject-method}
\alias{plotSpatialComorbidity}
\title{Map Comorbidity Interaction Scores onto Spatial Transcriptomics Data}
\description{
  Projects the gene-level comorbidity interaction signatures (synergistic,
  antagonistic, or additive gene sets) onto a spatial transcriptomics
  dataset, producing per-spot interaction scores that reveal the tissue
  architecture of comorbid pathology.

  Two scoring methods are supported: simple mean expression of the
  interaction gene set (\code{"mean_expr"}) and single-sample gene set
  enrichment (ssGSEA, \code{"ssgsea"}).
}
\usage{
\S4method{mapSpatialComorbidity}{CoMorTraxObject}(object, spatialData,
    method = c("mean_expr","ssgsea"),
    geneSet = c("synergistic","antagonistic","all"),
    cellType = NULL, verbose = TRUE)

plotSpatialComorbidity(object, geneSet = "synergistic",
    title = NULL, lowColor = "#f7f7f7",
    highColor = "#E63946", pointSize = 1.5)
}
\arguments{
  \item{object}{A \code{\link{CoMorTraxObject}} with interaction scores
    computed.}
  \item{spatialData}{A \code{SpatialExperiment} object containing spatial
    transcriptomics data (spots × genes).}
  \item{method}{Scoring method. \code{"mean_expr"} (default) computes the
    mean log-normalised expression of the interaction gene set per spot.
    \code{"ssgsea"} uses single-sample GSEA (requires \pkg{GSVA}).}
  \item{geneSet}{Which interaction gene set to score. One of
    \code{"synergistic"}, \code{"antagonistic"}, or \code{"all"}.}
  \item{cellType}{Character. Restrict to interaction genes from a specific
    cell type. \code{NULL} uses all cell types.}
  \item{verbose}{Logical. Default \code{TRUE}.}
  \item{title}{Character. Plot title.}
  \item{lowColor, highColor}{Hex colour strings for the score gradient.}
  \item{pointSize}{Numeric. Size of spot points in the tissue plot.}
}
\value{
  \code{mapSpatialComorbidity()} returns the updated
  \code{\link{CoMorTraxObject}} with \code{object@spatialResults} populated:
  a named list of per-spot score vectors, one per \code{geneSet}.

  \code{plotSpatialComorbidity()} returns a \pkg{ggplot2} object showing
  spots on tissue coordinates coloured by comorbidity score.
}
\author{Ahmed Mokhtar Ramzy Salem}
\seealso{
  \code{\link{quantifyInteractionEffects}},
  \code{\link{plotCoMorTrax}}
}
