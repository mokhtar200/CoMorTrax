\name{plotCoMorTrax}
\docType{methods}
\alias{plotCoMorTrax}
\alias{plotCoMorTrax,CoMorTraxObject-method}
\title{Unified Visualisation Interface for CoMorTraxObject}
\description{
  Single entry-point for all CoMorTrax visualisations. The \code{type}
  argument selects the plot type; additional arguments are passed to the
  corresponding plotting function.
}
\usage{
\S4method{plotCoMorTrax}{CoMorTraxObject}(object,
    type = c("umap","interactionHeatmap","volcanoInteraction",
             "classificationBar","pathwayDotplot","vulnerabilityRank",
             "networkHubs","classifierROC","featureImportance","overview"),
    cellType = NULL, topN = 50, labelTopN = 10,
    colorBy = "disease_state", title = NULL,
    synColor = "#E63946", antColor = "#457B9D",
    addColor = "#6c757d", ...)
}
\arguments{
  \item{object}{A \code{\link{CoMorTraxObject}}.}
  \item{type}{Plot type (see Details).}
  \item{cellType}{Character. Filter to a specific cell type where
    applicable.}
  \item{topN}{Integer. Number of top genes / pathways / hubs to display.}
  \item{labelTopN}{Integer. Number of gene labels on volcano plots.}
  \item{colorBy}{Character. \code{colData} column for UMAP colouring.}
  \item{title}{Character. Custom plot title. \code{NULL} uses a default.}
  \item{synColor}{Hex colour for synergistic elements. Default
    \code{"#E63946"}.}
  \item{antColor}{Hex colour for antagonistic elements. Default
    \code{"#457B9D"}.}
  \item{addColor}{Hex colour for additive elements. Default
    \code{"#6c757d"}.}
  \item{...}{Additional arguments forwarded to specific plot functions.}
}
\details{
  \describe{
    \item{\code{"umap"}}{UMAP embedding coloured by \code{colorBy}.}
    \item{\code{"interactionHeatmap"}}{Heatmap of top interaction genes
      × disease states.}
    \item{\code{"volcanoInteraction"}}{Volcano plot: \eqn{\Delta} vs.
      \eqn{-\log_{10}(p_\text{adj})}.}
    \item{\code{"classificationBar"}}{Bar chart of interaction-class
      gene counts per cell type.}
    \item{\code{"pathwayDotplot"}}{Dot plot of comorbidity-enriched
      pathways.}
    \item{\code{"vulnerabilityRank"}}{Ranked bar chart of cell-type
      vulnerability scores.}
    \item{\code{"networkHubs"}}{Bar chart of top hub genes by
      betweenness centrality.}
    \item{\code{"classifierROC"}}{ROC curve for the trained classifier.}
    \item{\code{"featureImportance"}}{Lollipop chart of feature
      importance scores.}
    \item{\code{"overview"}}{Composite summary figure (all key plots).}
  }
}
\value{
  A \pkg{ggplot2} object (or a \code{list} of ggplot objects for
  \code{type = "overview"}).
}
\author{Ahmed Mokhtar Ramzy Salem}
\seealso{
  \code{\link{summaryCoMorTrax}},
  \code{\link{quantifyInteractionEffects}}
}
