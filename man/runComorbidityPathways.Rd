\name{runComorbidityPathways}
\docType{methods}
\alias{runComorbidityPathways}
\alias{runComorbidityPathways,CoMorTraxObject-method}
\alias{comparePathways}
\title{Comorbidity-Stratified Pathway Enrichment Analysis}
\description{
  Runs gene-set enrichment analysis separately on the synergistic,
  antagonistic, and additive gene sets produced by
  \code{\link{quantifyInteractionEffects}}, using GO (BP/MF/CC), KEGG,
  and/or Reactome databases via \pkg{clusterProfiler}.
}
\usage{
\S4method{runComorbidityPathways}{CoMorTraxObject}(object,
    genesets = c("GO_BP","KEGG","Reactome"),
    organism = "hsa", minGSSize = 10, maxGSSize = 500,
    pCutoff = 0.05, verbose = TRUE)

comparePathways(object, cellType = NULL, topN = 20)
}
\arguments{
  \item{object}{A \code{\link{CoMorTraxObject}} with interaction scores
    computed (see \code{\link{quantifyInteractionEffects}}).}
  \item{genesets}{Character vector. Any combination of \code{"GO_BP"},
    \code{"GO_MF"}, \code{"GO_CC"}, \code{"KEGG"}, \code{"Reactome"}.}
  \item{organism}{KEGG organism code. \code{"hsa"} (human, default) or
    \code{"mmu"} (mouse).}
  \item{minGSSize}{Integer. Minimum gene-set size. Default 10.}
  \item{maxGSSize}{Integer. Maximum gene-set size. Default 500.}
  \item{pCutoff}{Numeric. BH-adjusted p-value cutoff. Default 0.05.}
  \item{verbose}{Logical. Default \code{TRUE}.}
  \item{cellType}{Character. Cell type to compare. \code{NULL} uses all.}
  \item{topN}{Integer. Number of top pathways per gene set. Default 20.}
}
\value{
  \code{runComorbidityPathways()} returns the updated
  \code{\link{CoMorTraxObject}} with \code{pathwayResults(object)}
  populated.

  \code{comparePathways()} returns a \code{data.frame} with enriched
  pathways side by side for synergistic and antagonistic gene sets.
}
\author{Ahmed Mokhtar Ramzy Salem}
\seealso{
  \code{\link{quantifyInteractionEffects}},
  \code{\link{scoreCellVulnerability}}
}
