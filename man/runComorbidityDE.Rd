\name{runComorbidityDE}
\docType{methods}
\alias{runComorbidityDE}
\alias{runComorbidityDE,CoMorTraxObject-method}
\alias{getDeTable}
\title{Comorbidity-Aware Differential Expression}
\description{
  Runs differential expression analysis across all comorbidity comparisons
  enumerated by \code{\link{encodeDiseaseLabels}}, stratified by cell type.
  Three statistical back-ends are supported: Wilcoxon rank-sum (fast,
  non-parametric), edgeR glmQLFTest (negative-binomial GLM), and MAST
  (hurdle model for single-cell data).
}
\usage{
\S4method{runComorbidityDE}{CoMorTraxObject}(object, cellTypeCol = "cell_type",
    comparisons = "auto", method = c("wilcoxon","edgeR","MAST"),
    minCells = 10, logFCcutoff = 0.5, adjPcutoff = 0.05,
    verbose = TRUE, ...)

getDeTable(object, comparison, cellType = NULL)
}
\arguments{
  \item{object}{A \code{\link{CoMorTraxObject}} with disease labels encoded
    (see \code{\link{encodeDiseaseLabels}}).}
  \item{cellTypeCol}{Name of the \code{colData} column containing cell-type
    annotations. Default \code{"cell_type"}.}
  \item{comparisons}{Either \code{"auto"} (run all comparisons returned by
    \code{\link{listComorbidityComparisons}}) or a character vector of
    specific comparison strings such as \code{"AD_T2D_vs_Control"}.}
  \item{method}{Statistical method. One of \code{"wilcoxon"} (default),
    \code{"edgeR"}, or \code{"MAST"}.}
  \item{minCells}{Integer. Minimum number of cells per group required for
    a comparison within a cell type to be attempted. Default 10.}
  \item{logFCcutoff}{Numeric. Minimum absolute log2 fold-change to call a
    gene significant. Default 0.5.}
  \item{adjPcutoff}{Numeric. Adjusted p-value threshold. Default 0.05.}
  \item{verbose}{Logical. Print comparison-level progress. Default
    \code{TRUE}.}
  \item{...}{Additional arguments passed to the chosen DE method.}
  \item{comparison}{Character. A comparison string (e.g.
    \code{"AD_T2D_vs_Control"}) for \code{getDeTable}.}
  \item{cellType}{Character. Cell-type name for \code{getDeTable}. If
    \code{NULL}, results across all cell types are concatenated.}
}
\value{
  \code{runComorbidityDE()} returns the updated
  \code{\link{CoMorTraxObject}} with \code{deResults(object)} populated.

  \code{getDeTable()} returns a \code{data.frame} with columns:
  \code{gene}, \code{logFC}, \code{pvalue}, \code{padj},
  \code{significant} (logical), \code{direction} (\code{"up"},
  \code{"down"}, or \code{"ns"}).
}
\author{Ahmed Mokhtar Ramzy Salem}
\seealso{
  \code{\link{encodeDiseaseLabels}},
  \code{\link{quantifyInteractionEffects}}
}
\examples{
\dontrun{
library(CoMorTrax)
set.seed(1)
sce <- simulateCoMorTrax(n_cells = 300, n_genes = 120,
                          diseases = c("AD","T2D"))
cmt <- createCoMorTraxObject(sce) |>
  encodeDiseaseLabels(diseaseCol = c("disease_AD","disease_T2D"),
                      verbose = FALSE) |>
  runComorbidityDE(method = "wilcoxon", minCells = 5)

de <- getDeTable(cmt, comparison = "AD_T2D_vs_Control",
                 cellType = "Neuron")
head(de)
}
}
