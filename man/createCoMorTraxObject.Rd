\name{createCoMorTraxObject}
\docType{methods}
\alias{createCoMorTraxObject}
\alias{createCoMorTraxObject,SingleCellExperiment-method}
\alias{createCoMorTraxObject,matrix-method}
\alias{createCoMorTraxObject,dgCMatrix-method}
\title{Create a CoMorTraxObject}
\description{
  Constructs a \code{\link{CoMorTraxObject}} from an existing
  \code{SingleCellExperiment}, a dense \code{matrix}, or a sparse
  \code{dgCMatrix} of raw counts.
}
\usage{
\S4method{createCoMorTraxObject}{SingleCellExperiment}(data, params = list(), ...)
\S4method{createCoMorTraxObject}{matrix}(data, colData = NULL, rowData = NULL,
    params = list(), ...)
\S4method{createCoMorTraxObject}{dgCMatrix}(data, colData = NULL, rowData = NULL,
    params = list(), ...)
}
\arguments{
  \item{data}{A \code{SingleCellExperiment}, \code{matrix} (genes × cells),
    or \code{dgCMatrix} of raw integer counts.}
  \item{colData}{Optional \code{DataFrame} or \code{data.frame} of cell
    metadata (rows = cells). Must include disease annotation columns for
    downstream encoding.}
  \item{rowData}{Optional \code{DataFrame} of gene metadata (rows = genes).}
  \item{params}{Named \code{list} of analysis parameters. Recognised keys:
    \code{minCells} (default 3), \code{minGenes} (default 200),
    \code{maxMitoPCT} (default 20).}
  \item{...}{Additional arguments passed to \code{SingleCellExperiment}.}
}
\value{
  A \code{\link{CoMorTraxObject}} with the raw count matrix stored in the
  \code{"counts"} assay slot and all CoMorTrax-specific slots initialised to
  empty defaults.
}
\author{Ahmed Mokhtar Ramzy Salem}
\seealso{
  \code{\link{preprocessCoMorTrax}},
  \code{\link{encodeDiseaseLabels}},
  \code{\link{simulateCoMorTrax}}
}
\examples{
library(CoMorTrax)
set.seed(1)
sce <- simulateCoMorTrax(n_cells = 100, n_genes = 80,
                          diseases = c("AD","T2D"))
cmt <- createCoMorTraxObject(sce)
cmt

# From a raw count matrix
mat <- matrix(rpois(80 * 100, lambda = 5), nrow = 80,
              dimnames = list(paste0("Gene", seq_len(80)),
                              paste0("Cell", seq_len(100))))
cmt2 <- createCoMorTraxObject(mat)
}
