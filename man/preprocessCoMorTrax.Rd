\name{preprocessCoMorTrax}
\docType{methods}
\alias{preprocessCoMorTrax}
\alias{preprocessCoMorTrax,CoMorTraxObject-method}
\title{Preprocess a CoMorTraxObject}
\description{
  Runs a complete, reproducible preprocessing pipeline on a
  \code{\link{CoMorTraxObject}}: QC metric computation, mitochondrial
  gene filtering, normalisation, highly variable gene selection, PCA,
  optional batch correction, and UMAP embedding.
}
\usage{
\S4method{preprocessCoMorTrax}{CoMorTraxObject}(object, normMethod = "logNorm",
    nHVG = 2000, nPCs = 30, batchCorrect = FALSE,
    batchMethod = "harmony", batchCol = "batch",
    runUMAP = TRUE, mitoPattern = "^MT-", verbose = TRUE, ...)
}
\arguments{
  \item{object}{A \code{\link{CoMorTraxObject}}.}
  \item{normMethod}{Normalisation method. One of \code{"logNorm"}
    (fast log-normalisation), \code{"scran"} (pooling-based size-factor
    estimation via \pkg{scran}), or \code{"scTransform"} (regularised
    negative-binomial regression, requires \pkg{sctransform}).}
  \item{nHVG}{Integer. Number of highly variable genes to select.
    Default 2000.}
  \item{nPCs}{Integer. Number of principal components to compute.
    Default 30.}
  \item{batchCorrect}{Logical. Whether to apply batch correction.
    Requires a \code{batchCol} column in \code{colData}. Default
    \code{FALSE}.}
  \item{batchMethod}{Batch correction method. One of \code{"harmony"}
    (requires \pkg{harmony}) or \code{"fastMNN"} (requires \pkg{batchelor}).}
  \item{batchCol}{Name of the \code{colData} column containing batch
    labels. Default \code{"batch"}.}
  \item{runUMAP}{Logical. Compute UMAP embedding after PCA. Default
    \code{TRUE}.}
  \item{mitoPattern}{Regular expression for identifying mitochondrial
    genes in \code{rownames}. Default \code{"^MT-"} (human). Use
    \code{"^mt-"} for mouse.}
  \item{verbose}{Logical. Print progress messages. Default \code{TRUE}.}
  \item{...}{Additional arguments passed to internal normalisation
    functions.}
}
\value{
  The input \code{\link{CoMorTraxObject}} updated with:
  \itemize{
    \item \code{logcounts} assay (normalised expression).
    \item \code{reducedDim(object, "PCA")} — principal component scores.
    \item \code{reducedDim(object, "UMAP")} — UMAP coordinates (if
      \code{runUMAP = TRUE}).
    \item Updated \code{colData} with QC metrics (\code{total_counts},
      \code{detected}, \code{subsets_mito_percent}) and cluster
      assignments (\code{leiden_cluster}).
  }
}
\author{Ahmed Mokhtar Ramzy Salem}
\seealso{
  \code{\link{createCoMorTraxObject}},
  \code{\link{encodeDiseaseLabels}}
}
\examples{
\dontrun{
library(CoMorTrax)
set.seed(1)
sce <- simulateCoMorTrax(n_cells = 200, n_genes = 150,
                          diseases = c("AD","T2D"))
cmt <- createCoMorTraxObject(sce) |>
  preprocessCoMorTrax(normMethod = "logNorm", nHVG = 100,
                       runUMAP = FALSE, verbose = FALSE)
}
}
