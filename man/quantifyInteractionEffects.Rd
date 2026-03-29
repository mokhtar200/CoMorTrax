\name{quantifyInteractionEffects}
\docType{methods}
\alias{quantifyInteractionEffects}
\alias{quantifyInteractionEffects,CoMorTraxObject-method}
\alias{getSynergisticGenes}
\alias{getAntagonisticGenes}
\alias{getAdditiveGenes}
\title{Quantify Non-Additive Interaction Effects in Comorbidity}
\description{
  Computes, for every gene (and, optionally, every cell type), an
  \emph{interaction delta} that captures the degree to which the
  transcriptional response under comorbidity deviates from the sum of
  the two single-disease responses:

  \deqn{\Delta = \text{FC}(D_1 \cap D_2\,\text{vs. Ctrl}) -
    [\text{FC}(D_1\,\text{vs. Ctrl}) + \text{FC}(D_2\,\text{vs. Ctrl})]}

  Genes with \eqn{\Delta > \text{synThreshold}} are classified
  \emph{synergistic}; those with \eqn{\Delta < \text{antThreshold}} are
  \emph{antagonistic}; the remainder are \emph{additive}. Permutation-based
  p-values are computed via \pkg{BiocParallel}.

  This is the analytical core of \pkg{CoMorTrax} and the primary
  contribution absent from existing single-cell tools.
}
\usage{
\S4method{quantifyInteractionEffects}{CoMorTraxObject}(object, disease1, disease2,
    cellTypeCol = "cell_type", synThreshold = 0.25,
    antThreshold = -0.25, nPermutations = 1000,
    adjPcutoff = 0.05, BPPARAM = BiocParallel::SerialParam(),
    verbose = TRUE)

getSynergisticGenes(object, cellType = NULL, padj = 0.05)
getAntagonisticGenes(object, cellType = NULL, padj = 0.05)
getAdditiveGenes(object, cellType = NULL, padj = 0.05)
}
\arguments{
  \item{object}{A \code{\link{CoMorTraxObject}} with DE results populated
    (see \code{\link{runComorbidityDE}}).}
  \item{disease1}{Character. Name of the first disease (must match a column
    in \code{diseaseMatrix(object)}).}
  \item{disease2}{Character. Name of the second disease.}
  \item{cellTypeCol}{Name of the \code{colData} column containing cell-type
    labels. Default \code{"cell_type"}.}
  \item{synThreshold}{Numeric. Minimum \eqn{\Delta} to classify a gene as
    synergistic. Default \code{0.5}.}
  \item{antThreshold}{Numeric. Maximum \eqn{\Delta} to classify a gene as
    antagonistic. Default \code{-0.5}.}
  \item{nPermutations}{Integer. Number of label-permutation iterations for
    null-distribution estimation. Default \code{1000}.}
  \item{adjPcutoff}{Numeric. BH-adjusted p-value threshold. Default
    \code{0.05}.}
  \item{BPPARAM}{A \code{BiocParallelParam} instance controlling
    parallelisation. Default \code{SerialParam()} (single-threaded).}
  \item{verbose}{Logical. Default \code{TRUE}.}
  \item{cellType}{Character. Filter results to a specific cell type.
    \code{NULL} returns all cell types.}
  \item{padj}{Numeric. Adjusted p-value filter applied in the accessor
    functions. Default \code{0.05}.}
}
\value{
  \code{quantifyInteractionEffects()} returns the updated
  \code{\link{CoMorTraxObject}} with \code{interactionScores(object)}
  set to a \code{DataFrame} with columns:
  \describe{
    \item{\code{gene}}{Gene identifier.}
    \item{\code{cellType}}{Cell-type label (or \code{"all"}).}
    \item{\code{observed}}{Observed log2FC in comorbid vs. control.}
    \item{\code{expected}}{Sum of single-disease log2FCs.}
    \item{\code{interactionDelta}}{\eqn{\Delta = \text{observed} -
      \text{expected}}.}
    \item{\code{permPvalue}}{Permutation-derived p-value.}
    \item{\code{permPadj}}{BH-adjusted permutation p-value.}
    \item{\code{interactionClass}}{One of \code{"synergistic"},
      \code{"antagonistic"}, \code{"additive"}.}
  }

  The accessor functions return character vectors of gene names.
}
\author{Ahmed Mokhtar Ramzy Salem}
\references{
  Bliss, C.I. (1939). The toxicity of poisons applied jointly.
  \emph{Annals of Applied Biology}, 26, 585–615.

  Loewe, S. (1953). The problem of synergism and antagonism of combined
  drugs. \emph{Arzneimittelforschung}, 3, 285–290.
}
\seealso{
  \code{\link{runComorbidityDE}},
  \code{\link{runComorbidityPathways}},
  \code{\link{scoreCellVulnerability}}
}
\examples{
\dontrun{
library(CoMorTrax)
set.seed(42)
sce <- simulateCoMorTrax(n_cells = 400, n_genes = 200,
                          n_syn = 20, n_ant = 15,
                          diseases = c("AD","T2D"))
cmt <- createCoMorTraxObject(sce) |>
  encodeDiseaseLabels(diseaseCol = c("disease_AD","disease_T2D")) |>
  runComorbidityDE(method = "wilcoxon", minCells = 5) |>
  quantifyInteractionEffects(disease1 = "AD", disease2 = "T2D",
                              nPermutations = 200)

head(interactionScores(cmt))
getSynergisticGenes(cmt, cellType = "Neuron")
}
}
