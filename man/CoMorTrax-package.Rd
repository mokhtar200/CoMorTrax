\name{CoMorTrax-package}
\alias{CoMorTrax-package}
\docType{package}
\title{CoMorTrax: Comorbidity-Aware Multi-Omics Disease Modeling Engine}
\description{
  \pkg{CoMorTrax} is a Bioconductor framework for explicitly modeling
  multi-disease molecular interactions at single-cell resolution.

  Unlike conventional single-cell tools that treat each disease in
  isolation, \pkg{CoMorTrax} quantifies the \emph{non-additive}
  transcriptional crosstalk that emerges when diseases co-exist in the
  same cell — synergistic amplification, antagonistic suppression, and
  pathway-level comorbidity burden — enabling the discovery of
  comorbidity-specific biomarkers invisible to standard analyses.
}
\details{
  \strong{Core workflow:}
  \enumerate{
    \item \code{\link{createCoMorTraxObject}} — construct the container.
    \item \code{\link{preprocessCoMorTrax}} — QC, normalisation, HVG, PCA, UMAP.
    \item \code{\link{encodeDiseaseLabels}} — binary multi-label disease matrix.
    \item \code{\link{runComorbidityDE}} — DE across all comorbid comparisons.
    \item \code{\link{quantifyInteractionEffects}} — synergy / antagonism scoring.
    \item \code{\link{runComorbidityPathways}} — comorbidity-enriched pathways.
    \item \code{\link{scoreCellVulnerability}} — composite vulnerability ranking.
    \item \code{\link{buildComorbidityNetworks}} — GRN construction and rewiring.
    \item \code{\link{trainComorbidityClassifier}} — ML comorbidity state prediction.
    \item \code{\link{mapSpatialComorbidity}} — spatial projection (optional).
    \item \code{\link{plotCoMorTrax}} / \code{\link{summaryCoMorTrax}} — output.
    \item \code{\link{exportCoMorTrax}} — structured result export.
  }

  \strong{Testing and benchmarking:}
  \code{\link{simulateCoMorTrax}} generates synthetic data with injected
  synergistic / antagonistic signals of known identity, enabling
  truth-recoverable benchmarking of every analysis module.
  \code{\link{validateCoMorTrax}} provides biological, robustness, and
  cross-dataset validation layers.
}
\author{
  Ahmed Mokhtar Ramzy Salem \email{ahmed.salem@bioinformatics.org}

  Maintainer: Ahmed Mokhtar Ramzy Salem \email{ahmed.salem@bioinformatics.org}
}
\references{
  Salem, A.M.R. (2024). CoMorTrax: A comorbidity-aware multi-omics
  framework for modeling disease interactions at single-cell resolution.
  \emph{bioRxiv}.

  Amezquita, R.A., et al. (2020). Orchestrating single-cell analysis
  with Bioconductor. \emph{Nature Methods}, 17, 137–145.
}
\seealso{
  \code{\link{createCoMorTraxObject}},
  \code{\link{simulateCoMorTrax}},
  \code{browseVignettes("CoMorTrax")}
}
\examples{
library(CoMorTrax)
set.seed(42)
sce <- simulateCoMorTrax(
  n_cells  = 200,
  n_genes  = 100,
  n_syn    = 10,
  n_ant    = 8,
  diseases = c("AD","T2D")
)
cmt <- createCoMorTraxObject(sce)
summaryCoMorTrax(cmt)
}
\keyword{package}
