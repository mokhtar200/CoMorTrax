\name{validateCoMorTrax}
\alias{validateCoMorTrax}
\title{Multi-Layer Validation of CoMorTrax Results}
\description{
  Performs three layers of validation on a completed
  \code{\link{CoMorTraxObject}}:
  \enumerate{
    \item \strong{Biological validity}: hypergeometric enrichment test
      of discovered interaction genes against a user-supplied list of
      known disease genes.
    \item \strong{Robustness}: bootstrap subsampling of cells; the
      proportion of interaction genes recovered across \code{nBootstrap}
      subsamples forms a stability score.
    \item \strong{Cross-dataset concordance}: if a second
      \code{\link{CoMorTraxObject}} is supplied, computes Jaccard
      similarity of interaction gene sets between datasets.
  }
}
\usage{
validateCoMorTrax(object, knownGenes = NULL, nBootstrap = 100,
    subsampleFrac = 0.8, referenceObject = NULL,
    verbose = TRUE)
}
\arguments{
  \item{object}{A \code{\link{CoMorTraxObject}} with interaction scores
    computed.}
  \item{knownGenes}{Optional character vector of known disease-associated
    genes for hypergeometric testing.}
  \item{nBootstrap}{Integer. Number of bootstrap iterations for
    robustness estimation. Default 100.}
  \item{subsampleFrac}{Numeric in (0, 1). Fraction of cells to sample per
    bootstrap iteration. Default 0.8.}
  \item{referenceObject}{Optional second \code{\link{CoMorTraxObject}} for
    cross-dataset concordance.}
  \item{verbose}{Logical. Default \code{TRUE}.}
}
\value{
  A named \code{list} with components:
  \describe{
    \item{\code{biological}}{Hypergeometric test results (p-value, fold
      enrichment, overlap genes).}
    \item{\code{robustness}}{Numeric stability score in [0, 1]; the mean
      Jaccard similarity of interaction gene sets across bootstrap
      replicates.}
    \item{\code{concordance}}{Jaccard similarity to the reference object
      (if supplied); \code{NA} otherwise.}
  }
}
\author{Ahmed Mokhtar Ramzy Salem}
\seealso{
  \code{\link{quantifyInteractionEffects}},
  \code{\link{simulateCoMorTrax}}
}
\examples{
\dontrun{
library(CoMorTrax)
# (After full pipeline)
known <- c("APOE","CLU","TREM2","TCF7L2","PPARG")
val <- validateCoMorTrax(cmt, knownGenes = known, nBootstrap = 50)
print(val)
}
}
