\name{summaryCoMorTrax}
\docType{methods}
\alias{summaryCoMorTrax}
\alias{summaryCoMorTrax,CoMorTraxObject-method}
\title{Print a Comprehensive Summary of a CoMorTraxObject}
\description{
  Displays a formatted summary of all analysis modules completed in the
  \code{\link{CoMorTraxObject}}, including dataset dimensions, disease
  group sizes, number of significant DE genes per comparison, interaction
  gene counts by class, pathway counts, cell-type vulnerability rankings,
  and classifier performance.
}
\usage{
\S4method{summaryCoMorTrax}{CoMorTraxObject}(object, ...)
}
\arguments{
  \item{object}{A \code{\link{CoMorTraxObject}}.}
  \item{...}{Currently unused.}
}
\value{
  Invisibly returns a named \code{list} of summary statistics.
  Called for its side-effect of printing.
}
\author{Ahmed Mokhtar Ramzy Salem}
\seealso{\code{\link{plotCoMorTrax}}, \code{\link{exportCoMorTrax}}}
\examples{
library(CoMorTrax)
set.seed(1)
sce <- simulateCoMorTrax(n_cells = 100, n_genes = 80,
                          diseases = c("AD","T2D"))
cmt <- createCoMorTraxObject(sce)
summaryCoMorTrax(cmt)
}
