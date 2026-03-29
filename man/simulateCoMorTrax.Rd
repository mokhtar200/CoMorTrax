\name{simulateCoMorTrax}
\alias{simulateCoMorTrax}
\title{Simulate a Comorbidity Single-Cell Dataset}
\description{
  Generates a synthetic \code{SingleCellExperiment} with biologically
  realistic count distributions and \emph{injected} synergistic and
  antagonistic signals. The simulation is truth-recoverable: the identities
  of all synergistic and antagonistic genes are known, enabling rigorous
  benchmarking and unit testing of the full CoMorTrax pipeline.

  Counts are sampled from a negative-binomial model with cell-type-specific
  baseline expression. Disease effects are added as multiplicative fold
  changes. Synergistic genes receive an additional amplification in the
  double-disease group beyond the sum of single-disease effects.
  Antagonistic genes receive a corresponding attenuation.
}
\usage{
simulateCoMorTrax(n_cells = 500, n_genes = 300,
    n_syn = 20, n_ant = 15, n_cell_types = 4,
    cell_types = NULL,
    diseases = c("AD","T2D","Dementia"),
    effect_size = 1.5, syn_boost = 2.0, ant_boost = -1.5,
    dispersion = 0.1, seed = 42, verbose = FALSE)
}
\arguments{
  \item{n_cells}{Integer. Number of cells to simulate. Default 500.}
  \item{n_genes}{Integer. Number of genes. Default 300.}
  \item{n_syn}{Integer. Number of injected synergistic genes. Default 20.}
  \item{n_ant}{Integer. Number of injected antagonistic genes. Default 15.}
  \item{n_cell_types}{Integer. Number of cell types (ignored if
    \code{cell_types} is supplied).}
  \item{cell_types}{Optional character vector of cell-type names.
    Default \code{c("Neuron","Microglia","Astrocyte","OPC")}.}
  \item{diseases}{Character vector of disease names to simulate. A
    multi-label comorbidity structure is automatically generated with
    all single, pairwise, and higher-order combinations. Default
    \code{c("AD","T2D","Dementia")}.}
  \item{effect_size}{Numeric. Mean log2FC magnitude for each disease
    effect. Default 1.5.}
  \item{syn_boost}{Numeric. Additional log2FC for synergistic genes in the
    comorbid group (on top of individual effects). Default 2.0.}
  \item{ant_boost}{Numeric. Log2FC adjustment for antagonistic genes in the
    comorbid group (negative = suppression). Default -1.5.}
  \item{dispersion}{Numeric. Negative-binomial dispersion parameter.
    Default 0.1.}
  \item{seed}{Integer. Random seed. Default 42.}
  \item{verbose}{Logical. Print simulation parameters. Default
    \code{FALSE}.}
}
\value{
  A \code{SingleCellExperiment} with:
  \itemize{
    \item \code{assay(sce, "counts")} — raw integer count matrix
      (genes × cells).
    \item \code{colData} columns: \code{cell_type},
      one binary column per disease (\code{disease_<name>}),
      \code{is_control}, \code{disease_state}.
    \item \code{metadata(sce)$sim_params} — list of all simulation
      parameters and the indices of true synergistic / antagonistic genes.
  }
}
\author{Ahmed Mokhtar Ramzy Salem}
\seealso{
  \code{\link{createCoMorTraxObject}},
  \code{\link{validateCoMorTrax}}
}
\examples{
library(CoMorTrax)
set.seed(42)
sce <- simulateCoMorTrax(
  n_cells    = 300,
  n_genes    = 150,
  n_syn      = 15,
  n_ant      = 10,
  diseases   = c("AD","T2D"),
  cell_types = c("Neuron","Microglia","Astrocyte")
)
sce

# True synergistic genes are stored in metadata
true_syn <- metadata(sce)$sim_params$syn_genes
cat("True synergistic genes:", paste(head(true_syn, 5), collapse=", "), "\n")
}
