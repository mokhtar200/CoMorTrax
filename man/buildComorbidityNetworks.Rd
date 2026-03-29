\name{buildComorbidityNetworks}
\docType{methods}
\alias{buildComorbidityNetworks}
\alias{buildComorbidityNetworks,CoMorTraxObject-method}
\alias{computeRewiring}
\alias{getNetworkHubs}
\title{Build and Compare Comorbidity-Stratified Gene Regulatory Networks}
\description{
  Constructs gene co-expression or regulatory networks separately for
  each disease state (e.g. \code{"Control"}, \code{"AD"}, \code{"T2D"},
  \code{"AD_T2D"}) and annotates each network with degree centrality,
  betweenness centrality, and Louvain community membership. The
  \code{computeRewiring()} function then detects structural changes
  (gained / lost edges, hub rewiring) between any two states.
}
\usage{
\S4method{buildComorbidityNetworks}{CoMorTraxObject}(object,
    method = c("correlation","WGCNA","GENIE3"),
    states = "auto", corMethod = "pearson",
    corCutoff = 0.6, pCutoff = 0.01, nTopGenes = 500,
    verbose = TRUE)

computeRewiring(object, state1, state2, topHubs = 20)
getNetworkHubs(object, state, top = 20)
}
\arguments{
  \item{object}{A \code{\link{CoMorTraxObject}} with interaction scores
    (for gene selection) and logcounts assay.}
  \item{method}{Network inference method. One of \code{"correlation"}
    (Pearson/Spearman co-expression), \code{"WGCNA"} (requires
    \pkg{WGCNA}), or \code{"GENIE3"} (requires \pkg{GENIE3}).}
  \item{states}{Character vector of disease states for which to build
    networks. \code{"auto"} uses all states in
    \code{diseaseGroups(object)}.}
  \item{corMethod}{Correlation method: \code{"pearson"} (default) or
    \code{"spearman"}.}
  \item{corCutoff}{Numeric. Minimum absolute correlation to retain an
    edge. Default 0.6.}
  \item{pCutoff}{Numeric. Significance threshold for edge retention
    (after FDR correction). Default 0.01.}
  \item{nTopGenes}{Integer. Number of top interaction genes to include
    in the network (selected by \eqn{|\Delta|}). Default 500.}
  \item{verbose}{Logical. Default \code{TRUE}.}
  \item{state1, state2}{Character. Disease state names for rewiring
    comparison.}
  \item{topHubs}{Integer. Number of hub genes to highlight in rewiring
    output. Default 20.}
  \item{state}{Character. Disease state from which to retrieve hubs.}
  \item{top}{Integer. Number of top hub genes to return. Default 20.}
}
\value{
  \code{buildComorbidityNetworks()} returns the updated
  \code{\link{CoMorTraxObject}} with \code{networks(object)} containing
  per-state \pkg{igraph} objects annotated with vertex attributes
  (\code{degree}, \code{betweenness}, \code{community}).

  \code{computeRewiring()} returns a \code{list} with components:
  \code{edge_jaccard}, \code{gained_edges}, \code{lost_edges},
  \code{hub_overlap}, \code{rewired_hubs}.

  \code{getNetworkHubs()} returns a named numeric vector of betweenness
  centrality scores for the top hub genes in the specified state.
}
\author{Ahmed Mokhtar Ramzy Salem}
\seealso{
  \code{\link{buildComorbidityNetworks}},
  \code{\link{plotCoMorTrax}}
}
