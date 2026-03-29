## =============================================================================
## network_modeling.R - Gene Regulatory Network Modeling (Module 6)
## Author: Ahmed Mokhtar Ramzy Salem
## =============================================================================

#' @describeIn buildComorbidityNetworks Build GRNs for each disease state.
#'
#' @description
#' \strong{Module 6: Gene Regulatory Network Modeling}
#'
#' Constructs gene co-expression / regulatory networks for each disease group
#' and identifies network rewiring under comorbidity - hub genes gained or
#' lost, edge weight changes, and community structure reorganization.
#'
#' Network construction methods:
#' \describe{
#'   \item{\code{"correlation"}}{Pearson/Spearman correlation-based network.
#'     Edges = significant gene-gene correlations.}
#'   \item{\code{"WGCNA"}}{Weighted Gene Co-expression Network Analysis
#'     (requires \pkg{WGCNA}). Produces module-aware networks.}
#'   \item{\code{"GENIE3"}}{Tree-based regression network inference
#'     (requires \pkg{GENIE3}).}
#' }
#'
#' @param object A \code{\link{CoMorTraxObject}} with encoded disease labels
#'   and (optionally) interaction scores.
#' @param method Character; network construction method. One of
#'   \code{"correlation"} (default), \code{"WGCNA"}, \code{"GENIE3"}.
#' @param corrMethod Character; for \code{method = "correlation"}, the
#'   correlation measure: \code{"pearson"} or \code{"spearman"}.
#' @param corrThreshold Numeric; minimum absolute correlation to retain an
#'   edge. Default \code{0.3}.
#' @param pThreshold Numeric; p-value threshold for correlation edges.
#'   Default \code{0.05}.
#' @param geneSubset Character vector or \code{NULL}; restrict to a set of
#'   genes (e.g., interaction genes) for computational efficiency.
#'   If \code{NULL}, top HVGs are used (up to \code{maxGenes}).
#' @param maxGenes Integer; maximum genes for network computation.
#'   Default \code{1000L}.
#' @param states Character vector or \code{NULL}; disease states for which
#'   to build networks. If \code{NULL}, all groups are used.
#' @param minCellsNetwork Integer; minimum cells per state. Default \code{20L}.
#' @param computeRewiring Logical; compute inter-network rewiring statistics.
#'   Default \code{TRUE}.
#' @param BPPARAM A \code{\link[BiocParallel]{BiocParallelParam}} object.
#' @param verbose Logical. Default \code{TRUE}.
#' @param ... Ignored.
#'
#' @return An updated \code{\link{CoMorTraxObject}} with:
#'   \describe{
#'     \item{\code{networks}}{Named list of \code{igraph} objects, one per state.}
#'     \item{\code{networkDiff}}{List of rewiring statistics between states.}
#'   }
#'
#' @author Ahmed Mokhtar Ramzy Salem
#' @export
#' @importFrom igraph graph_from_adjacency_matrix V E degree betweenness
#'   cluster_louvain modularity edge_density graph.difference
#'   induced_subgraph delete.vertices
#' @examples
#' \dontrun{
#' obj <- buildComorbidityNetworks(obj,
#'                                 method = "correlation",
#'                                 corrThreshold = 0.3,
#'                                 maxGenes = 500)
#' }
setMethod("buildComorbidityNetworks", "CoMorTraxObject",
    function(object,
             method          = "correlation",
             corrMethod      = "pearson",
             corrThreshold   = 0.3,
             pThreshold      = 0.05,
             geneSubset      = NULL,
             maxGenes        = 1000L,
             states          = NULL,
             minCellsNetwork = 5L,
             computeRewiring = TRUE,
             BPPARAM         = BiocParallel::SerialParam(),
             verbose         = TRUE,
             ...) {

        .check_disease_encoded(object)
        if (verbose)
            .msg("Module 6: Gene Regulatory Network Modeling ...", "step")

        grps  <- diseaseGroups(object)
        if (is.null(states)) states <- names(grps)

        # Gene subset
        if (is.null(geneSubset)) {
            rd <- SummarizedExperiment::rowData(object)
            if ("isHVG" %in% colnames(rd))
                geneSubset <- rownames(object)[rd$isHVG]
            else
                geneSubset <- rownames(object)
        }
        if (length(geneSubset) > maxGenes) {
            if (verbose)
                .msg(sprintf("  Subsetting to top %d genes.", maxGenes), "info")
            # Prefer interaction genes if available
            if (nrow(interactionScores(object)) > 0) {
                int_genes <- unique(as.character(
                    interactionScores(object)$gene[
                        interactionScores(object)$classification != "additive"]))
                priority <- intersect(int_genes, geneSubset)
                rest     <- setdiff(geneSubset, priority)
                geneSubset <- c(priority, rest)[seq_len(maxGenes)]
            } else {
                geneSubset <- geneSubset[seq_len(maxGenes)]
            }
        }

        expr_mat <- .get_norm_matrix(object)[geneSubset, , drop = FALSE]

        # Build one network per state
        net_list <- list()

        for (state in states) {
            cell_idx <- grps[[state]]
            if (length(cell_idx) < minCellsNetwork) {
                if (verbose)
                    message(sprintf("  Skipping '%s': < %d cells.",
                                    state, minCellsNetwork))
                next
            }

            sub_mat <- expr_mat[, cell_idx, drop = FALSE]

            if (verbose)
                .msg(sprintf("  Building %s network (%d cells) ...",
                             state, length(cell_idx)), "info")

            adj_mat <- switch(method,
                correlation = .build_corr_network(sub_mat, corrMethod,
                                                  corrThreshold, pThreshold),
                WGCNA       = .build_wgcna_network(sub_mat, verbose),
                GENIE3      = .build_genie3_network(sub_mat, BPPARAM),
                stop("method must be 'correlation', 'WGCNA', or 'GENIE3'.")
            )

            g <- igraph::graph_from_adjacency_matrix(
                adj_mat, mode = "undirected",
                weighted = TRUE, diag = FALSE)

            # Replace any NaN/NA edge weights with 0
            ew <- igraph::E(g)$weight
            ew[is.nan(ew) | is.na(ew)] <- 0
            igraph::E(g)$weight <- ew

            # Node attributes ??? use weights=NA to ignore weights if any are 0
            igraph::V(g)$degree      <- igraph::degree(g)
            igraph::V(g)$betweenness <- igraph::betweenness(g,
                weights = NA_real_, normalized = TRUE)

            # Community detection (requires positive weights)
            w_pos <- pmax(igraph::E(g)$weight, 1e-9)
            comm <- igraph::cluster_louvain(g, weights = w_pos)
            igraph::V(g)$community   <- igraph::membership(comm)
            igraph::V(g)$module_size <- as.integer(
                table(igraph::membership(comm))[
                    as.character(igraph::membership(comm))])

            # Graph-level attributes
            igraph::graph_attr(g, "state")      <- state
            igraph::graph_attr(g, "n_cells")    <- length(cell_idx)
            igraph::graph_attr(g, "modularity") <- igraph::modularity(comm)
            igraph::graph_attr(g, "density")    <- igraph::edge_density(g)

            net_list[[state]] <- g
            if (verbose)
                .msg(sprintf("    %d nodes, %d edges, modularity=%.3f",
                             igraph::vcount(g), igraph::ecount(g),
                             igraph::modularity(comm)), "info")
        }

        networks(object) <- net_list

        # Rewiring analysis
        if (computeRewiring && length(net_list) >= 2) {
            if (verbose) .msg("  Computing network rewiring ...", "info")
            object@networkDiff <- .compute_rewiring(net_list, verbose)
        }

        if (verbose) .msg("Network modeling complete.", "success")
        object
    }
)

## --- Network constructors ----------------------------------------------------

#' @keywords internal
.build_corr_network <- function(mat, method, threshold, p_thresh) {
    ct <- t(as.matrix(mat))
    n_cells <- nrow(ct)

    # Compute correlation
    cr <- tryCatch(
        stats::cor(ct, method = method, use = "pairwise.complete.obs"),
        error = function(e) stats::cor(ct, method = "spearman",
                                       use = "pairwise.complete.obs")
    )
    cr[is.na(cr)] <- 0

    # Adaptive threshold: if threshold too strict, relax it
    abs_cr  <- abs(cr)
    diag(abs_cr) <- 0
    q80 <- stats::quantile(abs_cr[abs_cr > 0], 0.80, na.rm = TRUE)
    eff_thresh <- min(threshold, max(q80, 0.1))

    # Significance via t-statistic
    df_val <- max(n_cells - 2, 1)
    t_stat <- cr * sqrt(df_val) / sqrt(pmax(1 - cr^2, 1e-10))
    p_mat  <- 2 * stats::pt(-abs(t_stat), df = df_val)
    p_mat[is.na(p_mat)] <- 1

    # Build adjacency: keep edges above threshold and significant
    adj <- abs(cr)   # use absolute correlation as weight
    adj[abs(cr) < eff_thresh] <- 0
    adj[p_mat > p_thresh]     <- 0
    diag(adj) <- 0
    adj
}

#' @keywords internal
.build_wgcna_network <- function(mat, verbose) {
    .require_pkg("WGCNA")
    WGCNA::allowWGCNAThreads(nThreads = 1)
    expr_t   <- t(as.matrix(mat))
    powers   <- c(1:10, seq(12, 20, by = 2))
    sft      <- WGCNA::pickSoftThreshold(expr_t, powerVector = powers,
                                         verbose = 0)
    power    <- sft$powerEstimate
    if (is.na(power)) power <- 6L  # fallback
    adjacency <- WGCNA::adjacency(expr_t, power = power, type = "signed")
    adjacency
}

#' @keywords internal
.build_genie3_network <- function(mat, BPPARAM) {
    .require_pkg("GENIE3")
    weight_mat <- GENIE3::GENIE3(as.matrix(mat),
                                 nTrees = 500, nCores = 1)
    # Symmetrize
    (weight_mat + t(weight_mat)) / 2
}

## --- Rewiring analysis -------------------------------------------------------

#' @keywords internal
.compute_rewiring <- function(net_list, verbose) {
    state_names <- names(net_list)
    rewiring    <- list()

    for (i in seq_len(length(state_names) - 1)) {
        for (j in seq(i + 1, length(state_names))) {
            s1 <- state_names[i]
            s2 <- state_names[j]
            g1 <- net_list[[s1]]
            g2 <- net_list[[s2]]

            key <- paste0(s1, "_vs_", s2)

            # Common genes
            common_genes <- intersect(igraph::V(g1)$name, igraph::V(g2)$name)
            g1_sub <- igraph::induced_subgraph(g1, common_genes)
            g2_sub <- igraph::induced_subgraph(g2, common_genes)

            # Edge overlap
            e1 <- igraph::as_edgelist(g1_sub)
            e2 <- igraph::as_edgelist(g2_sub)
            e1_key <- apply(e1, 1, function(r) paste(sort(r), collapse = "--"))
            e2_key <- apply(e2, 1, function(r) paste(sort(r), collapse = "--"))

            shared_edges  <- length(intersect(e1_key, e2_key))
            gained_edges  <- length(setdiff(e2_key, e1_key))
            lost_edges    <- length(setdiff(e1_key, e2_key))
            jaccard_edges <- shared_edges / max(1,
                length(union(e1_key, e2_key)))

            # Hub rewiring (top 10% by degree)
            deg1 <- sort(igraph::degree(g1_sub), decreasing = TRUE)
            deg2 <- sort(igraph::degree(g2_sub), decreasing = TRUE)
            n_hub <- max(1, round(length(common_genes) * 0.1))
            hubs1 <- names(deg1)[seq_len(n_hub)]
            hubs2 <- names(deg2)[seq_len(n_hub)]
            hub_jaccard <- length(intersect(hubs1, hubs2)) /
                max(1, length(union(hubs1, hubs2)))

            rewiring[[key]] <- list(
                state1       = s1,
                state2       = s2,
                commonGenes  = length(common_genes),
                sharedEdges  = shared_edges,
                gainedEdges  = gained_edges,
                lostEdges    = lost_edges,
                jaccardEdges = jaccard_edges,
                hubJaccard   = hub_jaccard,
                gainedHubs   = setdiff(hubs2, hubs1),
                lostHubs     = setdiff(hubs1, hubs2)
            )

            if (verbose)
                .msg(sprintf("    %s: +%d/-%d edges, hub J=%.2f",
                             key, gained_edges, lost_edges, hub_jaccard),
                     "info")
        }
    }
    rewiring
}

## --- Accessors ---------------------------------------------------------------

#' @title Get Network Hub Genes
#' @description Extract top hub genes (by betweenness centrality) for a
#'   given disease state network.
#' @param object A \code{\link{CoMorTraxObject}} with networks built.
#' @param state Character; disease state name.
#' @param topN Integer; number of hub genes to return.
#' @return A named numeric vector of betweenness scores.
#' @export
getNetworkHubs <- function(object, state, topN = 20L) {
    nets <- networks(object)
    if (length(nets) == 0) stop("Run buildComorbidityNetworks() first.")
    if (!state %in% names(nets))
        stop("State '", state, "' not found.")
    g <- nets[[state]]
    bet <- igraph::betweenness(g, weights = NA_real_, normalized = TRUE)
    sort(bet, decreasing = TRUE)[seq_len(min(topN, length(bet)))]
}

## --- computeRewiring (exported accessor) -------------------------------------

#' @title Retrieve Network Rewiring Statistics Between Two Disease States
#' @description Extracts the pre-computed rewiring statistics (gained/lost
#'   edges, edge Jaccard, hub rewiring) for a specific pair of disease states.
#'   Rewiring is computed inside \code{\link{buildComorbidityNetworks}} when
#'   \code{computeRewiring = TRUE}.
#' @param object A \code{\link{CoMorTraxObject}} with networks built.
#' @param state1 Character; name of the first disease state.
#' @param state2 Character; name of the second disease state.
#' @return A named list with components: \code{state1}, \code{state2},
#'   \code{commonGenes}, \code{sharedEdges}, \code{gainedEdges},
#'   \code{lostEdges}, \code{jaccardEdges}, \code{hubJaccard},
#'   \code{gainedHubs}, \code{lostHubs}.
#' @export
computeRewiring <- function(object, state1, state2) {
    nd <- object@networkDiff
    if (length(nd) == 0)
        stop("No rewiring data found. Run buildComorbidityNetworks() ",
             "with computeRewiring = TRUE first.")
    key  <- paste(state1, state2, sep = "_vs_")
    key2 <- paste(state2, state1, sep = "_vs_")
    if (key %in% names(nd))  return(nd[[key]])
    if (key2 %in% names(nd)) return(nd[[key2]])
    stop("No rewiring comparison found for '", state1, "' vs '", state2, "'.")
}
