## =============================================================================
## pathway_analysis.R - Pathway-Level Comorbidity Analysis (Module 4)
## Author: Ahmed Mokhtar Ramzy Salem
## =============================================================================

#' @describeIn runComorbidityPathways Pathway enrichment for comorbidity gene sets.
#'
#' @description
#' \strong{Module 4: Pathway-Level Interaction Analysis}
#'
#' Performs gene ontology (GO) and KEGG pathway enrichment \emph{separately}
#' for synergistic, antagonistic, and additive gene sets, revealing pathways
#' that are specifically activated, suppressed, or unaffected under comorbidity.
#'
#' The key biological insight is that pathways enriched only in synergistic
#' genes represent \strong{comorbidity-emergent biology} - processes that
#' neither disease drives alone.
#'
#' @param object A \code{\link{CoMorTraxObject}} with interaction scores.
#' @param geneSets Character vector; which gene sets to analyze. Any
#'   combination of \code{"synergistic"}, \code{"antagonistic"},
#'   \code{"additive"}. Default: all three.
#' @param databases Character vector; enrichment databases. Any combination
#'   of \code{"GO_BP"}, \code{"GO_MF"}, \code{"GO_CC"}, \code{"KEGG"},
#'   \code{"Reactome"}. Default: \code{c("GO_BP", "KEGG")}.
#' @param organism Character; organism code for \pkg{clusterProfiler}.
#'   Default \code{"hsa"} (human). Use \code{"mmu"} for mouse.
#' @param pvalueCutoff Numeric; p-value cutoff for enrichment. Default \code{0.05}.
#' @param qvalueCutoff Numeric; q-value cutoff. Default \code{0.2}.
#' @param minGSSize Integer; minimum gene set size. Default \code{10L}.
#' @param maxGSSize Integer; maximum gene set size. Default \code{500L}.
#' @param cellType Character or \code{NULL}; restrict to one cell type's
#'   interaction scores.
#' @param backgroundGenes Character vector or \code{NULL}; background gene
#'   universe. If \code{NULL}, all expressed genes are used.
#' @param verbose Logical. Default \code{TRUE}.
#' @param ... Ignored.
#'
#' @return An updated \code{\link{CoMorTraxObject}} with \code{pathwayResults}
#'   populated. Structure: \code{list(geneSet -> list(database -> enrichResult))}.
#'
#' @author Ahmed Mokhtar Ramzy Salem
#' @export
#' @importFrom clusterProfiler enrichGO enrichKEGG bitr
#' @importFrom org.Hs.eg.db org.Hs.eg.db
#' @examples
#' \dontrun{
#' obj <- runComorbidityPathways(obj,
#'                               geneSets = c("synergistic", "antagonistic"),
#'                               databases = c("GO_BP", "KEGG"),
#'                               organism = "hsa")
#' }
setMethod("runComorbidityPathways", "CoMorTraxObject",
    function(object,
             geneSets       = c("synergistic", "antagonistic", "additive"),
             databases      = c("GO_BP", "KEGG"),
             organism       = "hsa",
             pvalueCutoff   = 0.05,
             qvalueCutoff   = 0.2,
             minGSSize      = 10L,
             maxGSSize      = 500L,
             cellType       = NULL,
             backgroundGenes = NULL,
             verbose        = TRUE,
             ...) {

        if (nrow(interactionScores(object)) == 0)
            stop("Run quantifyInteractionEffects() first.")

        if (verbose)
            .msg("Module 4: Pathway-Level Comorbidity Analysis ...", "step")

        geneSets <- match.arg(geneSets,
                              c("synergistic", "antagonistic", "additive"),
                              several.ok = TRUE)
        databases <- match.arg(databases,
                               c("GO_BP", "GO_MF", "GO_CC", "KEGG", "Reactome"),
                               several.ok = TRUE)

        # Background universe
        if (is.null(backgroundGenes))
            backgroundGenes <- rownames(object)

        # Map gene symbols to Entrez IDs
        if (verbose) .msg("  Converting gene symbols to Entrez IDs ...", "info")

        org_db <- .get_org_db(organism)
        bg_entrez <- .symbol_to_entrez(backgroundGenes, org_db)

        results <- list()

        for (gs in geneSets) {
            if (verbose)
                .msg(sprintf("  Processing %s genes ...", gs), "info")

            gene_set <- switch(gs,
                synergistic  = getSynergisticGenes(object, cellType),
                antagonistic = getAntagonisticGenes(object, cellType),
                additive     = getAdditiveGenes(object, cellType)
            )

            if (length(gene_set) < 5) {
                warning(sprintf("Gene set '%s' has < 5 genes; skipping.", gs))
                next
            }

            gs_entrez <- .symbol_to_entrez(gene_set, org_db)
            results[[gs]] <- list()

            for (db in databases) {
                if (verbose)
                    .msg(sprintf("    %s | %s ...", gs, db), "info")

                enrich_res <- tryCatch(
                    .run_enrichment(gs_entrez, bg_entrez, db, organism,
                                    pvalueCutoff, qvalueCutoff,
                                    minGSSize, maxGSSize, org_db),
                    error = function(e) {
                        warning(sprintf("Enrichment failed for %s|%s: %s",
                                        gs, db, conditionMessage(e)))
                        NULL
                    }
                )

                if (!is.null(enrich_res)) {
                    results[[gs]][[db]] <- enrich_res
                    n_terms <- nrow(as.data.frame(enrich_res))
                    if (verbose)
                        .msg(sprintf("      %d enriched terms found.", n_terms),
                             "info")
                }
            }
        }

        pathwayResults(object) <- results

        if (verbose)
            .msg("Pathway analysis complete.", "success")
        object
    }
)

## --- Internal enrichment runners ---------------------------------------------

#' @keywords internal
.run_enrichment <- function(genes, universe, db, organism,
                            pCutoff, qCutoff, minGS, maxGS, org_db) {

    if (db %in% c("GO_BP", "GO_MF", "GO_CC")) {
        ont <- sub("GO_", "", db)
        clusterProfiler::enrichGO(
            gene          = genes,
            universe      = universe,
            OrgDb         = org_db,
            ont           = ont,
            keyType       = "ENTREZID",
            pAdjustMethod = "BH",
            pvalueCutoff  = pCutoff,
            qvalueCutoff  = qCutoff,
            minGSSize     = minGS,
            maxGSSize     = maxGS,
            readable      = TRUE
        )
    } else if (db == "KEGG") {
        clusterProfiler::enrichKEGG(
            gene          = genes,
            universe      = universe,
            organism      = organism,
            keyType       = "ncbi-geneid",
            pAdjustMethod = "BH",
            pvalueCutoff  = pCutoff,
            qvalueCutoff  = qCutoff,
            minGSSize     = minGS,
            maxGSSize     = maxGS
        )
    } else if (db == "Reactome") {
        .require_pkg("ReactomePA")
        ReactomePA::enrichPathway(
            gene          = genes,
            universe      = universe,
            organism      = .reactome_organism(organism),
            pAdjustMethod = "BH",
            pvalueCutoff  = pCutoff,
            qvalueCutoff  = qCutoff,
            minGSSize     = minGS,
            maxGSSize     = maxGS,
            readable      = TRUE
        )
    }
}

#' @keywords internal
.symbol_to_entrez <- function(genes, org_db) {
    tryCatch({
        mapping <- clusterProfiler::bitr(
            genes, fromType = "SYMBOL", toType = "ENTREZID",
            OrgDb = org_db)
        unique(mapping$ENTREZID)
    }, error = function(e) {
        warning("Gene ID conversion failed: ", conditionMessage(e))
        character(0)
    })
}

#' @keywords internal
.get_org_db <- function(organism) {
    switch(organism,
        hsa = {
            if (requireNamespace("org.Hs.eg.db", quietly = TRUE))
                org.Hs.eg.db::org.Hs.eg.db
            else stop("Install org.Hs.eg.db for human enrichment analysis.")
        },
        mmu = {
            .require_pkg("org.Mm.eg.db")
            org.Mm.eg.db::org.Mm.eg.db
        },
        rno = {
            .require_pkg("org.Rn.eg.db")
            org.Rn.eg.db::org.Rn.eg.db
        },
        stop("Organism '", organism, "' not supported yet.")
    )
}

#' @keywords internal
.reactome_organism <- function(org_code) {
    switch(org_code,
        hsa = "human",
        mmu = "mouse",
        rno = "rat",
        stop("Unsupported Reactome organism: ", org_code)
    )
}

## --- Convenience: compare pathways across gene sets --------------------------

#' @title Compare Pathways Across Comorbidity Gene Sets
#' @description
#' Create a side-by-side comparison of enriched pathways across synergistic,
#' antagonistic, and additive gene sets for a given database.
#'
#' @param object A \code{\link{CoMorTraxObject}} with pathway results.
#' @param database Character; which database to compare.
#' @param topN Integer; top N terms per gene set.
#' @return A \code{data.frame} with columns \code{Term}, \code{geneSet},
#'   \code{pvalue}, \code{qvalue}, \code{geneRatio}.
#' @export
comparePathways <- function(object, database = "GO_BP", topN = 20L) {
    pr <- pathwayResults(object)
    if (length(pr) == 0)
        stop("Run runComorbidityPathways() first.")

    out <- lapply(names(pr), function(gs) {
        db_res <- pr[[gs]][[database]]
        if (is.null(db_res)) return(NULL)
        df <- as.data.frame(db_res)
        if (nrow(df) == 0) return(NULL)
        df <- head(df[order(df$pvalue), ], topN)
        data.frame(
            Term      = df$Description,
            geneSet   = gs,
            pvalue    = df$pvalue,
            qvalue    = df$qvalue,
            geneRatio = df$GeneRatio,
            stringsAsFactors = FALSE
        )
    })
    do.call(rbind, Filter(Negate(is.null), out))
}
