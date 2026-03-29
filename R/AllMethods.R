## =============================================================================
## AllMethods.R - S4 Method Implementations for CoMorTrax
## Author: Ahmed Mokhtar Ramzy Salem
## =============================================================================

## --- show method -------------------------------------------------------------

#' @describeIn CoMorTraxObject Display a concise summary of the object.
#' @return Invisibly returns the object. Called for its side-effect of printing.
#' @importFrom methods show
#' @export
setMethod("show", "CoMorTraxObject", function(object) {
    cat("class: CoMorTraxObject\n")
    cat("dim:", nrow(object), "genes x", ncol(object), "cells\n")

    # Disease info
    dm <- object@diseaseMatrix
    if (nrow(dm) > 0) {
        cat("diseases:", paste(colnames(dm), collapse = ", "), "\n")
        grps <- object@diseaseGroups
        if (length(grps) > 0) {
            cat("disease groups (n cells):\n")
            for (g in names(grps)) {
                cat("  ", g, ":", length(grps[[g]]), "\n")
            }
        }
    } else {
        cat("diseases: <not encoded yet - run encodeDiseaseLabels()>\n")
    }

    # Module completion status
    cat("modules completed:\n")
    cat("  [", ifelse(nrow(dm) > 0, "x", " "), "] Disease encoding\n", sep="")
    cat("  [", ifelse(length(object@deResults) > 0, "x", " "),
        "] Differential expression (", length(object@deResults),
        " comparisons)\n", sep="")
    cat("  [", ifelse(nrow(object@interactionScores) > 0, "x", " "),
        "] Interaction effects\n", sep="")
    cat("  [", ifelse(length(object@pathwayResults) > 0, "x", " "),
        "] Pathway analysis\n", sep="")
    cat("  [", ifelse(nrow(object@cellVulnerability) > 0, "x", " "),
        "] Cell vulnerability\n", sep="")
    cat("  [", ifelse(length(object@networks) > 0, "x", " "),
        "] Network modeling (", length(object@networks), " states)\n", sep="")
    cat("  [", ifelse(length(object@classifier) > 0, "x", " "),
        "] Classifier\n", sep="")
    cat("  [", ifelse(length(object@spatialResults) > 0, "x", " "),
        "] Spatial mapping\n", sep="")
    cat("CoMorTrax version:", object@coMorTraxVersion, "\n")
    invisible(object)
})

## --- Accessor methods --------------------------------------------------------

setMethod("diseaseMatrix", "CoMorTraxObject", function(x) x@diseaseMatrix)
setReplaceMethod("diseaseMatrix", "CoMorTraxObject", function(x, value) {
    x@diseaseMatrix <- value
    validObject(x)
    x
})

setMethod("diseaseGroups", "CoMorTraxObject", function(x) x@diseaseGroups)
setReplaceMethod("diseaseGroups", "CoMorTraxObject", function(x, value) {
    x@diseaseGroups <- value
    x
})

setMethod("deResults", "CoMorTraxObject",
    function(x, comparison = NULL) {
        res <- x@deResults
        if (!is.null(comparison)) {
            if (!comparison %in% names(res))
                stop("Comparison '", comparison, "' not found. ",
                     "Available: ", paste(names(res), collapse = ", "))
            return(res[[comparison]])
        }
        res
    }
)
setReplaceMethod("deResults", "CoMorTraxObject", function(x, value) {
    x@deResults <- value
    x
})

setMethod("interactionScores", "CoMorTraxObject",
    function(x) x@interactionScores)
setReplaceMethod("interactionScores", "CoMorTraxObject", function(x, value) {
    x@interactionScores <- value
    validObject(x)
    x
})

setMethod("pathwayResults", "CoMorTraxObject",
    function(x, geneSet = NULL) {
        res <- x@pathwayResults
        if (!is.null(geneSet)) {
            if (!geneSet %in% names(res))
                stop("Gene set '", geneSet, "' not found. ",
                     "Available: ", paste(names(res), collapse = ", "))
            return(res[[geneSet]])
        }
        res
    }
)
setReplaceMethod("pathwayResults", "CoMorTraxObject", function(x, value) {
    x@pathwayResults <- value
    x
})

setMethod("cellVulnerability", "CoMorTraxObject",
    function(x) x@cellVulnerability)
setReplaceMethod("cellVulnerability", "CoMorTraxObject", function(x, value) {
    x@cellVulnerability <- value
    x
})

setMethod("networks", "CoMorTraxObject",
    function(x, state = NULL) {
        nets <- x@networks
        if (!is.null(state)) {
            if (!state %in% names(nets))
                stop("State '", state, "' not found.")
            return(nets[[state]])
        }
        nets
    }
)
setReplaceMethod("networks", "CoMorTraxObject", function(x, value) {
    x@networks <- value
    x
})

setMethod("classifier", "CoMorTraxObject", function(x) x@classifier)
setReplaceMethod("classifier", "CoMorTraxObject", function(x, value) {
    x@classifier <- value
    x
})

## --- summaryCoMorTrax --------------------------------------------------------

#' @export
setMethod("summaryCoMorTrax", "CoMorTraxObject", function(object, ...) {
    cat("======================================================\n")
    cat("  CoMorTrax Analysis Summary\n")
    cat("  Author: Ahmed Mokhtar Ramzy Salem\n")
    cat("======================================================\n\n")

    cat("Dataset:\n")
    cat("  Genes :", nrow(object), "\n")
    cat("  Cells :", ncol(object), "\n")

    dm <- object@diseaseMatrix
    if (nrow(dm) > 0) {
        cat("\nDiseases modeled:\n")
        for (d in colnames(dm)) {
            n_pos <- sum(dm[[d]] == 1, na.rm = TRUE)
            cat(sprintf("  %-20s : %d cells positive\n", d, n_pos))
        }
        cat("\nDisease combinations:\n")
        for (g in names(object@diseaseGroups)) {
            cat(sprintf("  %-30s : %d cells\n", g,
                        length(object@diseaseGroups[[g]])))
        }
    }

    de <- object@deResults
    if (length(de) > 0) {
        cat("\nDE comparisons (", length(de), "):\n", sep = "")
        for (comp in names(de)) {
            n_sig <- if (!is.null(de[[comp]]$padj))
                sum(de[[comp]]$padj < 0.05, na.rm = TRUE)
            else NA
            cat(sprintf("  %-40s : %s sig. genes\n", comp,
                        ifelse(is.na(n_sig), "?", n_sig)))
        }
    }

    is_df <- object@interactionScores
    if (nrow(is_df) > 0) {
        cat("\nInteraction effects:\n")
        tbl <- table(is_df$classification)
        for (cls in names(tbl)) {
            cat(sprintf("  %-15s : %d genes\n", cls, tbl[[cls]]))
        }
    }

    cv <- object@cellVulnerability
    if (nrow(cv) > 0) {
        cat("\nCell-type vulnerability (top 5):\n")
        top5 <- head(cv[order(-cv$vulnerabilityScore), ], 5)
        for (i in seq_len(nrow(top5))) {
            cat(sprintf("  %d. %-25s score = %.3f\n", i,
                        top5$cellType[i], top5$vulnerabilityScore[i]))
        }
    }

    nets <- object@networks
    if (length(nets) > 0) {
        cat("\nNetworks built:\n")
        for (s in names(nets)) {
            g <- nets[[s]]
            cat(sprintf("  %-25s : %d nodes, %d edges\n", s,
                        igraph::vcount(g), igraph::ecount(g)))
        }
    }

    clf <- object@classifier
    if (length(clf) > 0 && !is.null(clf$performance)) {
        cat("\nClassifier performance:\n")
        perf <- clf$performance
        if (!is.null(perf$auc))
            cat(sprintf("  AUC: %.3f\n", perf$auc))
        if (!is.null(perf$accuracy))
            cat(sprintf("  Accuracy: %.3f\n", perf$accuracy))
    }

    cat("\n======================================================\n")
    invisible(object)
})
