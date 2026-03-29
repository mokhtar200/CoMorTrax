## =============================================================================
## visualization.R - Comprehensive Visualization Module
## Author: Ahmed Mokhtar Ramzy Salem
## =============================================================================

#' @describeIn plotCoMorTrax Unified CoMorTrax plotting interface.
#'
#' @description
#' Produces publication-quality figures for all CoMorTrax analysis modules.
#'
#' @param object A \code{\link{CoMorTraxObject}}.
#' @param type Character; plot type. One of:
#'   \describe{
#'     \item{\code{"umap"}}{UMAP embedding colored by disease group.}
#'     \item{\code{"interactionHeatmap"}}{Heatmap of interaction delta scores
#'       across genes and cell types.}
#'     \item{\code{"volcanoInteraction"}}{Volcano plot of interaction effects:
#'       x = interactionDelta, y = -log10(padj).}
#'     \item{\code{"classificationBar"}}{Bar chart of synergistic/antagonistic/
#'       additive gene counts per cell type.}
#'     \item{\code{"pathwayDotplot"}}{Dot plot comparing pathway enrichment
#'       across gene sets.}
#'     \item{\code{"vulnerabilityRank"}}{Ranked bar chart of cell-type
#'       vulnerability scores.}
#'     \item{\code{"networkHubs"}}{Bar chart of top network hub genes.}
#'     \item{\code{"classifierROC"}}{ROC curve for the trained classifier.}
#'     \item{\code{"featureImportance"}}{Top feature importance bar chart.}
#'     \item{\code{"overview"}}{Framework overview diagram.}
#'   }
#' @param ... Additional arguments passed to the specific plot function:
#'   \describe{
#'     \item{\code{cellType}}{Character; restrict to one cell type.}
#'     \item{\code{topN}}{Integer; show top N features/genes.}
#'     \item{\code{state}}{Character; disease state for network hub plot.}
#'     \item{\code{database}}{Character; pathway database for dot plot.}
#'     \item{\code{colors}}{Named character vector of colors.}
#'     \item{\code{title}}{Character; plot title override.}
#'   }
#' @return A \code{ggplot} object.
#'
#' @author Ahmed Mokhtar Ramzy Salem
#' @export
#' @importFrom ggplot2 ggplot aes aes_string geom_point geom_bar geom_tile
#'   geom_col geom_hline geom_vline geom_text geom_label scale_color_manual
#'   scale_fill_manual scale_color_viridis_d scale_fill_viridis_c
#'   scale_x_continuous scale_y_continuous labs theme theme_classic
#'   theme_bw element_text element_blank element_line coord_flip
#'   facet_wrap facet_grid guides guide_legend arrow unit
#' @importFrom viridis scale_color_viridis scale_fill_viridis viridis
#' @examples
#' \dontrun{
#' # UMAP with disease group coloring
#' p1 <- plotCoMorTrax(obj, type = "umap")
#'
#' # Interaction volcano
#' p2 <- plotCoMorTrax(obj, type = "volcanoInteraction", cellType = "Neuron")
#'
#' # Vulnerability ranking
#' p3 <- plotCoMorTrax(obj, type = "vulnerabilityRank")
#' }
setMethod("plotCoMorTrax", "CoMorTraxObject",
    function(object, type = "umap", ...) {
        args <- list(...)
        switch(type,
            umap               = .plot_umap(object, args),
            interactionHeatmap = .plot_interaction_heatmap(object, args),
            volcanoInteraction = .plot_volcano_interaction(object, args),
            classificationBar  = .plot_classification_bar(object, args),
            pathwayDotplot     = .plot_pathway_dotplot(object, args),
            vulnerabilityRank  = .plot_vulnerability_rank(object, args),
            networkHubs        = .plot_network_hubs(object, args),
            classifierROC      = .plot_classifier_roc(object, args),
            featureImportance  = .plot_feature_importance(object, args),
            overview           = .plot_overview(object, args),
            stop("Unknown plot type: '", type, "'. See ?plotCoMorTrax.")
        )
    }
)

## --- UMAP ---------------------------------------------------------------

#' @keywords internal
.plot_umap <- function(object, args) {
    umap_df <- as.data.frame(
        SingleCellExperiment::reducedDim(object, "UMAP"))
    colnames(umap_df) <- c("UMAP1", "UMAP2")
    cd <- SummarizedExperiment::colData(object)
    umap_df$diseaseGroup <- if ("diseaseGroup" %in% colnames(cd))
        as.character(cd$diseaseGroup)
    else
        rep("unknown", nrow(umap_df))

    n_grps <- length(unique(umap_df$diseaseGroup))
    cols   <- if (!is.null(args$colors)) args$colors
              else setNames(viridis::viridis(n_grps),
                            unique(umap_df$diseaseGroup))
    ttl <- .arg(args, "title", "CoMorTrax: Disease Groups (UMAP)")

    ggplot2::ggplot(umap_df,
        ggplot2::aes(x = UMAP1, y = UMAP2, color = diseaseGroup)) +
    ggplot2::geom_point(size = 0.5, alpha = 0.7) +
    ggplot2::scale_color_manual(values = cols, name = "Disease Group") +
    ggplot2::labs(title = ttl, x = "UMAP 1", y = "UMAP 2") +
    .comortrax_theme()
}

## --- Interaction Heatmap -------------------------------------------------

#' @keywords internal
.plot_interaction_heatmap <- function(object, args) {
    df <- as.data.frame(interactionScores(object))
    if (nrow(df) == 0) stop("Run quantifyInteractionEffects() first.")

    top_n <- .arg(args, "topN", 50L)

    # Select top genes by |delta|
    df <- df[order(-abs(df$interactionDelta)), ]
    top_genes <- unique(df$gene)[seq_len(min(top_n, length(unique(df$gene))))]
    sub <- df[df$gene %in% top_genes, ]

    if (!"cellType" %in% colnames(sub)) sub$cellType <- "all_cells"

    p <- ggplot2::ggplot(sub,
             ggplot2::aes(x = cellType, y = gene,
                          fill = interactionDelta)) +
         ggplot2::geom_tile(color = "grey90", linewidth = 0.2) +
         ggplot2::scale_fill_gradient2(
             low = "#2166AC", mid = "white", high = "#D6604D",
             midpoint = 0, name = "Interaction\nDelta") +
         ggplot2::labs(
             title = "Interaction Effect Heatmap",
             subtitle = "Red = Synergistic | Blue = Antagonistic",
             x = "Cell Type", y = "Gene") +
         ggplot2::theme_bw() +
         ggplot2::theme(
             axis.text.x = ggplot2::element_text(angle = 45, hjust = 1),
             axis.text.y = ggplot2::element_text(size = 6))
    p
}

## --- Volcano: Interaction effects ----------------------------------------

#' @keywords internal
.plot_volcano_interaction <- function(object, args) {
    df <- as.data.frame(interactionScores(object))
    if (nrow(df) == 0) stop("Run quantifyInteractionEffects() first.")

    ct <- .arg(args, "cellType", NULL)
    if (!is.null(ct) && "cellType" %in% colnames(df))
        df <- df[df$cellType == ct, ]

    df$negLog10P <- -log10(df$padj + 1e-300)
    df$col <- df$classification
    colors <- c(synergistic  = "#D6604D",
                antagonistic = "#2166AC",
                additive     = "grey70")

    top_n   <- .arg(args, "topN", 15L)
    top_syn <- head(df[df$classification == "synergistic",
                       c("gene","interactionDelta")],
                   top_n)
    top_ant <- head(df[df$classification == "antagonistic",
                       c("gene","interactionDelta")],
                   top_n)
    top_lab <- rbind(top_syn, top_ant)

    p <- ggplot2::ggplot(df,
             ggplot2::aes(x = interactionDelta, y = negLog10P,
                          color = col)) +
         ggplot2::geom_point(size = 1.2, alpha = 0.6) +
         ggplot2::scale_color_manual(values = colors,
                                     name = "Classification") +
         ggplot2::geom_vline(xintercept = c(-0.5, 0.5),
                             linetype = "dashed", color = "grey40") +
         ggplot2::geom_hline(yintercept = -log10(0.05),
                             linetype = "dotted", color = "grey40") +
         ggplot2::labs(
             title = "Interaction Effect Volcano Plot",
             subtitle = if (!is.null(ct)) paste("Cell type:", ct) else NULL,
             x = "Interaction Delta (Observed - Expected log2FC)",
             y = "-log10(adjusted p-value)") +
         .comortrax_theme()
    p
}

## --- Classification bar chart --------------------------------------------

#' @keywords internal
.plot_classification_bar <- function(object, args) {
    df <- as.data.frame(interactionScores(object))
    if (nrow(df) == 0) stop("Run quantifyInteractionEffects() first.")
    if (!"cellType" %in% colnames(df)) df$cellType <- "all_cells"

    counts <- as.data.frame(table(df$cellType, df$classification))
    colnames(counts) <- c("cellType", "classification", "count")

    colors <- c(synergistic  = "#D6604D",
                antagonistic = "#2166AC",
                additive     = "#4DAC26")

    ggplot2::ggplot(counts,
        ggplot2::aes(x = cellType, y = count, fill = classification)) +
    ggplot2::geom_col(position = "stack") +
    ggplot2::scale_fill_manual(values = colors, name = "Gene Class") +
    ggplot2::coord_flip() +
    ggplot2::labs(
        title = "Interaction Gene Classification by Cell Type",
        x = "Cell Type", y = "Number of Genes") +
    .comortrax_theme()
}

## --- Pathway dot plot ----------------------------------------------------

#' @keywords internal
.plot_pathway_dotplot <- function(object, args) {
    database <- .arg(args, "database", "GO_BP")
    topN     <- .arg(args, "topN", 15L)
    comp_df  <- comparePathways(object, database = database, topN = topN)
    if (nrow(comp_df) == 0) stop("No pathway results to plot.")

    comp_df$geneRatioNum <- vapply(comp_df$geneRatio, function(r) {
        parts <- strsplit(r, "/")[[1]]
        as.numeric(parts[1]) / as.numeric(parts[2])
    }, numeric(1))

    colors <- c(synergistic  = "#D6604D",
                antagonistic = "#2166AC",
                additive     = "#4DAC26")

    ggplot2::ggplot(comp_df,
        ggplot2::aes(x = geneSet, y = Term,
                     size = geneRatioNum, color = geneSet)) +
    ggplot2::geom_point() +
    ggplot2::scale_color_manual(values = colors, name = "Gene Set") +
    ggplot2::scale_size_continuous(name = "Gene Ratio", range = c(2, 8)) +
    ggplot2::labs(
        title = paste("Comorbidity Pathway Enrichment:", database),
        subtitle = "Pathways active under specific interaction types",
        x = "Interaction Gene Set", y = "Pathway") +
    .comortrax_theme() +
    ggplot2::theme(axis.text.y = ggplot2::element_text(size = 8))
}

## --- Vulnerability rank bar chart ----------------------------------------

#' @keywords internal
.plot_vulnerability_rank <- function(object, args) {
    cv <- as.data.frame(cellVulnerability(object))
    if (nrow(cv) == 0) stop("Run scoreCellVulnerability() first.")

    cv$cellType <- factor(cv$cellType,
                          levels = rev(cv$cellType[order(cv$vulnerabilityScore)]))

    ggplot2::ggplot(cv,
        ggplot2::aes(x = cellType, y = vulnerabilityScore,
                     fill = vulnerabilityScore)) +
    ggplot2::geom_col(width = 0.7) +
    ggplot2::scale_fill_viridis_c(option = "plasma",
                                   name = "Vulnerability\nScore") +
    ggplot2::geom_text(ggplot2::aes(label = sprintf("%.2f", vulnerabilityScore)),
                       hjust = -0.1, size = 3.5) +
    ggplot2::coord_flip() +
    ggplot2::labs(
        title = "Cell-Type Vulnerability to Comorbidity",
        subtitle = paste("Composite score: synergistic + antagonistic",
                         "genes + pathway disruption"),
        x = "Cell Type", y = "Vulnerability Score") +
    ggplot2::ylim(0, max(cv$vulnerabilityScore) * 1.15) +
    .comortrax_theme()
}

## --- Network hub bar chart -----------------------------------------------

#' @keywords internal
.plot_network_hubs <- function(object, args) {
    state <- .arg(args, "state", names(networks(object))[1])
    topN  <- .arg(args, "topN", 20L)
    hubs  <- getNetworkHubs(object, state = state, topN = topN)
    df    <- data.frame(gene = names(hubs), betweenness = hubs,
                        stringsAsFactors = FALSE)
    df$gene <- factor(df$gene, levels = rev(df$gene))

    ggplot2::ggplot(df, ggplot2::aes(x = gene, y = betweenness,
                                      fill = betweenness)) +
    ggplot2::geom_col(width = 0.7) +
    ggplot2::scale_fill_viridis_c(option = "magma", name = "Betweenness") +
    ggplot2::coord_flip() +
    ggplot2::labs(
        title = paste("Network Hub Genes:", state),
        subtitle = "Normalized betweenness centrality",
        x = "Gene", y = "Betweenness Centrality") +
    .comortrax_theme()
}

## --- ROC curve -----------------------------------------------------------

#' @keywords internal
.plot_classifier_roc <- function(object, args) {
    clf <- classifier(object)
    if (length(clf) == 0) stop("Run trainComorbidityClassifier() first.")

    preds <- clf$predictions
    cd    <- SummarizedExperiment::colData(object)
    y     <- factor(preds$actual)

    roc_data <- lapply(levels(y), function(lv) {
        bin_y <- as.numeric(y == lv)
        if (!is.null(clf$predictions) && lv %in% colnames(clf$model$votes)) {
            probs <- clf$model$votes[, lv]
            roc_obj <- pROC::roc(bin_y, probs, quiet = TRUE)
            data.frame(
                fpr   = 1 - roc_obj$specificities,
                tpr   = roc_obj$sensitivities,
                class = lv,
                auc   = as.numeric(pROC::auc(roc_obj))
            )
        } else NULL
    })
    roc_data <- do.call(rbind, Filter(Negate(is.null), roc_data))

    if (nrow(roc_data) == 0) {
        message("ROC data unavailable; returning empty plot.")
        return(ggplot2::ggplot() +
               ggplot2::labs(title = "ROC - no data") +
               .comortrax_theme())
    }

    auc_labels <- unique(roc_data[, c("class", "auc")])
    auc_labels$label <- sprintf("%s (AUC=%.2f)",
                                auc_labels$class, auc_labels$auc)

    ggplot2::ggplot(roc_data,
        ggplot2::aes(x = fpr, y = tpr, color = class)) +
    ggplot2::geom_line(linewidth = 1) +
    ggplot2::geom_abline(slope = 1, intercept = 0,
                         linetype = "dashed", color = "grey60") +
    ggplot2::scale_color_viridis_d(
        name = "Disease Group",
        labels = setNames(auc_labels$label, auc_labels$class)) +
    ggplot2::labs(
        title = "Comorbidity Classifier ROC Curves",
        subtitle = paste("Method:", clf$method),
        x = "False Positive Rate", y = "True Positive Rate") +
    ggplot2::scale_x_continuous(limits = c(0, 1)) +
    ggplot2::scale_y_continuous(limits = c(0, 1)) +
    .comortrax_theme()
}

## --- Feature importance --------------------------------------------------

#' @keywords internal
.plot_feature_importance <- function(object, args) {
    clf <- classifier(object)
    if (length(clf) == 0) stop("Run trainComorbidityClassifier() first.")
    if (is.null(clf$featureImportance))
        stop("Feature importance only available for randomForest.")

    topN <- .arg(args, "topN", 30L)
    fi   <- clf$featureImportance[seq_len(min(topN,
                                               length(clf$featureImportance)))]
    df   <- data.frame(gene = names(fi), importance = as.numeric(fi),
                       stringsAsFactors = FALSE)
    df$gene <- factor(df$gene, levels = rev(df$gene))

    ggplot2::ggplot(df,
        ggplot2::aes(x = gene, y = importance, fill = importance)) +
    ggplot2::geom_col(width = 0.7) +
    ggplot2::scale_fill_viridis_c(option = "viridis", name = "Importance") +
    ggplot2::coord_flip() +
    ggplot2::labs(
        title = "Top Classifier Feature Importances",
        subtitle = "Mean Decrease Accuracy (MDA)",
        x = "Gene", y = "Feature Importance (MDA)") +
    .comortrax_theme()
}

## --- Overview diagram (text-based) ---------------------------------------

#' @keywords internal
.plot_overview <- function(object, args) {
    steps <- data.frame(
        x     = rep(1, 8),
        y     = seq(8, 1, by = -1),
        label = c(
            "1. Multi-Label Disease Encoding",
            "2. Comorbidity-Aware Differential Expression",
            "3. Interaction Effect Quantification\n   (Synergistic | Antagonistic | Additive)",
            "4. Pathway-Level Comorbidity Analysis",
            "5. Cell-Type Vulnerability Scoring",
            "6. Gene Regulatory Network Modeling",
            "7. Machine Learning Comorbidity Prediction",
            "8. Spatial Comorbidity Mapping (optional)"
        ),
        done = c(
            nrow(diseaseMatrix(object))    > 0,
            length(deResults(object))      > 0,
            nrow(interactionScores(object)) > 0,
            length(pathwayResults(object)) > 0,
            nrow(cellVulnerability(object)) > 0,
            length(networks(object))       > 0,
            length(classifier(object))     > 0,
            length(object@spatialResults) > 0
        ),
        stringsAsFactors = FALSE
    )
    steps$status <- ifelse(steps$done, "? Complete", "? Pending")
    steps$color  <- ifelse(steps$done, "#2CA02C", "#AEC7E8")

    ggplot2::ggplot(steps,
        ggplot2::aes(x = x, y = y, label = paste(status, label),
                     color = color)) +
    ggplot2::geom_text(hjust = 0, size = 4.5, fontface = "bold") +
    ggplot2::scale_color_identity() +
    ggplot2::xlim(0.9, 3) + ggplot2::ylim(0, 9) +
    ggplot2::labs(
        title = "CoMorTrax Analysis Pipeline",
        subtitle = "Author: Ahmed Mokhtar Ramzy Salem") +
    ggplot2::theme_void() +
    ggplot2::theme(
        plot.title    = ggplot2::element_text(size = 16, face = "bold",
                                               hjust = 0.5),
        plot.subtitle = ggplot2::element_text(size = 11, hjust = 0.5,
                                               color = "grey40"))
}

## --- Theme helper --------------------------------------------------------

#' @keywords internal
.comortrax_theme <- function() {
    ggplot2::theme_classic() +
    ggplot2::theme(
        plot.title    = ggplot2::element_text(size = 14, face = "bold"),
        plot.subtitle = ggplot2::element_text(size = 10, color = "grey40"),
        axis.title    = ggplot2::element_text(size = 11),
        axis.text     = ggplot2::element_text(size = 9),
        legend.title  = ggplot2::element_text(size = 10),
        legend.text   = ggplot2::element_text(size = 9),
        strip.text    = ggplot2::element_text(size = 10, face = "bold")
    )
}

#' @keywords internal
.arg <- function(args, name, default) {
    if (!is.null(args[[name]])) args[[name]] else default
}
