## =============================================================================
## disease_encoding.R - Multi-Label Disease Encoding Engine (Module 1)
## Author: Ahmed Mokhtar Ramzy Salem
## =============================================================================

#' @describeIn encodeDiseaseLabels Encode disease labels for a CoMorTraxObject.
#'
#' @description
#' \strong{Module 1: Multi-Label Disease Encoding Engine}
#'
#' Transforms raw disease metadata into a structured binary disease matrix,
#' enumerates all observed disease combinations, and maps each cell to its
#' comorbidity group. This step is \strong{required} before any downstream
#' analysis.
#'
#' @param object A \code{\link{CoMorTraxObject}}.
#' @param diseaseCol Character; name of column(s) in \code{colData(object)}
#'   to use as disease labels. Can be:
#'   \itemize{
#'     \item A single column with multi-label encoding (comma-separated strings,
#'       e.g. "AD,DM") - set \code{multiLabelSep}.
#'     \item A vector of binary (0/1) column names, one per disease.
#'   }
#' @param controlLabel Character; value in \code{diseaseCol} indicating a
#'   healthy control cell (e.g., \code{"control"}, \code{"healthy"}).
#'   Used when \code{diseaseCol} is a single categorical column.
#' @param diseases Character vector; explicit disease names to model. If
#'   \code{NULL}, inferred automatically from the data.
#' @param multiLabelSep Character; separator between disease labels in a
#'   multi-label string column. Default \code{","}.
#' @param minCellsPerGroup Integer; minimum cells required in a disease
#'   group for it to be retained. Groups below this threshold are flagged
#'   but not removed. Default \code{10L}.
#' @param verbose Logical. Default \code{TRUE}.
#' @param ... Ignored.
#'
#' @return An updated \code{\link{CoMorTraxObject}} with:
#' \describe{
#'   \item{\code{diseaseMatrix}}{Binary \code{DataFrame}: rows = cells,
#'     columns = diseases. 1 = disease present, 0 = absent.}
#'   \item{\code{diseaseGroups}}{Named list mapping group IDs
#'     (e.g., \code{"AD_only"}, \code{"AD_DM"}, \code{"control"}) to
#'     cell column indices.}
#'   \item{\code{colData(object)$diseaseGroup}}{Per-cell group annotation.}
#' }
#'
#' @section Disease Matrix Layout:
#' \preformatted{
#'   Cell    AD  DM  Dementia
#'   Cell1    1   0     0      ? "AD_only"
#'   Cell2    1   1     0      ? "AD_DM"
#'   Cell3    0   0     0      ? "control"
#' }
#'
#' @author Ahmed Mokhtar Ramzy Salem
#' @export
#' @importFrom SummarizedExperiment colData colData<-
#' @importFrom S4Vectors DataFrame
#' @examples
#' library(SingleCellExperiment)
#' n_cells <- 60
#' counts_mat <- matrix(rpois(n_cells * 100, 5), nrow = 100, ncol = n_cells)
#' rownames(counts_mat) <- paste0("Gene", seq_len(100))
#' colnames(counts_mat) <- paste0("Cell", seq_len(n_cells))
#' disease_vec <- rep(c("AD", "DM", "AD,DM", "control"), each = 15)
#' sce <- SingleCellExperiment(
#'     assays = list(counts = counts_mat),
#'     colData = DataFrame(disease = disease_vec))
#' obj <- createCoMorTraxObject(sce)
#' obj <- encodeDiseaseLabels(obj, diseaseCol = "disease",
#'                            controlLabel = "control",
#'                            multiLabelSep = ",")
#' diseaseMatrix(obj)
setMethod("encodeDiseaseLabels", "CoMorTraxObject",
    function(object,
             diseaseCol      = "disease",
             controlLabel    = "control",
             diseases        = NULL,
             multiLabelSep   = ",",
             minCellsPerGroup = 10L,
             verbose          = TRUE,
             ...) {

        if (verbose) .msg("Module 1: Multi-Label Disease Encoding ...", "step")

        cd <- SummarizedExperiment::colData(object)
        n_cells <- ncol(object)

        # ------------------------------------------------------------------
        # Case A: one or more binary (0/1) indicator columns
        # ------------------------------------------------------------------
        # Detect if this is a binary column (single col with 0/1 values only)
        .is_binary_col <- function(x) {
            v <- suppressWarnings(as.integer(as.character(x)))
            # Must have at least some non-NA values AND all non-NA must be 0/1
            non_na <- v[!is.na(v)]
            length(non_na) > 0 && all(non_na %in% c(0L, 1L))
        }
        single_binary <- (length(diseaseCol) == 1 &&
                          diseaseCol %in% colnames(cd) &&
                          .is_binary_col(cd[[diseaseCol[1]]]))

        if (length(diseaseCol) > 1 || single_binary) {
            missing_cols <- setdiff(diseaseCol, colnames(cd))
            if (length(missing_cols) > 0)
                stop("diseaseCol column(s) not found in colData: ",
                     paste(missing_cols, collapse = ", "))
            # Strip common prefixes (disease_, Disease_, dx_) from column names
            # to get clean disease labels (e.g. "disease_AD" -> "AD")
            clean_names <- sub("^(?:disease_|Disease_|dx_|dis_)", "",
                               diseaseCol, ignore.case = FALSE)
            disease_names <- clean_names
            dm_mat <- as.data.frame(cd[, diseaseCol, drop = FALSE])
            dm_mat <- lapply(dm_mat, function(x) as.integer(x == 1))
            dm_mat <- as.data.frame(dm_mat)
            # Rename columns to clean disease names
            colnames(dm_mat) <- disease_names
        } else {
            # ------------------------------------------------------------------
            # Case B: single categorical / multi-label string column
            # ------------------------------------------------------------------
            if (!diseaseCol %in% colnames(cd))
                stop("Column '", diseaseCol, "' not found in colData.")

            raw_labels <- as.character(cd[[diseaseCol]])

            # Split multi-label entries
            label_lists <- strsplit(raw_labels, multiLabelSep, fixed = TRUE)
            label_lists <- lapply(label_lists, function(x) trimws(x))

            # Infer disease names (exclude control label)
            if (is.null(diseases)) {
                all_labels <- unique(unlist(label_lists))
                diseases   <- sort(setdiff(all_labels, controlLabel))
            }
            disease_names <- diseases

            # Build binary matrix
            dm_mat <- as.data.frame(
                matrix(0L, nrow = n_cells, ncol = length(disease_names),
                       dimnames = list(colnames(object), disease_names))
            )
            for (i in seq_len(n_cells)) {
                for (d in label_lists[[i]]) {
                    if (d %in% disease_names)
                        dm_mat[i, d] <- 1L
                }
            }
        }

        # ------------------------------------------------------------------
        # Build group labels
        # ------------------------------------------------------------------
        group_labels <- apply(dm_mat, 1, function(row) {
            active <- disease_names[row == 1]
            if (length(active) == 0) return(controlLabel)
            paste(sort(active), collapse = "_")
        })

        # Map groups to cell indices
        unique_groups <- unique(group_labels)
        group_list <- lapply(unique_groups, function(g) {
            which(group_labels == g)
        })
        names(group_list) <- unique_groups

        # Sort groups: control first, then single diseases, then combos
        n_diseases_in_group <- vapply(names(group_list), function(g) {
            if (g == controlLabel) -1L
            else as.integer(length(strsplit(g, "_")[[1]]))
        }, integer(1))
        group_list <- group_list[order(n_diseases_in_group)]

        # ------------------------------------------------------------------
        # QC: warn about small groups
        # ------------------------------------------------------------------
        for (g in names(group_list)) {
            n_g <- length(group_list[[g]])
            if (n_g < minCellsPerGroup && verbose) {
                warning(sprintf(
                    "Group '%s' has only %d cells (< minCellsPerGroup=%d). ",
                    g, n_g, minCellsPerGroup,
                    "Downstream analyses may be unreliable for this group."),
                    call. = FALSE)
            }
        }

        # ------------------------------------------------------------------
        # Store results
        # ------------------------------------------------------------------
        dm_df <- S4Vectors::DataFrame(dm_mat)
        rownames(dm_df) <- colnames(object)

        diseaseMatrix(object) <- dm_df
        diseaseGroups(object) <- group_list

        # Annotate each cell
        SummarizedExperiment::colData(object)$diseaseGroup <- group_labels

        if (verbose) {
            .msg(sprintf("  Diseases modeled: %s",
                         paste(disease_names, collapse = ", ")), "info")
            .msg(sprintf("  Disease groups identified: %d",
                         length(group_list)), "info")
            for (g in names(group_list))
                .msg(sprintf("    %-30s %d cells",
                             g, length(group_list[[g]])), "info")
            .msg("Disease encoding complete.", "success")
        }

        object
    }
)

## --- Helper: enumerate all theoretically possible comparisons ----------------

#' @title List All Comorbidity Comparisons
#' @description
#' Enumerate every pairwise and multi-group comparison supported by the
#' current disease encoding. Useful for inspecting which DE comparisons
#' \code{\link{runComorbidityDE}} will perform.
#'
#' @param object A \code{\link{CoMorTraxObject}} with encoded disease labels.
#' @param controlLabel Character; label for the control group.
#' @return A \code{data.frame} with columns \code{comparison},
#'   \code{numerator}, and \code{denominator}.
#' @export
#' @examples
#' \dontrun{
#' listComorbidityComparisons(obj)
#' }
listComorbidityComparisons <- function(object, controlLabel = "control") {
    .check_disease_encoded(object)
    grps <- names(diseaseGroups(object))
    grps_no_ctrl <- grps[grps != controlLabel]

    rows <- list()

    # 1. Each disease group vs control
    for (g in grps_no_ctrl) {
        rows[[length(rows) + 1]] <- data.frame(
            comparison  = paste0(g, "_vs_", controlLabel),
            numerator   = g,
            denominator = controlLabel,
            stringsAsFactors = FALSE
        )
    }

    # 2. Comorbid groups vs each constituent single-disease group
    multi_groups <- grps_no_ctrl[vapply(grps_no_ctrl, function(g)
        length(strsplit(g, "_")[[1]]) > 1, logical(1))]

    for (mg in multi_groups) {
        diseases_in <- strsplit(mg, "_")[[1]]
        for (d in diseases_in) {
            single_grp <- grps[vapply(grps, function(g) {
                g != controlLabel && g != mg &&
                identical(strsplit(g, "_")[[1]], d)
            }, logical(1))]
            if (length(single_grp) == 1) {
                rows[[length(rows) + 1]] <- data.frame(
                    comparison  = paste0(mg, "_vs_", single_grp),
                    numerator   = mg,
                    denominator = single_grp,
                    stringsAsFactors = FALSE
                )
            }
        }
    }

    if (length(rows) == 0)
        return(data.frame(comparison=character(), numerator=character(),
                          denominator=character(), stringsAsFactors=FALSE))
    do.call(rbind, rows)
}
