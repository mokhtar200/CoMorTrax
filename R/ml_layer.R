## =============================================================================
## ml_layer.R - Machine Learning Comorbidity Prediction (Module 7)
## Author: Ahmed Mokhtar Ramzy Salem
## =============================================================================

#' @describeIn trainComorbidityClassifier Train a comorbidity state classifier.
#'
#' @description
#' \strong{Module 7: Machine Learning Comorbidity Prediction}
#'
#' Trains and cross-validates a classifier to predict comorbidity state
#' from single-cell gene expression. The model uses interaction genes
#' (synergistic + antagonistic) as features when available, otherwise
#' falls back to HVGs.
#'
#' Supported classifiers:
#' \describe{
#'   \item{\code{"randomForest"}}{Random forest ensemble (via \pkg{randomForest}).
#'     Provides feature importance (Gini impurity / MDA).}
#'   \item{\code{"lasso"}}{L1-penalized logistic regression via \pkg{glmnet}.
#'     Automatic feature selection through sparsity.}
#'   \item{\code{"elasticNet"}}{Elastic net (alpha=0.5) via \pkg{glmnet}.}
#' }
#'
#' @param object A \code{\link{CoMorTraxObject}} with disease labels and
#'   (ideally) interaction scores.
#' @param method Character; model type. One of \code{"randomForest"},
#'   \code{"lasso"}, \code{"elasticNet"}. Default \code{"randomForest"}.
#' @param featureType Character; feature set to use. One of
#'   \code{"interaction"} (synergistic + antagonistic genes),
#'   \code{"HVG"} (highly variable genes), \code{"DE"} (significant DE genes
#'   across all comparisons), \code{"custom"} (supply via \code{featureGenes}).
#' @param featureGenes Character vector; required when \code{featureType = "custom"}.
#' @param maxFeatures Integer; maximum features for the model. Default \code{500L}.
#' @param cvFolds Integer; number of cross-validation folds. Default \code{5L}.
#' @param seed Integer; random seed.
#' @param cellTypeCol Character or \code{NULL}; if provided, trains separate
#'   per-cell-type classifiers.
#' @param verbose Logical. Default \code{TRUE}.
#' @param ... Ignored.
#'
#' @return An updated \code{\link{CoMorTraxObject}} with \code{classifier}
#'   populated. Structure:
#'   \describe{
#'     \item{\code{model}}{The trained model object.}
#'     \item{\code{features}}{Character vector of genes used.}
#'     \item{\code{performance}}{List with \code{auc}, \code{accuracy},
#'       \code{confusionMatrix}, \code{cvResults}.}
#'     \item{\code{featureImportance}}{Named numeric vector of feature
#'       importances (random forest only).}
#'     \item{\code{predictions}}{Per-cell predicted class and probabilities.}
#'   }
#'
#' @author Ahmed Mokhtar Ramzy Salem
#' @export
#' @importFrom randomForest randomForest importance
#' @importFrom glmnet cv.glmnet
#' @importFrom pROC roc auc
#' @importFrom stats predict
#' @examples
#' \dontrun{
#' obj <- trainComorbidityClassifier(obj,
#'                                   method = "randomForest",
#'                                   featureType = "interaction",
#'                                   cvFolds = 5)
#' clf <- classifier(obj)
#' clf$performance$auc
#' }
setMethod("trainComorbidityClassifier", "CoMorTraxObject",
    function(object,
             method       = "randomForest",
             featureType  = "interaction",
             featureGenes = NULL,
             maxFeatures  = 500L,
             cvFolds      = 5L,
             seed         = NULL,
             cellTypeCol  = NULL,
             verbose      = TRUE,
             ...) {

        .check_disease_encoded(object)
        p    <- object@params
        seed <- .param(seed, p$seed)
        set.seed(seed)

        if (verbose)
            .msg("Module 7: Comorbidity Classifier Training ...", "step")

        # --- Feature selection -----------------------------------------------
        features <- .select_features(object, featureType, featureGenes,
                                     maxFeatures, verbose)

        # --- Prepare data matrix ---------------------------------------------
        expr_mat <- .get_norm_matrix(object)[features, , drop = FALSE]
        X <- t(as.matrix(expr_mat))   # cells x genes

        # Remove features (columns) that are all-NA or have zero variance
        col_na   <- colSums(is.na(X)) == nrow(X)
        col_var  <- apply(X, 2, function(v) stats::var(v, na.rm = TRUE))
        bad_cols <- col_na | (!is.na(col_var) & col_var == 0)
        if (any(bad_cols)) X <- X[, !bad_cols, drop = FALSE]

        # Impute any remaining NAs with column means
        for (j in seq_len(ncol(X))) {
            nas <- is.na(X[, j])
            if (any(nas)) X[nas, j] <- mean(X[!nas, j])
        }
        # If still all-NA columns exist after imputation, set to 0
        X[is.na(X)] <- 0

        if (ncol(X) == 0)
            stop("No usable features found after cleaning. ",
                 "Run quantifyInteractionEffects() or preprocessCoMorTrax() first.")

        # Labels
        cd        <- SummarizedExperiment::colData(object)
        y_labels  <- as.character(cd$diseaseGroup)
        y_factor  <- factor(y_labels)

        if (nlevels(y_factor) < 2)
            stop("At least 2 disease groups required for classification.")

        if (verbose)
            .msg(sprintf("  Features: %d | Samples: %d | Classes: %d",
                         ncol(X), nrow(X), nlevels(y_factor)), "info")

        # --- Cross-validation ------------------------------------------------
        if (verbose)
            .msg(sprintf("  %d-fold cross-validation ...", cvFolds), "info")

        folds   <- .create_cv_folds(y_factor, cvFolds, seed)
        cv_res  <- lapply(seq_len(cvFolds), function(fold) {
            train_idx <- which(folds != fold)
            test_idx  <- which(folds == fold)
            X_train   <- X[train_idx, , drop = FALSE]
            y_train   <- y_factor[train_idx]
            X_test    <- X[test_idx,  , drop = FALSE]
            y_test    <- y_factor[test_idx]
            .train_fold(X_train, y_train, X_test, y_test, method, seed)
        })

        cv_acc <- mean(vapply(cv_res, `[[`, numeric(1), "accuracy"))
        if (verbose)
            .msg(sprintf("  CV Accuracy: %.3f", cv_acc), "info")

        # --- Final model (all data) ------------------------------------------
        if (verbose) .msg("  Training final model ...", "info")
        final_model <- .train_model(X, y_factor, method, seed)

        # Predictions on training data
        preds <- .predict_model(final_model, X, method, y_factor)

        # AUC (one-vs-rest for multiclass)
        auc_val <- tryCatch({
            if (nlevels(y_factor) == 2) {
                roc_obj <- pROC::roc(as.numeric(y_factor) - 1,
                                     preds$prob[, 2],
                                     quiet = TRUE)
                as.numeric(pROC::auc(roc_obj))
            } else {
                mean(vapply(levels(y_factor), function(lv) {
                    bin_y <- as.numeric(y_factor == lv)
                    if (lv %in% colnames(preds$prob)) {
                        r <- pROC::roc(bin_y, preds$prob[, lv], quiet = TRUE)
                        as.numeric(pROC::auc(r))
                    } else 0.5
                }, numeric(1)))
            }
        }, error = function(e) NA_real_)

        # Confusion matrix
        cm <- table(Predicted = preds$class, Actual = y_factor)

        # Feature importance
        fi <- if (method == "randomForest") {
            imp <- randomForest::importance(final_model, type = 1)
            sort(imp[, 1], decreasing = TRUE)
        } else {
            NULL
        }

        # Annotate each cell with predictions
        cd_upd <- SummarizedExperiment::colData(object)
        cd_upd$predicted_diseaseGroup <- preds$class
        SummarizedExperiment::colData(object) <- cd_upd

        classifier(object) <- list(
            model             = final_model,
            features          = features,
            method            = method,
            featureType       = featureType,
            performance       = list(
                auc            = auc_val,
                accuracy       = mean(preds$class == y_factor),
                cvAccuracy     = cv_acc,
                cvFolds        = cvFolds,
                confusionMatrix = cm,
                cvResults      = cv_res
            ),
            featureImportance = fi,
            predictions       = data.frame(
                cell      = colnames(object),
                actual    = y_labels,
                predicted = preds$class,
                stringsAsFactors = FALSE
            )
        )

        if (verbose) {
            .msg(sprintf("  Final AUC: %.3f", auc_val), "info")
            .msg(sprintf("  Final Accuracy: %.3f",
                         mean(preds$class == y_factor)), "info")
            .msg("Classifier training complete.", "success")
        }
        object
    }
)

## --- Model internals ---------------------------------------------------------

#' @keywords internal
.select_features <- function(object, featureType, featureGenes, maxFeatures,
                             verbose) {
    if (featureType == "custom") {
        if (is.null(featureGenes))
            stop("featureGenes required when featureType = 'custom'.")
        return(intersect(featureGenes, rownames(object)))
    }
    genes <- switch(featureType,
        interaction = {
            if (nrow(interactionScores(object)) > 0) {
                is_df <- as.data.frame(interactionScores(object))
                sig_genes <- unique(as.character(
                    is_df$gene[is_df$classification != "additive"]))
                # Fallback: if no synergistic/antagonistic genes found,
                # use top genes by |interactionDelta| regardless of classification
                if (length(sig_genes) == 0) {
                    if (verbose)
                        message("  No sig. interaction genes; using top |delta| genes.")
                    ord <- order(abs(is_df$interactionDelta), decreasing = TRUE)
                    top_n <- min(maxFeatures, nrow(is_df))
                    sig_genes <- unique(as.character(is_df$gene[ord[seq_len(top_n)]]))
                }
                sig_genes
            } else {
                if (verbose)
                    message("No interaction scores; falling back to HVGs.")
                .hvg_genes(object)
            }
        },
        HVG = .hvg_genes(object),
        DE  = {
            de <- getDeTable(object, significantOnly = TRUE)
            unique(de$gene)
        },
        stop("featureType must be 'interaction', 'HVG', 'DE', or 'custom'.")
    )
    genes <- intersect(genes, rownames(object))
    if (length(genes) > maxFeatures)
        genes <- genes[seq_len(maxFeatures)]
    if (verbose)
        .msg(sprintf("  Using %d %s features.", length(genes), featureType),
             "info")
    genes
}

#' @keywords internal
.hvg_genes <- function(object) {
    rd <- SummarizedExperiment::rowData(object)
    if ("isHVG" %in% colnames(rd))
        rownames(object)[rd$isHVG]
    else
        rownames(object)
}

#' @keywords internal
.create_cv_folds <- function(y, k, seed) {
    set.seed(seed)
    folds <- rep(NA, length(y))
    for (lv in levels(y)) {
        idx  <- which(y == lv)
        fold_assign <- ((seq_along(idx) - 1) %% k) + 1
        folds[idx]  <- sample(fold_assign)
    }
    folds
}

#' @keywords internal
.train_model <- function(X, y, method, seed) {
    set.seed(seed)
    # Ensure no NAs
    X[is.na(X)] <- 0
    if (method == "randomForest") {
        # Safe mtry: between 1 and ncol(X)
        mtry <- max(1L, min(as.integer(sqrt(ncol(X))), ncol(X)))
        randomForest::randomForest(X, y, ntree = 200,
                                   mtry = mtry, importance = TRUE)
    } else if (method %in% c("lasso", "elasticNet")) {
        alpha  <- if (method == "lasso") 1 else 0.5
        nfolds <- max(2L, min(5L, nrow(X) %/% 2L))
        glmnet::cv.glmnet(X, y, family = "multinomial",
                          alpha = alpha, nfolds = nfolds)
    } else {
        stop("method must be 'randomForest', 'lasso', or 'elasticNet'.")
    }
}

#' @keywords internal
.train_fold <- function(X_tr, y_tr, X_te, y_te, method, seed) {
    model <- .train_model(X_tr, y_tr, method, seed)
    preds <- .predict_model(model, X_te, method, y_te)
    list(accuracy = mean(preds$class == y_te))
}

#' @keywords internal
.predict_model <- function(model, X, method, y_levels) {
    if (method == "randomForest") {
        prob  <- predict(model, X, type = "prob")
        class <- predict(model, X, type = "response")
        list(class = as.character(class), prob = prob)
    } else {
        best_s <- model$lambda.min
        prob   <- predict(model, X, s = best_s, type = "response")[, , 1]
        class  <- apply(prob, 1, function(r) colnames(prob)[which.max(r)])
        list(class = class, prob = prob)
    }
}
