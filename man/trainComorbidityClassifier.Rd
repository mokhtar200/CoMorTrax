\name{trainComorbidityClassifier}
\docType{methods}
\alias{trainComorbidityClassifier}
\alias{trainComorbidityClassifier,CoMorTraxObject-method}
\title{Train a Machine-Learning Comorbidity State Classifier}
\description{
  Trains a supervised classifier to predict a cell's comorbidity disease
  state from its gene expression profile. Features are drawn from the
  interaction gene signature, HVGs, or significant DE genes. Three
  model families are supported: Random Forest, LASSO-regularised
  logistic regression, and Elastic Net. Model performance is estimated
  via stratified k-fold cross-validation with AUC and balanced accuracy.
}
\usage{
\S4method{trainComorbidityClassifier}{CoMorTraxObject}(object,
    method = c("randomForest","lasso","elasticNet"),
    featureSource = c("interaction","hvg","de"),
    nFeatures = 200, nTrees = 500,
    alpha = 1, kFold = 5, seed = 42, verbose = TRUE)
}
\arguments{
  \item{object}{A \code{\link{CoMorTraxObject}} with disease labels
    encoded and at least one upstream analysis completed.}
  \item{method}{Classification algorithm. One of \code{"randomForest"},
    \code{"lasso"}, or \code{"elasticNet"}.}
  \item{featureSource}{Gene selection strategy. \code{"interaction"} uses
    the top interaction genes (by \eqn{|\Delta|}); \code{"hvg"} uses the
    highly variable genes selected during preprocessing; \code{"de"} uses
    the union of significant DE genes across all comparisons.}
  \item{nFeatures}{Integer. Maximum number of features to include.
    Default 200.}
  \item{nTrees}{Integer. Number of trees for \code{"randomForest"}.
    Default 500.}
  \item{alpha}{Numeric in [0, 1]. Elastic-net mixing parameter
    (\code{alpha = 1}: LASSO; \code{alpha = 0}: Ridge). Only used when
    \code{method = "elasticNet"}. Default 1.}
  \item{kFold}{Integer. Number of cross-validation folds. Default 5.}
  \item{seed}{Integer. Random seed for reproducibility. Default 42.}
  \item{verbose}{Logical. Default \code{TRUE}.}
}
\value{
  The updated \code{\link{CoMorTraxObject}} with \code{classifier(object)}
  set to a \code{list} containing:
  \describe{
    \item{\code{model}}{The fitted model object.}
    \item{\code{features}}{Character vector of selected feature genes.}
    \item{\code{cv_auc}}{Mean cross-validated AUC.}
    \item{\code{balanced_accuracy}}{Mean balanced accuracy across folds.}
    \item{\code{confusion_matrix}}{Confusion matrix from held-out
      predictions.}
    \item{\code{feature_importance}}{\code{data.frame} of genes ranked
      by importance.}
    \item{\code{predictions}}{Per-cell predicted class labels (also added
      to \code{colData} as \code{cmt_predicted_state}).}
  }
}
\author{Ahmed Mokhtar Ramzy Salem}
\seealso{
  \code{\link{quantifyInteractionEffects}},
  \code{\link{plotCoMorTrax}}
}
\examples{
\dontrun{
library(CoMorTrax)
# (After full upstream pipeline)
cmt <- trainComorbidityClassifier(cmt,
    method = "randomForest", featureSource = "interaction",
    nTrees = 200, kFold = 3)
cat("AUC:", classifier(cmt)$cv_auc, "\n")
}
}
