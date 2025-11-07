#' Generic predict function for LogisticEstimators package
#'
#' This function dispatches to the appropriate `predict.<estimator>` method
#' depending on the class of the fitted model.
#'
#' @param object A fitted model object (e.g., of class "ridge_logit", "liu_logit", etc.)
#' @param ... Additional arguments passed to specific methods.
#' @export
predict <- function(object, ...) {
  UseMethod("predict")
}

#' Generic summary function
#' @export
summary <- function(object, ...) {
  UseMethod("summary")
}
