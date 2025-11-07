#' Logistic Dorugade Estimator (LDE)
#'
#' Fits the Logistic Dorugade Estimator (LDE) using the correct formula:
#' \deqn{\hat{\beta}_{LDE} = (X'WX + K d I)^{-1}(X'WX)\hat{\beta}_{MLE}}
#'
#' Introduces both ridge-type (K) and Liu-type (d) biasing parameters to
#' handle multicollinearity in logistic regression.
#'
#' @param X Matrix or data frame of predictors (without intercept)
#' @param y Binary response vector (0 or 1)
#' @param K Ridge-type parameter (K >= 0)
#' @param d Liu-type biasing parameter (0 <= d <= 1)
#' @return An object of class "lde_logit"
#' @export
lde_logit <- function(X, y, K = 1, d = 0.5) {
  if (d < 0 || d > 1) stop("d must be between 0 and 1")
  if (K < 0) stop("K must be non-negative")

  # Ensure intercept column
  if (!all(X[, 1] == 1)) {
    X <- cbind(1, X)
  }
  colnames(X) <- c("(Intercept)", paste0("X", 1:(ncol(X) - 1)))
  n <- nrow(X)
  p <- ncol(X)

  # Step 1: Fit standard logistic MLE
  glm_fit <- glm(y ~ X[, -1], family = binomial)
  beta_mle <- coef(glm_fit)
  W <- diag(as.vector(glm_fit$fit * (1 - glm_fit$fit)))

  # Step 2: Compute LDE coefficients
  XWX <- t(X) %*% W %*% X
  I <- diag(p)
  I[1, 1] <- 0  # no penalty on intercept

  beta_lde <- solve(XWX + K * d * I, XWX %*% beta_mle)

  # Step 3: Compute fitted values and model statistics
  eta <- X %*% beta_lde
  mu <- 1 / (1 + exp(-eta))
  eps <- 1e-15
  deviance <- -2 * sum(y * log(mu + eps) + (1 - y) * log(1 - mu + eps))
  null_mu <- mean(y)
  null_dev <- -2 * sum(y * log(null_mu + eps) + (1 - y) * log(1 - null_mu + eps))
  aic <- deviance + 2 * p

  # Step 4: Standard errors (approx)
  fisher_info <- solve(XWX + K * d * I)
  se <- sqrt(diag(fisher_info))
  z_val <- beta_lde / se
  p_val <- 2 * (1 - pnorm(abs(z_val)))

  coef_table <- data.frame(
    Estimate = beta_lde,
    Std.Error = se,
    z.value = z_val,
    Pr...z.. = p_val
  )
  rownames(coef_table) <- colnames(X)

  # Step 5: Return structured object
  result <- list(
    coefficients = coef_table,
    fitted.values = mu,
    deviance = deviance,
    null.deviance = null_dev,
    aic = aic,
    d = d,
    K = K,
    df.residual = n - p,
    df.total = n,
    family = "binomial(logit)",
    converged = glm_fit$converged,
    call = match.call()
  )
  class(result) <- "lde_logit"
  return(result)
}

# --- Summary Method ---
#' @export
summary.lde_logit <- function(object, ...) {
  cat("\nCall:\n")
  print(object$call)

  cat("\nCoefficients:\n")
  printCoefmat(object$coefficients, P.values = TRUE, has.Pvalue = TRUE)

  cat("\nNull Deviance:", round(object$null.deviance, 4),
      "\nResidual Deviance:", round(object$deviance, 4),
      "\nAIC:", round(object$aic, 4),
      "\nParameters: K =", object$K, ", d =", object$d,
      "\nConverged:", object$converged, "\n")
}

# --- Predict Method ---
#' @export
predict.lde_logit <- function(object, newdata,
                              type = c("link", "response", "class"),
                              threshold = 0.5, ...) {
  type <- match.arg(type)
  X <- as.matrix(newdata)
  if (!all(X[, 1] == 1)) {
    X <- cbind(1, X)
  }
  beta <- object$coefficients$Estimate
  eta <- X %*% beta
  mu <- 1 / (1 + exp(-eta))

  if (type == "link") return(drop(eta))
  if (type == "response") return(drop(mu))
  if (type == "class") return(ifelse(mu >= threshold, 1, 0))
}

#' Tune Logistic Dorugade Estimator (LDE)
#'
#' Performs a grid search over (K, d) parameters to find the best
#' combination minimizing AIC or maximizing accuracy.
#'
#' @param X Matrix or data frame of predictors (without intercept)
#' @param y Binary response vector (0 or 1)
#' @param K_seq Sequence of ridge-type values to test (default = seq(0, 2, 0.2))
#' @param d_seq Sequence of Liu-type values to test (default = seq(0, 1, 0.1))
#' @param criterion Criterion for tuning: "AIC" or "accuracy"
#' @param folds Number of folds for cross-validation (if criterion = "accuracy")
#' @return A list containing:
#' \itemize{
#'   \item `best_K` - best ridge-like shrinkage parameter
#'   \item `best_d` - best Liu-type biasing parameter
#'   \item `best_value` - best (lowest AIC or highest accuracy)
#'   \item `results` - full tuning grid
#' }
#' @export
tune_lde <- function(X, y,
                     K_seq = seq(0, 2, 0.2),
                     d_seq = seq(0, 1, 0.1),
                     criterion = c("AIC", "accuracy"),
                     folds = 5) {
  criterion <- match.arg(criterion)
  X <- as.matrix(X)
  n <- length(y)
  grid_results <- data.frame()

  # --- AIC-based tuning ---
  if (criterion == "AIC") {
    for (K in K_seq) {
      for (d in d_seq) {
        fit <- tryCatch(lde_logit(X, y, K = K, d = d),
                        error = function(e) NULL)
        if (!is.null(fit)) {
          grid_results <- rbind(grid_results,
                                data.frame(K = K, d = d, AIC = fit$aic))
        }
      }
    }
    best <- grid_results[which.min(grid_results$AIC), ]
    return(list(best_K = best$K, best_d = best$d,
                best_value = best$AIC,
                criterion = "AIC",
                results = grid_results))
  }

  # --- Accuracy-based tuning ---
  if (criterion == "accuracy") {
    folds_idx <- sample(rep(1:folds, length.out = n))
    for (K in K_seq) {
      for (d in d_seq) {
        accs <- c()
        for (f in 1:folds) {
          X_train <- X[folds_idx != f, ]; y_train <- y[folds_idx != f]
          X_test  <- X[folds_idx == f, ]; y_test  <- y[folds_idx == f]

          fit <- tryCatch(lde_logit(X_train, y_train, K = K, d = d),
                          error = function(e) NULL)
          if (!is.null(fit)) {
            preds <- predict(fit, X_test, type = "class")
            accs <- c(accs, mean(preds == y_test))
          }
        }
        grid_results <- rbind(grid_results,
                              data.frame(K = K, d = d,
                                         accuracy = mean(accs, na.rm = TRUE)))
      }
    }
    best <- grid_results[which.max(grid_results$accuracy), ]
    return(list(best_K = best$K, best_d = best$d,
                best_value = best$accuracy,
                criterion = "accuracy",
                results = grid_results))
  }
}
