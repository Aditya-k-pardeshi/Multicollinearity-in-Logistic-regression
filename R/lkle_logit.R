#' Logistic Kibria and Lukman Estimator (LKLE)
#'
#' Fits the Logistic Kibria and Lukman Estimator using:
#' \deqn{\hat{\beta}_{LKLE} = (X'WX + kI)^{-1}(X'WX + KI)\hat{\beta}_{MLE}}
#'
#' This is a two-parameter modification of the logistic MLE, designed to handle
#' multicollinearity by combining ridge-type (k) and Liu-type (K) shrinkage.
#' The intercept is not penalized.
#'
#' @param X Predictor matrix or data frame (without intercept)
#' @param y Binary response vector (0 or 1)
#' @param k Ridge shrinkage parameter (k ≥ 0)
#' @param K Liu adjustment parameter (0 ≤ K ≤ 1)
#' @param tol Convergence tolerance (unused, for API consistency)
#' @param max_iter Iterations (unused, for API consistency)
#' @return An object of class "lkle_logit"
#' @export
lkle_logit <- function(X, y, k = 0.1, K = 0.3, tol = 1e-6, max_iter = 50) {
  # --- Intercept handling ---
  if (!all(X[, 1] == 1)) {
    X <- cbind(1, X)
  }
  colnames(X) <- c("(Intercept)", paste0("X", 1:(ncol(X) - 1)))
  n <- nrow(X)
  p <- ncol(X)

  # --- Step 1: Fit standard logistic MLE ---
  glm_fit <- glm(y ~ X[, -1], family = binomial)
  beta_mle <- coef(glm_fit)
  W <- diag(as.vector(glm_fit$fit * (1 - glm_fit$fit)))

  # --- Step 2: Compute LKLE coefficients ---
  XWX <- t(X) %*% W %*% X
  I <- diag(p)
  I[1, 1] <- 0  # no penalty on intercept

  beta_lkle <- solve(XWX + k * I, (XWX + K * I) %*% beta_mle)

  # --- Step 3: Fitted values and deviance ---
  eta <- X %*% beta_lkle
  mu <- 1 / (1 + exp(-eta))
  eps <- 1e-15
  deviance <- -2 * sum(y * log(mu + eps) + (1 - y) * log(1 - mu + eps))
  null_mu <- mean(y)
  null_dev <- -2 * sum(y * log(null_mu + eps) + (1 - y) * log(1 - null_mu + eps))
  aic <- deviance + 2 * p

  # --- Step 4: Standard errors (approx.) ---
  fisher_info <- solve(XWX + k * I)
  se <- sqrt(diag(fisher_info))
  z_val <- beta_lkle / se
  p_val <- 2 * (1 - pnorm(abs(z_val)))

  coef_table <- data.frame(
    Estimate = beta_lkle,
    Std.Error = se,
    z.value = z_val,
    Pr...z.. = p_val
  )
  rownames(coef_table) <- colnames(X)

  # --- Step 5: Return object ---
  result <- list(
    coefficients = coef_table,
    fitted.values = mu,
    deviance = deviance,
    null.deviance = null_dev,
    df.residual = n - p,
    df.total = n,
    aic = aic,
    k = k,
    K = K,
    family = "binomial(logit)",
    converged = glm_fit$converged,
    call = match.call()
  )
  class(result) <- "lkle_logit"
  return(result)
}

# --- Summary Method ---
#' @export
summary.lkle_logit <- function(object, ...) {
  cat("\nCall:\n")
  print(object$call)

  cat("\nCoefficients:\n")
  printCoefmat(object$coefficients, P.values = TRUE, has.Pvalue = TRUE)

  cat("\nNull Deviance:", round(object$null.deviance, 4),
      "\nResidual Deviance:", round(object$deviance, 4),
      "\nAIC:", round(object$aic, 4),
      "\nRidge (k):", object$k,
      "\nLiu (K):", object$K,
      "\nConverged:", object$converged, "\n")
}

# --- Predict Method ---
#' @export
predict.lkle_logit <- function(object, newdata,
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

#' Tune the Logistic Kibria and Lukman Estimator (LKLE)
#'
#' Performs a grid search over k and K parameters to find the best
#' (k, K) combination minimizing AIC or maximizing prediction accuracy.
#'
#' @param X Matrix or data frame of predictors (without intercept)
#' @param y Binary response vector (0/1)
#' @param k_seq Sequence of ridge-type shrinkage values to test (default = seq(0, 2, 0.2))
#' @param K_seq Sequence of Liu-type adjustment values to test (default = seq(0, 1, 0.1))
#' @param criterion Criterion for tuning: "AIC" or "accuracy"
#' @param folds Number of folds for cross-validation if criterion = "accuracy" (default = 5)
#' @return A list containing:
#' \itemize{
#'   \item `best_k` - best ridge shrinkage parameter
#'   \item `best_K` - best Liu adjustment parameter
#'   \item `best_value` - best AIC or accuracy achieved
#'   \item `results` - full tuning grid results
#' }
#' @export
tune_lkle <- function(X, y,
                      k_seq = seq(0, 2, 0.2),
                      K_seq = seq(0, 1, 0.1),
                      criterion = c("AIC", "accuracy"),
                      folds = 5) {
  criterion <- match.arg(criterion)
  X <- as.matrix(X)
  n <- length(y)
  grid_results <- data.frame()

  if (criterion == "AIC") {
    for (k in k_seq) {
      for (K in K_seq) {
        fit <- tryCatch(lkle_logit(X, y, k = k, K = K),
                        error = function(e) NULL)
        if (!is.null(fit)) {
          grid_results <- rbind(grid_results,
                                data.frame(k = k, K = K, AIC = fit$aic))
        }
      }
    }
    best <- grid_results[which.min(grid_results$AIC), ]
    return(list(best_k = best$k, best_K = best$K,
                best_value = best$AIC,
                results = grid_results))
  }

  if (criterion == "accuracy") {
    folds_idx <- sample(rep(1:folds, length.out = n))
    for (k in k_seq) {
      for (K in K_seq) {
        accs <- c()
        for (f in 1:folds) {
          X_train <- X[folds_idx != f, ]
          y_train <- y[folds_idx != f]
          X_test <- X[folds_idx == f, ]
          y_test <- y[folds_idx == f]
          fit <- tryCatch(lkle_logit(X_train, y_train, k = k, K = K),
                          error = function(e) NULL)
          if (!is.null(fit)) {
            preds <- predict(fit, X_test, type = "class")
            accs <- c(accs, mean(preds == y_test))
          }
        }
        grid_results <- rbind(grid_results,
                              data.frame(k = k, K = K,
                                         accuracy = mean(accs, na.rm = TRUE)))
      }
    }
    best <- grid_results[which.max(grid_results$accuracy), ]
    return(list(best_k = best$k, best_K = best$K,
                best_value = best$accuracy,
                results = grid_results))
  }
}

