#' Logistic Modified Ridge Estimator (LMRT)
#'
#' Fits the Logistic Modified Ridge Estimator (LMRT) using:
#' \deqn{\hat{\beta}_{MRT} = (X'WX + k(1 + d)I)^{-1}(X'WX)\hat{\beta}_{MLE}}
#'
#' @param X Matrix or data frame of predictors (without intercept)
#' @param y Binary response vector (0 or 1)
#' @param k Ridge-type biasing parameter (k >= 0)
#' @param d Additional bias parameter (0 <= d <= 1)
#' @return An object of class "lmrt_logit"
#' @export
lmrt_logit <- function(X, y, k = 1, d = 0.5) {
  # --- Check parameters ---
  if (k < 0) stop("k must be non-negative")
  if (d < 0 || d > 1) stop("d must be between 0 and 1")

  # --- Ensure intercept ---
  if (!all(X[, 1] == 1)) {
    X <- cbind(1, X)
  }
  colnames(X) <- c("(Intercept)", paste0("X", 1:(ncol(X) - 1)))
  n <- nrow(X)
  p <- ncol(X)

  # --- Step 1: Fit standard MLE ---
  glm_fit <- glm(y ~ X[, -1], family = binomial)
  beta_mle <- coef(glm_fit)
  W <- diag(as.vector(glm_fit$fit * (1 - glm_fit$fit)))

  # --- Step 2: Compute LMRT coefficients ---
  XWX <- t(X) %*% W %*% X
  I <- diag(p)
  I[1, 1] <- 0  # No penalty on intercept

  beta_lmrt <- solve(XWX + k * (1 + d) * I, XWX %*% beta_mle)

  # --- Step 3: Fitted values and deviance ---
  eta <- X %*% beta_lmrt
  mu <- 1 / (1 + exp(-eta))
  eps <- 1e-15
  deviance <- -2 * sum(y * log(mu + eps) + (1 - y) * log(1 - mu + eps))
  null_mu <- mean(y)
  null_dev <- -2 * sum(y * log(null_mu + eps) + (1 - y) * log(1 - null_mu + eps))
  aic <- deviance + 2 * p

  # --- Step 4: Standard errors ---
  fisher_info <- solve(XWX + k * (1 + d) * I)
  se <- sqrt(diag(fisher_info))
  z_val <- beta_lmrt / se
  p_val <- 2 * (1 - pnorm(abs(z_val)))

  coef_table <- data.frame(
    Estimate = beta_lmrt,
    Std.Error = se,
    z.value = z_val,
    Pr...z.. = p_val
  )
  rownames(coef_table) <- colnames(X)

  # --- Step 5: Return structured result ---
  result <- list(
    coefficients = coef_table,
    fitted.values = mu,
    deviance = deviance,
    null.deviance = null_dev,
    aic = aic,
    k = k,
    d = d,
    df.residual = n - p,
    df.total = n,
    family = "binomial(logit)",
    converged = glm_fit$converged,
    call = match.call()
  )
  class(result) <- "lmrt_logit"
  return(result)
}

# --- Summary method ---
#' @export
summary.lmrt_logit <- function(object, ...) {
  cat("\nCall:\n")
  print(object$call)
  cat("\nCoefficients:\n")
  printCoefmat(object$coefficients, P.values = TRUE, has.Pvalue = TRUE)
  cat("\nNull Deviance:", round(object$null.deviance, 4),
      "\nResidual Deviance:", round(object$deviance, 4),
      "\nAIC:", round(object$aic, 4),
      "\nk:", object$k,
      "\nd:", object$d,
      "\nConverged:", object$converged, "\n")
}

# --- Predict method ---
#' @export
predict.lmrt_logit <- function(object, newdata,
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

#' Tune Parameters for Logistic Modified Ridge Estimator (LMRT)
#'
#' Automatically searches over grids of k and d values to find the combination
#' that minimizes AIC (or maximizes classification accuracy).
#'
#' @param X Matrix or data frame of predictors.
#' @param y Binary response vector (0/1).
#' @param k_grid Sequence of k values to test (default: seq(0, 5, 0.5)).
#' @param d_grid Sequence of d values to test (default: seq(0, 1, 0.1)).
#' @param metric Optimization metric: "AIC" (default) or "accuracy".
#' @param verbose Logical; if TRUE, prints progress.
#' @return A list with best parameters, tuned model, and performance table.
#' @export
tune_lmrt <- function(X, y,
                      k_grid = seq(0, 5, 0.5),
                      d_grid = seq(0, 1, 0.1),
                      metric = c("AIC", "accuracy"),
                      verbose = TRUE) {
  metric <- match.arg(metric)
  results <- data.frame()

  best_score <- ifelse(metric == "AIC", Inf, -Inf)
  best_model <- NULL
  best_k <- NA
  best_d <- NA

  # Loop over all parameter combinations
  for (k in k_grid) {
    for (d in d_grid) {
      fit <- tryCatch(
        lmrt_logit(X, y, k = k, d = d),
        error = function(e) NULL
      )
      if (!is.null(fit)) {
        if (metric == "AIC") {
          score <- fit$aic
          better <- score < best_score
        } else {
          pred <- predict(fit, X, type = "class")
          score <- mean(pred == y)
          better <- score > best_score
        }

        results <- rbind(results, data.frame(k = k, d = d, score = score))

        if (better) {
          best_score <- score
          best_model <- fit
          best_k <- k
          best_d <- d
        }
      }
      if (verbose) message(sprintf("Checked: k=%.2f, d=%.2f, %s=%.4f", k, d, metric, score))
    }
  }

  best <- list(
    best_k = best_k,
    best_d = best_d,
    best_score = best_score,
    best_model = best_model,
    results = results,
    metric = metric
  )
  class(best) <- "tuned_lmrt"
  return(best)
}

# Summary method for tuned LMRT
#' @export
summary.tuned_lmrt <- function(object, ...) {
  cat("\nTuning Summary for Logistic Modified Ridge Estimator (LMRT)\n")
  cat("------------------------------------------------------------\n")
  cat("Optimization Metric:", object$metric, "\n")
  cat("Best k:", object$best_k, "\n")
  cat("Best d:", object$best_d, "\n")
  cat("Best Score:", round(object$best_score, 4), "\n")
  cat("------------------------------------------------------------\n")
  cat("\nFinal Model Summary:\n")
  summary(object$best_model)
}
