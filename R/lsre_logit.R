#' Logistic Stein Ridge Estimator (LSRE)
#'
#' Implements the Logistic Stein Ridge Estimator (LSRE):
#' \deqn{β̂_LSRE = c_r A β̂_MLE}
#' where
#' \deqn{A = (X'WX + kI)^{-1} X'WX}
#' and
#' \deqn{c_r = (β̂_MLE'A'Aβ̂_MLE) / (β̂_MLE'A'Aβ̂_MLE + tr(MSE * (A(X'WX)^{-1}A')))}
#'
#' @param X Matrix or data frame of predictors (without intercept).
#' @param y Binary response vector (0/1).
#' @param k Ridge-type shrinkage parameter (k >= 0).
#' @return An object of class "lsre_logit" similar to glm.
#' @export
lsre_logit <- function(X, y, k = 1) {
  # --- Ensure data shape ---
  X <- as.matrix(X)
  y <- as.numeric(y)
  if (!all(X[, 1] == 1)) {
    X <- cbind(1, X)
  }
  colnames(X) <- c("(Intercept)", paste0("X", 1:(ncol(X) - 1)))
  n <- nrow(X)
  p <- ncol(X)

  # --- Step 1: Fit MLE ---
  glm_fit <- glm(y ~ X[, -1], family = binomial)
  beta_mle <- coef(glm_fit)
  mu <- glm_fit$fitted.values
  W <- diag(as.vector(mu * (1 - mu)))
  XWX <- t(X) %*% W %*% X
  I <- diag(p)
  I[1, 1] <- 0  # no penalty on intercept

  # --- Step 2: Compute A matrix ---
  A <- solve(XWX + k * I) %*% XWX

  # --- Step 3: Estimate MSE from MLE residuals ---
  eta <- X %*% beta_mle
  mu <- 1 / (1 + exp(-eta))
  residuals <- y - mu
  MSE <- sum(residuals^2) / (n - p)

  # --- Step 4: Compute shrinkage factor c_r ---
  num <- as.numeric(t(beta_mle) %*% t(A) %*% A %*% beta_mle)
  den <- num + sum(diag(MSE * (A %*% solve(XWX) %*% t(A))))
  c_r <- num / den
  c_r <- max(min(c_r, 1), 0)  # ensure 0 < c_r < 1

  # --- Step 5: Compute LSRE coefficients ---
  beta_lsre <- c_r * (A %*% beta_mle)

  # --- Step 6: Fitted values & deviance ---
  eta <- X %*% beta_lsre
  mu <- 1 / (1 + exp(-eta))
  eps <- 1e-15
  deviance <- -2 * sum(y * log(mu + eps) + (1 - y) * log(1 - mu + eps))
  null_mu <- mean(y)
  null_dev <- -2 * sum(y * log(null_mu + eps) + (1 - y) * log(1 - null_mu + eps))
  aic <- deviance + 2 * p

  # --- Step 7: Standard errors (approximate) ---
  fisher_info <- solve(XWX + k * I)
  se <- sqrt(diag(fisher_info))
  z_values <- beta_lsre / se
  p_values <- 2 * (1 - pnorm(abs(z_values)))

  coef_table <- data.frame(
    Estimate = beta_lsre,
    Std.Error = se,
    z.value = z_values,
    Pr...z.. = p_values
  )
  rownames(coef_table) <- colnames(X)

  # --- Step 8: Return structured result ---
  result <- list(
    coefficients = coef_table,
    fitted.values = mu,
    deviance = deviance,
    null.deviance = null_dev,
    aic = aic,
    k = k,
    c_r = c_r,
    df.residual = n - p,
    df.total = n,
    family = "binomial(logit)",
    converged = glm_fit$converged,
    call = match.call()
  )
  class(result) <- "lsre_logit"
  return(result)
}

# Summary method
#' @export
summary.lsre_logit <- function(object, ...) {
  cat("\nCall:\n")
  print(object$call)
  cat("\nShrinkage factor (c_r):", round(object$c_r, 6), "\n")
  cat("\nCoefficients:\n")
  printCoefmat(object$coefficients, P.values = TRUE, has.Pvalue = TRUE)
  cat("\nNull Deviance:", round(object$null.deviance, 4),
      "\nResidual Deviance:", round(object$deviance, 4),
      "\nAIC:", round(object$aic, 4),
      "\nLambda (k):", object$k, "\n")
}

# Predict method
#' @export
predict.lsre_logit <- function(object, newdata, type = c("link", "response", "class"), threshold = 0.5, ...) {
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
#' Tune Ridge Parameter for Logistic Stein Ridge Estimator (LSRE)
#'
#' Searches for the optimal ridge parameter k that minimizes AIC
#' or maximizes prediction accuracy for the Logistic Stein Ridge Estimator.
#'
#' @param X Matrix or data frame of predictors.
#' @param y Binary response vector (0/1).
#' @param k_grid Sequence of k values to test (default: seq(0, 5, 0.5)).
#' @param metric Optimization metric: "AIC" (default) or "accuracy".
#' @param verbose Logical; if TRUE, prints progress for each k.
#' @return A list containing the best k, best model, best score, and all results.
#' @export
tune_lsre <- function(X, y,
                      k_grid = seq(0, 5, 0.5),
                      metric = c("AIC", "accuracy"),
                      verbose = TRUE) {
  metric <- match.arg(metric)
  results <- data.frame()

  best_score <- ifelse(metric == "AIC", Inf, -Inf)
  best_model <- NULL
  best_k <- NA

  # Iterate over grid of k values
  for (k in k_grid) {
    fit <- tryCatch(lsre_logit(X, y, k = k), error = function(e) NULL)
    if (!is.null(fit)) {
      if (metric == "AIC") {
        score <- fit$aic
        better <- score < best_score
      } else {
        pred <- predict(fit, X, type = "class")
        score <- mean(pred == y)
        better <- score > best_score
      }

      results <- rbind(results, data.frame(k = k, score = score))

      if (better) {
        best_score <- score
        best_model <- fit
        best_k <- k
      }

      if (verbose) message(sprintf("Checked: k=%.2f, %s=%.4f", k, metric, score))
    }
  }

  out <- list(
    best_k = best_k,
    best_score = best_score,
    best_model = best_model,
    results = results,
    metric = metric
  )
  class(out) <- "tuned_lsre"
  return(out)
}

#' Summary Method for tuned LSRE
#' @export
summary.tuned_lsre <- function(object, ...) {
  cat("\nTuning Summary for Logistic Stein Ridge Estimator (LSRE)\n")
  cat("--------------------------------------------------------\n")
  cat("Optimization Metric:", object$metric, "\n")
  cat("Best k:", object$best_k, "\n")
  cat("Best Score:", round(object$best_score, 4), "\n")
  cat("--------------------------------------------------------\n")
  cat("\nFinal Model Summary:\n")
  summary(object$best_model)
}
