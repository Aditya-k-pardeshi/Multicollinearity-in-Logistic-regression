#' Logistic Stein K-Liu Estimator (LSKE)
#'
#' Implements the Logistic Stein K-Liu Estimator (LSKE) as defined in:
#' \deqn{β̂_LSKE = c_k K β̂_MLE}
#' where
#' \deqn{K = (X'WX + kI)^{-1}(X'WX - kI)}
#' and
#' \deqn{c_k = (β̂_MLE'K'β̂_MLE) / (β̂_MLE'K'Kβ̂_MLE + tr(MSE * (K(X'WX)^{-1}K')))}
#'
#' @param X Matrix or data frame of predictors (without intercept).
#' @param y Binary response vector (0/1).
#' @param k Biasing parameter (k >= 0).
#' @return An object of class "lske_logit" similar to glm output.
#' @export
lske_logit <- function(X, y, k = 1) {
  # --- Ensure X has intercept ---
  X <- as.matrix(X)
  y <- as.numeric(y)
  if (!all(X[, 1] == 1)) X <- cbind(1, X)
  colnames(X) <- c("(Intercept)", paste0("X", 1:(ncol(X) - 1)))
  n <- nrow(X)
  p <- ncol(X)

  # --- Step 1: Fit standard MLE ---
  glm_fit <- glm(y ~ X[, -1], family = binomial)
  beta_mle <- coef(glm_fit)
  mu <- glm_fit$fitted.values
  W <- diag(as.vector(mu * (1 - mu)))
  XWX <- t(X) %*% W %*% X
  I <- diag(p)
  I[1, 1] <- 0  # no penalty on intercept

  # --- Step 2: Compute K matrix ---
  K_mat <- solve(XWX + k * I) %*% (XWX - k * I)

  # --- Step 3: Estimate MSE from MLE residuals ---
  eta <- X %*% beta_mle
  mu <- 1 / (1 + exp(-eta))
  residuals <- y - mu
  MSE <- sum(residuals^2) / (n - p)

  # --- Step 4: Compute shrinkage factor c_k ---
  num <- as.numeric(t(beta_mle) %*% t(K_mat) %*% beta_mle)
  den <- as.numeric(t(beta_mle) %*% t(K_mat) %*% K_mat %*% beta_mle) +
    sum(diag(MSE * (K_mat %*% solve(XWX) %*% t(K_mat))))
  c_k <- num / den
  c_k <- max(min(c_k, 1), 0)  # constrain to [0,1]

  # --- Step 5: Compute LSKE coefficients ---
  beta_lske <- c_k * (K_mat %*% beta_mle)

  # --- Step 6: Compute fitted values & deviance ---
  eta <- X %*% beta_lske
  mu <- 1 / (1 + exp(-eta))
  eps <- 1e-15
  deviance <- -2 * sum(y * log(mu + eps) + (1 - y) * log(1 - mu + eps))
  null_mu <- mean(y)
  null_dev <- -2 * sum(y * log(null_mu + eps) + (1 - y) * log(1 - null_mu + eps))
  aic <- deviance + 2 * p

  # --- Step 7: Standard errors (approximate) ---
  fisher_info <- solve(XWX + k * I)
  se <- sqrt(diag(fisher_info))
  z_values <- beta_lske / se
  p_values <- 2 * (1 - pnorm(abs(z_values)))

  coef_table <- data.frame(
    Estimate = beta_lske,
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
    c_k = c_k,
    df.residual = n - p,
    df.total = n,
    family = "binomial(logit)",
    converged = glm_fit$converged,
    call = match.call()
  )
  class(result) <- "lske_logit"
  return(result)
}

# --- Summary Method ---
#' @export
summary.lske_logit <- function(object, ...) {
  cat("\nCall:\n")
  print(object$call)
  cat("\nShrinkage factor (c_k):", round(object$c_k, 6), "\n")
  cat("\nCoefficients:\n")
  printCoefmat(object$coefficients, P.values = TRUE, has.Pvalue = TRUE)
  cat("\nNull Deviance:", round(object$null.deviance, 4),
      "\nResidual Deviance:", round(object$deviance, 4),
      "\nAIC:", round(object$aic, 4),
      "\nLambda (k):", object$k,
      "\nConverged:", object$converged, "\n")
}

# --- Predict Method ---
#' @export
predict.lske_logit <- function(object, newdata, type = c("link", "response", "class"), threshold = 0.5, ...) {
  type <- match.arg(type)
  X <- as.matrix(newdata)
  if (!all(X[, 1] == 1)) X <- cbind(1, X)
  beta <- object$coefficients$Estimate
  eta <- X %*% beta
  mu <- 1 / (1 + exp(-eta))
  if (type == "link") return(drop(eta))
  if (type == "response") return(drop(mu))
  if (type == "class") return(ifelse(mu >= threshold, 1, 0))
}

#' Tune Logistic Stein K-Liu Estimator (LSKE)
#'
#' Tunes the biasing parameter `k` in the Logistic Stein K-Liu Estimator (LSKE)
#' by searching over a grid and selecting the model with minimum AIC (default)
#' or highest accuracy if response `y` is given.
#'
#' @param X Matrix or data frame of predictors (without intercept)
#' @param y Binary response vector (0/1)
#' @param k_seq Sequence of k values to try (default: seq(0, 5, by = 0.1))
#' @param criterion Model selection criterion: "AIC" (default) or "accuracy"
#' @param verbose Logical; if TRUE, prints progress for each k.
#' @return A list containing:
#'   - `best_model`: LSKE model with optimal k
#'   - `best_k`: selected k
#'   - `results`: data.frame of k and performance metric
#' @export
tune_lske <- function(X, y, k_seq = seq(0, 5, by = 0.1),
                      criterion = c("AIC", "accuracy"),
                      verbose = TRUE) {
  criterion <- match.arg(criterion)
  results <- data.frame(k = numeric(0), metric = numeric(0))

  for (k in k_seq) {
    fit <- tryCatch(lske_logit(X, y, k = k), error = function(e) NULL)
    if (is.null(fit)) next

    if (criterion == "AIC") {
      metric <- fit$aic
    } else {
      preds <- predict.lske_logit(fit, X, type = "class")
      metric <- mean(preds == y)
    }

    if (verbose)
      cat("k =", k, "→", criterion, "=", round(metric, 4), "\n")

    results <- rbind(results, data.frame(k = k, metric = metric))
  }

  if (criterion == "AIC") {
    best_k <- results$k[which.min(results$metric)]
  } else {
    best_k <- results$k[which.max(results$metric)]
  }

  best_model <- lske_logit(X, y, k = best_k)

  return(list(best_model = best_model, best_k = best_k, results = results))
}
