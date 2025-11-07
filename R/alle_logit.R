#' Almost Unbiased Liu Logistic Estimator (ALLE)
#'
#' Fits an Almost Unbiased Liu Logistic Regression model using IRLS.
#' This estimator reduces the bias of the standard Liu logistic estimator.
#'
#' The estimator is defined as:
#' \deqn{\hat{\beta}_{ALLE} = [I - (1 + d)(X'WX + I)^{-1}] \hat{\beta}_{MLE}}
#'
#' @param X Matrix or data frame of predictors.
#' @param y Binary response vector (0/1).
#' @param d Liu shrinkage parameter, numeric in [0, 1].
#' @param tol Convergence tolerance for IRLS (default = 1e-6).
#' @param max_iter Maximum number of IRLS iterations (default = 50).
#' @return An object of class "alle_logit" similar to glm output.
#' @export
alle_logit <- function(X, y, d = 0.5,  standardize = FALSE, tol = 1e-6, max_iter = 50) {
  if (d < 0 || d > 1) stop("Parameter d must satisfy 0 <= d <= 1")

  # Ensure intercept
  if (!all(X[, 1] == 1)) {
    X <- cbind(1, X)
  }
  colnames(X) <- c("(Intercept)", paste0("X", 1:(ncol(X) - 1)))
  n <- nrow(X)
  p <- ncol(X)

  # Optional standardization
  if (standardize) {
    X_mean <- apply(X[, -1], 2, mean)
    X_sd <- apply(X[, -1], 2, sd)
    X[, -1] <- scale(X[, -1])
  }

  # Step 1: Fit MLE logistic model
  glm_fit <- glm(y ~ X[, -1], family = binomial)
  beta_mle <- coef(glm_fit)

  # Step 2: Compute Fisher information matrix X'WX
  eta <- X %*% beta_mle
  mu <- 1 / (1 + exp(-eta))
  W <- diag(as.vector(mu * (1 - mu)))
  XWX <- t(X) %*% W %*% X

  # Step 3: Apply ALLE formula
  beta_alle <- (diag(p) - (1 + d) * solve(XWX + diag(p))) %*% beta_mle

  # Step 4: Compute fitted values and deviance
  eta <- X %*% beta_alle
  mu <- 1 / (1 + exp(-eta))
  eps <- 1e-15
  deviance <- -2 * sum(y * log(mu + eps) + (1 - y) * log(1 - mu + eps))
  null_mu <- mean(y)
  null_deviance <- -2 * sum(y * log(null_mu + eps) + (1 - y) * log(1 - null_mu + eps))

  # Step 5: Standard errors and inference
  fisher_info <- solve(XWX)
  se <- sqrt(diag(fisher_info))
  z_values <- beta_alle / se
  p_values <- 2 * (1 - pnorm(abs(z_values)))

  coef_table <- data.frame(
    Estimate = beta_alle,
    Std.Error = se,
    z.value = z_values,
    Pr...z.. = p_values
  )
  rownames(coef_table) <- colnames(X)

  # Step 6: Return results
  result <- list(
    coefficients = coef_table,
    fitted.values = mu,
    deviance = deviance,
    null.deviance = null_deviance,
    df.residual = n - p,
    df.total = n,
    d = d,
    family = "binomial(logit)",
    call = match.call()
  )
  class(result) <- "alle_logit"
  return(result)
}

# Summary method
#' @export
summary.alle_logit <- function(object, ...) {
  cat("\nCall:\n")
  print(object$call)

  cat("\nCoefficients:\n")
  printCoefmat(object$coefficients, P.values = TRUE, has.Pvalue = TRUE)

  cat("\nNull Deviance:", round(object$null.deviance, 4),
      "\nResidual Deviance:", round(object$deviance, 4),
      "\nShrinkage d:", object$d, "\n")
}

#' Predict method for Almost Unbiased Liu Logistic Estimator (ALLE)
#'
#' @param object An object of class "alle_logit" from alle_logit().
#' @param newdata A matrix or data frame of new predictor values.
#' @param type Type of prediction: "link" (linear predictor),
#'        "response" (probabilities), or "class" (0/1).
#' @param threshold Classification cutoff for type = "class" (default = 0.5).
#' @param ... Additional arguments (ignored).
#'
#' @return A numeric vector of predictions.
#' @export
predict.alle_logit <- function(object, newdata,
                               type = c("link", "response", "class"),
                               threshold = 0.5, ...) {
  type <- match.arg(type)

  # Ensure input is matrix
  X <- as.matrix(newdata)

  # Add intercept if not present
  if (!all(X[, 1] == 1)) {
    X <- cbind(1, X)
  }

  # Compute linear predictor
  beta <- object$coefficients$Estimate
  eta <- X %*% beta
  mu <- 1 / (1 + exp(-eta))

  # Return based on type
  if (type == "link") return(drop(eta))
  if (type == "response") return(drop(mu))
  if (type == "class") return(ifelse(mu >= threshold, 1, 0))
}

#' Tune Almost Unbiased Liu Logistic Estimator (ALLE)
#'
#' Tunes the biasing parameter `d` in the ALLE model by evaluating
#' a sequence of `d` values and selecting the model with minimum AIC
#' or highest accuracy.
#'
#' @param X Matrix or data frame of predictors (without intercept)
#' @param y Binary response vector (0/1)
#' @param d_seq Sequence of d values to try (default: seq(0, 1, by = 0.05))
#' @param criterion Model selection criterion: "AIC" (default) or "accuracy"
#' @param verbose Logical; if TRUE, prints progress
#' @return A list containing:
#'   - `best_model`: ALLE model with optimal d
#'   - `best_d`: selected d value
#'   - `results`: data.frame of d values and corresponding metric
#' @export
tune_alle <- function(X, y,
                      d_seq = seq(0, 1, by = 0.05),
                      criterion = c("AIC", "accuracy"),
                      verbose = TRUE) {

  criterion <- match.arg(criterion)
  results <- data.frame(d = numeric(0), metric = numeric(0))

  for (d in d_seq) {
    # Try fitting the model
    fit <- tryCatch(alle_logit(X, y, d = d), error = function(e) NULL)

    # Skip failed fits
    if (is.null(fit)) {
      if (verbose) cat("d =", round(d, 3), "→ Model failed. Skipping.\n")
      next
    }

    # Safely compute metric
    metric <- tryCatch({
      if (criterion == "AIC") {
        if (is.null(fit$aic)) NA else as.numeric(fit$aic)
      } else {
        preds <- tryCatch(predict.alle_logit(fit, X, type = "class"), error = function(e) NULL)
        if (is.null(preds)) NA else mean(preds == y)
      }
    }, error = function(e) NA)

    # Skip invalid metrics (NULL, NA, or non-numeric)
    if (is.null(metric) || is.na(metric) || !is.numeric(metric)) {
      if (verbose) cat("d =", round(d, 3), "→ Invalid metric. Skipping.\n")
      next
    }

    if (verbose)
      cat("d =", round(d, 3), "→", criterion, "=", round(metric, 4), "\n")

    results <- rbind(results, data.frame(d = d, metric = metric))
  }

  if (nrow(results) == 0)
    stop("No valid models were fitted. Please check your alle_logit() implementation.")

  # Choose best d based on criterion
  if (criterion == "AIC") {
    best_d <- results$d[which.min(results$metric)]
  } else {
    best_d <- results$d[which.max(results$metric)]
  }

  best_model <- alle_logit(X, y, d = best_d)

  return(list(best_model = best_model, best_d = best_d, results = results))
}


