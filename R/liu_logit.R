#' Liu Logistic Regression
#'
#' Fits a Liu-penalized logistic regression model.
#' Returns an object of class "liu_logit" similar to glm output.
#'
#' @param X Matrix or data frame of predictors (including intercept or not).
#' @param y Binary response vector (0/1).
#' @param d Liu biasing parameter (0 <= d < 1).
#' @param tol Convergence tolerance (default = 1e-6).
#' @param max_iter Maximum number of iterations (default = 50).
#' @return Object of class "liu_logit" containing coefficients, deviance, etc.
#' @export
liu_logit <- function(X, y, d = 0, standardize = FALSE, tol = 1e-6, max_iter = 50) {
  # Check d is within [0,1]
  if (d < 0 || d > 1) {
    stop("Error: Liu biasing parameter 'd' must satisfy 0 <= d <= 1")
  }

  # Ensure X has intercept
  if (!all(X[, 1] == 1)) {
    X <- cbind(1, X)
  }
  colnames(X) <- c("(Intercept)", paste0("X", 1:(ncol(X) - 1)))
  n <- nrow(X)
  p <- ncol(X)

  if (standardize) {
    X_mean <- apply(X[, -1], 2, mean)
    X_sd <- apply(X[, -1], 2, sd)
    X[, -1] <- scale(X[, -1])
  }

  # Case 1: d == 0 → use standard glm MLE
  if (d == 0) {
    glm_fit <- glm(y ~ X[, -1], family = binomial)
    coef_table <- summary(glm_fit)$coefficients

    out <- list(
      coefficients = coef_table,
      fitted.values = fitted(glm_fit),
      null.deviance = glm_fit$null.deviance,
      residual.deviance = glm_fit$deviance,
      aic = glm_fit$aic,
      df.residual = glm_fit$df.residual,
      df.total = length(y),
      d = 0,
      family = "binomial",
      converged = glm_fit$converged,
      X = X,
      y = y
    )
    class(out) <- "liu_logit"
    return(out)
  }

  # Step 1: Obtain MLE via glm
  glm_fit <- glm(y ~ X[, -1], family = binomial)
  beta_mle <- coef(glm_fit)
  eta <- glm_fit$linear.predictors
  mu <- glm_fit$fitted.values
  W <- diag(as.vector(mu * (1 - mu)))

  # Step 2: Compute Liu estimator
  XWX <- t(X) %*% W %*% X
  I <- diag(p)
  I[1, 1] <- 0  # don't penalize intercept

  beta_liu <- solve(XWX + I) %*% (XWX + d * I) %*% beta_mle

  # Step 3: Compute fitted values & deviance
  eta_liu <- X %*% beta_liu
  mu_liu <- 1 / (1 + exp(-eta_liu))
  eps <- 1e-15
  deviance <- -2 * sum(y * log(mu_liu + eps) + (1 - y) * log(1 - mu_liu + eps))
  null_mu <- mean(y)
  null_deviance <- -2 * sum(y * log(null_mu + eps) + (1 - y) * log(1 - null_mu + eps))

  # Step 4: Standard errors (approximate)
  fisher_info <- solve(XWX + I) %*% (XWX + d * I) %*% solve(XWX + I)
  se <- sqrt(diag(fisher_info))
  z_values <- beta_liu / se
  p_values <- 2 * (1 - pnorm(abs(z_values)))

  coef_table <- data.frame(
    Estimate = drop(beta_liu),
    Std.Error = se,
    z.value = z_values,
    Pr...z.. = p_values
  )
  rownames(coef_table) <- colnames(X)

  # Step 5: Collect output
  out <- list(
    coefficients = coef_table,
    fitted.values = mu_liu,
    null.deviance = null_deviance,
    deviance = deviance,
    aic = deviance + 2 * p,
    df.residual = n - p,
    df.total = n,
    d = d,
    family = "binomial(logit)",
    converged = TRUE,
    call = match.call()
  )

  class(out) <- "liu_logit"
  return(out)
}

#' Summary Method for Liu Logistic Regression
#' @export
summary.liu_logit <- function(object, ...) {
  cat("\nCall:\n")
  print(object$call)
  cat("\nCoefficients:\n")
  printCoefmat(object$coefficients, P.values = TRUE, has.Pvalue = TRUE)

  cat("\n(Dispersion parameter for binomial family taken as 1)\n")
  cat("\nNull Deviance:", round(object$null.deviance, 4),
      "\nResidual Deviance:", round(object$deviance, 4),
      "\nAIC:", round(object$aic, 4),
      "\nBiasing Parameter (d):", object$d,
      "\nConverged:", object$converged, "\n")
}
#' @export
predict.liu_logit <- function(object, newdata, type = c("link", "response", "class")) {
  type <- match.arg(type)
  X <- as.matrix(newdata)
  if (!all(X[, 1] == 1)) {
    X <- cbind(1, X)
  }

  eta <- X %*% object$coefficients$Estimate
  mu <- 1 / (1 + exp(-eta))

  if (type == "link") return(eta)
  if (type == "response") return(mu)
  if (type == "class") return(ifelse(mu >= 0.5, 1, 0))
}
