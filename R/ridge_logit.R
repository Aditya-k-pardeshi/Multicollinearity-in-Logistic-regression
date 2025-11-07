#' Ridge Logistic Regression
#'
#' Fits a ridge-penalized logistic regression model using Iteratively Reweighted Least Squares (IRLS).
#' Returns an object of class "ridge_logit" similar to glm output.
#'
#' @param X A matrix or data frame of predictors.
#' @param y A binary response vector (0 or 1).
#' @param lambda The ridge penalty parameter (default = 1).
#' @param tol Convergence tolerance (default = 1e-6).
#' @param maxit Maximum number of iterations (default = 50).
#' @return An object of class "ridge_logit" similar to glm output.
#' @export
ridge_logit <- function(X, y, lambda = 0,standardize = FALSE, tol = 1e-6, max_iter = 50) {
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


  # Initialize beta
  beta <- rep(0, p)
  iter <- 0
  converged <- FALSE

  # Case 1: lambda == 0 → use standard glm MLE
  if (lambda == 0) {
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
      lambda = 0,
      family = "binomial",
      converged = glm_fit$converged,
      X = X,
      y = y
    )
    class(out) <- "ridge_logit"
    return(out)
  }


  # IRLS loop (for logistic regression with ridge penalty)
  while (iter < max_iter) {
    eta <- X %*% beta
    mu <- 1 / (1 + exp(-eta))
    W <- diag(as.vector(mu * (1 - mu)))
    z <- eta + (y - mu) / (mu * (1 - mu))

    # Penalized weighted least squares step
    XWX <- t(X) %*% W %*% X
    XWz <- t(X) %*% W %*% z

    # Ridge penalty applied to non-intercept terms only
    P <- diag(p)
    P[1, 1] <- 0

    # Solve system
    beta_new <- solve(XWX + lambda * P, XWz)

    # Check convergence
    if (max(abs(beta_new - beta)) < tol) {
      converged <- TRUE
      beta <- beta_new
      break
    }
    beta <- beta_new
    iter <- iter + 1
  }

  # Compute fitted values and deviance
  eta <- X %*% beta
  mu <- 1 / (1 + exp(-eta))
  eps <- 1e-15
  deviance <- -2 * sum(y * log(mu + eps) + (1 - y) * log(1 - mu + eps))
  null_mu <- mean(y)
  null_deviance <- -2 * sum(y * log(null_mu + eps) + (1 - y) * log(1 - null_mu + eps))

  # Standard errors (approximate from penalized Fisher info)
  fisher_info <- solve(XWX + lambda * P)
  se <- sqrt(diag(fisher_info))

  z_values <- beta / se
  p_values <- 2 * (1 - pnorm(abs(z_values)))

  coef_table <- data.frame(
    Estimate = beta,
    Std.Error = se,
    z.value = z_values,
    Pr...z.. = p_values
  )
  rownames(coef_table) <- colnames(X)

  # AIC and DF
  df_total <- n
  df_residual <- n - p
  aic <- deviance + 2 * p

  # Return structured output
  result <- list(
    coefficients = coef_table,
    fitted.values = mu,
    deviance = deviance,
    null.deviance = null_deviance,
    df.residual = df_residual,
    df.total = df_total,
    lambda = lambda,
    family = "binomial(logit)",
    converged = converged,
    iterations = iter,
    call = match.call()
  )
  class(result) <- "ridge_logit"
  return(result)
}

# Summary method
#' @export
summary.ridge_logit <- function(object, ...) {
  cat("\nCall:\n")
  print(object$call)

  cat("\nCoefficients:\n")
  printCoefmat(object$coefficients, P.values = TRUE, has.Pvalue = TRUE)

  cat("\n(Dispersion parameter for binomial family taken as 1)\n")
  cat("\nNull Deviance:", round(object$null.deviance, 4),
      "\nResidual Deviance:", round(object$deviance, 4),
      "\nAIC:", round(object$deviance + 2 * nrow(object$coefficients), 4),
      "\nLambda:", object$lambda,
      "\nConverged:", object$converged, "\n")
}
#' @export
predict.ridge_logit <- function(object, newdata, type = c("link", "response", "class")) {
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

