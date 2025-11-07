#' Almost Unbiased Logistic Estimator (AULE)
#'
#' Fits an AULE penalized logistic regression model using Iteratively Reweighted Least Squares (IRLS).
#' The estimator is defined as:
#' \deqn{\hat{\beta}_{AULE} = (I - d^2 (X'WX + I)^{-2}) \hat{\beta}_{MLE}},
#' where \eqn{0 \le d \le 1} is the shrinkage parameter.
#'
#' @param X Matrix or data frame of predictors (including intercept column if desired)
#' @param y Binary response vector (0/1)
#' @param d Shrinkage parameter, numeric in [0, 1]
#' @param tol Convergence tolerance for IRLS (default 1e-6)
#' @param max_iter Maximum number of IRLS iterations (default 50)
#' @return An object of class "aule_logit", similar to `glm` output
#' @export
aule_logit <- function(X, y, d = 0.5, standardize = FALSE, tol = 1e-6, max_iter = 50) {
  # Check bounds of d
  if (d < 0 || d > 1) stop("Parameter 'd' must be between 0 and 1.")

  # Ensure intercept
  if (!all(X[, 1] == 1)) {
    X <- cbind(1, X)
  }
  colnames(X) <- c("(Intercept)", paste0("X", 1:(ncol(X) - 1)))
  n <- nrow(X)
  p <- ncol(X)

  # Standardize (optional)
  if (standardize) {
    X_mean <- apply(X[, -1], 2, mean)
    X_sd <- apply(X[, -1], 2, sd)
    X[, -1] <- scale(X[, -1])
  }

  # Initialize
  beta <- rep(0, p)
  iter <- 0
  converged <- FALSE

  # Step 1: get MLE estimate (for Liu-type adjustment)
  glm_fit <- glm(y ~ X[, -1], family = binomial)
  beta_mle <- coef(glm_fit)

  # IRLS loop
  while (iter < max_iter) {
    eta <- X %*% beta
    mu <- 1 / (1 + exp(-eta))
    W <- diag(as.vector(mu * (1 - mu)))
    z <- eta + (y - mu) / (mu * (1 - mu))

    XWX <- t(X) %*% W %*% X
    XWz <- t(X) %*% W %*% z
    I_p <- diag(p)
    I_p[1, 1] <- 0  # don't penalize intercept

    # AULE formula (Liu-type adjusted)
    beta_new <- solve(XWX + I_p, (XWX + d * I_p) %*% beta_mle)

    if (max(abs(beta_new - beta)) < tol) {
      converged <- TRUE
      beta <- beta_new
      break
    }

    beta <- beta_new
    iter <- iter + 1
  }

  # Fitted values & Deviances
  eta <- X %*% beta
  mu <- 1 / (1 + exp(-eta))
  eps <- 1e-15
  deviance <- -2 * sum(y * log(mu + eps) + (1 - y) * log(1 - mu + eps))
  null_mu <- mean(y)
  null_deviance <- -2 * sum(y * log(null_mu + eps) + (1 - y) * log(1 - null_mu + eps))

  # Standard errors (approx)
  fisher_info <- solve(XWX + I_p)
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

  aic <- deviance + 2 * p

  # Return
  result <- list(
    coefficients = coef_table,
    fitted.values = mu,
    deviance = deviance,
    null.deviance = null_deviance,
    df.residual = n - p,
    df.total = n,
    d = d,
    family = "binomial(logit)",
    converged = converged,
    iterations = iter,
    call = match.call()
  )
  class(result) <- "aule_logit"
  return(result)
}
#' @export
predict.aule_logit <- function(object, newdata, type = c("link", "response", "class")) {
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
