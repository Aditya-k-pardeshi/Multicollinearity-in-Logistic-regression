#' Almost Unbiased Ridge Logistic Estimator (AURE)
#'
#' Fits an AURE penalized logistic regression model using the
#' Almost Unbiased Ridge Estimation approach. This estimator reduces
#' the bias introduced by standard Ridge Logistic Regression while
#' maintaining stability under multicollinearity.
#'
#' The model is estimated using the formula:
#' \deqn{\beta_{AURE} = (I - (1 + \lambda)(X'X + \lambda I)^{-1}) \beta_{MLE}}
#'
#' @param X Matrix or data frame of predictors.
#' @param y Binary response vector.
#' @param lambda Ridge shrinkage parameter (0 <= lambda <= 1).
#' @param tol Convergence tolerance for IRLS algorithm (default = 1e-6).
#' @param max_iter Maximum number of IRLS iterations (default = 100).
#'
#' @return A list of class "aure_logit" containing:
#' \item{coefficients}{Estimated regression coefficients.}
#' \item{fitted.values}{Predicted probabilities.}
#' \item{lambda}{Shrinkage parameter used.}
#'
#' @details
#' The Almost Unbiased Ridge Estimator (AURE) adjusts the traditional Ridge
#' estimator to achieve reduced bias while preserving variance stability.
#' It performs well in the presence of multicollinearity.
#'
#' @references
#' Kibria, B. M. G. (2003). Performance of some new ridge regression estimators.
#' *Communications in Statistics – Simulation and Computation*, 32(2), 419–435.
#'
#' @examples
#' set.seed(123)
#' X <- matrix(rnorm(100*4), ncol=4)
#' y <- rbinom(100, 1, 0.5)
#' model <- aure_logit(X, y, lambda = 0.5)
#' preds <- predict(model, X, type = "response")
#' mean(ifelse(preds > 0.5, 1, 0) == y)
#'
#' @export
aure_logit <- function(X, y, lambda = 1, standardize = FALSE, tol = 1e-6, max_iter = 50) {
  # Ensure valid lambda
  if (lambda < 0) stop("Lambda must be >= 0.")

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

  # Initialize
  beta <- rep(0, p)
  iter <- 0
  converged <- FALSE

  # IRLS loop (Ridge first)
  while (iter < max_iter) {
    eta <- X %*% beta
    mu <- 1 / (1 + exp(-eta))
    W <- diag(as.vector(mu * (1 - mu)))
    z <- eta + (y - mu) / (mu * (1 - mu))

    XWX <- t(X) %*% W %*% X
    XWz <- t(X) %*% W %*% z
    I_p <- diag(p)
    I_p[1, 1] <- 0  # don't penalize intercept

    # Ridge estimate
    beta_ridge <- solve(XWX + lambda * I_p, XWz)

    # Apply AURE correction
    beta_new <- solve(diag(p) + lambda * I_p) %*% beta_ridge

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
  fisher_info <- solve(XWX + lambda * I_p)
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
  class(result) <- "aure_logit"
  return(result)
}
#' @export
predict.aure_logit <- function(object, newdata, type = c("link", "response", "class")) {
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
