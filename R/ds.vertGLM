ds.vertGLM <- function(Xa, Xb, y, family = gaussian(), max_iter = 100, tol = 1e-6, scale_data = TRUE, partition_size = NULL) {
  n <- length(y)
  p_a <- ncol(Xa)
  p_b <- ncol(Xb)

  if (scale_data) {
    Xa <- scale(Xa)
    Xb <- scale(Xb)
    if (family$family == "gaussian") {
      y <- scale(y)
    }
  }

  beta_a <- rep(0, p_a)
  beta_b <- rep(0, p_b)

  #eta <- rep(0, n)
  eta_a <- rep(0, p_a)
  eta_b <- rep(0, p_b)

  for (iter in 1:max_iter) {
    beta_a_old <- beta_a
    beta_b_old <- beta_b

    #update beta_a
    results_a <- vertGLMDS(Xa, y, eta_b, family, beta_a, regularization = 1e-4)
    beta_a <- results_a$beta

    eta_a <- results_a$eta

    #update beta_b
    results_b <- vertGLMDS(Xb, y, eta_a, family, beta_b, regularization = 1e-4)
    beta_b <- results_b$beta

    eta_b <- results_b$eta

    #check convergence criteria
    diff_a <- sum(abs(beta_a - beta_a_old))
    diff_b <- sum(abs(beta_b - beta_b_old))
    if (diff_a < tol && diff_b < tol) {
      message(paste("Converged after", iter, "iterations"))
      break
    }

    #sometimes I would get increasingly divergent results with ballooning estimates
    if (any(abs(beta_a) > 1e8) || any(abs(beta_b) > 1e8)) {
      message("Estimates became too large. Stopping early.")
      break
    }
  }

  return(list(beta_a = beta_a, beta_b = beta_b, iterations = iter))
}
