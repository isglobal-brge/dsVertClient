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

vertGLMDS <- function(X_partition, y, eta_partition, family, beta_partition, regularization = 1e-4) {
  n <- length(y)

  eta <- eta_partition + X_partition %*% beta_partition

  #change mean response based on family
  if (family$family == "gaussian") {
    mu <- eta
    W <- diag(1, n)
    z <- y
  } else if (family$family == "binomial") {
    mu <- 1 / (1 + exp(-eta))
    W <- diag(as.vector(mu * (1 - mu)))
    z <- eta + (y - mu) / (mu * (1 - mu))
  } else if (family$family == "poisson") {
    mu <- exp(eta)
    W <- diag(as.vector(mu))
    z <- eta + (y - mu) / mu
  } else {
    stop("Unsupported family")
  }

  #get new beta estimates regarding regularization
  XtWX <- t(X_partition) %*% W %*% X_partition
  XtWX <- XtWX + diag(regularization, nrow(XtWX))

  beta_update <- solve(XtWX, t(X_partition) %*% (W %*% (z - eta)))

  #again ran into large updates, so Im trying to prevent this. Maybe it is bad for certain data
  if (any(abs(beta_update) > 1e2)) {
    beta_update <- beta_update / max(abs(beta_update)) * 1e2
    message("Large beta update detected.")
  }

  beta_partition <- beta_partition + beta_update

  eta <- X_partition %*% beta_partition

  return(list(beta = beta_partition, eta = eta))
}


simulate_data <- function(n, p_a, p_b, family = "gaussian") {
  Xa <- matrix(rnorm(n * p_a), n, p_a)
  Xb <- matrix(rnorm(n * p_b), n, p_b)
  beta_true <- c(rep(1, p_a), rep(1, p_b))

  if (family == "gaussian") {
    y <- cbind(Xa, Xb) %*% beta_true + rnorm(n)
  } else if (family == "binomial") {
    eta <- cbind(Xa, Xb) %*% beta_true
    prob <- 1 / (1 + exp(-eta))
    y <- rbinom(n, 1, prob)
  } else if (family == "poisson") {
    eta <- cbind(Xa, Xb) %*% beta_true
    lambda <- exp(eta)
    y <- rpois(n, lambda)
  } else {
    stop("Unsupported family")
  }

  return(list(Xa = Xa, Xb = Xb, y = y))
}

set.seed(2024)
n <- 10000
p_a <- 5
p_b <- 5

#test for the linear regression case
data <- simulate_data(n, p_a, p_b, "gaussian")
results_gaussian <- ds.vertGLM(data$Xa, data$Xb, data$y, family = gaussian(), scale_data = TRUE)

X_full <- scale(cbind(data$Xa, data$Xb))
y_scaled <- scale(data$y)
glm_fit_gaussian <- glm(y_scaled ~ X_full - 1)

beta_gaussian <- c(results_gaussian$beta_a, results_gaussian$beta_b)
names(beta_gaussian) <- NULL
beta_glm_gaussian <- coef(glm_fit_gaussian)
names(beta_glm_gaussian) <- NULL



#test the binomial family of GLMs
data_binomial <- simulate_data(n, p_a, p_b, family = "binomial")

results_binomial <- ds.vertGLM(data_binomial$Xa, data_binomial$Xb, data_binomial$y, family = binomial(), scale_data = TRUE)

X_full_binomial <- scale(cbind(data_binomial$Xa, data_binomial$Xb))
glm_fit_binomial <- glm(data_binomial$y ~ X_full_binomial - 1, family = binomial())

beta_binomial <- c(results_binomial$beta_a, results_binomial$beta_b)
names(beta_binomial) <- NULL
beta_glm_binomial <- coef(glm_fit_binomial)
names(beta_glm_binomial) <- NULL


print(all.equal(beta_gaussian, beta_glm_gaussian, tolerance = 1e-6))
print(all.equal(beta_binomial, beta_glm_binomial, tolerance = 1e-6))

