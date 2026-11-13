# 3. Define the functions

# Estimate a loading-space basis using PCA
loading_basis <- function(X, r) {
  svd(X, nu = 0, nv = r)$v[, 1:r, drop = FALSE]
}

# Select the number of factors and estimate their values
estimate_factors <- function(Y, rmax = 8) {
  
  n <- nrow(Y)
  d <- ncol(Y)
  
  # Keep the candidate count below the matrix dimensions
  rmax <- min(rmax, n - 1, d - 1)
  
  # Calculate the PCA decomposition once
  fit <- svd(Y, nu = rmax, nv = 0)
  
  # Compare models containing 0, 1, ..., rmax factors
  candidates <- 0:rmax
  IC <- numeric(length(candidates))
  
  # Bai-Ng ICp2 penalty for each additional factor
  penalty <- ((n + d) / (n * d)) * log(min(n, d))
  
  for (k in candidates) {
    
    # Average squared error remaining after retaining k factors
    residual_variance <- (
      sum(fit$d^2) - sum(fit$d[seq_len(k)]^2)
    ) / (n * d)
    
    # Prevent taking the logarithm of zero
    residual_variance <- max(
      residual_variance,
      .Machine$double.eps
    )
    
    IC[k + 1] <- log(residual_variance) + k * penalty
  }
  
  # Choose the factor count with the smallest criterion
  selected_r <- candidates[which.min(IC)]
  
  # Normalize factors so that t(F) %*% F / n = I
  F <- sqrt(n) *
    fit$u[, seq_len(selected_r), drop = FALSE]
  
  list(r = selected_r, F = F)
}


# Calculate the Wald statistic from estimated factors
wald_test <- function(F, bandwidth = ceiling(nrow(F)^(1 / 3))) {
  
  n <- nrow(F)
  r <- ncol(F)
  
  stopifnot(n %% 2 == 0, r > 0)
  
  # Split the observations at the midpoint
  first <- 1:(n / 2)
  second <- (n / 2 + 1):n
  
  # Select the lower triangle, including the diagonal
  lower <- lower.tri(matrix(0, r, r), diag = TRUE)
  
  products <- do.call(
    rbind,
    lapply(1:n, function(t) {
      tcrossprod(F[t, ])[lower]
    })
  )
  
  # Compare factor second moments between the two halves
  V <- (
    colSums(products[first, , drop = FALSE]) -
      colSums(products[second, , drop = FALSE])
  ) / sqrt(n)
  
  # Subtract the identity entries
  Z <- sweep(products, 2, diag(r)[lower], "-")
  
  # Start with the lag-zero covariance.
  Omega <- crossprod(Z) / n
  
  # Add lagged covariances using Bartlett weights
  bandwidth <- max(1, min(bandwidth, n - 1))
  
  if (bandwidth > 1) {
    for (lag in 1:(bandwidth - 1)) {
      
      Gamma <- crossprod(
        Z[(lag + 1):n, , drop = FALSE],
        Z[1:(n - lag), , drop = FALSE]
      ) / n
      
      weight <- 1 - lag / bandwidth
      Omega <- Omega + weight * (Gamma + t(Gamma))
    }
  }
  
  # Stop rather than report an unreliable inverse.
  if (rcond(Omega) < 1e-10) {
    stop("The estimated covariance is nearly singular.")
  }
  
  # Calculate the statistic and its chi-square p-value.
  W <- as.numeric(t(V) %*% solve(Omega, V))
  df <- r * (r + 1) / 2
  p_value <- pchisq(W, df = df, lower.tail = FALSE)
  
  c(W = W, p_value = p_value, bandwidth = bandwidth)
}


# 4. State the hypothesis

cat(
  "H0: Lambda_k = Lambda Phi_k, for k = 1, 2.\n",
  "H1: The loading spaces differ.\n"
)

# Compare the TRUE loading spaces from the simulation.
space_difference <- numeric(length(delta_values))

for (j in seq_along(delta_values)) {
  space_difference[j] <- sqrt(
    sum(projection_differences[[j]]^2)
  )
}

kable(
  data.frame(
    delta = delta_values,
    projection_difference = space_difference
  ),
  digits = 3,
  caption = "Differences between the true loading spaces"
)


# 5. Transform the data and run the test
# we should be able to cut T in half
stopifnot(T %% 2 == 0)

# Estimate loading space.
Q1 <- loading_basis(X1, r)

first <- 1:(T / 2)
second <- (T / 2 + 1):T

results <- list()

for (j in seq_along(delta_values)) {
  
  # Select Xs for this delta.
  X <- X2[[j]]
  
  # Estimate X2 loading space.
  Q2 <- loading_basis(X, r)
  
  # Transform the first half using X1's projection.
  Y <- X
  
  Y[first, ] <- (
    X[first, ] + X[first, ] %*% Q1 %*% t(Q1)
  ) / 2
  
  # Transform the second half using X2's projection.
  Y[second, ] <- (
    X[second, ] + X[second, ] %*% Q2 %*% t(Q2)
  ) / 2
  
  # Select and estimate factors from the entire transformed series.
  fit <- estimate_factors(Y, rmax = 8)
  
  W <- NA_real_
  p_value <- NA_real_
  bandwidth <- NA_real_
  
  # Run the test if at least one factor was selected.
  if (fit$r > 0) {
    test <- wald_test(fit$F)
    
    W <- unname(test["W"])
    p_value <- unname(test["p_value"])
    bandwidth <- unname(test["bandwidth"])
  }
  
  results[[j]] <- data.frame(
    delta = delta_values[j],
    selected_factors = fit$r,
    bandwidth = bandwidth,
    df = fit$r * (fit$r + 1) / 2,
    W = W,
    p_value = p_value
  )
}


# 6. Display the conclusions
test_results <- do.call(rbind, results)

test_results$conclusion <- ifelse(
  is.na(test_results$p_value),
  "Test unavailable",
  ifelse(
    test_results$p_value < alpha,
    "Reject H0",
    "Do not reject H0"
  )
)

kable(
  test_results,
  digits = 3,
  caption = "Simulation results for each delta"
)
