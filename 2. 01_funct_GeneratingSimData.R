# This code is for the simulated data based off some factor models as seen in Duker and Pipiras (2024)

# set seed
set.seed(293) # Used class number as the seed

# Define the factor models

# Number of variables, observations, and factors
d <- 200 
T <- 200
r <- 3

# Significance level
alpha <- 0.05

# Magnitudes of departure from the null
delta_values <- c(0, 0.25, 0.50, 1)

# Use delta = 0.25 for the example
example_j <- which(delta_values == 0.25)

# First loading matrix and a random perturbation matrix
Lambda1 <- matrix(rnorm(d * r), nrow = d, ncol = r)
B <- matrix(rnorm(d * r), nrow = d, ncol = r)

# Nonsingular transformation: rescales the three loading columns
Phi <- diag(c(0.8, 1, 1.2))

# Generate independent standard normal factors
F1 <- matrix(rnorm(T * r), nrow = T, ncol = r)
F2 <- matrix(rnorm(T * r), nrow = T, ncol = r)

# Generate independent standard normal errors
E1 <- matrix(rnorm(T * d), nrow = T, ncol = d)
E2 <- matrix(rnorm(T * d), nrow = T, ncol = d)

# Generate the first dataset and subtract each column's mean
X1 <- F1 %*% t(Lambda1) + E1
X1 <- scale(X1, center = TRUE, scale = FALSE)

# Store one loading matrix and dataset for each delta
Lambda2 <- vector("list", length(delta_values))
X2 <- vector("list", length(delta_values))

for (j in seq_along(delta_values)) {
  delta <- delta_values[j]
  
  # delta = 0 satisfies the null; positive delta adds a perturbation
  Lambda2[[j]] <- Lambda1 %*% Phi + delta * B
  
  # Generate and center the second dataset
  X2[[j]] <- F2 %*% t(Lambda2[[j]]) + E2
  X2[[j]] <- scale(X2[[j]], center = TRUE, scale = FALSE)
}


# Visualization of Data

# Calculate the projection matrix for the first loading space.
Q1 <- qr.Q(qr(Lambda1))
P1 <- Q1 %*% t(Q1)

# Calculate P2 - P1 for each delta.
projection_differences <- vector("list", length(delta_values))

for (j in seq_along(delta_values)) {
  Q2 <- qr.Q(qr(Lambda2[[j]]))
  P2 <- Q2 %*% t(Q2)
  
  projection_differences[[j]] <- P2 - P1
}

# The first 3 factors
series_colors <- c("steelblue", "firebrick", "darkgreen")

# Use the same vertical limits across all series plots.
series_limits <- range(
  unlist(lapply(X2, function(X) X[, 1:3]))
)

# Display the first 30 variables in each projection heatmap.
n_show <- min(30, d)

# Use the same color scale across all projection plots.
# Blue = negative difference; white = zero; red = positive difference.
projection_colors <- colorRampPalette(
  c("blue", "white", "red")
)(101)

color_limit <- max(abs(unlist(projection_differences)))
color_limit <- max(color_limit, .Machine$double.eps)


# Function for one observed-series plot

# j identifies the delta scenario to display.
plot_series <- function(j) {
  matplot(
    1:T,
    X2[[j]][, 1:3],
    type = "l",
    lty = 1,
    col = series_colors,
    ylim = series_limits,
    xlab = "Time",
    ylab = "Observed value",
    main = bquote(delta == .(delta_values[j]))
  )
  
  legend(
    "topright",
    legend = paste("Series", 1:3),
    col = series_colors,
    lty = 1,
    bty = "n",
    cex = 0.7
  )
}


# Function for one projection-difference plot

plot_projection <- function(j) {
  image(
    1:n_show,
    1:n_show,
    projection_differences[[j]][1:n_show, 1:n_show],
    col = projection_colors,
    zlim = c(-color_limit, color_limit),
    xlab = "Variable",
    ylab = "Variable",
    main = bquote(P[2] - P[1] ~ ":" ~ delta == .(delta_values[j]))
  )
}


# Display all layouts

show_figures <- function(plot_function) {
  
  old_par <- par(no.readonly = TRUE)
  on.exit(par(old_par))
  
  # Find the scenario corresponding to delta = 0.
  zero <- which(delta_values == 0)
  
  # Each delta as a separate figure.
  par(mfrow = c(1, 1), mar = c(4, 4, 3, 1))
  
  for (j in seq_along(delta_values)) {
    plot_function(j)
  }
  
  # All four deltas in a 2 x 2 figure.
  par(mfrow = c(2, 2), mar = c(4, 4, 3, 1))
  
  for (j in seq_along(delta_values)) {
    plot_function(j)
  }
  
  # Compare delta = 0 with each positive delta (Pairs)
  for (j in which(delta_values > 0)) {
    par(mfrow = c(1, 2), mar = c(4, 4, 3, 1))
    
    plot_function(zero)
    plot_function(j)
  }
}


#Display data figures
show_figures(plot_series)


# Display the heatmap figures
show_figures(plot_projection)

