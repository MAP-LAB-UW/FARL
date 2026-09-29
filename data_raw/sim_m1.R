library(MASS)
library(dplyr)

set.seed(3)

n <- 3000
rho <- 1
sigma_e <- 1
scaling <- 10
K <- 2
p <- 60
percent <- 20
D <- 5
J <- 6                       # number of items per dimension
J_total <- J * D

# -------------------------------------------------------------------------
# Generate common factors
# -------------------------------------------------------------------------
Sigma <- diag(K)
Sigma[upper.tri(Sigma)] <- runif(K * (K - 1) / 2, -0.2, 0.2)
Sigma[lower.tri(Sigma)] <- t(Sigma)[lower.tri(Sigma)]

F <- mvrnorm(
  n = n,
  mu = rep(0, K),
  Sigma = Sigma
)
F <- scale(F, center = TRUE, scale = TRUE)

B1 <- matrix(
  data = c(
    rep(c(0.5, 1, 1.5), p / (3 * K)),
    rep(0, p / K)
  ),
  byrow = FALSE,
  nrow = p / K,
  ncol = K
)

B2 <- matrix(
  data = c(
    rep(0, p / K),
    rep(c(0.5, 1, 1.5), p / (3 * K))
  ),
  byrow = FALSE,
  nrow = p / K,
  ncol = K
)

B <- cbind(t(B1), t(B2))

index0.5 <- seq(1, p, by = 3)
index1 <- seq(2, p, by = 3)
index1.5 <- seq(3, p, by = 3)

# -------------------------------------------------------------------------
# Generate X = F B + E
# -------------------------------------------------------------------------
sigma_E <- diag(rho, p)
E <- mvrnorm(n, rep(0, p), sigma_E)

FB <- F %*% B
var_FB <- sd(FB)^2
target_var_E <- (0.1 / 0.9) * var_FB
current_var_E <- sd(E)^2
scaling_factor <- sqrt(target_var_E / current_var_E)
E <- E * scaling_factor

X <- F %*% B + E
colnames(X) <- paste0("X", seq_len(p))

# -------------------------------------------------------------------------
# Generate five dimension-specific latent regressions
# -------------------------------------------------------------------------
n_nonzero <- as.integer(p * percent / 100)
cat("num of nonzero Xs per dimension = ", n_nonzero, "\n")

beta <- matrix(
  runif(K * D, 0.75, 1.25),
  nrow = K,
  ncol = D,
  dimnames = list(
    paste0("F", seq_len(K)),
    paste0("Dimension", seq_len(D))
  )
)

v <- matrix(
  runif(n_nonzero * D, 0.75, 1.25),
  nrow = n_nonzero,
  ncol = D
)

veta <- matrix(
  0,
  nrow = p,
  ncol = D,
  dimnames = list(
    colnames(X),
    paste0("Dimension", seq_len(D))
  )
)

bin <- c(
  index0.5[1],
  index1[1],
  index1[10],
  index1.5[5],
  index1.5[15]
)
main_vars <- paste0("X", bin)

for (b_index in bin) {
  X[, b_index] <- ifelse(X[, b_index] < 0, 0, 1)
}
veta[bin, ] <- v[seq_along(bin), , drop = FALSE]

binary_index <- sample(
  setdiff(seq_len(p), bin),
  size = p / 2 - length(bin),
  replace = FALSE
)

for (b_index in binary_index) {
  X[, b_index] <- ifelse(X[, b_index] < 0, 0, 1)
}

X <- X[complete.cases(X), , drop = FALSE]
X <- X[, apply(X, 2, sd, na.rm = TRUE) != 0, drop = FALSE]

X_discrete <- apply(
  X[, bin, drop = FALSE],
  2,
  function(column) ifelse(column > mean(column, na.rm = TRUE), 1, 0)
)

# Exclude every third variable from the additional active-variable candidates.
excluded_variables <- index1.5
set_of_important_Xs <- setdiff(seq_len(p), excluded_variables)
candidates <- setdiff(set_of_important_Xs, bin)
active_index <- vector("list", D)

for (dimension in seq_len(D)) {
  active_index[[dimension]] <- sample(
    candidates,
    size = n_nonzero - length(bin),
    replace = FALSE
  )

  veta[active_index[[dimension]], dimension] <-
    v[(length(bin) + 1):n_nonzero, dimension]
}

names(active_index) <- paste0("Dimension", seq_len(D))

main <- bin
s <- length(bin)
combinations <- expand.grid(replicate(s, 0:1, simplify = FALSE))
colnames(combinations) <- main_vars

Y_temp <- F %*% beta + E %*% veta

# Generate correlated residual errors and scale each dimension to SNR = 10:1.
Sigma_error <- matrix(0.75, nrow = D, ncol = D)
diag(Sigma_error) <- 1

error <- mvrnorm(
  n = n,
  mu = rep(0, D),
  Sigma = Sigma_error
)

for (dimension in seq_len(D)) {
  current_var_e <- sd(error[, dimension])^2
  var_Y <- sd(Y_temp[, dimension])^2
  target_var_e <- var_Y / (scaling - 1)
  error[, dimension] <-
    error[, dimension] * sqrt(target_var_e / current_var_e)
}

Y <- Y_temp + error
theta_center <- colMeans(Y)
theta_scale <- apply(Y, 2, sd)

theta <- sweep(Y, 2, theta_center, "-")
theta <- sweep(theta, 2, theta_scale, "/")
colnames(theta) <- paste0("Theta", seq_len(D))

# True latent-regression mean and residual covariance on the standardized
# theta scale used to generate the item responses.
true_prior_mean <- sweep(Y_temp, 2, theta_center, "-")
true_prior_mean <- sweep(true_prior_mean, 2, theta_scale, "/")
colnames(true_prior_mean) <- colnames(theta)

standardized_error <- sweep(error, 2, theta_scale, "/")
true_sigma <- stats::cov(standardized_error)
dimnames(true_sigma) <- list(colnames(theta), colnames(theta))

# -------------------------------------------------------------------------
# Generate 10 simple-structure 2PL items per dimension
# -------------------------------------------------------------------------
a_value <- rlnorm(J_total, 0, 0.25)
b <- runif(J_total, -2, 2)
d <- -a_value * b

a <- matrix(
  0,
  nrow = J_total,
  ncol = D,
  dimnames = list(
    paste0("item", seq_len(J_total)),
    paste0("Dimension", seq_len(D))
  )
)

item_index <- vector("list", D)

for (dimension in seq_len(D)) {
  items <-
    ((dimension - 1) * J + 1):
    (dimension * J)

  item_index[[dimension]] <- items
  a[items, dimension] <- a_value[items]
}

names(item_index) <- paste0("Dimension", seq_len(D))

resp <- mirt::simdata(
  a = a,
  d = d,
  itemtype = rep("dich", J_total),
  Theta = theta
)

colnames(resp) <- paste0("i", sprintf("%03d", seq_len(J_total)))

item_dimension <- rep(seq_len(D), each = J)
item_params <- data.frame(
  item = seq_len(J_total),
  dimension = item_dimension,
  b = b,
  a = a_value,
  d = d,
  c = 0
)

subject <- factor(seq_len(n))
resp0 <- data.frame(resp, subject = subject)

itemNames <- paste0("i", sprintf("%03d", seq_len(J_total)))
stuItems <- reshape(
  data = resp0,
  varying = itemNames,
  idvar = "subject",
  direction = "long",
  v.names = "score",
  times = itemNames,
  timevar = "key"
)

new_itemNames <- paste0("item", seq_len(J_total))
stuItems$key <- rep(new_itemNames, each = n)

parTab <- item_params %>%
  mutate(
    ItemID = new_itemNames,
    test = "mcomp",
    subtest = paste0("dimension", dimension),
    slope = a,
    difficulty = b,
    guessing = c,
    D = 1
  )

sim_m1 <- list(
  X = X,
  Y = resp,
  theta = theta,
  a = a,
  b = b,
  d = d,
  parTab = parTab,
  item_index = item_index,
  main = main,
  main_vars = main_vars,
  X_discrete = X_discrete,
  beta = beta,
  veta = veta,
  active_index = active_index,
  Sigma_error = Sigma_error,
  true_prior_mean = true_prior_mean,
  true_sigma = true_sigma,
  theta_center = theta_center,
  theta_scale = theta_scale
)

stopifnot(
  all(dim(sim_m1$X) == c(n, p)),
  all(dim(sim_m1$Y) == c(n, J_total)),
  all(dim(sim_m1$theta) == c(n, D)),
  all(dim(sim_m1$a) == c(J_total, D)),
  all(dim(sim_m1$true_prior_mean) == c(n, D)),
  all(dim(sim_m1$true_sigma) == c(D, D)),
  length(sim_m1$item_index) == D,
  all(colSums(sim_m1$a != 0) == J)
)

usethis::use_data(sim_m1, overwrite = TRUE, compress = "xz")
