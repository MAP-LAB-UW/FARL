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
v1 <- 0.5
v2 <- 1 # for 10:1
J <- 10

Sigma <- diag(K)
Sigma[Sigma == 0] <- runif(K * K - K, -0.2, 0.2)
F <- mvrnorm(n = n, rep(0, K), Sigma)
F <- scale(F, center = TRUE, scale = TRUE)

B1 <- matrix(
  data = c(rep(c(0.5, 1, 1.5), p / (3 * K)), rep(0, p / K)),
  byrow = FALSE,
  nrow = p / K,
  ncol = K
)
B2 <- matrix(
  data = c(rep(0, p / K), rep(c(0.5, 1, 1.5), p / (3 * K))),
  byrow = FALSE,
  nrow = p / K,
  ncol = K
)
B <- cbind(t(B1), t(B2))

index0.5 <- seq(1, p, by = 3)
index1 <- seq(2, p, by = 3)
index1.5 <- seq(3, p, by = 3)

## Generate X = F B + E, E ~ N(0, rho * I_p)
sigma_E <- diag(rho, p)
E <- mvrnorm(n, rep(0, p), sigma_E)

FB <- F %*% B
var_FB <- sd(FB)^2
target_var_E <- (0.1 / 0.9) * var_FB
current_var_E <- sd(E)^2
scaling_factor <- sqrt(target_var_E / current_var_E)
E <- E * scaling_factor

X <- F %*% B + E

## Generate the latent outcome
beta <- runif(K, 0.75, 1.25)
v <- runif(p * percent / 100, 0.75, 1.25)
cat("num of nonXs = ", p * percent / 100, "\n")
e <- rnorm(n, 0, sigma_e)
current_var_e <- sd(e)^2

q <- 1
veta <- rep(0, p)
bin <- c(index0.5[1], index1[1], index1[10], index1.5[5], index1.5[15])
main_vars <- paste0("X", bin)

for (b_index in bin) {
  X[, b_index] <- ifelse(X[, b_index] < 0, 0, 1)
  veta[b_index] <- v[q]
  q <- q + 1
}

binary_index <- sample(
  setdiff(1:p, bin),
  size = p / 2 - 5,
  replace = FALSE
)
for (b_index in binary_index) {
  X[, b_index] <- ifelse(X[, b_index] < 0, 0, 1)
}

X <- X[complete.cases(X), ]
X <- X[, apply(X, 2, sd, na.rm = TRUE) != 0]
coviates <- X

X_discrete <- apply(X[, bin], 2, function(col) {
  ifelse(col > mean(col, na.rm = TRUE), 1, 0)
})

examp <- if (any(c(v1, v2) == 0.5)) index1.5 else index0.5
cat("examp = ", examp, "\n")
setofimportantXs <- setdiff(1:p, examp)
candidates <- setdiff(setofimportantXs, bin)
index <- sample(candidates, size = length(v) - 5, replace = FALSE)

for (j in seq_len(length(v) - 5)) {
  veta[index[j]] <- v[j + 5]
}

main <- bin
s <- length(bin)
combinations <- expand.grid(replicate(s, 0:1, simplify = FALSE))
colnames(combinations) <- main[1:s]

Y_temp <- F %*% beta + E %*% veta
var_Y <- sd(Y_temp)^2
target_var_e <- (1 / (10 - 1)) * var_Y
scaling <- sqrt(target_var_e / current_var_e)
e <- e * scaling
Y <- Y_temp + e
norm_c <- sd(Y)
theta <- Y / norm_c
mean_c <- mean(theta)
theta <- theta - mean_c

## Generate mixed-format item responses
# Change this vector to choose the item type for each of the J items.
# Supported here: "2PL", "3PL", and three-category "gpcm".
itemtype <- rep(c("3PL", "gpcm"), length.out = J)

a <- rlnorm(J, 0, 0.25)

three_pl_items <- which(itemtype == "3PL")
two_pl_items <- which(itemtype == "2PL")
gpcm_items <- which(itemtype == "gpcm")
dich_items <- c(two_pl_items, three_pl_items)

# The usual b and c parameters apply to dichotomous items.
b <- rep(NA_real_, J)
b[dich_items] <- runif(length(dich_items), -2, 2)

c_param <- rep(0, J)
c_param[three_pl_items] <- runif(length(three_pl_items), 0.10, 0.25)

# GPCM items have categories 0, 1, and 2, hence two step difficulties.
b1 <- rep(NA_real_, J)
b2 <- rep(NA_real_, J)
if (length(gpcm_items) > 0) {
  step_center <- runif(length(gpcm_items), -1, 1)
  step_gap <- runif(length(gpcm_items), 0.5, 1.5)
  b1[gpcm_items] <- step_center - step_gap / 2
  b2[gpcm_items] <- step_center + step_gap / 2
}

# Convert the familiar IRT parameters to the intercept form used by
# mirt::simdata(). For a three-category GPCM:
#   d0 = 0, d1 = -a*b1, d2 = -a*(b1 + b2).
d <- matrix(
  NA_real_,
  nrow = J,
  ncol = 3,
  dimnames = list(NULL, c("d0", "d1", "d2"))
)
d[dich_items, 1] <- -a[dich_items] * b[dich_items]
d[gpcm_items, 1] <- 0
d[gpcm_items, 2] <- -a[gpcm_items] * b1[gpcm_items]
d[gpcm_items, 3] <- -a[gpcm_items] * (b1[gpcm_items] + b2[gpcm_items])

resp <- mirt::simdata(
  a = a,
  d = d,
  itemtype = itemtype,
  Theta = theta,
  guess = c_param
)

item_params <- data.frame(
  item = 1:J,
  b = b,
  a = a,
  c = c_param,
  b1 = b1,
  b2 = b2,
  itemtype = itemtype,
  d0 = d[, 1],
  d1 = d[, 2],
  d2 = d[, 3]
)

colnames(resp) <- paste0("i", sprintf("%03d", 1:J))
subject <- factor(1:n)
resp0 <- data.frame(cbind(resp, subject))

itemNames <- paste0("i", sprintf("%03d", 1:J))
stuItems <- reshape(
  data = resp0,
  varying = itemNames,
  idvar = "subject",
  direction = "long",
  v.names = "score",
  times = itemNames,
  timevar = "key"
)

new_itemNames <- paste0("item", 1:J)
stuItems$key <- rep(new_itemNames, each = n)

parTab <- item_params %>%
  mutate(
    ItemID = new_itemNames,
    test = "comp",
    subtest = rep("main", each = J),
    slope = a,
    difficulty = b,
    guessing = c,
    D = 1
  )

sim_a2 <- list(
  X = X,
  Y = resp,
  a = a,
  b = b,
  c = c_param,
  b1 = b1,
  b2 = b2,
  d = d,
  itemtype = itemtype,
  parTab = parTab
)

# In an R package, this saves data/sim_a2.rda.
usethis::use_data(sim_a2, overwrite = TRUE)
