library(MASS)
library(dplyr)
set.seed(3)
n=3000; rho = 1; sigma_e=1; scaling <- 10; K = 2; p = 60; percent = 20; v1 = 0.5; v2 = 1# for 10:1
J <- 15
Sigma <-diag(K)
Sigma[Sigma == 0] <- runif(K*K - K, -0.2,0.2)
F <- mvrnorm(n = n, rep(0,K), Sigma)

F <- scale(F, center = TRUE, scale = TRUE)

B1 <- matrix(data = c(rep(c(0.5,1,1.5), (p/(3*K))), rep(0, p/K)), byrow = FALSE,
             nrow = p/K, ncol = K)
B2 <- matrix(data = c(rep(0, p/K), rep(c(0.5,1,1.5), (p/(3*K)))), byrow = FALSE,
             nrow = p/K, ncol = K)
B <- cbind(t(B1), t(B2))
index0.5 <- seq(1, p, by = 3)
index1 <- seq(2, p, by = 3)
index1.5 <- seq(3, p, by = 3)

## Generate X = W^T U + E, E ~ N( 0, rho * I_p)
sigma_E <- diag(rho,p)
E <- mvrnorm(n, rep(0,p), sigma_E)

FB <- F %*% B
var_FB <- sd(FB)^2
target_var_E <- (0.1 / 0.9) * var_FB
current_var_E <- sd(E)^2
scaling_factor <- sqrt(target_var_E / current_var_E)
E <- E * scaling_factor

X <- F %*% B + E
#prob = 1/(1+exp(D))
#X <- matrix(rbinom(length(prob), size = 1, prob = as.vector(prob)), nrow = nrow(prob), ncol = ncol(prob))
#X <- scale(X, center = TRUE, scale = TRUE)

## Generate Y
beta <- runif(K, 0.75,1.25)
v <- runif((p*percent/100), 0.75,1.25)
cat("num of nonXs = ", p*percent/100, "\n")
#v <- sample(c(0.5,1),0.2*p, replace = T)
e <- rnorm(n, 0, sigma_e)
current_var_e <- sd(e)^2

q <- 1
veta <- rep(0,p)
bin <- c(index0.5[1],index1[1], index1[10], index1.5[5], index1.5[15])
main_vars <- paste0("X", bin, sept = "")
#bin <- c(1,32)
for (b in bin) {
  X[,b] <- ifelse(X[,b] <0, 0, 1)
  veta[b] <- v[q]
  q <- q + 1
}
binary_index <- sample( (setdiff(c(1:p), bin)), size = (p/2 - 5),replace = FALSE)
for (b in binary_index) {
  X[,b] <- ifelse(X[,b] <0, 0, 1)
}
X <- X[complete.cases(X),]
X <- X[, apply(X, 2, sd, na.rm = TRUE) != 0] # remove column with sd = 0
coviates <- X
X_discrete <- apply(X[,bin], 2, function(col) {
  ifelse(col > mean(col, na.rm = TRUE), 1, 0)
})
examp <- if(any(c(v1, v2)==0.5)) (index1.5) else (index0.5)
cat("examp = ", examp, "\n")
setofimportantXs <- setdiff(1:p, examp)
candidates <- setdiff(setofimportantXs, bin)
index <- sample(candidates, size = (length(v) - 5), replace = FALSE)
for (j in c(1:(length(v)-5))) {
  veta[index[j]] <- v[j + 5]
}
main <- bin
s <- length(bin)
combinations <- expand.grid(replicate(s, 0:1, simplify = FALSE))
colnames(combinations) <- main[1:s]  # Assuming 'main' contains names of your 5 binary variables

Y_temp <- F %*% (beta) + (E) %*% veta
var_Y <- sd(Y_temp)^2
target_var_e <- (1/(10-1)) * var_Y
scaling <- sqrt(target_var_e / current_var_e)
e <- e * scaling
Y <- Y_temp + e
norm_c <- sd(Y)
theta <- Y/norm_c
mean_c <- mean(theta)
theta <- theta - mean_c

## Generate
a <- rlnorm(J,0,0.25)
b <- runif(J,-2,2)
d <- -a*b

resp <- mirt::simdata(a=a,d=d,itemtype = "dich",Theta = theta)
item_params <- data.frame(cbind(
  item = c(1:J),
  b = b,
  a = a,
  c = 0
))
colnames(resp) <- paste0("i", sprintf("%03d", 1:J))
subject <- factor(c(1:n))
resp0 <- data.frame(cbind(resp, subject))

itemNames <- paste0("i", sprintf("%03d", 1:J))
stuItems <- reshape(data=resp0, varying=itemNames, idvar="subject",
                    direction="long", v.names="score",
                    times=itemNames, timevar="key")
new_itemNames <- (paste0("item",1:J))
stuItems$key <- rep(new_itemNames,each=n)

parTab <- item_params %>%
  mutate(ItemID = new_itemNames,
         test = "comp",
         subtest = rep("main",each=J),
         slope = a,
         difficulty = b,
         guessing = c,
         D = 1)
sim_a1 <- list(X =X, Y = resp, a = a, b = b, d = d, parTab = parTab)
usethis::use_data(sim_a1, overwrite = TRUE)
