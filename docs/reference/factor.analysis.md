# Factor analysis

The main function for factor analysis with potentially high dimensional
variables. Here we implement some recent algorithms that is optimized
for the high dimensional problem where the number of samples n is less
than the number of variables p.

## Usage

``` r
factor.analysis(Y, r, method = c("ml", "pc", "esa"))
```

## Arguments

- Y:

  data matrix, a n\*p matrix

- r:

  number of factors

- method:

  algorithm to be used

## Value

a list of objects

- Gamma:

  estimated factor loadings

- Z:

  estimated latent factors

- Sigma:

  estimated noise variance matrix

## Details

The three methods are quasi-maximum likelihood (ml), principal component
analysis (pc), and factor analysis using an early stopping criterion
(esa).

The ml is iteratively solved the Expectation-Maximization algorithm
using the PCA solution as the initial value. See Bai and Li (2012) and
for more details. For the esa method, see Owen and Wang (2015) for more
details.

## References

Bai, J. and Li, K. (2012). Statistical analysis of factor models of high
dimension. *The Annals of Statistics 40*, 436-465. Owen, A. B. and Wang,
J. (2015). Bi-cross-validation for factor analysis. *arXiv:1503.03515*.

## See also

`fa.pc`, `fa.em`, `ESA`

## Examples

``` r
## a factor model
n <- 100
p <- 1000
r <- 5
Z <- matrix(rnorm(n * r), n, r)
Gamma <- matrix(rnorm(p * r), p, r)
Y <- Z %*% t(Gamma) + rnorm(n * p)

## to check the results, verify the true factors are in the linear span of the estimated factors.
pc.results <- factor.analysis(Y, r = 5, "pc")
#> Error in factor.analysis(Y, r = 5, "pc"): could not find function "factor.analysis"
sapply(summary(lm(Z ~ pc.results$Z)), function(x) x$r.squared)
#> Error in eval(predvars, data, env): object 'pc.results' not found

ml.results <- factor.analysis(Y, r = 5, "ml")
#> Error in factor.analysis(Y, r = 5, "ml"): could not find function "factor.analysis"
sapply(summary(lm(Z ~ ml.results$Z)), function(x) x$r.squared)
#> Error in eval(predvars, data, env): object 'ml.results' not found

esa.results <- factor.analysis(Y, r = 5, "esa")
#> Error in factor.analysis(Y, r = 5, "esa"): could not find function "factor.analysis"
sapply(summary(lm(Z ~ esa.results$Z)), function(x) x$r.squared)
#> Error in eval(predvars, data, env): object 'esa.results' not found
```
