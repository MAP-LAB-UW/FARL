# Simulated Dataset: 1D FARLR-Style Dichotomous Item Responses

A simulated one-dimensional dichotomous response dataset generated under
a FARLR-style latent regression formulation. The latent trait \\\theta\\
is constructed from a linear component plus scaled noise to achieve a
target signal-to-noise ratio (SNR), and is then centered and
standardized. Item responses are subsequently generated using
[`mirt::simdata()`](https://philchalmers.github.io/mirt/reference/simdata.html)
with dichotomous (2PL) items.

## Usage

``` r
data(sim_a1)
```

## Format

A list with the following components:

|  |  |
|----|----|
| `X` | Covariate/design matrix used for generating \\\theta\\. |
| `Y` | Simulated dichotomous item response matrix (`N` by `J`). |
| `a` | True item discrimination parameters (slopes), length `J`. |
| `b` | True item difficulty parameters, length `J`. |
| `d` | True item intercept parameters used by [`mirt::simdata()`](https://philchalmers.github.io/mirt/reference/simdata.html), computed as `d = -a*b`. |
| `parTab` | Parameter table used internally for model fitting. |

## Details

The latent trait is generated from a linear predictor and additive
noise:

- \\\theta = F \beta + E \nu + e\\, where \\e\\ is scaled to match a
  target SNR.

Item parameters are generated as:

- \\a_j \sim \mathrm{Lognormal}(0, 0.25)\\

- \\b_j \sim \mathrm{Uniform}(-2, 2)\\

- \\d_j = -a_j b_j\\

Responses are simulated via `mirt::simdata(itemtype = "dich")`.

## Author

Yijun Cheng \<chengxb@uw.edu\>
