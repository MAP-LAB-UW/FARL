# Simulated Dataset: 1D FARLR-Style Mixed-Format Item Responses

A simulated one-dimensional mixed-format response dataset generated
under a FARLR-style latent regression formulation. The latent trait
\\\theta\\ is constructed from a linear component plus scaled noise to
achieve a target signal-to-noise ratio (SNR), and is then centered and
standardized. Item responses are subsequently generated using
[`mirt::simdata()`](https://philchalmers.github.io/mirt/reference/simdata.html)
with a mixture of three-parameter logistic (3PL) and generalized partial
credit model (GPCM) items.

## Usage

``` r
data(sim_a2)
```

## Format

A list with the following components:

|  |  |
|----|----|
| `X` | Covariate/design matrix used for generating \\\theta\\. |
| `Y` | Simulated mixed-format item response matrix (`N` by `J`). The 3PL items have scores 0/1, while the GPCM items have scores 0/1/2. |
| `a` | True item discrimination parameters (slopes), length `J`. |
| `b` | True difficulty parameters for 2PL/3PL items. Entries corresponding to GPCM items are `NA`. |
| `c` | True lower-asymptote (guessing) parameters for 3PL items. Entries for non-3PL items are zero. |
| `b1` | First GPCM step-difficulty parameter. Entries for dichotomous items are `NA`. |
| `b2` | Second GPCM step-difficulty parameter. Entries for dichotomous items are `NA`. |
| `d` | Item intercept matrix used by [`mirt::simdata()`](https://philchalmers.github.io/mirt/reference/simdata.html), with columns `d0`, `d1`, and `d2`. Unused entries are `NA`. |
| `itemtype` | Character vector identifying each item as `"3PL"`, `"2PL"`, or `"gpcm"`. |
| `parTab` | Parameter table containing item types, slopes, difficulties, guessing parameters, GPCM step difficulties, and mirt intercepts. |

## Details

The latent trait is generated from a linear predictor and additive
noise:

- \\\theta = F \beta + E \nu + e\\, where \\e\\ is scaled to match a
  target SNR.

The default dataset contains 15 items: eight 3PL items and seven
three-category GPCM items. Item slopes are generated as \\a_j \sim
\mathrm{Lognormal}(0, 0.25)\\.

For a 3PL item, the probability of a correct response is \$\$P(Y\_{ij}=1
\mid \theta_i) = c_j + (1-c_j)
\mathrm{logit}^{-1}\\a_j(\theta_i-b_j)\\,\$\$ and the mirt intercept is
computed as \\d_j=-a_jb_j\\.

Each GPCM item has three response categories, 0, 1, and 2, with two step
difficulties, \\b\_{1j}\\ and \\b\_{2j}\\. These are converted to the
mirt intercept parameterization as \$\$d\_{0j}=0,\quad
d\_{1j}=-a_jb\_{1j},\quad d\_{2j}=-a_j(b\_{1j}+b\_{2j}).\$\$

Responses are simulated using the item-specific values in `itemtype`,
with the 3PL lower-asymptote parameters supplied through the `guess`
argument to
[`mirt::simdata()`](https://philchalmers.github.io/mirt/reference/simdata.html).

## Author

Yijun Cheng \<chengxb@uw.edu\>
