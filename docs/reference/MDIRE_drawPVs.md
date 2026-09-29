# Draw Plausible Values from a Multidimensional DIRE Model

Draws plausible values from an object fitted by \[MDIRE_mml()\]. The
function identifies the PV columns for every subtest without assuming
five dimensions or ten plausible values. Optional binary grouping
variables can be supplied to reproduce the subgroup summaries from the
original \`highDimDire()\` code.

## Usage

``` r
MDIRE_drawPVs(
  object,
  npv = 10L,
  pvVariableNameSuffix = "_dire",
  group_vars = NULL,
  ...
)
```

## Arguments

- object:

  An object returned by \[MDIRE_mml()\].

- npv:

  Positive integer. Number of plausible values per respondent.

- pvVariableNameSuffix:

  Suffix passed to \[Dire::drawPVs()\].

- group_vars:

  Optional character names or numeric indices of binary columns in
  \`object\$X\` used to form subgroup combinations.

- ...:

  Additional arguments passed to \[Dire::drawPVs()\].

## Value

A list containing raw PV data, dimension-specific PV matrices, optional
subgroup PV matrices and summaries, and the original DIRE PV result.
