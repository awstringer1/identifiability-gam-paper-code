# Identifiability Constraints in Generalized Additive Models

This repository contains code for reproducing the analysis in the paper *Identifiability Constraints in Generalized Additive Models*, [published](https://onlinelibrary.wiley.com/doi/10.1002/cjs.11786) in the *Canadian Journal of Statistics* (open-access):

Stringer, A. (2023). Identifiability Constraints in Generalized Additive Models. *Canadian Journal of Statistics* (to appear).

## Preparation

The code requires the following packages:
```
install.packages(
  "tidyverse",
  "here",
  "mgcv",
  "TMB",
  "Matrix",
  "parallel",
  "gamair"
)
```

Note: the `here` package is only required for portability of relative file paths, to make it a little easier to execute these scripts yourself.
It isn't required to reproduce the results in the paper.

## Code

The analyses are reproduced by the following scripts:

- **Section 4.1**: `01-gaussian-simulations.R`
- **Section 4.2**: `02-poisson-simulations.R`, `03-bernoulli-simulations.R`
- **Section 5**: `04-chicago.R`
