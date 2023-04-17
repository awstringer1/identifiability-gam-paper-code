# Identifiability Constraints in Generalized Additive Models

This repository contains code for reproducing the analysis in the paper *Identifiability Constraints in Generalized Additive Models*.

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
