# Identifiability Constraints in Generalized Additive Models

This repository contains code for reproducing the analysis in the paper *Identifiability Constraints in Generalized Additive Models*.

## Preparation

The code requires the following packages:
```
install.packages(
  "tidyverse",
  "mgcv",
  "TMB",
  "Matrix",
  "parallel",
  "gamair"
)
```

To run scripts `02` and `03`, currently you have to first manually set your working directory to be the root directory of this repository (where the scripts are), or otherwise change `tmbpath` at the top of each script. This will be fixed when the `github` repository eventually becomes public.

## Code

The analyses are reproduced by the following scripts:

- **Section 4.1**: `01-gaussian-simulations.R`
- **Section 4.2**: `02-poisson-simulations.R`, `03-bernoulli-simulations.R`
- **Section 5**: `04-chicago.R`
