# R codes for Infection risk estimation and Modeling the antibody response

[![DOI](https://zenodo.org/badge/898850069.svg)](https://doi.org/10.5281/zenodo.14906348)
## Summary
This repository contains computational codes used in [Miyamoto et al., Commun Med, 2025](https://doi.org/10.1038/s43856-025-00894-8).

## Contents
- script
  - Anti_N antibody_response.R : Modeling the antibody response main script
    - nmodel01.stan : Modeling the antibody response stan model file
  - infection_risk_estimation.R


## Dependencies
- R version 4.4.1 (2024-06-14)
- R packages:
  - rstan_2.32.6
  - dplyr_1.1.4
  - ggplot2_3.5.1
  - brms_2.21.0
  - gratia_0.9.2
  - tidyr_1.3.1
  - viridis_0.6.5
  - scales_1.3.0
  - cmdstanr_0.8.1
