# antiN_CoP
# A Bayesian hierarchical model to estimate maximum cross-neutralizing titers against SARS-CoV-2 variants

## Summary
This repository contains computational codes used in [Miyamoto et al., iScience, 2023](https://t.co/I24zCMfYdN).

For each exposure group, we estimated the neutralization titer and time of vaccination using a Bayesian hierarchical model. The log10 neutralization titer (NT) after breakthrough infection or booster vaccination was described using a three-parameter logistic model for each time interval between the second vaccination and the third exposure (vaccination or breakthrough infection). We inferred population means (μv) separately for neutralization titers against the ancestral strain, BA.1, BA.2, BA.2.75, BA.5, and BQ.1.1. We used a hierarchical structure to describe the distribution of µhv for each exposure group. Arrays in the model index over one or more indices: H=3 exposure history h; N=108 participants n; V=6 target viruses v. The model was as follows:

- NTnvt ~ Normal (µhv / (1 + αv exp(- βv tn)), σ_NTv)
- µhv ~ Normal (µv, σ_µv) [0, 5]
- µv ~ Normal (2.5, 1) [0, 5]
- αv ~ Normal (2.5, 1) [0, 5]
- βv ~ Student_t (4, 0.05, 0.1) [0, 1]
- σ_µv ~ Student_t (4, 0, 0.5) [0, ∞]

The values in square brackets denote the truncation bounds of the distributions. The explanatory variable was time, tn, and the outcome variable was NTnvt, which represented the neutralization titers against the target virus v in participant n at time t. A non-informative prior was set for the standard distribution σ_NTv. The parameters αv and βv controlled the intercept and the steepness of the logistic function, respectively. The mean parameter for neutralization titers against target virus v according to the exposure history h, µhv, was generated from a normal distribution with hyperparameters of the mean, µv, and standard deviation, σ_µv. For the distribution generating βv and σ_µv, we used a Student’s t distribution with four degrees of freedom, instead of a normal distribution, to reduce the effects of outlier values of βv and σ_µv.
The time interval of days to 90% maximum neutralization titers against each virus (tMNT90v) was calculated according to the parameters αv and βv as follows:

- tMNT90v = log(9αv) / βv

Parameter estimation was performed via the Markov chain Monte Carlo (MCMC) approach implemented in rstan. Four independent MCMC chains were run with 5,000 steps in the warm-up and sampling iterations, with subsampling every five iterations.

## Contents
- script
  - MaxNT_Interval_model.R :main script
  - interval_model.stan :stan model file
- input
  - Table S3.xlsx : Raw data available from [Miyamoto et al., iScience, 2023](https://t.co/I24zCMfYdN).
- output

## Dependencies
- R version 4.1.2 (2021-11-01)
- R packages:
  - rstan_2.26.13
  - dplyr_1.0.9
  - openxlsx_4.2.5.1
