
<!-- README.md is generated from README.Rmd. Please edit that file -->

# Integrated Population Model for winter steelhead on the WA coast

## Summary

This repository contains code to fit a hierarchical Integrated
Population Model (IPM) to multiple populations of wild winter steelhead
(*Oncorhynchus mykiss*) on the Washington coast. The IPM is a
statistical population dynamics model that integrates information on
spawner abundances, total harvest, and spawner age structure into a
combined run-reconstruction and spawner-recruitment model. It accounts
for iteroparity by distinguishing maiden from repeat spawners and
estimates time-varying kelt survival rates. The model also estimates
time-varying recruitment residuals and population parameters such as
productivity and capacity.

## Repository Structure

The repository contains an `R` folder with the following directories:

| Directory | Description                                           |
|-----------|-------------------------------------------------------|
| `code`    | Contains code to run the analyses and produce figures |
| `data`    | Contains fish and covariate data used in the analysis |
| `output`  | Used to save model output and data after fitting IPMs |

## Analyses

The `code` directory contains the following scripts:

| File                   | Description                              |
|------------------------|------------------------------------------|
| `IPM_sthd_multipop.R`  | Fits integrated population models (IPMs) |
| `posthoc_analysis.R`   | Posthoc analysis for covariate selection |
| `manuscript_figures.R` | Produces the main text figures           |

To produce the main results and figures, run code as described below:

1.  `IPM_sthd_multipop.R` using covar_effects=FALSE
2.  `posthoc_analysis.R` to identify covariates to be included in the
    IPM
3.  `IPM_sthd_multipop.R` using covar_effects=TRUE
4.  `manuscript_figures.R` to produce the main text figures

## Dependencies

The model fitting relies on the R package ‘salmonIPM’. This package is
available on GitHub at <https://github.com/ebuhle>.

## Output

The code described above produces the main results figures shown in
Ohlberger et al. (unpublished).
