# flu-model


### Summary

This repository contains code for analysis of influenza antibody dynamics on multiple timescales.

### Guide to files

`main_model.R` Main model fitting and plotting code - calls following 3 source files:

> `load_data_[LOCATION].R` Loads data and reshapes to get in format efficient for simulation/inference

> `sero_functions.R` Functions for model and inference

> `posterior_analysis_flu.R` Reads in MCMC output and plots posteriors and model estimates for main figures

### Reference:

Kucharski AJ, Lessler J, Cummings DATC, Riley S (2018) [Timescales of influenza A/H3N2 antibody dynamics](https://journals.plos.org/plosbiology/article?id=10.1371/journal.pbio.2004974). PLOS Biology