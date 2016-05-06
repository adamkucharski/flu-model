# flu-model


### Summary

This repository contains code for extended annual version of following model:

Kucharski AJ, Lessler J, Read JM, Zhu H, Jiang CQ et al. (2015) Estimating the life course of influenza A(H3N2) antibody responses from cross-sectional data. PLOS Biol.

Model includes boosting and cross-reactivity, with Poisson observation process. Also include antigenic seniority (boosting + suppression), short-term waning and uniform measurement error.

### Guide to files

`main_model.R` Main model fitting code - calls following 2 source files.

> `load_data.R` Loads data from HaNam and reshapes to get in format efficient for simulation/inference

> `sero_functions.R` Functions for model and inference

`simulation_diagnostics.R` Reads in MCMC output and plots posteriors and model estimates

`plot_data.R` Gives basic summary plots of data, and some animations
