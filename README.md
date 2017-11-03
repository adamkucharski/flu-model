# flu-model


### Summary

This repository contains code to accompany the paper "Timescales of influenza A/H3N2 antibody dynamics".

### Guide to files

`main_model.R` Main model fitting and plotting code - calls following 3 source files:

> `load_data_[LOCATION].R` Loads data and reshapes to get in format efficient for simulation/inference

> `sero_functions.R` Functions for model and inference

> `posterior_analysis_flu.R` Reads in MCMC output and plots posteriors and model estimates for main figures
