# flu-model


### Summary

This repository contains code to accompany the paper "Timescales of influenza A/H3N2 antibody dynamics".

### Guide to files

`main_model.R` Main model fitting and plotting code - calls following 2 source files:

> `load_data.R` Loads data from Vietnam and China and reshapes to get in format efficient for simulation/inference

> `sero_functions.R` Functions for model and inference

`simulation_diagnostics.R` Reads in MCMC output and plots posteriors and model estimates. Not required for final fitting and plotting.
