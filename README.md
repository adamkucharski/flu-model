# flu-model

# SR Test of ability to push

### Summary

This repository contains code for annual version of following model:

Kucharski AJ, Lessler J, Read JM, Zhu H, Jiang CQ et al. (2015) Estimating the life course of influenza A(H3N2) antibody responses from cross-sectional data. PLOS Biol.

Currently, model just includes boosting and cross-reactivity, with Poisson observation process. Need to work out efficient way of including antigenic seniority and other mechanisms.

### Guide to files

`main_model.R` Main model fitting code - calls following 2 source files.

> `load_data.R` Loads data from HaNam and reshapes to get in format efficient for simulation/inference

> `sero_functions.R` Functions for model and inference

`plot_data.R` Gives basic summary plots of data, and some animations
