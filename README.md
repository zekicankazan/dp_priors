# Priors with Differential Privacy

This repository contains code to perform the experiments in "Prior Distributions for Bayesian Inference for Gaussian Model Parameters Under Differential Privacy Ensured by Enforcing Bounds on Data Values" (Kazan and Reiter, 2025). The file `functions.R` contains a function to draw samples from the TGM distribution and a function to run our Gibbs sampler. The `Figures_scripts` folder contains scripts to produce each figure in the manuscript. The figures are saved to the `Figures` folder.

The file `simulation_study.R` runs the simulation study used to produce Figures 2 and 3 (which we ran on a shared compute cluster). The results are saved in the compressed `coverage_analysis` folder, and can be read in via the file `simulation_study_saved.R`.

The folder `JAGR_comparison` contains the material necessary to create Figures 6-8, which compare our proposal to a constrained likelihood alternative, implemented via the Gibbs sampler of Ju et al. (2022). The file `JAGR_function` contains a function to run the Gibbs sampler of Ju et al. (2022). The file `simulation_study_JAGR.R` runs the simulation study used to produce Figure 8 (which we ran on a shared compute cluster). The results are saved in the compressed `coverage_analysis_JAGR` subfolder, and code to read them in is provided in `Figures_scripts/Figure8.R`.

This research was supported by NSF grant SES-2217456.
