# kNNDM: k-fold Nearest Neighbour Distance Matching Cross-Validation for map accuracy estimation

This repository contains the R-scripts used in the paper "kNNDM: k-fold Nearest Neighbour Distance Matching Cross-Validation for map accuracy estimation" by Jan Linnenbrink, Carles Milà, Marvin Ludwig and Hanna Meyer. The manuscript has been submitted to the journal *Geoscientific Model Development*.

The figures are stored in [figures](figures/) and can be reproduced using the following functions:

* [figures_simulation1_virtualSpecies.Rmd](figures_simulation1_virtualSpecies.Rmd): R-Markdown file to reproduce the figures of simulation 1, which is presented in the main manuscript. [figures_simulation1_virtualSpecies.pdf](figures_simulation1_virtualSpecies.pdf) shows the output as one pdf document.

* [figures_simulation2_AGB.qmd](figures_simulation2_AGB.qmd): Quarto file to reproduce the figures of simulation 2, which is presented in the supplementary. [figures_simulation2_AGB.pdf](figures_simulation2_AGB.pdf) shows the output as one pdf document.



The two simulation studies [simulation1_virtualSpecies](simulation1_virtualSpecies/) and [simulation2_AGB](simulation2_AGB/) are organized in two directories with similar structures:

* [code] contains the code to reproduce the simulations.
* [data] contains the data used as inputs for the simulations. For the first simulation, these are the predictor stack ([species_stack.grd]) and the sampling area ([species_vdata.gpkg]). For the second simulation, the directory contains the prediction points ([ppoints.gpkg]), along with the 100 realizations of each sampling design.
* [results] contains the results of the simulations.
