# kNNDM_paper
# Nearest Neighbour Distance Matching k-fold Cross-Validation for map validation

This repository contains the R-scripts used in the paper "k-fold Nearest Neighbour Distance Matching Cross-Validation for map validation in the CAST R package" by Jan Linnenbrink, Carles Milà, Marvin Ludwig and Hanna Meyer. The manuscript has been submitted to the journal *Methods in Ecology and Evolution*.

The repository is structured as follows:

* [code](code/): contains the R-code to reproduce all analyses.
** [simulation](code/simulation/): contains the R-code to reproduce the simulation.
*** [sim_analysis.R](code/simulation/sim_analysis.R): R-Script to run the simulation.
*** [sim_functions.R](code/simulation/sim_functions.R): contains the functions called from [sim_analysis.R](code/simulation/sim_analysis.R).
*** [sim_utils.R](code/simulation/sim_utils.R): contains additional helper functions used in [sim_analysis.R](code/simulation/sim_analysis.R) and [sim_analysis_W.R](code/simulation/sim_analysis_W.R).
*** [sim_landscape.R](code/simulation/sim_landscape.R): R-Script to create the landscape dataset used in [sim_analysis.R](code/simulation/sim_analysis.R) and [sim_analysis_W.R](code/simulation/sim_analysis_W.R).
*** [sim_analysis_W.R](code/simulation/sim_analysis_W.R): R-Script to reproduce the simulation while keeping all outcomes of kNNDM and establishing an link between (CV estimated - True) ~ W.
*** [sim_functions_W.R](code/simulation/sim_functions_W.R): contains the functions called from [sim_analysis_W.R](code/simulation/sim_analysis_W.R).
*** [knndm_W.R](code/simulation/knndm_W.R): contains the modified kNNDM-function used in [sim_analysis_W.R](code/simulation/sim_analysis_W.R).
** [fig_utils.R](code/figures_utils.R): helper functions used in [figures.Rmd](figures.Rmd).

** [case_study](code/case_study/): contains the R-code to reproduce the case study.
*** [reproduce_sla.R](code/case_study/reproduce_sla.R): R-Script to reproduce the case study.
*** [reproduce_utils.R](code/case_study/reproduce_utils.R): contains additional helper functions used in [reproduce_sla.R](code/case_study/reproduce_sla.R).

* [data](data/): contains the data used in the analyses.
** [simulation](data/simulation): contains the data used in the simulation.
** [case_study](data/case_study): contains the data used in the case study.

* [results](results/): contains the generated results
** [sim_res.csv](results/sim_res.csv): the results from the simulation.
** [sim_res_W.csv](results/sim_res_W.csv): the results from the modified simulation.
** [sla_mod.Rdata](results/sla_mod.Rdata): the reproduced model and CV-split.

* [figures](figures/): contains the figures.

* [figures.Rmd](figures.Rmd): R-Markdown file to reproduce the figures.
