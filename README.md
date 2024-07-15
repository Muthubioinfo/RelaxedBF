# Bayesian clock model selection
This repository contains the files used to recreate all the analysis in Bayesian clock/relaxed clock model selection paper.

The user is required to have R/Rstudio and install EVOLVER, BASEML and MCMCTree programs available in PAML package (see [http://abacus.gene.ucl.ac.uk/software/paml.html](http://abacus.gene.ucl.ac.uk/software/paml.html). Use bash/LINUX commands to run PAML programs.

(1) The "simulate_trees.R" is used to simulate phylogenies with branch lengths under the different rate models - Strict clock (STR), Independent-log normal seriously violated or clock-like (ILN-SV or ILN-CL) and geometric brownian motion serious violated or clock-like (GBM-SV or GBM-CL). The "tree_fig1_sim.tree" is called as an input to "simulate_tree.R".

The "simulation_files" contains the files such as the '.ctl' files used for MCMCTree analysis (see ./ctl_files_mcmctree), and the "simulation_files/calibrated_trees" is used for analysing under different calibration setting (no fossil calibrations, slightly misspecified calibrations, badly misspecified). The "simulation_files/simulation_results" contains the result tables for each calibration setting. 

"simulation_files/simulation_results/Table_simulation_analysis_full_headings.xlsx" contains the percentage of simulations that are identifiable under the true model, among the total number of 1000 simulations. 

(2) Real dataset: The small_dataset, primate_dataset and plant_dataset.zip contains the input files used for marginal likelihood estimation using approximation for faster analysis.



                                                  
