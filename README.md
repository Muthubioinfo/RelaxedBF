# Bayesian clock model selection
This repository contains the files used to recreate all the analysis from Bayesian Selection of Relaxed-clock model selection paper.

The users are required to have R/Rstudio and install EVOLVER, BASEML and MCMCTree programs available in PAML package, see [http://abacus.gene.ucl.ac.uk/software/paml.html](http://abacus.gene.ucl.ac.uk/software/paml.html). Use bash/LINUX commands to run PAML programs.

#  Simulating Phylogenetic alignments under the clock and Relaxed-clock models: 

Use "simulate_trees.R" to simulate phylogenies with branch lengths under the different rate model settings (shown below), for a given tree topology called "tree_fig1_sim.tree".
This input phylogeny was obtained from figure 3 of the Rannala and Yang (2007) paper. 

(1) Strict clock (STR), 
(2) Independent-log normal seriously violated rates (ILN-SV),
(3) Independent-log normal with clock-like variation (ILN-CL),
(4) Geometric brownian motion serious violated (GBM-SV), 
(5) Geometric brownian motion clock-like (GBM-CL). 


The directory ["simulation_files"](https://github.com/Muthubioinfo/RelaxedBF/tree/main/simulation_files) contains the input files used for MCMCTree analysis. This includes mcmctree.ctl](https://github.com/Muthubioinfo/RelaxedBF/tree/main/simulation_files/ctl_files_mcmctree/iln_sv) files, and the '[calibrated.tree](https://github.com/Muthubioinfo/RelaxedBF/tree/main/simulation_files/calibrated_trees)' is used for analysing under different calibration setting ( [slightly misspecified calibrations](https://github.com/Muthubioinfo/RelaxedBF/blob/main/simulation_files/calibrated_trees/slightly_misp_cal.tree), [badly misspecified calibrations](https://github.com/Muthubioinfo/RelaxedBF/blob/main/simulation_files/calibrated_trees/bad_misp_cal.tree), and under [no calibration](https://github.com/Muthubioinfo/RelaxedBF/blob/main/simulation_files/calibrated_trees/no_fossil_cal.tree) ). The "simulation_files/simulation_results" contains the result tables for each calibration setting. 

#  Simulation results 

The [Table_simulation_analysis_full_headings.xlsx](https://github.com/Muthubioinfo/RelaxedBF/blob/main/simulation_files/simulation_results/Table_simulation_analysis_full_headings.xlsx) contains the percentage of simulations that are identifiable under the true model, among the total number of 1000 simulations. 

(2) Real dataset: The small_dataset, primate_dataset and plant_dataset.zip contains the input files used for marginal likelihood estimation using approximation for faster analysis.



                                                  
