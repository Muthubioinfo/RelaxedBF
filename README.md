# Bayesian clock model selection
Files for Bayesian clock/relaxed clock model selection paper

(1) The "simulate_trees.R" is used to simulate phylogenies with branch lengths under the different rate models - Strict clock, Independent-log normal seriously violated or clock-like (ILN-SV and ILN-CL) and geometric brownian motion serious violated or clock-like (GBM-SV and GBM-CL). An input tree "tree_fig1_sim.tree" is called as an input to "simulate_tree.R".

The simulation_files contains the files such as the '.ctl' files used for MCMCTree analysis (see ./ctl_files_mcmctree) and the 'calibrated_trees' for the fossil calibration settings is used for each rate model (slightly misspecified, badly misspecified, and no calibration). In "simulation_files/simulation_results", the tables with the results is provided for each calibration setting. 
The "Table_simulation_analysis_full_headings.xlsx" contains the percentage of simulations that are identifiable under the true model analysis, among the total number of simulations (i.e., among 1000 simulations). 

(2) Real dataset: The small_dataset, primate_dataset and plant_dataset.zip contains the input files used for marginal likelihood estimation using approximation for faster analysis.



                                                  
