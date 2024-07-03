This folder contains the input files used for simulating rates along the phylogenetic trees and 
generate sequence alignments under the different rate models - Strict clock, Independent log-normal and geometric brownian motion.
The results for this analysis can be found in "simulation_results.zip".

The "tree_fig1_sim.trees" is obtained from figure 3 of Yang and Rannala (2006), "Bayesian Estimation of Species Divergence 
Times Under a Molecular Clock Using Multiple Fossil Calibrations with Soft Bounds". This "tree_fig1_sim.trees" is used as input to
to the R scripts "simulate_trees_..._loci.R" to simulate rates in the phylogeny and generate 
phylogenetic trees informed with calculated branch lengths under different models. 

The simulated phylogenetic trees can then be used in "MCbase.dat" 
to generate sequence alignments using the evolver program in PAML package.

The ctl_files contains the control files for 'MCMCTree' program. 
This is used to carry out Bayesian model selection using marginal likelihood estimation.
See tutorial of marginal likelihood estimation in "https://dosreislab.github.io/2017/10/24/marginal-likelihood-mcmc3r.html"

For each rate model setting (i.e., str_nn, iln_sv, iln_cl, gbm_sv and gbm_cl), there are 1-8 stepping-stones with one control file 
("mcmctree.ctl") for Bayesian MCMCTree analysis from PAML.
