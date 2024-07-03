# Simulation under rate models


This folder contains the input files used for simulating rates along the phylogenetic trees and 
generate sequence alignments under the different rate models - Strict clock (STR), Independent log-normal (ILN) and geometric brownian motion (GBM).
The results for this analysis can be found in "simulation_results.zip".

The "tree_fig1_sim.trees" is obtained from figure 3 of Yang and Rannala (2006), "Bayesian Estimation of Species Divergence 
Times Under a Molecular Clock Using Multiple Fossil Calibrations with Soft Bounds". This "tree_fig1_sim.trees" is used as input to
to "simulate_trees_..._loci.R" to simulate rates in the phylogeny under different loci settings, loci = 1, 2, and 5, and generate phylogenetic trees informed with branch lengths under different rate models. These R scripts uses 'simclock' package available in 
https://github.com/dosreislab/simclock.

The simulated phylogenetic trees can then be used in "MCbase.dat" to generate molecular sequence alignments using the evolver program in PAML package. Then, these alignments can be used to estimate the marginal likelihood under each rate model (STR, ILN and GBM). 
See tutorial of marginal likelihood estimation in https://dosreislab.github.io/2017/10/24/marginal-likelihood-mcmc3r.html

The ctl_files contains the control files for 'MCMCTree' program to carry out the marginal likelihood estimation.  
For each rate model setting (i.e., str_nn, iln_sv, iln_cl, gbm_sv and gbm_cl), there are 1-8 stepping-stones with one control file 
("mcmctree.ctl") for Bayesian MCMCTree analysis from PAML.
