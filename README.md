# Bayesian clock model selection
This repository contains the files used to recreate all the analysis from Bayesian Selection of Relaxed-clock model selection paper.

The users are required to have R/Rstudio and install ```EVOLVER```, ```BASEML``` and ```MCMCTree``` programs available in ```PAML``` package, see [http://abacus.gene.ucl.ac.uk/software/paml.html](http://abacus.gene.ucl.ac.uk/software/paml.html). Use bash/LINUX commands to run PAML programs.

<br/>

## Simulating Phylogenetic alignments under the clock and Relaxed-clock models

Use "simulate_trees.R" to simulate phylogenies with branch lengths under the different rate model settings (shown below), for a given tree topology called "tree_fig1_sim.tree".
This input phylogeny was obtained from figure 3 of the Rannala and Yang (2007) paper. 

* Strict clock (```STR```)
* Independent-log normal seriously violated rates (```ILN-SV```)
* Independent-log normal with clock-like variation (```ILN-CL```)
* Geometric brownian motion serious violated (```GBM-SV```)
* Geometric brownian motion clock-like (```GBM-CL```)

Simulate trees under five rate model configurations (STR, ILN-CL, ILNSV, GBM-CL, GBM-SV) and three locus configurations _L_ = 1, 2 and 5.

Then, simulate nucleotide alignments for the 15,000 simulated trees using the ```EVOLVER``` program from ```PAML``` with the [MCbase.dat](https://github.com/Muthubioinfo/RelaxedBF/blob/main/simulation_files/MCbase.dat) control file. The nucleotide substitution model is JC69 model (Jukes and Cantor, 1969), assuming uniform codon frequencies. The output of ```EVOLVER``` generates nucleotide sequence alignments corresponding for each of the 1,000 trees under each rate model setting.


<br/>


## Bayesian rate model selection and rate prior
The directory ["simulation_files"](https://github.com/Muthubioinfo/RelaxedBF/tree/main/simulation_files) contains the input files used for MCMCTree analysis. This includes [mcmctree.ctl](https://github.com/Muthubioinfo/RelaxedBF/tree/main/simulation_files/ctl_files_mcmctree/iln_sv) and the [calibrated.tree](https://github.com/Muthubioinfo/RelaxedBF/tree/main/simulation_files/calibrated_trees) files are used for analysing under different calibration setting ( [slightly misspecified calibrations](https://github.com/Muthubioinfo/RelaxedBF/blob/main/simulation_files/calibrated_trees/slightly_misp_cal.tree), [badly misspecified calibrations](https://github.com/Muthubioinfo/RelaxedBF/blob/main/simulation_files/calibrated_trees/bad_misp_cal.tree), and under [no calibration](https://github.com/Muthubioinfo/RelaxedBF/blob/main/simulation_files/calibrated_trees/no_fossil_cal.tree). 

<br/>

## Simulation results 

All the supplementary tables are in [simulation_files](https://github.com/Muthubioinfo/RelaxedBF/tree/main/simulation_files/simulation_results).

The [Table_simulation_analysis_full_headings.xlsx](https://github.com/Muthubioinfo/RelaxedBF/blob/main/simulation_files/simulation_results/Table_simulation_analysis_full_headings.xlsx) contains the percentage of simulations that are identifiable under the true model, among the total number of 1,000 simulations. The [no_calibrations](https://github.com/Muthubioinfo/RelaxedBF/tree/main/simulation_files/simulation_results/no_calibrations), [misspecied_calibration](https://github.com/Muthubioinfo/RelaxedBF/tree/main/simulation_files/simulation_results/misspecied_calibration), [bad_calibrations](https://github.com/Muthubioinfo/RelaxedBF/tree/main/simulation_files/simulation_results/bad_calibration) contains results of posterior model probabilities and marginal likelihoods (minimum, median and maximum value among 1,000 simulations) analysed under all the Relaxed-clock model settings. 

<br/>

## Real dataset

The small_dataset, primate_dataset and plant_dataset.zip contains the input files used for approximating marginal likelihood under ```ARCSIN```, ```SQRT```, ```LOG``` and ```NT``` for faster Bayesian selection analysis. Best approaches are ```ARCSIN``` and ```SQRT```, which provide the best approximations. 

                                                  
