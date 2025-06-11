# Bayesian rate model selection
This repository provides the pipeline to study **Bayesian selection of relaxed-clock models**.
Use ```pipeline_bayesian-model-selection.md``` to run Bayesian model selection for your sequence alignment.


# 📦 Requirements

Please install the following tools:
- [R / RStudio](https://cran.r-project.org/)
- [PAML package (EVOLVER, BASEML, MCMCTree)](http://abacus.gene.ucl.ac.uk/software/paml.html)
- UNIX/Linux environment for bash scripting

<br/>

## Simulating Phylogenetic alignments under the clock and Relaxed-clock models

## 🧬 Step 1: Simulating phylogenetic trees under relaxed-clock models:

In the first step, sample the parameters rates, ```r``` and rate drift, ``` $\sigma^2$ ``` from gamma distribution for _L_ number of Loci.
Then, for a given input tree topology called "tree_fig1_sim.tree" (obtained from figure 3 of the [Rannala and Yang (2007)](https://academic.oup.com/sysbio/article/56/3/453/1657118)), simulate phylogenies with branch lengths under different rate model settings: 

* Strict clock (```STR```)
* Independent-log normal seriously violated rates (```ILN-SV```)
* Independent-log normal with clock-like variation (```ILN-CL```)
* Geometric brownian motion serious violated (```GBM-SV```)
* Geometric brownian motion clock-like (```GBM-CL```)

Run the Rscript ```simulate_trees.R```, to generate output of phylogenies,

```
Rscript simulate_trees.R

```

Three different Loci, _L_ = 1, 2 and 5 are simulated each under five rate model configurations (STR, ILN-CL, ILNSV, GBM-CL, GBM-SV), we get 
3 x 5 = 15 trees. As in (paper)[https://academic.oup.com/sysbio/article-abstract/74/2/323/7906181?redirectedFrom=fulltext&login=false] For repeating simulation of _N_ = 1000, we get a total of 15,000 trees as output.

## 🧬 Step 2: Simulating Phylogenetic Alignments
Then, simulate nucleotide alignment for each of the 15,000 simulated trees using the ```EVOLVER``` program from ```PAML``` using the [MCbase.dat](https://github.com/Muthubioinfo/RelaxedBF/blob/main/simulation_files/MCbase.dat) control file. The nucleotide substitution model is JC69 model (Jukes and Cantor, 1969), assuming uniform codon frequencies. The output of ```EVOLVER``` generates nucleotide sequence alignments corresponding for each of the 1,000 trees under each rate model setting.

<br/>

## Step 3: Bayesian rate model selection
For Bayesian model selection of rate models, estimate the marginal likelihood of each true rate model (STR, ILN and GBM)

For example, if the simulations are under ILN-SV, then marginal likelihoods should be estimated under STR, ILN-SV and GBM-SV. 
And if the simulations are under ILN-SV, then marginal likelihoods should be estimated under STR, ILN-SV and GBM-SV.
The directory ["simulation_files"](https://github.com/Muthubioinfo/RelaxedBF/tree/main/simulation_files) contains all the input files used for MCMCTree analysis. This includes [mcmctree.ctl](https://github.com/Muthubioinfo/RelaxedBF/tree/main/simulation_files/ctl_files_mcmctree),
Use the appropriate [control files](https://github.com/Muthubioinfo/RelaxedBF/tree/main/simulation_files/ctl_files_mcmctree) for each rate model setting. The [calibrated.tree](https://github.com/Muthubioinfo/RelaxedBF/tree/main/simulation_files/calibrated_trees) files are used for analysing simulations under different the following calibrations, (see figure 1 in [paper](https://academic.oup.com/sysbio/article-abstract/74/2/323/7906181?redirectedFrom=fulltext&login=false) :

- ( [slightly misspecified calibrations](https://github.com/Muthubioinfo/RelaxedBF/blob/main/simulation_files/calibrated_trees/slightly_misp_cal.tree),
- [badly misspecified calibrations](https://github.com/Muthubioinfo/RelaxedBF/blob/main/simulation_files/calibrated_trees/bad_misp_cal.tree),
- [no calibration](https://github.com/Muthubioinfo/RelaxedBF/blob/main/simulation_files/calibrated_trees/no_fossil_cal.tree). 

<br/>

### Simulation results 

The supplementary results from the simulation analysis are in [simulation_files](https://github.com/Muthubioinfo/RelaxedBF/tree/main/simulation_files/simulation_results).

The [Table_simulation_analysis_full_headings.xlsx](https://github.com/Muthubioinfo/RelaxedBF/blob/main/simulation_files/simulation_results/Table_simulation_analysis_full_headings.xlsx) contains the percentage of simulations that are identifiable under the true model, among the total number of 1,000 simulations. The [no_calibrations](https://github.com/Muthubioinfo/RelaxedBF/tree/main/simulation_files/simulation_results/no_calibrations), [misspecied_calibration](https://github.com/Muthubioinfo/RelaxedBF/tree/main/simulation_files/simulation_results/misspecied_calibration), [bad_calibrations](https://github.com/Muthubioinfo/RelaxedBF/tree/main/simulation_files/simulation_results/bad_calibration) contains results of posterior model probabilities and marginal likelihoods (minimum, median and maximum value among 1,000 simulations) analysed under all the Relaxed-clock model settings. 

<br/>

## Approximating marginal likelihood using real datasets

Likelihood approximation methods can provide a more faster and efficient way to estimate marginal likelihood for Bayesian model selection. There are four different methods to approximate likelihood based on the parameter transformations using ```ARCSIN```, ```SQRT```, ```LOG``` and ```NT```. Best approaches are ```ARCSIN``` and ```SQRT```, which provide the best approximations. See [dos Reis et al. 2011](https://academic.oup.com/mbe/article/28/7/2161/1051613) for more details on the applications and implementation of these approximations in MCMCTree program.

Three empirical datasets were tested here -

### Small dataset 

The "[small_dataset](https://github.com/Muthubioinfo/RelaxedBF/tree/main/small_dataset)" directory contains 13 alignment files used to test Bayesian model selection to compare exact likelihood and approximate likelihood methods. The first three alignments align1.phy, align2.phy and align3.phy in "Mito_primate_1-3.zip" is obtained from [Yang and Rannala (2006)](https://academic.oup.com/mbe/article/23/1/212/1193630). These are three partitioned alignments of seven ape mitochondrial genomes. The remaining ten alignments are align4.phy, align5.phy, ... and align13.phy are the the ten protein-coding (first and second codon positions) gene alignments of 72 mammalian genomes obtained from [Álvarez-Carretero et al. 2022](https://www.nature.com/articles/s41586-021-04341-1). The gene names of the protein-coding genes from reference human genome are provided below along with their ensemble IDs. 

Alignments obtained from [Yang and Rannala (2006)](https://academic.oup.com/mbe/article/23/1/212/1193630) consisting of seven ape mitochondrial genomes partitioned based on the three codon positions (c.p.) 

| Alignment (.phy) |  Name/Ensembl ID and codon position (c.p.) | 
| ---------------  | ------------------------------------------ |      
| align1.phy       |	Mitochondrial genome (1st c.p.)           |   
| align2.phy       |	Mitochondrial genome (2nd c.p.)           | 
| align3.phy       |	Mitochondrial genome (3rd c.p.)           |  

Dataset of ten protein-coding gene alignments (1st and 2nd codon positions) obtained from [Álvarez-Carretero et al. (2022)](https://www.nature.com/articles/s41586-021-04341-1) that consists of 72 mammal genomes. 


| Alignment (.phy) |  Name/Ensembl ID and codon position (c.p.) | Human reference gene names |
| ---------------- | ------------------------------------------ | -------------------------- |
| align4.phy 	     |  ENSG00000059588 (1st + 2nd c.p.)          | 		      TARBP1           | 
| align5.phy 	     |  ENSG00000103534 (1st + 2nd c.p.)          |		        TMC5             |  
| align6.phy 	     |  ENSG00000112818 (1st + 2nd c.p.)          |		        MEP1A            | 
| align7.phy 	     |  ENSG00000119684 (1st + 2nd c.p.)          |		        MLH3             |  
| align8.phy 	     |  ENSG00000130413 (1st + 2nd c.p.)          |		        STK33            |  
| align9.phy 	     |  ENSG00000134222 (1st + 2nd c.p.)          |		        PSRC1            |  
| align10.phy      |	ENSG00000136634 (1st + 2nd c.p.)          |		        IL10             | 
| align11.phy      |	ENSG00000162994 (1st + 2nd c.p.)          |		        CLHC1            |  
| align12.phy      |	ENSG00000165392 (1st + 2nd c.p.)          |		         WRN             |  
| align13.phy      |	ENSG00000182010 (1st + 2nd c.p.)          |		        RTKN2            |  


For the MCMCTree analysis for marginal likelihood estimation, the ```.trees``` and ```.ctl``` files is in ```Mito_primate_1-3.zip``` and ```mammal_dataset_4-13.zip``` files. 

### Primate_dataset 

The primate dataset are obtained from [dos Reis et al. 2018](https://academic.oup.com/sysbio/article/67/4/594/4802240) that consisting of six alignment partitions each with around 372 primates 

The given dataset has six alignment partitions. 

| Partitions (.phy)  |   Number of species  | Dataset                                      |
| ------------------ | -------------------- | -------------------------------------------- |
| Partition 1        |        330           | Mitochondrial genes 1st + 2nd codon position |
| Partition 2        |        330           | Mitochondrial genes 3rd codon position       |
| Partition 3        |        220           | Mitochondrial RNA                            |
| Partition 4        |        239           | Nuclear 1st + 2nd codon positions            |
| Partition 5        |        239           | Nuclear 3rd codon positions                  |
| Partition 6        |        220           | Nuclear noncoding segments (UTR and introns) |

For Bayesian model selection analysis, users should analyse each alignment (```.phy```) with the corresponding ```.trees``` calibrated with ```B(0.999,1.001)``` at the root branch. These files are in the "primate_dataset" directory. The control files for MCMCTree is also found (```.ctl```) for each likelihood approximation transformation (```ARCSIN```, ```SQRT```, ```LOG``` and ```NT```).


### Plant dataset

The [plant dataset](https://github.com/Muthubioinfo/RelaxedBF/tree/main/plant_dataset) is obtained from [Barba-Montoya et al. 2018](https://nph.onlinelibrary.wiley.com/doi/10.1111/nph.15011). This is the largest dataset that has been tested for Bayesian model selection with 644 flowering plants including tracheophytes and Angiosperms.


From the given dataset, there are three alignment partitions.

| Partitions (.phy)  | Number of species  | Dataset                                      | 
| ------------------ | ------------------ | -------------------------------------------- |
| Partition 1        |        643         | Plastid 1st and 2nd codon position           |
| Partition 2        |        515         | Mitochondrial 1st and 2nd codon position     |
| Partition 3        |        540         | Nuclear RNA                                  |


The three partitioned alignments are ```partition1.phy```, ```partition2.phy``` and ```partition3.phy``` files with ```.trees``` calibrated with B(0.999,1.001) at the root. Use the ```.ctl``` to run MCMCTree program to analyse for Bayesian model selection.

## References
                                                  
Barba‐Montoya, J., Dos Reis, M., Schneider, H., Donoghue, P.C. and Yang, Z., 2018. Constraining uncertainty in the timescale of angiosperm evolution and the veracity of a Cretaceous Terrestrial Revolution. New Phytologist, 218(2), pp.819-834.

Rannala, B. and Yang, Z., 2007. Inferring speciation times under an episodic molecular clock. Systematic biology, 56(3), pp.453-466.

Reis, M.D., Gunnell, G.F., Barba-Montoya, J., Wilkins, A., Yang, Z. and Yoder, A.D., 2018. Using phylogenomic data to explore the effects of relaxed clocks and calibration strategies on divergence time estimation: primates as a test case. Systematic Biology, 67(4), pp.594-615.

Yang, Z. and Rannala, B., 2006. Bayesian estimation of species divergence times under a molecular clock using multiple fossil calibrations with soft bounds. Molecular biology and evolution, 23(1), pp.212-226.

