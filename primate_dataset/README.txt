#############################################
####PRIMATE DATASET


The dataset of 372 primates obtained from dos Reis et al. 2018. 

Reis, M.D., Gunnell, G.F., Barba-Montoya, J., Wilkins, A., Yang, Z. and Yoder, A.D., 2018. Using phylogenomic data to explore the effects of relaxed clocks and calibration strategies on divergence time estimation: primates as a test case. Systematic Biology, 67(4), pp.594-615.

The given dataset has six alignment partitions.

Partition 1 - Mitochondrial genes 1st + 2nd codon position
Partition 2 - Mitochondrial genes 3rd codon position
Partition 3 - Mitochondrial RNA
Partition 4 - Nuclear 1st + 2nd codon sites
Partition 5 - Nuclear 3rd codon site
Partition 6 - Nuclear noncoding segments (UTR and introns)


The ".phy" file represents the alignment files, and ".tree" files are rerooted trees for each partition when no calibrations ( t1 ~ 'B(0.999,1.001)' ) are considered in the analysis. And "_gam1000.trees" is used for analysis with narrow density calibration ( t1 ~ G(1000,1000) ) . Also have enclosed the ".ctl" files for each likelihood approximate method.

