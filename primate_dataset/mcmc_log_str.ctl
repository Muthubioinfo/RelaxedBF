seed = -1
seqfile = align.phy
treefile = tree.trees
outfile = out

ndata = 1      * number of partitions
usedata = 2 ./in.BV 2  * Use in.BV (0: NOTRANSFORM, 1: SQRT, 2: LOG, 3: ARCSIN) 
clock = 1  * 1: global clock; 2: independent rates; 3: correlated rates

model = 4    * 0:JC69, 1:K80, 2:F81, 3:F84, 4:HKY85
alpha = .5   * alpha for gamma rates at sites
ncatG = 5    * No. categories in discrete gamma

cleandata = 0    * remove sites with ambiguity data (1:yes, 0:no)?

BDparas = 1 1 0      * birth, death, sampling
kappa_gamma = 2 .2   * gamma prior for kappa
alpha_gamma = 2 4    * gamma prior for alpha

rgene_gamma = 2 20   * gamma prior for mean rates for genes
sigma2_gamma = 1 10  * gamma prior for sigma^2 (for clock=2 or 3)

print = 1
burnin = 4000
sampfreq = 5
nsample = 20000
