seed = -1
seqfile = sim.phy
treefile = treecal.tree
outfile = out

ndata = 2      * number of partitions
usedata = 1    * 0: no data; 1:seq like; 2:use in.BV; 3: out.BV
clock = 1      * 1: global clock; 2: independent rates; 3: correlated rates

model = 0    * 0:JC69, 1:K80, 2:F81, 3:F84, 4:HKY85
alpha = 0   * alpha for gamma rates at sites
ncatG = 0    * No. categories in discrete gamma

BDparas = 1 1 0.1      * birth, death, sampling
*kappa_gamma = 6 2   * gamma prior for kappa
*alpha_gamma = 1 1    * gamma prior for alpha

rgene_gamma = 2 2   * gamma prior for mean rates for genes
sigma2_gamma = 2 2  * gamma prior for sigma^2 (for clock=2 or 3)

print = 1
burnin = 10000
sampfreq = 5
nsample = 18000
BayesFactorBeta = 1e-300
