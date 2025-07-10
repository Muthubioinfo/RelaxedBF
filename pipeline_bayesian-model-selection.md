# Pipeline for Bayesian relaxed clock model selection
This is a pipeline to determine appropriate clock model (strict clock, independent log-normal and geometric brownian motion) suitable to your sequence alignment for Bayesian divergence time estimation. 

## Requirements
1) Sequence alignment, e.g., `sequence_alignment.phy`
2) Tree topology obtained from output of IQTREE, RAxML or other softwares used for phylogeny reconstruction. Add prior time calibrations to your tree file.
3) MCMCTree control files. You need to create three control files, with each file setting to `clock =   1` or `2` or `3`. These control files can be named as `mcmctree_str.ctl`, `mcmctree_iln.ctl` and `mcmctree_gbm.ctl` respectively.
 You will need to create another control file called `mcmctree_outbv.ctl`, where `usedata = 2`. This   file is option is used command MCMCTree to generate `out.BV` with hessian and gradient matrix (See     step-1).
5) [mcmc3r R package](https://github.com/dosreislab/mcmc3r)
  Install in R by typing `devtools::install_github ("dosreislab/mcmc3r")`
6) R script called `bvalstep64.R`, see Step-2.
7) [PAML programs - BASEMLand MCMCTree)](http://abacus.gene.ucl.ac.uk/software/paml.html)
8) UNIX/Linux environment for bash scripting and High-performance cluster computing


## Step-1: Estimate parameters at the maximum likelihood and generate into hessian and gradient matrix

Create a main directory.

```
main=/LOCATION-OF-MAIN-DIRECTORY
```

Create a subdirectory called 'files' within the main directory. Here, keep all the input files required for this analysis. This includes ```sequence_alignment.phy```, ```calibrated_tree.trees```, 
and ```mcmctree.ctl```. Also create `mcmctree_str.ctl`, `mcmctree_iln.ctl`, `mcmctree_gbm.ctl`. 

```
source=$main/files
```

Specify location of Rfile to creates n directories (for n beta points) and pastes 'mcmctree.ctl' files respectively. Each control file has beta-prior specified at the end line.

Specify the source files

```
cd $main
mkdir inBV

seq=$source/SEQUENCE-ALIGNMENT-TO-TEST-FOR-MODEL-SELECTION.TXT
tree=$source/LOCATION-OF-TREE-FILE.trees
ctl=$source/mcmctree.ctl

###Executable file - mcmctree
mcmctree=/data/home/btx709/PAML_programs/mcmctree

###Executable file - mcmctree
baseml=/data/home/btx709/PAML_programs/baseml

###Path to target directory
path=$main/inBV

###Copy files
ln -s $seq $path
ln -s $tree $path
ln -s $ctl $path
cp $mcmctree $path
cp $baseml $path
```

Execute command to estimate lnmax and convert parameter estimates into a hessian and gradient matrix

```
cd $path
export PATH=$PATH:$path
mcmctree mcmctree.ctl

#####The out.BV is contains the hessian and gradient matrix of parameter estimates
cp out.BV in.BV

###Print out the lnmax estimate in .txt
string=$( tail -n 1 $path/out )
value=${str#* = }
echo $value >> $path/lnmax.txt

printf "done"
```

### Step-2: Prepare files for marginal likelihood estimation using approximate likelihood methods

The R script ```bvalstep64.R``` helps to create ```n``` directories (for n beta points) to run MCMCTree across ```n``` stepping stones. In all the ```n``` directories, ```mcmctree.ctl``` file is created specified with a prior for each ```$\beta``` point. 

```<\r>
library('mcmc3r')
b <- mcmc3r::make.beta(n=64, a=5, method="step-stones")

#Create n directories containing .ctl files with the Bayesfactors b
mcmc3r::make.bfctlf(b, ctlf="mcmctree.ctl", betaf="beta.txt") 

#This will tell the MCMCTree to sample from the power posterior with b value
#Note that MCMCTree currently cannot sample log-likelihoods using b = 0, and 
#so b = 10e–300 (a tiny number) is used instead in 1/mcmctree.ctl
```

Run the following script 

```<\bash>
bvalueR=$source/bvalstep64.R

cd $main
mkdir arcsin

cd $main/arcsin
mkdir iln gbm

seq=$source/SEQUENCE-FILEtxt
tree=$source/TREE-FILE.tree
inbv=$source/in.BV
ctl_iln_arcsin=$source/mcmctree_iln_arcsin.ctl
ctl_gbm_arcsin=$source/mcmctree_gbm_arcsin.ctl

###Executable file - mcmctree
mcmctree=/data/home/btx709/PAML_programs/mcmctree

###The working directory to run 
pathiln_arcsin=$main/arcsin/iln
pathgbm_arcsin=$main/arcsin/gbm

###Step 3 - Copy control files
cp $ctl_iln_arcsin $pathiln_arcsin/mcmctree.ctl
cp $ctl_gbm_arcsin $pathgbm_arcsin/mcmctree.ctl

ln -s $bvalueR $pathiln_arcsin
ln -s $bvalueR $pathgbm_arcsin

cd $pathiln_arcsin
Rscript bvalstep64.R

cd $pathgbm_arcsin
Rscript bvalstep64.R
```

Now copy or link files to prepare files for marginal likelihood estimation

```<\bash>
for i in {1..64}
do
ln -s $seq $pathiln_arcsin/$i
ln -s $tree $pathiln_arcsin/$i
ln -s $seq $pathgbm_arcsin/$i
ln -s $tree $pathgbm_arcsin/$i
ln -s $inbv $pathiln_arcsin/$i
ln -s $inbv $pathgbm_arcsin/$i
cp $mcmctree $pathiln_arcsin/$i
cp $mcmctree $pathgbm_arcsin/$i

done

printf "done"
```

### Step-3: Run ```MCMCTree``` to estimate marginal likelihood under the relaxed clock models - strict clock (STR), independent log-normal (ILN) and geometric brownian motion (GBM). 
Note that likelihood approximations are less to estimate marginal likelihood under strick clock models. Hence, use approximate likelihood for relaxed clock models only.

Use array jobs to run ```MCMCTree``` to efficiently save computational time.

```<\bash>
#########
cd $pathiln_arcsin/${SGE_TASK_ID}
./mcmctree

cd $pathgbm_arcsin/${SGE_TASK_ID}
./mcmctree

```
Finally, run the following Rscripts in `$pathiln_arcsin` and `cd $pathgbm_arcsin` to compute the log-marginal likelihood and standard deviation for each clock/relaxed model. 

In `$pathiln_arcsin`, run

```<\r>
iln <- mcmc3r::stepping.stones()
iln$logml; iln$se
```

In `$pathgbm_arcsin`, run
```<\r>
gbm <- mcmc3r::stepping.stones()
gbm$logml; gbm$se
```

Then calculate Bayes factor and Posterior model probabilties.

```<\r>
BF_iln <- exp(  iln$logml - max( c(iln$logml, gbm$logml) )
BF_gbm <- exp(  gbm$logml - max( c(iln$logml, gbm$logml) )

posterior_prob_iln <- BF_iln / sum( c(iln$logml, gbm$logml) )
posterior_prob_gbm <- BF_gbm / sum( c(iln$logml, gbm$logml) )
```

The model with largest posterior probability is considered to be the chosen model for the given observed data.





