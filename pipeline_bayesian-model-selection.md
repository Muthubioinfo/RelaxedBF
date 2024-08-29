## Pipeline for Bayesian model selection


## Step-1: Estimate parameters at the maximum likelihood and generate into hessian and gradient matrix


Create a main directory.

```
main=/LOCATION-OF-MAIN-DIRECTORY
```

Create a subdirectory called 'files' within the main directory. Here, keep all the input files required for this analysis. This includes ```sequence_alignment.phy```, ```calibrated_tree.trees```, ```mcmctree.ctl``` and ```mcmctree_outbv.ctl``` . 

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
ctl=$source/mcmctree_outbv.ctl

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
mcmctree mcmctree_outbv.ctl


#####The out.BV is contains the hessian and gradient matrix of parameter estimates
cp out.BV in.BV

###Print out the lnmax estimate in .txt
string=$( tail -n 1 $path/out )
value=${str#* = }
echo $value >> $path/lnmax.txt

printf "done"
```

Step-2: Prepare files for marginal likelihood estimation using approximate likelihood methods - ```ARCSIN``` ```SQRT``` ```LOG``` and ```NT```.

Specify location of ```bvalstep64.R```. This Rfile helps create ```n``` directories (for n beta points) to run mcmctree to estimate marginal likelihood at ```n``` stepping stones. In all the ```n``` directories, ```mcmctree.ctl``` file is created specified with a prior for each ```$\beta``` point. 

```
bvalueR=$source/bvalstep64.R

cd $main
mkdir arcsin sqrt log nt

cd $main/arcsin
mkdir iln gbm
cd $main/sqrt
mkdir iln gbm
cd $main/log
mkdir iln gbm
cd $main/nt
mkdir iln gbm

seq=$source/SEQUENCE-FILEtxt
tree=$source/TREE-FILE.tre
inbv=$source/in.BV
ctl_iln_arcsin=$source/mcmctree_iln_arcsin.ctl
ctl_iln_sqrt=$source/mcmctree_iln_sqrt.ctl
ctl_iln_log=$source/mcmctree_iln_log.ctl
ctl_iln_nt=$source/mcmctree_iln_nt.ctl

ctl_gbm_arcsin=$source/mcmctree_gbm_arcsin.ctl
ctl_gbm_sqrt=$source/mcmctree_gbm_sqrt.ctl
ctl_gbm_log=$source/mcmctree_gbm_log.ctl
ctl_gbm_nt=$source/mcmctree_gbm_nt.ctl

###Executable file - mcmctree
mcmctree=/data/home/btx709/PAML_programs/mcmctree

###The working directory to run 
pathiln_arcsin=$main/arcsin/iln
pathiln_sqrt=$main/sqrt/iln
pathiln_log=$main/log/iln
pathiln_nt=$main/nt/iln

pathgbm_arcsin=$main/arcsin/gbm
pathgbm_sqrt=$main/sqrt/gbm
pathgbm_log=$main/log/gbm
pathgbm_nt=$main/nt/gbm

###Step 3 - Copy control files
cp $ctl_iln_arcsin $pathiln_arcsin/mcmctree.ctl
cp $ctl_iln_sqrt $pathiln_sqrt/mcmctree.ctl
cp $ctl_iln_log $pathiln_log/mcmctree.ctl
cp $ctl_iln_nt $pathiln_nt/mcmctree.ctl

cp $ctl_gbm_arcsin $pathgbm_arcsin/mcmctree.ctl
cp $ctl_gbm_sqrt $pathgbm_sqrt/mcmctree.ctl
cp $ctl_gbm_log $pathgbm_log/mcmctree.ctl
cp $ctl_gbm_nt $pathgbm_nt/mcmctree.ctl

ln -s $bvalueR $pathiln_arcsin
ln -s $bvalueR $pathiln_sqrt
ln -s $bvalueR $pathiln_log
ln -s $bvalueR $pathiln_nt
ln -s $bvalueR $pathgbm_arcsin
ln -s $bvalueR $pathgbm_sqrt
ln -s $bvalueR $pathgbm_log
ln -s $bvalueR $pathgbm_nt


cd $pathiln_arcsin
Rscript bvalstep64.R
cd $pathiln_sqrt
Rscript	bvalstep64.R
cd $pathiln_log
Rscript bvalstep64.R
cd $pathiln_nt
Rscript bvalstep64.R

cd $pathgbm_arcsin
Rscript bvalstep64.R
cd $pathgbm_sqrt
Rscript bvalstep64.R
cd $pathgbm_log
Rscript bvalstep64.R
cd $pathgbm_nt
Rscript bvalstep64.R
```

Now copy or link files to prepare files for marginal likelihood estimation

```
for i in {1..64}
do
ln -s $seq $pathiln_arcsin/$i
ln -s $seq $pathiln_sqrt/$i
ln -s $seq $pathiln_log/$i
ln -s $seq $pathiln_nt/$i

ln -s $tree $pathiln_arcsin/$i
ln -s $tree $pathiln_sqrt/$i
ln -s $tree $pathiln_log/$i
ln -s $tree $pathiln_nt/$i

ln -s $seq $pathgbm_arcsin/$i
ln -s $seq $pathgbm_sqrt/$i
ln -s $seq $pathgbm_log/$i
ln -s $seq $pathgbm_nt/$i

ln -s $tree $pathgbm_arcsin/$i
ln -s $tree $pathgbm_sqrt/$i
ln -s $tree $pathgbm_log/$i
ln -s $tree $pathgbm_nt/$i

ln -s $inbv $pathiln_arcsin/$i
ln -s $inbv $pathiln_sqrt/$i
ln -s $inbv $pathiln_log/$i
ln -s $inbv $pathiln_nt/$i

ln -s $inbv $pathgbm_arcsin/$i
ln -s $inbv $pathgbm_sqrt/$i
ln -s $inbv $pathgbm_log/$i
ln -s $inbv $pathgbm_nt/$i

cp $mcmctree $pathiln_arcsin/$i
cp $mcmctree $pathiln_sqrt/$i
cp $mcmctree $pathiln_log/$i
cp $mcmctree $pathiln_nt/$i

cp $mcmctree $pathgbm_arcsin/$i
cp $mcmctree $pathgbm_sqrt/$i
cp $mcmctree $pathgbm_log/$i
cp $mcmctree $pathgbm_nt/$i
done

printf "done"

```

Step-3: Run ```MCMCTree``` to estimate marginal likelihood under the relaxed clock models - strict clock (STR), independent log-normal (ILN) and geometric brownian motion (GBM). 
Note that likelihood approximations are less to estimate marginal likelihood under strick clock models. Hence, use approximate likelihood for relaxed clock models only.

Use array jobs to run ```MCMCTree``` to efficiently save computational time.

```
#########
cd $pathstr_arcsin/${SGE_TASK_ID}
./mcmctree
cd $pathstr_sqrt/${SGE_TASK_ID}
./mcmctree
cd $pathstr_log/${SGE_TASK_ID}
./mcmctree
cd $pathstr_nt/${SGE_TASK_ID}
./mcmctree

cd $pathiln_arcsin/${SGE_TASK_ID}
./mcmctree
cd $pathiln_sqrt/${SGE_TASK_ID}
./mcmctree
cd $pathiln_log/${SGE_TASK_ID}
./mcmctree
cd $pathiln_nt/${SGE_TASK_ID}
./mcmctree


cd $pathgbm_arcsin/${SGE_TASK_ID}
./mcmctree
cd $pathgbn_sqrt/${SGE_TASK_ID}
./mcmctree
cd $pathgbm_log/${SGE_TASK_ID}
./mcmctree
cd $pathiln_nt/${SGE_TASK_ID}
./mcmctree
```
Finally, run the following Rscript to compute the log-marginal likelihood and standard deviation for each clock/relaxed model. Then calculate bayes factor and posterior model probabilties.

```
cd $pathiln_arcsin
iln <- mcmc3r::stepping.stones()
iln$logml; iln$se

cd $pathgbm_arcsin
gbm <- mcmc3r::stepping.stones()
gbm$logml; gbm$se

BF_iln <- exp(  iln$logml - max( c(iln$logml, gbm$logml) )
BF_gbm <- exp(  gbm$logml - max( c(iln$logml, gbm$logml) )

posterior_prob_iln <- BF_iln / sum( c(iln$logml, gbm$logml) )
posterior_prob_gbm <- BF_gbm / sum( c(iln$logml, gbm$logml) )

```
The model with largest posterior probability is consideredn to be the chosen model for the given observed data.





