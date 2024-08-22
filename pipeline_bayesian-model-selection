## Pipeline for Bayesian model selection


```
#Create source directory with the required files such as sequence alignment, newick tree topology with fossil calibrations (if any), inbv file
#output from previous step, the control files for the models and the mcmctree executable file in the compiled location.
main=/data/scratch/btx709/viral_phy/TipDate.HIV2

##Directory where all the input files are created and stored
source=$main/files

####Prepare files in directory
##Specify location of Rfile that creates n directories (for n beta points) and pastes 'mcmctree.ctl' files respectively. 
#Each control file has beta-prior specified at the end line.

cd $main
mkdir inBV

####Specify the source and file names
seq=$source/HIV2ge.txt
tree=$source/HIV2ge.tre
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

###Execute command - estimate lnmax and convert parameter estimates into a hessian and gradient matrix
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




