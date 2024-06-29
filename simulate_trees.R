
#' Rscript to simulate rates to simulate rates on a phylogeny under 
#' the rate models -
#' 
#' Strict clock - STR
#' Independent log-normal - ILN
#' Geometric brownian motion - GBM
#' 
#' 
#' @include tree_fig1_sim.tree
#' 
#' @param tree an object of class phylo
#' @param N numeric, the number of nucleotide sequence length
#' @param L numeric, the number of loci/partition
#' @param model character, the relaxed clock model
#' @param rate numeric, the mean rate in substitutions per site
#' 
#' @param sigma_sv numeric, the rate "diffusion" parameter for the relaxed 
#' clocks sampled from gamma(2/L,2/L)
#' 
#' @param sigma_cl numeric, the rate "diffusion" parameter for the relaxed 
#' clocks sampled from gamma(2/L,20/L)
#' 

#'Install 'simclock' package from github: https://github.com/dosreislab/simclock
devtools::install_github ("dosreislab/simclock")

#' Load 'ape' and 'simclock' libraries #
library('simclock')
library('ape')

#' Read main tree topology. Tree source obtained from Figure 3 of Rannala & Yang 2007. ##
tree_fig1 <- ape::read.tree("tree_fig1_sim.tree")

#' Define the number of nucleotide sequence length and the number of loci to be simulated.
N <- 1000 #Number of nucleotides
L <- 5 #Number of Loci 

#' Define data.frame for rate, sigma_sv and sigma_cl
rate <- data.frame(matrix(ncol = L, nrow = N))
colnames(rate) <- paste("L",1:L,sep = "")
sigma_sv <- data.frame(matrix(ncol = L, nrow = N))
colnames(sigma_sv) <- paste("L",1:L,sep = "")
sigma_cl <- data.frame(matrix(ncol = L, nrow = N)) 
colnames(sigma_cl) <- paste("L",1:L,sep = "")

#' Sample rate and sigma^2 from gamma distribution five times for each Loci and
#' simulate trees under different rate model setting -
#' STR: Strict clock model
#' ILN-SV: Independent log-normal distribution with serious violation
#' ILN-CL: Independent log-normal distribution with clock-like
#' GBM-SV: Geometric Brownian motion with serious violation
#' GBM-CL: Geometric Brownian motion with clock-like
for (i in 1:N) {
  rate[i] <- rgamma(N, shape = 2/L, rate = 2/L)
  sigma_sv[i] <- rgamma(N, shape = 2/L, rate = 2/L) # Seriously-violated
  sigma_cl[i] <- rgamma(N, shape = 2/L, rate = 20/L) # Clock-like
  
  for (j in 1:L) {
    str_clock <- simclock::relaxed.tree(tree_fig1, model = "clk", r = rate[[j]][i])
    iln_clock_sv <- simclock::relaxed.tree(tree_fig1, model = "iln", r = rate[[j]][i], s2 = sigma_sv[[j]][i])
    iln_clock_cl <- simclock::relaxed.tree(tree_fig1, model = "iln", r = rate[[j]][i], s2 = sigma_cl[[j]][i])
    gbm_clock_sv <- simclock::relaxed.tree(tree_fig1, model = "gbm", r = rate[[j]][i], s2 = sigma_sv[[j]][i])
    gbm_clock_cl <- simclock::relaxed.tree(tree_fig1, model = "gbm", r = rate[[j]][i], s2 = sigma_cl[[j]][i])
    write.tree(str_clock, file = paste0("str_clock_nn_",i,"_",j,".tree"))
    write.tree(iln_clock_sv, file = paste0("iln_clock_sv_",i,"_",j,".tree"))
    write.tree(iln_clock_cl, file = paste0("iln_clock_cl_",i,"_",j,".tree"))
    write.tree(gbm_clock_sv, file = paste0("gbm_clock_sv_",i,"_",j,".tree"))
    write.tree(gbm_clock_cl, file = paste0("gbm_clock_cl_",i,"_",j,".tree"))
  }
}



