
#' Simulating evolutionary rates and branch lengths along a phylogeny 
#' under strict clock and relaxed clock models - 
#' 
#' The \code{model} have three model settings such as "clk", "iln", "gbm_RY07",
#' for the strict clock (constant rates along the phylogeny), 
#' the independent log-normal rates (independent rates along the phylogeny), and
#' geometric-Brownian motion rates (autocorrelated rates along the phylogeny)
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
#' @details
#' To simulate mean rate, one can use 'rgamma()' function to generate mean
#' The 'simclock' R package can be used to generate phylogenetic trees 
#' under the respective model. The \code{model == clk} does not require 
#' sigma^2 or log-variance parameter, while both log-variance and mean rates
#' are simulated, when \code{model == "iln"} or \code{model == "gbm"} is used.
#' See the implementation of 'simclock' R package published in 
#' https://github.com/dosreislab/simclock. 


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

#' Sample rate and sigma^2 from gamma distribution for each Loci
#' and generate trees under different five clock model settings -
#' STR - Strict clock model
#' ILN-SV - Independent log-normal distribution with seriously violated rates
#' ILN-CL - Independent log-normal distribution with clock-like rates
#' GBM-SV - Geometric Brownian motion with serious violation
#' GBM-CL - Geometric Brownian motion with clock-like
for (i in 1:N) {
  rate[i] <- rgamma(N, shape = 2/L, rate = 2/L)
  sigma_sv[i] <- rgamma(N, shape = 2/L, rate = 2/L) # Seriously-violated
  sigma_cl[i] <- rgamma(N, shape = 2/L, rate = 20/L) # Clock-like
  
  for (j in 1:L) {
    str_clock <- simclock::relaxed.tree(tree_fig1, model = "clk", r = rate[[j]][i])
    iln_clock_sv <- simclock::relaxed.tree(tree_fig1, model = "iln", r = rate[[j]][i], s2 = sigma_sv[[j]][i])
    iln_clock_cl <- simclock::relaxed.tree(tree_fig1, model = "iln", r = rate[[j]][i], s2 = sigma_cl[[j]][i])
    gbm_clock_sv <- simclock::relaxed.tree(tree_fig1, model = "gbm_RY07", r = rate[[j]][i], s2 = sigma_sv[[j]][i])
    gbm_clock_cl <- simclock::relaxed.tree(tree_fig1, model = "gbm_RY07", r = rate[[j]][i], s2 = sigma_cl[[j]][i])
    write.tree(str_clock, file = paste0("str_clock_nn_",i,"_",j,".tree"))
    write.tree(iln_clock_sv, file = paste0("iln_clock_sv_",i,"_",j,".tree"))
    write.tree(iln_clock_cl, file = paste0("iln_clock_cl_",i,"_",j,".tree"))
    write.tree(gbm_clock_sv, file = paste0("gbm_clock_sv_",i,"_",j,".tree"))
    write.tree(gbm_clock_cl, file = paste0("gbm_clock_cl_",i,"_",j,".tree"))
  }
}

