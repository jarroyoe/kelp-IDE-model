#load the required packages
closeAllConnections()
rm(list = ls())
library(parallel)
library(EasyABC)
library(dplyr)
library(here)

#load the information
source(paste(here(),'/Code/ABC/globalVariables.R',sep=''))
source(paste(here(),'/Code/ABC/model.R',sep=''))

#train the model using the data from 2007 and 2008 (columns 2 and 3 of urchin, columns 23 and 24 of algae)
n <- 300
tol <- 0.1
stats_target <- summaryStatistics(posterior_A,posterior_A)

results <- ABC_mcmc(
  method = 'Marjoram',
  model = model,
  prior = priors,
  summary_stat_target = stats_target,
  n_rec=n,
  #tol=0.1,
  tolerance_quantile=tol,
  #tolerance_tab = tol,
  #nb_simul = n,
  verbose=TRUE,
  use_seed = TRUE,
  n_cluster = detectCores()
)
save.image(file=paste('Results',currentRegion,'.RData',sep=''))

