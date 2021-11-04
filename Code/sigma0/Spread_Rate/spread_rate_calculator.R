library(parallel)
library(dplyr)
library(here)
load(paste(sep='',here::here(),'/Results/ABC/final_distributions.RData'))
load(paste(sep = '', here::here(), '/Results/sigma0/Results_MCMC_1.RData'))
source(paste(
  sep = '',
  here::here(),
  '/Code/sigma0/Spread_Rate/model_with_spread_rate.R'
))

#define model parameters
regionToExplore <- 1
best_RMSE <- results$param[which(results$stats[,1]==min(results$stats[,1])),]

delta_A <- best_RMSE[1]
gamma_A <- best_RMSE[2]
sigma_A <- 0
gamma_S <- best_RMSE[3]
beta <- unlist(allEstimates %>% filter(Region==regionToExplore,Parameter=='beta') %>% select(Value))
a_A <- best_RMSE[4]
R <- best_RMSE[5]
delta_U <- best_RMSE[6]
gamma_G <- best_RMSE[7]
a_U <- best_RMSE[8]

#define control parameters
eta <- 1
Ac <- 0.5*beta
strong <- 0.7
mid <- 0.4
weak <- 0.1

#define the cluster
uthreshold<- 80
# controlstrategies <-
#   list(
#     c(0, 0, 0),
#     c(weak, 0, 0),
#     c(mid, 0, 0),
#     c(strong, 0, 0),
#     c(0, weak, 0),
#     c(0, mid, 0),
#     c(0, strong, 0),
#     c(0, 0, weak),
#     c(0, 0, mid),
#     c(0, 0, strong)
#   )
controlStrategies <- expand.grid(rep(list(c(0,weak,mid,strong)),3))
controlStrategies <-  setNames(split(as.matrix(controlStrategies),
                                     row(controlStrategies)), paste0("Row",1:3))
cores <- detectCores()
cl <- makeCluster(cores)
clusterExport(cl, ls())
clusterEvalQ(cl, as.vector(lsf.str(.GlobalEnv)))

#test all combinations
spreadRates <- parSapply(cl,X=controlStrategies,FUN=function(x){
  calculateSpreadRate(u_0=uthreshold,Tmax=12,delta_A=delta_A,
             gamma_A=gamma_A,
             sigma_A=sigma_A,
             gamma_S=gamma_S,
             beta=beta,
             a_A=a_A,
             R=R,
             delta_U=delta_U,
             gamma_G=gamma_G,
             a_U=a_U,control=x,eta=eta,Ac=Ac)})
save.image(file = 'spread_rate_control_strategies.Rdata')
