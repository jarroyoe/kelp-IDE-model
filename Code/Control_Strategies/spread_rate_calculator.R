closeAllConnections()
rm(list=ls())
library(parallel)
library(dplyr)
library(here)
source(paste(
  sep = '',
  here::here(),
  '/Code/Control_Strategies/model_with_spread_rate.R'
))
load(paste(sep = '', here::here(), '/Results/ABC/final_distributions.RData'))

#define model parameters
regionToExplore <- 2
parameters <-
  allEstimates %>% filter(Region == regionToExplore) %>% group_by(Parameter) %>% dplyr::summarize(MeanVal =
                                                                                                    mean(Value))
delta_A <- best_RMSE_2[1]
gamma_A <- best_RMSE_2[2]
sigma_A <- best_RMSE_2[3]
gamma_S <- best_RMSE_2[4]
beta <- unlist(parameters[parameters$Parameter == 'beta', 2])
a_A <- best_RMSE_2[5]
R <- best_RMSE_2[6]
delta_U <- best_RMSE_2[7]
gamma_G <- best_RMSE_2[8]
a_U <- best_RMSE_2[9]
epsilon <- best_RMSE_2[10]

#define control parameters
eta <- 1
Ac <- 0.5*beta
strong <- 0.7
mid <- 0.4
weak <- 0.1

#define the cluster
uthreshold <- 82
uthresholds<- seq(0,100,by=10)
kelpDensities <- seq(1,10)/beta
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

#test all combinations of control strategies
spreadRates <- parSapply(cl,X=controlStrategies,FUN=function(x){
  calculateSpreadRate(a_0=1/beta,u_0=uthreshold,Tmax=12,delta_A=delta_A,
             gamma_A=gamma_A,
             sigma_A=sigma_A,
             gamma_S=gamma_S,
             beta=beta,
             a_A=a_A,
             R=R,
             delta_U=delta_U,
             gamma_G=gamma_G,
             a_U=a_U,
             epsilon=epsilon,control=x,eta=eta,Ac=Ac)})
save.image(file = 'spread_rate_calculator_eta_1.Rdata')

spreadRates <- parSapply(cl,X=controlStrategies,FUN=function(x){
  calculateSpreadRate(a_0=1/beta,u_0=uthreshold,Tmax=12,delta_A=delta_A,
                      gamma_A=gamma_A,
                      sigma_A=sigma_A,
                      gamma_S=gamma_S,
                      beta=beta,
                      a_A=a_A,
                      R=R,
                      delta_U=delta_U,
                      gamma_G=gamma_G,
                      a_U=a_U,
                      epsilon=epsilon,control=x,eta=10*eta,Ac=Ac)})
save.image(file = 'spread_rate_calculator_eta_10.Rdata')

#test all initial urchin threshold densities
spreadRates <- parSapply(cl,X=uthresholds,FUN=function(x){
  calculateSpreadRate(a_0=1/beta,u_0=x,Tmax=12,delta_A=delta_A,
                      gamma_A=gamma_A,
                      sigma_A=sigma_A,
                      gamma_S=gamma_S,
                      beta=beta,
                      a_A=a_A,
                      R=R,
                      delta_U=delta_U,
                      gamma_G=gamma_G,
                      a_U=a_U,
                      epsilon=epsilon,control=c(0,strong,0),eta=eta,Ac=Ac)})
save.image(file = 'spread_rate_urchin_densities_strong.Rdata')

#test all initial kelp threshold densities
spreadRates <- parSapply(cl,X=kelpDensities,FUN=function(x){
  calculateSpreadRate(a_0=x,u_0=uthreshold,Tmax=12,delta_A=delta_A,
                      gamma_A=gamma_A,
                      sigma_A=sigma_A,
                      gamma_S=gamma_S,
                      beta=beta,
                      a_A=a_A,
                      R=R,
                      delta_U=delta_U,
                      gamma_G=gamma_G,
                      a_U=a_U,
                      epsilon=epsilon,control=c(0,strong,0),eta=eta,Ac=Ac)})
save.image(file = 'spread_rate_kelp_densities_strong.Rdata')

stopCluster(cl)
