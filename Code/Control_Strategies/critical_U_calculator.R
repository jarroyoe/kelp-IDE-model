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
# delta_A <- unlist(parameters[parameters$Parameter == 'delta_A', 2])
# gamma_A <- unlist(parameters[parameters$Parameter == 'gamma_A', 2])
# sigma_A <- unlist(parameters[parameters$Parameter == 'sigma_A', 2])
# gamma_S <- unlist(parameters[parameters$Parameter == 'gamma_S', 2])
# beta <- unlist(parameters[parameters$Parameter == 'beta', 2])
# a_A <- unlist(parameters[parameters$Parameter == 'a_A', 2])
# R <- unlist(parameters[parameters$Parameter == 'R', 2])
# delta_U <- unlist(parameters[parameters$Parameter == 'delta_U', 2])
# gamma_G <- unlist(parameters[parameters$Parameter == 'gamma_G', 2])
# a_U <- unlist(parameters[parameters$Parameter == 'a_U', 2])
delta_A <- best_RMSE[1]
gamma_A <- best_RMSE[2]
sigma_A <- best_RMSE[3]
gamma_S <- best_RMSE[4]
beta <- unlist(parameters[parameters$Parameter == 'beta', 2])
a_A <- best_RMSE[5]
R <- best_RMSE[6]
delta_U <- best_RMSE[7]
gamma_G <- best_RMSE[8]
a_U <- best_RMSE[9]
epsilon <- best_RMSE[10]

#define control parameters
eta <- 5
Ac <- 0.5*beta
strong <- 0.7
mid <- 0.4
weak <- 0.1

#define the cluster
maxU <- 1000
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
criticalUs <- parSapply(cl,X=controlStrategies,FUN=function(x){
  critical_U(a_0=1/beta,maxU=maxU,Tmax=12,delta_A=delta_A,
                        gamma_A=gamma_A,
                        sigma_A=sigma_A,
                        gamma_S=gamma_S,
                        beta=beta,
                        a_A=a_A,
                        R=R,
                        delta_U=delta_U,
                        gamma_G=gamma_G,
                        a_U=a_U,epsilon=epsilon,control=x,eta=eta,Ac=Ac,L=1)})
save.image(file = 'critical_U_control_strategies.Rdata')

criticalUs <- parSapply(cl,X=controlStrategies,FUN=function(x){
  critical_U(a_0=1/beta,maxU=maxU,Tmax=12,delta_A=delta_A,
             gamma_A=gamma_A,
             sigma_A=0.9*sigma_A,
             gamma_S=gamma_S,
             beta=beta,
             a_A=a_A,
             R=R,
             delta_U=delta_U,
             gamma_G=gamma_G,
             a_U=a_U,epsilon=epsilon,control=x,eta=eta,Ac=Ac,L=1)})
save.image(file = 'critical_U_control_strategies_0_9.Rdata')

criticalUs <- parSapply(cl,X=controlStrategies,FUN=function(x){
  critical_U(a_0=1/beta,maxU=maxU,Tmax=12,delta_A=delta_A,
             gamma_A=gamma_A,
             sigma_A=1.1*sigma_A,
             gamma_S=gamma_S,
             beta=beta,
             a_A=a_A,
             R=R,
             delta_U=delta_U,
             gamma_G=gamma_G,
             a_U=a_U,epsilon=epsilon,control=x,eta=eta,Ac=Ac,L=1)})
save.image(file = 'critical_U_control_strategies_1_1.Rdata')

stopCluster(cl)
