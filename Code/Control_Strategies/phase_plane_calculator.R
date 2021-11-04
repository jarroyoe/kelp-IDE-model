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

#Define grid
uthresholds<- 10^seq(-1,4,length = 20)
kelpDensities <- 10^seq(0,4,length = 20)/beta

#Define cluster
cores <- detectCores()
cl <- makeCluster(cores)
clusterExport(cl, ls())
clusterEvalQ(cl, as.vector(lsf.str(.GlobalEnv)))

#Test all combinations
nonSpatialVectors <- parApply(cl,X=expand.grid(uthresholds,kelpDensities),MARGIN=1,FUN=function(x){
  n <- 1000
  U0 <- rep(x[1], n)
  A0 <- rep(x[2], n)
  
  IDE_procedure_non_spatial(A_0=A0,
                            U_0=U0,
                            Tmax=50,
                            delta_A=delta_A,
                            gamma_A=gamma_A,
                            sigma_A=sigma_A,
                            gamma_S=gamma_S,
                            beta=beta,
                            a_A=a_A,
                            R=R,
                            delta_U=delta_U,
                            gamma_G=gamma_G,
                            a_U=a_U,
                            epsilon=epsilon,
                            control=c(0,0,0),
                            eta=0,
                            Ac=0)
})

kelp <- nonSpatialVectors[1:50,]
urchins <- nonSpatialVectors[51:100,]
kelp_urchin <- data.frame(Kelp=c(kelp),Urchin=c(urchins))
# kelp_urchin$firstKelp <- rep(0,10000)
# kelp_urchin$firstUrchin <- rep(0,10000)
# kelp_urchin$lastKelp <- rep(0,10000)
# kelp_urchin$lastUrchin <- rep(0,10000)
for(i in 1:20000){
  kelp_urchin$firstUrchin[i] <- kelp_urchin$Urchin[25*floor((i-1)/25)+1]
  kelp_urchin$firstKelp[i] <- kelp_urchin$Kelp[25*floor((i-1)/25)+1]
  kelp_urchin$lastKelp[i] <- kelp_urchin$Kelp[25*floor((i-1)/25)+25]
  kelp_urchin$lastUrchin[i] <- kelp_urchin$Urchin[25*floor((i-1)/25)+25]
  kelp_urchin$ID[i] <- floor((i-1)/25)+1
}
library(ggplot2)
library(scales)
kelp_urchin$time <- rep(1:50,400)
kelp_urchin %>% filter(time<=24) %>% group_by(ID) %>% mutate(Forest = lastKelp > 0.1) %>% ggplot(aes(
  colour = as.factor(ID),
  x = time,
  y = (Kelp + 1)
)) + geom_line(aes(linetype = Forest), alpha = 0.8, size = 0.5) + theme_bw() + theme(text = element_text(size = 14), legend.position = 'none') +
  scale_y_log10(
    breaks = trans_breaks("log10", function(x)
      10 ^ x),
    labels = trans_format("log10", math_format(10 ^ .x))
  ) + xlab("Time") + ylab("Total Kelp Density")+ ggtitle("a)")
kelp_urchin%>% filter(time<=24) %>% group_by(ID) %>% mutate(Forest = lastKelp > 0.1) %>% ggplot(aes(
  colour = as.factor(ID),
  x = time,
  y = (Urchin + 1)
)) + geom_line(aes(linetype = Forest), alpha = 0.8, size = 0.5) + theme_bw() + theme(text = element_text(size = 14), legend.position = 'none') +
  scale_y_log10(
    breaks = trans_breaks("log10", function(x)
      10 ^ x),
    labels = trans_format("log10", math_format(10 ^ .x))
  ) + xlab("Time") + ylab("Total Urchin Density") + ggtitle("b)")
save.image(file='nonSpatialVectors.RData')