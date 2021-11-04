library(parallel)
library(here)
library(dplyr)
library(randomForest)
load(paste(sep='',here::here(),'/Results/ABC/final_distributions.RData'))
load(paste(sep = '', here::here(), '/Results/sigma0/Results_MCMC_1.RData'))
source(paste(sep='',here::here(),  '/Code/sigma0/Spread_Rate/model_with_spread_rate.R'))

regionToExplore <- 1
########
#Generate empirical densities of parameters
pdf_delta_A <- density(results$param[,1],from=0,to=1)
pdf_beta <- density(unlist(allEstimates %>% filter(Region==regionToExplore,Parameter=='beta') %>% select(Value)),from=0,to=1)
#pdf_sigma_A <- density(unlist(allEstimates %>% filter(Region==regionToExplore,Parameter=='sigma_A') %>% select(Value)),from=0,to=100)
pdf_R <- density(results$param[,5],from=0,to=10)
pdf_gamma_S <- density(results$param[,3],from=0,to=100)
pdf_gamma_G <- density(results$param[,7],from=0,to=10)
pdf_gamma_A <- density(results$param[,2],from=0,to=1)
pdf_delta_U <- density(results$param[,6],from=0,to=1)
pdf_a_A <- density(results$param[,4],from=0,to=100)
pdf_a_U <- density(results$param[,8],from=0,to=100)

########
#Generate samples from empirical densities
sampleGenerator <- function(pdf, n) {
  samples <- approx(cumsum(pdf$y) / sum(pdf$y),
                    pdf$x,
                    runif(n))$y
  while (length(which(is.na(samples) > 0))) {
    m <- length(which(is.na(samples)))
    samples[which(is.na(samples))] <-
      approx(cumsum(pdf$y) / sum(pdf$y),
             pdf$x,
             runif(m))$y
  }
  samples
}
m <- 2000
delta_A <- sampleGenerator(pdf_delta_A,m)
beta <- sampleGenerator(pdf_beta,m)
#sigma_A <- sampleGenerator(pdf_sigma_A,m)
R <- sampleGenerator(pdf_R,m)
gamma_S <- sampleGenerator(pdf_gamma_S,m)
gamma_G <- sampleGenerator(pdf_gamma_G,m)
gamma_A <- sampleGenerator(pdf_gamma_A,m)
delta_U <- sampleGenerator(pdf_delta_U,m)
a_A <- sampleGenerator(pdf_a_A,m)
a_U <- sampleGenerator(pdf_a_U,m)
 mu_U <- runif(m)
 mu_S <- runif(m)
 mu_A <- runif(m)
 eta <- runif(m,max=100)
 Ac <- runif(m,max=1)
 L <- runif(m,max=10)
########
#Run Monte Carlo simulations
uThreshold=80
Tmax=12

cores <- detectCores()
cl <- makeCluster(cores)
clusterExport(cl, ls())
clusterEvalQ(cl, as.vector(lsf.str(.GlobalEnv)))
sampledRates <- parSapply(cl=cl,1:m,function(i){
  calculateSpreadRate(uThreshold,Tmax,delta_A[i],
             gamma_A[i],
             0,
             gamma_S[i],
             beta[i],
             a_A[i],
             R[i],
             delta_U[i],
             gamma_G[i],
             a_U[i],
             c(mu_U[i],mu_S[i],mu_A[i]),
             eta[i],
             Ac[i],
             L[i])
})
#turn any NA's (vertical lines) into 0s
sampledRates[which(is.na(sampledRates))] <- rep(0,length(which(is.na(sampledRates))))
save.image(file='monte_carlo_simulations.RData')
########
#Run the simulations over a Random Forest
RF_data <- data.frame(sampledRates=sampledRates,
                      delta_A=delta_A,
                      gamma_A=gamma_A,
                      #sigma_A=sigma_A,
                      gamma_S=gamma_S,
                      beta=beta,
                      a_A=a_A,
                      R=R,
                      delta_U=delta_U,
                      gamma_G=gamma_G,
                      a_U=a_U,
                      mu_U=mu_U,
                      mu_S=mu_S,
                      mu_A=mu_A,
                      eta=eta,
                      Ac=Ac,
                      L=L)
RF <- randomForest(sampledRates~.,data=RF_data,importance=TRUE,proximity=TRUE)
save.image(file='spread_rate_GSA.RData')
