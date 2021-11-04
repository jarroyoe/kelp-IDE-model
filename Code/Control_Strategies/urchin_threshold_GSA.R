library(parallel)
library(here)
library(dplyr)
library(randomForest)
source(paste(sep='',here::here(),'/Code/Control_Strategies/model_with_spread_rate.R'))
load(paste(sep='',here::here(),'/Results/ABC/final_distributions.RData'))

regionToExplore <- 2
########
#Generate empirical densities of parameters
pdf_delta_A <- density((allEstimates %>% filter(Region==regionToExplore,Parameter=='delta_A') %>% select(Value))$Value,from=0,to=1)
pdf_beta <- density((allEstimates %>% filter(Region==regionToExplore,Parameter=='beta') %>% select(Value))$Value,from=0,to=10)
pdf_sigma_A <- density((allEstimates %>% filter(Region==regionToExplore,Parameter=='sigma_A') %>% select(Value))$Value,from=0,to=100)
pdf_R <- density((allEstimates %>% filter(Region==regionToExplore,Parameter=='R') %>% select(Value))$Value,from=0,to=10)
pdf_gamma_S <- density((allEstimates %>% filter(Region==regionToExplore,Parameter=='gamma_S') %>% select(Value))$Value,from=0,to=100)
pdf_gamma_G <- density((allEstimates %>% filter(Region==regionToExplore,Parameter=='gamma_G') %>% select(Value))$Value,from=0,to=10)
pdf_gamma_A <- density((allEstimates %>% filter(Region==regionToExplore,Parameter=='gamma_A') %>% select(Value))$Value,from=0,to=1)
pdf_delta_U <- density((allEstimates %>% filter(Region==regionToExplore,Parameter=='delta_U') %>% select(Value))$Value,from=0,to=1)
pdf_a_A <- density((allEstimates %>% filter(Region==regionToExplore,Parameter=='a_A') %>% select(Value))$Value,from=0,to=100)
pdf_a_U <- density((allEstimates %>% filter(Region==regionToExplore,Parameter=='a_U') %>% select(Value))$Value,from=0,to=100)
pdf_epsilon <- density((allEstimates %>% filter(Region==regionToExplore,Parameter=='epsilon') %>% select(Value))$Value,from=0,to=10)


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
sigma_A <- sampleGenerator(pdf_sigma_A,m)
R <- sampleGenerator(pdf_R,m)
gamma_S <- sampleGenerator(pdf_gamma_S,m)
gamma_G <- sampleGenerator(pdf_gamma_G,m)
gamma_A <- sampleGenerator(pdf_gamma_A,m)
delta_U <- sampleGenerator(pdf_delta_U,m)
a_A <- sampleGenerator(pdf_a_A,m)
a_U <- sampleGenerator(pdf_a_U,m)
epsilon <- sampleGenerator(pdf_epsilon,m)
L <- runif(m,0,10)
a0 <- runif(m,0,5)
# mu_U <- runif(m)
# mu_S <- runif(m)
# mu_A <- runif(m)
# eta <- runif(m,max=100)
# Ac <- runif(m,max=1)
########
#Run Monte Carlo simulations
maxU=1000
Tmax=12

cores <- detectCores()
cl <- makeCluster(cores)
clusterExport(cl, ls())
clusterEvalQ(cl, as.vector(lsf.str(.GlobalEnv)))
sampledRates <- parSapply(cl=cl,1:m,function(i){
  critical_U(a0[i],maxU,Tmax,delta_A[i],
           gamma_A[i],
           sigma_A[i],
           gamma_S[i],
           beta[i],
           a_A[i],
           R[i],
           delta_U[i],
           gamma_G[i],
           a_U[i],
           epsilon[i],
           c(0,0,0),
           0,
           0,
           L[i])
})
#turn any NA's (vertical lines) into 0s
sampledRates[which(is.na(sampledRates))] <- rep(0,length(which(is.na(sampledRates))))
save.image(file='monte_carlo_simulations.RData')
########
#Run the simulations over a Random Forest
RF_data <- data.frame(sampledRates=sampledRates,
                      a0=a0,
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
                      L=L)
RF <- randomForest(sampledRates~.,data=RF_data,importance=TRUE,proximity=TRUE)
save.image(file='urchin_threshold_GSA.RData')

stopCluster(cl)