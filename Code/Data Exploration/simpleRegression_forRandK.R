library(here)
library(dplyr)
library(ggplot2)

#Load the data
load(paste(here::here(),'/Clean Data/Atx.rdata',sep=''))
regions <- matrix(c(c(38.5,38.6),c(39.27,39.31)),ncol=2)
currentRegion=1;
region <- which(between(x,regions[1,currentRegion],regions[2,currentRegion]))

Avalues <- A[region,20:25]
Ameans <- colMeans(Avalues)
R <- NULL
K <- NULL

for(i in 1:length(region)){
  X <- 1/Avalues[i,1:(length(20:25)-1)]
  Y <- NULL
  for(j in 1:(length(20:25)-1)){
    Y <- rbind(Y,1/Avalues[i,j+1])
  }
  if(length(which(is.infinite(X)))!=0||length(which(is.infinite(Y)))!=0){next}
  toRegress <- data.frame(x=X,y=Y)
  linearModel <- lm(y~x,toRegress)
  
  b <- linearModel$coefficients[1]
  m <- linearModel$coefficients[2]
  
  thisR <- 1/m
  thisK <- (thisR-1)/(b*thisR)
  R <- c(R,thisR)
  K <- c(K,thisK)
}
R <- R[R>-2]
K <- K[K<1]
hist(R)
hist(K)

save(R,K,currentRegion,file=paste('estimates_R_K_',currentRegion,'.RData',sep=''))