load(paste(here(),'Clean Data/Atx.rdata',sep='/'))
load(paste(here(),'Clean Data/Utx.rdata',sep='/'))
regions <- matrix(c(c(38.5,38.6),c(39.27,39.31)),ncol=2)

currentRegion=1
region <- which(between(x,regions[1,currentRegion],regions[2,currentRegion]))
xmin <- regions[1,currentRegion]
xmax <- regions[2,currentRegion]
prior_A <- A[region,23]
posterior_A <- A[region,24]
prior_U <- U[region,2]
posterior_U <- U[region,3]

#estimatedR <- c(4.841198,1.256952)
estimatedbeta <- c(0.4771063,0.4191067)
#estimatedDelta_A <- c(0.947,0.954)
