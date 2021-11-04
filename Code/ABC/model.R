#Definition of the model to use EasyABC procedure
#LIST OF PARAMETERS:
## 1: seed to use parallel ABC
## 2: delta_A unif 0,1
## 2: gamma_A unif 0 1
## 3: sigma_A unif 0 100
## 4: gamma_S unif 0 100
## 6: beta unif 0 2 OMMITED DUE TO DIFFERENT ESTIMATION
## 5: a_A unif 0 2
## 6: R unif 0 10
## 7: delta_U unif 0 1
## 8: gamma_G unif 0 10
## 9: a_U unif 0 2

#Drift algae is ommited from this iteration of the model
## 12: gamma_D unif 0 10
## 13: a_D unif 0 2

#load priors
priors <- list(c('unif',0,1),
  c('unif',0,1),c('unif',0,100),c('unif',0,100),#c('unif',0,100),
               c('unif',0,100),
               c('unif',0,10),c('unif',0,1),
               c('unif',0,10),c('unif',0,100),c('unif',0,10))#,
               #c('unif',0,10),c('unif',0,2))

#define functional forms inside the model
P_A <- function(A,U,delta_A,gamma_A,sigma_A,x){
  max(0,delta_A*A[x]*(1-gamma_A*A[x]*U[x]/(1+sigma_A*A[x]^2)))
}

S <- function(A,U,gamma_S,R,K,x){
  #exp(-gamma_S*U[x])
  max(1-gamma_S*U[x],0)/(1+(R-1)*A[x]/K)
}

R_A <- function(A,R,x){
  R*A[x]
}

P_U <- function(U,delta_U,x){
  delta_U*U[x]
}

k <- function(x,y,a){
  a*exp(-a*abs(x-y))/2
}

R_U <- function(A,U,x,gamma_G,sigma_A,epsilon){#,a_D,gamma_D){
  G <- gamma_G*U[x]*A[x]/(1+sigma_A*A[x]^2)+epsilon*A[x]*U[x]
  #D <- gamma_D*U[x]*wSimpson(w=A,f=k,xmin = xmin,xmax=xmax,y=xmin+(xmax-xmin)*(x-1)/(n-1),a=a_D)
  G#+D
}
####################################################################
#Integration procedure and model procedure
#Weighted Simpson's rule integration procedure
wSimpson <- function(w,f,xmin,xmax,...){
  odd <- 0
  even <- 0
  n <- length(w)-1
  step_size <- (xmax-xmin)/n
  xgrid <- seq(xmin,xmax,by=step_size)
  
  for(i in 1:n-1){
    fi <- f(xgrid[i+1],...)*w[i+1]
    if(i%%2==0){
      even <- even+fi
    }else{
      odd <- odd+fi
    }
  }
  
  (f(xgrid[1],...)*w[1]+f(xgrid[n+1],...)*w[n+1]+4*odd+2*even)*step_size/3
}

IDE_procedure <- function(A_0,U_0,Tmax=12,delta_A,gamma_A,sigma_A,gamma_S,K,a_A,R,delta_U,gamma_G,a_U,epsilon){
  A <- rbind(NULL, A_0)
  U <- rbind(NULL, U_0)
  n <- length(region)
  
  for(t in 1:Tmax){
    spores_A <- sapply(X=1:n,FUN=R_A,A=A[t,],R=R)
    A_hat <- sapply(X=1:n,FUN=P_A,A=A[t,],U=U[t,],delta_A=delta_A,gamma_A=gamma_A,sigma_A=sigma_A)+
      sapply(X=1:n,FUN=S,A=A[t,],U=U[t,],gamma_S=gamma_A,R=R,K=1/K)*
      sapply(X=(1:n)*(xmax-xmin)/n+xmin,FUN = wSimpson,w = spores_A,f = k,xmin = xmin,xmax = xmax,a=a_A)
    A <- rbind(A,A_hat)
    
    spores_U <- sapply(X=1:n,FUN=R_U,A=A[t,],U=U[t,],gamma_G=gamma_G,sigma_A=sigma_A,epsilon=epsilon)#,gamma_D=params[12],a_D=params[13])
    U_hat <- sapply(X=1:n,FUN=P_U,U=U[t,],delta_U=delta_U)+
      sapply(X=(1:n)*(xmax-xmin)/n+xmin,FUN = wSimpson,w = spores_A,f = k,xmin = xmin,xmax = xmax,a=a_U)
    U <- rbind(U,U_hat)
  }
  return(A)
}

####################################################################
#Summary statistics and model that enters into the ABC
summaryStatistics <- function(A,A_0){
 c(norm(A-A_0,"2")/norm(A_0,'2'))
}

model <- function(params){
  library(here)
  library(dplyr)
  set.seed(params[1])
  source(paste(here(),'/Code/ABC/globalVariables.R',sep=''))
  source(paste(here(),'/Code/ABC/model.R',sep=''))
  # source(paste(here(),'/Code/ABC/obtained_parameters.R',sep=''))
  
  # if(obtainedParameters==1){
    delta_A <- params[2]
    gamma_A <- params[3]
    sigma_A <- params[4]
    gamma_S <- params[5]
    K <- estimatedK[currentRegion]
    a_A <- params[6]
    R <- params[7]
    delta_U <- params[8]
    gamma_G <- params[9]
    a_U <- params[10]
    epsilon <- params[11]
  #   
  #   
  #   }else if(obtainedParameters==2){
  #   sigma_A <- params[2]
  #   gamma_S <- params[3]
  #   beta <- params[4]
  #   a_A <- params[5]
  #   R <- params[6]
  #   delta_U <- params[7]
  #   gamma_G <- params[8]
  #   a_U <- params[9]
  # 
  #   
  #   }else if(obtainedParameters==3){
  #   gamma_S <- params[2]
  #   beta <- params[3]
  #   a_A <- params[4]
  #   R <- params[5]
  #   delta_U <- params[6]
  #   gamma_G <- params[7]
  #   a_U <- params[8]
  # 
  #   
  #   }else if(obtainedParameters==4){
  #   beta <- params[2]
  #   a_A <- params[3]
  #   R <- params[4]
  #   delta_U <- params[5]
  #   gamma_G <- params[6]
  #   a_U <- params[7]
  # 
  #   
  #   }else if(obtainedParameters==5){
  #   a_A <- params[2]
  #   R <- params[3]
  #   delta_U <- params[4]
  #   gamma_G <- params[5]
  #   a_U <- params[6]
  # 
  #   
  #   }else if(obtainedParameters==6){
  #   R <- params[2]
  #   delta_U <- params[3]
  #   gamma_G <- params[4]
  #   a_U <- params[5]
  # 
  #   
  #   }else if(obtainedParameters==7){
  #   delta_U <- params[2]
  #   gamma_G <- params[3]
  #   a_U <- params[4]
  # 
  #   
  #   }else if(obtainedParameters==8){
  #   gamma_G <- params[2]
  #   a_U <- params[3]
  # 
  #   
  #   }else if(obtainedParameters==9){
  #   a_U <- params[2]
  # }
  # browser()
    At <-
    IDE_procedure(
      prior_A,
      prior_U,
      Tmax=12,
      delta_A,gamma_A,sigma_A,gamma_S,K,a_A,R,delta_U,gamma_G,a_U,epsilon
    )
  summaryStatistics(At[13,],posterior_A)
}

