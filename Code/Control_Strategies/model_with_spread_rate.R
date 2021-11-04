####################################################################
#define global variables
xmin <- -5
xmax <- 5
n <- 1000
####################################################################
#define functional forms inside the model
P_A <- function(A, U, delta_A, gamma_A, sigma_A, x,mu,eta,Ac) {
  delta <- sum(A[max(1,round(x-eta*n/(xmax-xmin))):min(n,round(x+eta*n/(xmax-xmin)))])>=Ac
  max(0, delta_A * A[x] * (1 +mu*delta - gamma_A * A[x] * U[x] / (1 + sigma_A * A[x] ^
                                                          2)))
}

S <- function(A, U, gamma_S, beta, x) {
  max(1 - gamma_S * U[x], 0) / (1 + beta * A[x])
}

R_A <- function(A, R, x,mu,eta,Ac) {
  delta <- sum(A[max(1,round(x-eta*n/(xmax-xmin))):min(n,round(x+eta*n/(xmax-xmin)))])>=Ac
  R*(1+mu*delta) * A[x]
}

P_U <- function(A,U, delta_U, x,mu,eta,Ac) {
  delta <- sum(A[max(1,round(x-eta*n/(xmax-xmin))):min(n,round(x+eta*n/(xmax-xmin)))])>=Ac
  delta_U*(1-delta*mu) * U[x]
}

k <- function(x, y, a) {
  a * exp(-a * abs(x - y)) / 2
}

R_U <- function(A, U, x, gamma_G, sigma_A,epsilon) {
  #,a_D,gamma_D){
  gamma_G * U[x] * A[x] / (1 + sigma_A * A[x] ^ 2)+epsilon*A[x]*U[x]
}
####################################################################
#Integration procedure and model procedure
#Weighted Simpson's rule integration procedure
wSimpson <- function(w, f, xmin, xmax, ...) {
  odd <- 0
  even <- 0
  n <- length(w) - 1
  step_size <- (xmax - xmin) / n
  xgrid <- seq(xmin, xmax, by = step_size)
  
  for (i in 1:n - 1) {
    fi <- f(xgrid[i + 1], ...) * w[i + 1]
    if (i %% 2 == 0) {
      even <- even + fi
    } else{
      odd <- odd + fi
    }
  }
  
  (f(xgrid[1], ...) * w[1] + f(xgrid[n + 1], ...) * w[n + 1] + 4 * odd +
      2 * even) * step_size / 3
}

IDE_procedure <-
  function(A_0,
           U_0,
           Tmax,
           delta_A,
           gamma_A,
           sigma_A,
           gamma_S,
           beta,
           a_A,
           R,
           delta_U,
           gamma_G,
           a_U,
           epsilon,
           control,
           eta,
           Ac) {
    A <- rbind(NULL, A_0)
    U <- rbind(NULL, U_0)

    
    for (t in 1:Tmax) {
      spores_A <-
        sapply(
          X = 1:n,
          FUN = R_A,
          A = A[t, ],
          R = R,
          mu = control[2],
          eta = eta,
          Ac = Ac
        )
      A_hat <-
        sapply(
          X = 1:n,
          FUN = P_A,
          A = A[t, ],
          U = U[t, ],
          delta_A = delta_A,
          gamma_A = gamma_A,
          sigma_A = sigma_A,
          mu = control[3],
          eta = eta,
          Ac = Ac
        ) +
        sapply(
          X = 1:n,
          FUN = S,
          A = A[t, ],
          U = U[t, ],
          gamma_S = gamma_A,
          beta = beta
        ) *
        sapply(
          X = (1:n)*(xmax-xmin)/n+xmin,
          FUN = wSimpson,
          w = spores_A,
          f = k,
          xmin = xmin,
          xmax = xmax,
          a = a_A
        )
      A <- rbind(A, A_hat)
      
      spores_U <-
        sapply(
          X = 1:n,
          FUN = R_U,
          A = A[t, ],
          U = U[t, ],
          gamma_G = gamma_G,
          sigma_A = sigma_A,
          epsilon=epsilon
        )#,gamma_D=params[12],a_D=params[13])
      U_hat <-
        sapply(
          X = 1:n,
          FUN = P_U,
          A = A[t, ],
          U = U[t, ],
          delta_U = delta_U,
          mu = control[1],
          eta = eta,
          Ac = Ac
        ) +
        sapply(
          X = (1:n)*(xmax-xmin)/n+xmin,
          FUN = wSimpson,
          w = spores_A,
          f = k,
          xmin = xmin,
          xmax = xmax,
          a = a_U
        )
      U <- rbind(U, U_hat)
    }
    return(A)
  }

IDE_procedure_non_spatial <-
  function(A_0,
           U_0,
           Tmax,
           delta_A,
           gamma_A,
           sigma_A,
           gamma_S,
           beta,
           a_A,
           R,
           delta_U,
           gamma_G,
           a_U,
           epsilon,
           control,
           eta,
           Ac) {
    A <- rbind(NULL, A_0)
    U <- rbind(NULL, U_0)
    
    
    for (t in 1:Tmax) {
      spores_A <-
        sapply(
          X = 1:n,
          FUN = R_A,
          A = A[t, ],
          R = R,
          mu = control[2],
          eta = eta,
          Ac = Ac
        )
      A_hat <-
        sapply(
          X = 1:n,
          FUN = P_A,
          A = A[t, ],
          U = U[t, ],
          delta_A = delta_A,
          gamma_A = gamma_A,
          sigma_A = sigma_A,
          mu = control[3],
          eta = eta,
          Ac = Ac
        ) +
        sapply(
          X = 1:n,
          FUN = S,
          A = A[t, ],
          U = U[t, ],
          gamma_S = gamma_A,
          beta = beta
        ) *
        sapply(
          X = (1:n)*(xmax-xmin)/n+xmin,
          FUN = wSimpson,
          w = spores_A,
          f = k,
          xmin = xmin,
          xmax = xmax,
          a = a_A
        )
      A <- rbind(A, A_hat)
      
      spores_U <-
        sapply(
          X = 1:n,
          FUN = R_U,
          A = A[t, ],
          U = U[t, ],
          gamma_G = gamma_G,
          sigma_A = sigma_A,
          epsilon=epsilon
        )#,gamma_D=params[12],a_D=params[13])
      U_hat <-
        sapply(
          X = 1:n,
          FUN = P_U,
          A = A[t, ],
          U = U[t, ],
          delta_U = delta_U,
          mu = control[1],
          eta = eta,
          Ac = Ac
        ) +
        sapply(
          X = (1:n)*(xmax-xmin)/n+xmin,
          FUN = wSimpson,
          w = spores_A,
          f = k,
          xmin = xmin,
          xmax = xmax,
          a = a_U
        )
      U <- rbind(U, U_hat)
    }
    sumA <- rowSums(A)*(xmax-xmin)/n
    sumU <- rowSums(U)*(xmax-xmin)/n
    return(cbind(sumA,sumU))
  }

####################################################################
#Calculation of spread rate
calculateSpreadRate <-
  function(a_0,u_0,
           Tmax,
           delta_A,
           gamma_A,
           sigma_A,
           gamma_S,
           beta,
           a_A,
           R,
           delta_U,
           gamma_G,
           a_U,
           epsilon,
           control = c(0, 0, 0),
           eta = 0,
           Ac = 0,
           L=0) {
    A0 <- rep(0,n)
    A0[round(n / 2-L*n/(xmax-xmin)):round(n / 2+L*n/(xmax-xmin))] <- a_0
    U0 <- rep(u_0, n)
    
    At <-
      IDE_procedure(
        A0,
        U0,
        Tmax,
        delta_A,
        gamma_A,
        sigma_A,
        gamma_S,
        beta,
        a_A,
        R,
        delta_U,
        gamma_G,
        a_U,
        epsilon,
        control,
        eta,
        Ac
      )
    invasion_width <-
      apply(At, function(x) {
        max(max(which(x > 0.01)),0) * (xmax-xmin)/n+xmin
      }, MARGIN = 1)
    # linearRegression <- lm(invasion_width ~ seq(0, Tmax))
    # return(linearRegression$coefficients[2])
    return(invasion_width[13]/13)
  }

critical_U <- function(a_0,maxU,Tmax,delta_A,
                       gamma_A,
                       sigma_A,
                       gamma_S,
                       beta,
                       a_A,
                       R,
                       delta_U,
                       gamma_G,
                       a_U,epsilon,control,eta,Ac,L){
  imax <- 10
  tol <- 1e-3
  a<-0
  b<-maxU
  i <- 0
  spreadRate <- 1
  while(i<imax && abs(spreadRate)>tol){
    spreadRate <- calculateSpreadRate(a_0,(a+b)/2,Tmax,
                                             delta_A,
                                             gamma_A,
                                             sigma_A,
                                             gamma_S,
                                             beta,
                                             a_A,
                                             R,
                                             delta_U,
                                             gamma_G,
                                             a_U,epsilon,control,eta,Ac,L)
    if(spreadRate<0){
      b <- (a+b)/2
    }else{
      a <- (a+b)/2
    }
    i <- i+1
  }
  return((a+b)/2)
}
