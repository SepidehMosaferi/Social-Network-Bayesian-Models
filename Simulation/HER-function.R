library(pgdraw)
library(MCMCpack)
library(progress)

###----------------------------------------------------###
###          Heterogeneous Error Model                 ###  
###----------------------------------------------------###
## Y: binary response for K-networks
## X: covariate matrix for each edge

## ID for edge should be implemented as follows: 
##  > ID_theta <- t(combn(1:N, 2))
##  > ID <- rbind(ID_theta, ID_theta[, c(2,1)])

HER <- function(Y, X, ID, Z=NULL, draw=5000, burn=2000){
  # settings 
  p <- dim(X)[2]
  n <- dim(X)[1]/2
  K <- dim(Y)[2]
  N <- max(ID)
  MC <- draw - burn 
  
  # check of ID ordering
  ID_theta <- t(combn(1:N, 2))
  ID_reference <- rbind(ID_theta, ID_theta[, c(2,1)])
  if(sum(abs(ID-ID_reference))>0){
    print("ID should be re-ordered")
    break()
  }
  
  # nonparametric term (Z) 
  qq <- 0
  if(is.null(Z)==F){
    qq <- dim(Z)[2]
    if(dim(Z)[1]!=(2*n)){
      print("Row dimensions of X and Z are not equal.")
      break()
    }
    lam_p <- lam_n <- 1
  }
  
  # design matrix
  XZ <- cbind(X, Z)
  
  # initial values 
  gam_p <- c(-2, rep(0, p-1))
  gam_n <- c(-2, rep(0, p-1))
  if(qq>0){
    gam_p <- c(gam_p, rep(0, qq))
    gam_n <- c(gam_n, rep(0, qq))
  }

  Theta <- matrix(NA, n, K)
  for(k in 1:K){
    Theta[,k] <- Y[1:n, k]
  }
  Phi <- apply(Y, 2, mean)
  
  # prior
  A0 <- rep(0, p+qq)
  invB0 <- 0.01*diag(p+qq)
  al <- be <- 1
  
  # objects for posterior samples 
  gam_p_pos <- gam_n_pos <- matrix(NA, MC, p+qq)
  Theta_pos <- array(NA, c(MC, n, K))
  Phi_pos <- matrix(NA, MC, K)
  
  # iteration 
  pb <- progress_bar$new(total=draw)   
  for(it in 1:draw){
    Theta_ex <- rbind(Theta, Theta)    # dimension is the same as that of X
    
    # omega_n
    Om <- matrix(NA, 2*n, K)
    psi_n <- c(XZ%*%gam_n)
    for(k in 1:K){
      Om[,k] <- Theta_ex[,k] * pgdraw(b=rep(1,2*n), c=psi_n)
    }
    
    # gam_n
    if(qq>0){ diag(invB0)[(p+1):(p+qq)] <- lam_n }
    mat <- invB0
    for(k in 1:K){
      mat <- mat + t(XZ)%*%(XZ*Om[,k])
    }
    B_n <- solve(mat)
    A_n <- B_n%*%( invB0%*%A0 + K*rowMeans(t(XZ)%*%((0.5-Y)*Theta_ex)) )
    gam_n <- mvrnorm(1, A_n, B_n)
    
    # lambda_n
    lam_n <- rgamma(1, 1+0.5*qq, 1+0.5*sum(gam_n[-(1:p)]^2))
    
    # omega_p
    Om <- matrix(NA, 2*n, K)
    psi_p <- c(XZ%*%gam_p)
    for(k in 1:K){
      Om[,k] <- (1-Theta_ex[,k])*pgdraw(b=rep(1,2*n), c=psi_p)
    }
    
    # gam_p
    if(qq>0){ diag(invB0)[(p+1):(p+qq)] <- lam_p }
    mat <- invB0
    for(k in 1:K){
      mat <- mat + t(XZ)%*%(XZ*Om[,k])
    }
    B_p <- solve(mat)
    A_p <- B_p%*%( invB0%*%A0 + K*rowMeans(t(XZ)%*%((Y-0.5)*(1-Theta_ex))) )
    gam_p <- mvrnorm(1, A_p, B_p)
    
    # lambda_p
    if(qq>0){ 
      lam_p <- rgamma(1, 1+0.5*qq, 1+0.5*sum(gam_p[-(1:p)]^2))
    }
    
    # Theta
    psi_n <- c(XZ%*%gam_n)
    psi_p <- c(XZ%*%gam_p)
    for(k in 1:K){
      p1_num <- Phi[k] * exp( psi_n[1:n]*(1-Y[1:n,k]) + psi_n[-(1:n)]*(1-Y[-(1:n),k]) ) 
      p1_den <- (1+exp(psi_n[1:n])) * (1+exp(psi_n[-(1:n)]))
      p0_num <- (1-Phi[k]) * exp( psi_p[1:n]*Y[1:n,k] + psi_p[-(1:n)]*Y[-(1:n),k] ) 
      p0_den <- (1+exp(psi_p[1:n])) * (1+exp(psi_p[-(1:n)]))
      prob <- 1 / (1 + p0_num*p1_den/(p0_den*p1_num))
      Theta[,k] <- rbinom(n, 1, prob)
    }
    
    # Phi
    ss <- apply(Theta, 2, sum)
    Phi <- rbeta(K, al+ss, be+(n-ss))
    
    # save 
    if(it>burn){
      cc <- it - burn
      gam_p_pos[cc,] <- gam_p
      gam_n_pos[cc,] <- gam_n
      Theta_pos[cc,,] <- Theta
      Phi_pos[cc,] <- Phi
    }
    
    pb$tick()
  }
  
  # results 
  Result <- list(gam_p=gam_p_pos, gam_n=gam_n_pos, Theta=Theta_pos, Phi=Phi_pos)
  return(Result)
}



