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

HER <- function(Y, X, ID, Z=NULL, Phi_est=F, draw=5000, burn=2000){
  # settings 
  p <- dim(X[[1]])[2]
  n <- unlist(lapply(Y, length))/2
  K <- length(Y)
  N <- unlist(lapply(ID, max))
  MC <- draw - burn 
  
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
  XZ <- X
  
  # initial values 
  gam_p <- c(-2, rep(0, p-1))    # coefficients for positive error prob
  gam_n <- c(-2, rep(0, p-1))    # coefficients for negative error prob
  if(qq>0){
    gam_p <- c(gam_p, rep(0, qq))
    gam_n <- c(gam_n, rep(0, qq))
  }

  Theta <- list()    # symmetric network 
  for(k in 1:K){
    Theta[[k]] <- Y[[k]][1:n[k]]
  }
  Phi <- unlist(lapply(Y, mean))
  
  # prior
  A0 <- rep(0, p+qq)
  invB0 <- 0.01*diag(p+qq)
  al <- be <- 1
  
  # objects for posterior samples 
  gam_p_pos <- gam_n_pos <- matrix(NA, MC, p+qq)
  Theta_pos <- list()
  for(k in 1:K){
    Theta_pos[[k]] <- matrix(NA, MC, n[k])
  }
  Phi_pos <- matrix(NA, MC, K)
  
  # iteration 
  pb <- progress_bar$new(total=draw)   
  for(it in 1:draw){
    # Theta (its dimension is the same as that of Y)
    Theta_ex <- list()
    for(k in 1:K){
      Theta_ex[[k]] <- c(Theta[[k]], Theta[[k]])     
    }
    
    # omega_n
    Om <- list()
    for(k in 1:K){
      Om[[k]] <- rep(NA, 2*n[k])
      psi_n <- c(XZ[[k]]%*%gam_n)
      Om[[k]] <- Theta_ex[[k]] * pgdraw(b=rep(1,2*n[k]), c=psi_n)
    }
    
    # gam_n
    if(qq>0){ diag(invB0)[(p+1):(p+qq)] <- lam_n }
    mat <- invB0
    vec <- invB0%*%A0
    for(k in 1:K){
      mat <- mat + t(XZ[[k]])%*%(XZ[[k]]*Om[[k]])
      vec <- vec + t(XZ[[k]])%*%((0.5-Y[[k]])*Theta_ex[[k]]) 
    }
    B_n <- solve(mat)
    A_n <- c(B_n%*%vec)
    gam_n <- mvrnorm(1, A_n, B_n)
    
    # lambda_n
    lam_n <- rgamma(1, 1+0.5*qq, 1+0.5*sum(gam_n[-(1:p)]^2))
    
    # omega_p
    Om <- list()
    for(k in 1:K){
      Om[[k]] <- rep(NA, 2*n[k])
      psi_p <- c(XZ[[k]]%*%gam_p)
      Om[[k]] <- (1-Theta_ex[[k]]) * pgdraw(b=rep(1,2*n[k]), c=psi_p)
    }
    
    # gam_p
    if(qq>0){ diag(invB0)[(p+1):(p+qq)] <- lam_p }
    mat <- invB0
    vec <- invB0%*%A0
    for(k in 1:K){
      mat <- mat + t(XZ[[k]])%*%(XZ[[k]]*Om[[k]])
      vec <- vec + t(XZ[[k]])%*%((Y[[k]]-0.5)*(1-Theta_ex[[k]])) 
    }
    B_p <- solve(mat)
    A_p <- c(B_p%*%vec)
    gam_p <- mvrnorm(1, A_p, B_p)
    
    # lambda_p
    if(qq>0){ 
      lam_p <- rgamma(1, 1+0.5*qq, 1+0.5*sum(gam_p[-(1:p)]^2))
    }
    
    # Theta
    for(k in 1:K){
      n_k <- n[k]
      psi_n <- c(XZ[[k]]%*%gam_n)
      psi_p <- c(XZ[[k]]%*%gam_p)
      p1_num <- Phi[k] * exp( psi_n[1:n_k]*(1-Y[[k]][1:n_k]) + psi_n[-(1:n_k)]*(1-Y[[k]][-(1:n_k)]) ) 
      p1_den <- (1+exp(psi_n[1:n_k])) * (1+exp(psi_n[-(1:n_k)]))
      p0_num <- (1-Phi[k]) * exp( psi_p[1:n_k]*Y[[k]][1:n_k] + psi_p[-(1:n_k)]*Y[[k]][-(1:n_k)] ) 
      p0_den <- (1+exp(psi_p[1:n_k])) * (1+exp(psi_p[-(1:n_k)]))
      prob <- 1 / (1 + p0_num*p1_den/(p0_den*p1_num))
      Theta[[k]] <- rbinom(n_k, 1, prob)
    }
    
    # Phi
    if(Phi_est){
      ss <- unlist( lapply(Theta, sum) )
      Phi <- rbeta(K, al+ss, be+(n-ss))
    }
    
    # save 
    if(it>burn){
      cc <- it - burn
      gam_p_pos[cc,] <- gam_p
      gam_n_pos[cc,] <- gam_n
      for(k in 1:K){
        Theta_pos[[k]][cc,] <- Theta[[k]]
      }
      Phi_pos[cc,] <- Phi
    }
    
    pb$tick()
  }
  
  # results 
  Result <- list(gam_p=gam_p_pos, gam_n=gam_n_pos, Theta=Theta_pos, Phi=Phi_pos)
  return(Result)
}



