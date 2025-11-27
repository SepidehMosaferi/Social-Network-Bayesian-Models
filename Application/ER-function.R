###----------------------------------------------------###
###          Homogeneous Error Model                   ###  
###----------------------------------------------------###

ER <- function(Y, ID, equal_prob=F, draw=5000, burn=2000){
  # settings 
  n <- unlist(lapply(Y, length))/2
  K <- length(Y)
  N <- unlist(lapply(ID, max))
  MC <- draw - burn 
  
  # initial values 
  e_p <- e_n <- 0.1
  Theta <- list()    # symmetric network 
  for(k in 1:K){
    Theta[[k]] <- Y[[k]][1:n[k]]
  }
  Phi <- unlist(lapply(Y, mean))
  alpha <- beta <- 1
  
  # prior
  A0 <- 0
  invB0 <- 0.01
  
  # objects for posterior samples 
  e_p_pos <- e_n_pos <- rep(NA, MC)
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
    
    # e_n
    a_n <- b_n <- 1
    for(k in 1:K){
      a_n <- a_n + sum((1-Y[[k]])*Theta_ex[[k]])
      b_n <- b_n + sum(Y[[k]]*Theta_ex[[k]])
    }
    e_n <- rbeta(1, a_n, b_n)
    
    # e_p
    a_p <- b_p <- 1
    for(k in 1:K){
      a_p <- a_p + sum(Y[[k]]*(1-Theta_ex[[k]]))
      b_p <- b_p + sum((1-Y[[k]])*(1-Theta_ex[[k]]))
    }
    e_p <- rbeta(1, a_p, b_p)
    
    # Theta
    psi_n <- log(e_n/(1-e_n))
    psi_p <- log(e_p/(1-e_p))
    for(k in 1:K){
      n_k <- n[k]
      p1 <- Phi[k]*exp( psi_n*(1-Y[[k]][1:n_k])+psi_n*(1-Y[[k]][-(1:n_k)]) ) / (1+exp(psi_n))^2
      p0 <- (1-Phi[k])*exp( psi_p*Y[[k]][1:n_k] + psi_p*Y[[k]][-(1:n_k)] ) / (1+exp(psi_p))^2
      prob <- p1 / (p1 + p0)
      Theta[[k]] <- rbinom(n_k, 1, prob)
    }
    
    # Phi
    ss <- unlist( lapply(Theta, sum) )
    Phi <- rbeta(K, alpha+ss, beta+(n-ss))
    if(equal_prob){
      Phi <- rep(rbeta(1, alpha+sum(ss), beta+sum(n)-sum(ss)), K)
    }
    
    # save 
    if(it>burn){
      cc <- it - burn
      e_p_pos[cc] <- e_p
      e_n_pos[cc] <- e_n
      for(k in 1:K){
        Theta_pos[[k]][cc,] <- Theta[[k]]
      }
      Phi_pos[cc,] <- Phi
    }
    
    pb$tick()
  }
  
  # results 
  Result <- list(e_p=e_p_pos, e_n=e_n_pos, Theta=Theta_pos, Phi=Phi_pos)
  return(Result)
}

