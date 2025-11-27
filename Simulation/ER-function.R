###----------------------------------------------------###
###          Homogeneous Error Model                   ###  
###----------------------------------------------------###

ER <- function(Y, ID, equal_prob=F, draw=5000, burn=2000){
  # settings 
  n <- dim(Y)[1]/2
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
  
  # initial values 
  e_p <- e_n <- 0.1
  Theta <- matrix(NA, n, K)
  for(k in 1:K){
    Theta[,k] <- Y[1:n, k]
  }
  Phi <- apply(Y, 2, mean)
  alpha <- beta <- 1
  
  # prior
  A0 <- 0
  invB0 <- 0.01
  
  # objects for posterior samples 
  e_p_pos <- e_n_pos <- rep(NA, MC)
  Theta_pos <- array(NA, c(MC, n, K))
  Phi_pos <- matrix(NA, MC, K)
  
  # iteration 
  pb <- progress_bar$new(total=draw)   
  for(it in 1:draw){
    Theta_ex <- rbind(Theta, Theta)    # same dimension as Y
    # e_n
    a_n <- 1 + sum((1-Y)*Theta_ex)
    b_n <- 1 + sum(Y*Theta_ex)
    e_n <- rbeta(1, a_n, b_n)
    
    # e_p
    a_p <- 1 + sum(Y*(1-Theta_ex))
    b_p <- 1 + sum((1-Y)*(1-Theta_ex))
    e_p <- rbeta(1, a_p, b_p)
    
    # Theta
    psi_n <- log(e_n/(1-e_n))
    psi_p <- log(e_p/(1-e_p))
    for(k in 1:K){
      p1 <- Phi[k]*exp( psi_n*(1-Y[1:n,k])+psi_n*(1-Y[-(1:n),k]) ) / (1+exp(psi_n))^2
      p0 <- (1-Phi[k])*exp( psi_p*Y[1:n,k] + psi_p*Y[-(1:n),k] ) / (1+exp(psi_p))^2
      prob <- p1 / (p1 + p0)
      Theta[,k] <- rbinom(n, 1, prob)
    }
    
    # Phi
    ss <- apply(Theta, 2, sum)
    Phi <- rbeta(K, alpha+ss, beta+(n-ss))
    if(equal_prob){
      Phi <- rep(rbeta(1, alpha+sum(ss), beta+n*K-sum(ss)), K)
    }
    
    # save 
    if(it>burn){
      cc <- it - burn
      e_p_pos[cc] <- e_p
      e_n_pos[cc] <- e_n
      Theta_pos[cc,,] <- Theta
      Phi_pos[cc,] <- Phi
    }
    
    pb$tick()
  }
  
  # results 
  Result <- list(e_p=e_p_pos, e_n=e_n_pos, Theta=Theta_pos, Phi=Phi_pos)
  return(Result)
}








###----------------------------------------------------###
###        Butts Model for symmetric network           ###  
###----------------------------------------------------###

Butts <- function(Y, ID, draw=5000, burn=2000){
  # settings 
  n <- dim(Y)[1]/2
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
  
  # initial values 
  e_p <- e_n <- rep(0.1, K)
  Theta <- matrix(NA, n, K)
  for(k in 1:K){
    Theta[,k] <- Y[1:n, k]
  }
  Phi <- mean(Y)
  alpha <- beta <- 1
  
  # prior
  A0 <- 0
  invB0 <- 0.01
  
  # objects for posterior samples 
  e_p_pos <- e_n_pos <- matrix(NA, MC, K)
  Theta_pos <- array(NA, c(MC, n, K))
  Phi_pos <- rep(NA, MC)
  
  # iteration 
  pb <- progress_bar$new(total=draw)   
  for(it in 1:draw){
    Theta_ex <- rbind(Theta, Theta)    # same dimension as Y
    
    # e_n
    a_n <- 1 + apply((1-Y)*Theta_ex, 2, sum)
    b_n <- 1 + apply(Y*Theta_ex, 2, sum)
    e_n <- rbeta(K, a_n, b_n)
    
    # e_p
    a_p <- 1 + apply(Y*(1-Theta_ex), 2, sum)
    b_p <- 1 + apply((1-Y)*(1-Theta_ex), 2, sum)
    e_p <- rbeta(K, a_p, b_p)
    
    # Theta
    psi_n <- log(e_n/(1-e_n))
    psi_p <- log(e_p/(1-e_p))
    for(k in 1:K){
      p1 <- Phi*exp( psi_n[k]*(1-Y[1:n,k])+psi_n[k]*(1-Y[-(1:n),k]) ) / (1+exp(psi_n[k]))^2
      p0 <- (1-Phi)*exp( psi_p[k]*Y[1:n,k] + psi_p[k]*Y[-(1:n),k] ) / (1+exp(psi_p[k]))^2
      prob <- p1 / (p1 + p0)
      Theta[,k] <- rbinom(n, 1, prob)
    }
    
    # Phi
    ss <- sum(Theta)
    Phi <- rbeta(1, alpha+ss, beta+(n*K-ss))
    
    # save 
    if(it>burn){
      cc <- it - burn
      e_p_pos[cc,] <- e_p
      e_n_pos[cc,] <- e_n
      Theta_pos[cc,,] <- Theta
      Phi_pos[cc] <- Phi
    }
    
    pb$tick()
  }
  
  # results 
  Result <- list(e_p=e_p_pos, e_n=e_n_pos, Theta=Theta_pos, Phi=Phi_pos)
  return(Result)
}

