###----------------------------------------------------###
###        Monte Carlo Simulation Study (K=6)          ###  
###----------------------------------------------------###
rm(list=ls())
set.seed(1)
source("HER-function.R")
source("ER-function.R")
library(PRROC)

# setting
sc <- 4      # 1-4
R <- 500     # Monte Carlo replications
K <- 6

# scenario
Phi_true <- seq(0.3, 0.8, length=K)
int <- (-1.5)

if(sc==1){
  gam_p_true <- c(int, 2, 2)
  gam_n_true <- c(int, -1.5, -1.5)
}

if(sc==2){
  gam_p_true <- c(int, 0, 0)
  gam_n_true <- c(int, 0, 0)
}

if(sc==3){
  gam_p_true <- c(int, 1, 1)
  gam_n_true <- c(int, -0.5, -0.5)
}

if(sc==4){
  gam_p_true <- c(int, 0, 0)
  gam_n_true <- c(int, 0, 0)
}

LGS <- function(x){ 1/(1+exp(-x)) }    # logistic function 


# settings 
N <- 50    
n <- N*(N-1)/2    # number of edges 
p <- 2      # number of covariates

ID_theta <- t(combn(1:N, 2))
ID <- rbind(ID_theta, ID_theta[,c(2,1)])


## Monte Carlo replications
Meth <- c("HER", "HER-np", "BER", "UNO", "ITS")
AB <- CP <- AL <- array(NA, c(R, K, 3))
dimnames(AB)[[3]] <- dimnames(CP)[[3]] <- dimnames(AL)[[3]] <- Meth[1:3]
AUC <- array(NA, c(R, K, 5))
dimnames(AUC)[[3]] <- Meth

for(r in 1:R){
  # covariate 
  X <- cbind(1, rnorm(2*n), rnorm(2*n))   # dyadic covariates
  
  # error prob
  e_p <- LGS( c(X%*%gam_p_true) )
  e_n <- LGS( c(X%*%gam_n_true) )
  if(sc==3 | sc==4){
    e_p <- LGS( c(X%*%gam_p_true) + 0.5*sin(X[,2]) + 0.5*X[,3]^2 )
    e_n <- LGS( c(X%*%gam_n_true) + 0.5*sin(X[,3]) + 0.5*X[,2]^2 )
  }
  
  # symmetric network (true)
  Theta_true <- matrix(NA, 2*n, K)   
  for(k in 1:K){
    rn <- rbinom(n, 1, Phi_true[k])
    Theta_true[,k] <- c(rn, rn)
  }
  
  # observed network 
  Y <- matrix(NA, 2*n, K)
  for(k in 1:K){
    prob <- Theta_true[,k]*(1-e_n) + (1-Theta_true[,k])*e_p
    Y[,k] <- rbinom(2*n, 1, prob)
  }
  
  # model fitting (HER)
  fit <- HER(Y, X, ID)   # proposed model
  Theta_HER <- apply(fit$Theta, c(2,3), mean)
  
  M <- 7
  knot1 <- quantile(X[,2], prob=seq(0.1, 0.9, length=M))
  knot2 <- quantile(X[,3], prob=seq(0.1, 0.9, length=M))
  Z <- cbind(X[,2]^2, X[,3]^2)
  Po <- function(xx){ xx*(xx>0) }
  for(k in 1:p){
    for(j in 1:M){
      Z <- cbind(Z, Po(X[,2]-knot1[j])^2)
    }
  }
  fit_np <- HER(Y, X, ID, Z=Z)   # proposed model with nonparametric regression
  Theta_HER_np <- apply(fit_np$Theta, c(2,3), mean)
  
  # model fitting (alternative methods)
  fit_BER <- ER(Y, ID)      # homogeneous model 
  Theta_BER <- apply(fit_BER$Theta, c(2,3), mean)
  
  Theta_UNO <- ifelse(Y[1:n,]+Y[-(1:n),]>0, 1, 0)    # Union
  Theta_ITS <- ifelse(Y[1:n,]*Y[-(1:n),]>0, 1, 0)    # Inersection
  
  # AUC
  true <- Theta_true[1:n,]
  for(sel in 1:K){
    ROC1 <- roc.curve(scores.class0=Theta_HER[true[,sel]==1,sel], scores.class1=Theta_HER[true[,sel]==0,sel])
    ROC2 <- roc.curve(scores.class0=Theta_HER_np[true[,sel]==1,sel], scores.class1=Theta_HER_np[true[,sel]==0,sel])
    ROC3 <- roc.curve(scores.class0=Theta_BER[true[,sel]==1,sel], scores.class1=Theta_BER[true[,sel]==0,sel])
    ROC4 <- roc.curve(scores.class0=Theta_UNO[true[,sel]==1,sel], scores.class1=Theta_UNO[true[,sel]==0,sel])
    ROC5 <- roc.curve(scores.class0=Theta_ITS[true[,sel]==1,sel], scores.class1=Theta_ITS[true[,sel]==0,sel])
    AUC[r,sel,] <- c(ROC1$auc, ROC2$auc, ROC3$auc, ROC4$auc, ROC5$auc)
  }
  
  # posterior mean of Phi
  hPhi_HER <- apply(fit$Phi, 2, mean)
  ad <- mean(Y) - mean(hPhi_HER)
  fit$Phi <- fit$Phi + ad
  hPhi_HER <- apply(fit$Phi, 2, mean)
  
  hPhi_HER_np <- apply(fit_np$Phi, 2, mean)
  hPhi_BER <- apply(fit_BER$Phi, 2, mean)
  
  # absolute bias
  AB[r,,1] <- hPhi_HER - Phi_true 
  AB[r,,2] <- hPhi_HER_np - Phi_true 
  AB[r,,3] <- hPhi_BER - Phi_true 
  
  # Credible interval 
  CI_HER <- apply(fit$Phi, 2, quantile, prob=c(0.025, 0.975))
  CI_HER_np <- apply(fit_np$Phi, 2, quantile, prob=c(0.025, 0.975))
  CI_BER <- apply(fit_BER$Phi, 2, quantile, prob=c(0.025, 0.975))
  
  # coverage and length 
  CP[r,,1] <- ifelse(CI_HER[1,]<Phi_true & CI_HER[2,]>Phi_true, 1, 0)
  CP[r,,2] <- ifelse(CI_HER_np[1,]<Phi_true & CI_HER_np[2,]>Phi_true, 1, 0)
  CP[r,,3] <- ifelse(CI_BER[1,]<Phi_true & CI_BER[2,]>Phi_true, 1, 0)
  AL[r,,1] <- CI_HER[2,] - CI_HER[1,]
  AL[r,,2] <- CI_HER_np[2,] - CI_HER_np[1,]
  AL[r,,3] <- CI_BER[2,] - CI_BER[1,]
  
  # print
  print(r)
  print( apply(CP[r,,], 2, mean) )
  print( apply(AUC[r,,], 2, mean) )
}



# save 
save(list=ls(), file=paste0("Result-s", sc, ".RData"))