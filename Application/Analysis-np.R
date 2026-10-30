rm(list=ls())

source("HER-function.R")
source("ER-function.R")

num <- 15    # 1, 11, 14, 15

## read dataset
load(paste0("Dataset", num, ".RData"))
lapply(Y, length)
#head(X)

## model fitting 
K <- length(X)
M <- 7     # number of knots 
Z <- list()
for(k in 1:K){
  XX <- X[[k]][,c(3,5)]
  knot1 <- quantile(XX[,1], prob=seq(0.1, 0.9, length=M))
  knot2 <- quantile(XX[,2], prob=seq(0.1, 0.9, length=M))
  ZZ <- cbind(XX[,1]^2, XX[,2]^2)
  Po <- function(xx){ xx*(xx>0) }
  for(sel in 1:2){
    for(j in 1:M){
      ZZ <- cbind(ZZ, Po(XX[,sel]-knot1[j])^2)
    }
  }
  Z[[k]] <- ZZ
}

fit_HER_np <- HER(Y, X, ID, Z=Z, Phi_est=F)   # proposed model with nonparametric terms 

## save 
save(list=ls(), file=paste0("Result", num, "-np.RData"))
