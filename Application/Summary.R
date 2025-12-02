rm(list=ls())

## Result 
num <- 15    # you can select from 1, 11, 14, 15
load(paste0("Result", num, "-np.RData"))
load(paste0("Result", num, ".RData"))
  

## results of intersection and union 
Theta_UNO <- Theta_ITS <- list()
for(k in 1:3){
  n_k <- length(Y[[k]])/2
  Theta_UNO[[k]] <- ifelse(Y[[k]][1:n_k]+Y[[k]][-(1:n_k)]>0, 1, 0)    # Union
  Theta_ITS[[k]] <- ifelse(Y[[k]][1:n_k]*Y[[k]][-(1:n_k)]>0, 1, 0)    # Intersection
}


## probability of link (three layers)
Link_prob <- matrix(NA, 6, 3)
dimnames(Link_prob)[[1]] <- c("HER-p", "HER-np", "Butts", "Obs", "Union", "Intersection")
dimnames(Link_prob)[[2]] <- paste0("Layer", 1:3)

for(k in 1:3){
  Link_prob[1,k] <- mean(apply(fit_HER$Theta[[k]]>0.5, 2, mean))      # proposed model (parametric)
  Link_prob[2,k] <- mean(apply(fit_HER_np$Theta[[k]]>0.5, 2, mean))   # proposed model (non-parametric)
  Link_prob[3,k] <- mean(apply(fit_BER$Theta[[k]]>0.5, 2, mean))      # Butts model 
  Link_prob[4,k] <- mean(Y[[k]])  
  Link_prob[5,k] <- mean(Theta_UNO[[k]])
  Link_prob[6,k] <- mean(Theta_ITS[[k]])
}

10^3*Link_prob




## error probability
# posterior 
LGS <- function(x){ 1/(1+exp(-x)) }
HER_ep_pos <- list()
HER_ep_pos[[1]] <- LGS( X[[1]]%*%t(fit_HER$gam_p) )
HER_ep_pos[[2]] <- LGS( X[[2]]%*%t(fit_HER$gam_p) )
HER_ep_pos[[3]] <- LGS( X[[3]]%*%t(fit_HER$gam_p) )

HER_en_pos <- list()
HER_en_pos[[1]] <- LGS( X[[1]]%*%t(fit_HER$gam_n) )
HER_en_pos[[2]] <- LGS( X[[2]]%*%t(fit_HER$gam_n) )
HER_en_pos[[3]] <- LGS( X[[3]]%*%t(fit_HER$gam_n) )

# posterior mean
HER_ep_pm <- HER_en_pm <- list()
for(k in 1:3){
  HER_ep_pm[[k]] <- apply(HER_ep_pos[[k]], 1, mean)
  HER_en_pm[[k]] <- apply(HER_en_pos[[k]], 1, mean)
}

plot(HER_ep_pm[[1]])
mean(fit_BER$e_p)


## non-linear effects of age
k <- 3
head(X[[k]])
XZ1 <- cbind(X[[k]][,3], Z[[k]][,c(1,3:(M+2))])
XZ2 <- cbind(X[[k]][,5], Z[[k]][,c(2,(M+3):(2*M+2))])
reg_pos1 <- fit_HER_np$gam_n[,c(3,7,8:(7+M))]%*%t(XZ1)
reg_pos2 <- fit_HER_np$gam_n[,c(5,9,(8+M):(7+2*M))]%*%t(XZ2)

plot(XZ1[,1], apply(reg_pos1, 2, mean))
plot(XZ2[,1], apply(reg_pos2, 2, mean))





## posterior summary of coefficients in error probability
# covariates are (intercept, gender1, age1, gender2, age2)
# 1 and 2 indicate the subject at the start and end of a directed arrow
apply(fit_HER$gam_n, 2, quantile, prob=c(0.05, 0.25, 0.5, 0.75, 0.95))
apply(fit_HER$gam_p, 2, quantile, prob=c(0.05, 0.25, 0.5, 0.75, 0.95))

apply(fit_HER_np$gam_n[,1:5], 2, quantile, prob=c(0.05, 0.25, 0.5, 0.75, 0.95))
apply(fit_HER_np$gam_p[,1:5], 2, quantile, prob=c(0.05, 0.25, 0.5, 0.75, 0.95))

## posterior densities of error rates
par(mfrow=c(1,2)) ## Give the number of village based on "num"
boxplot(c(HER_ep_pos[[1]]),main="nonrel (Village -) \n HER model",col="cyan",xlab="e+",
        cex.main=1.5,cex.lab=1.5,cex.axis=1.5)
boxplot(c(HER_en_pos[[1]]),main="nonrel (Village -) \n HER model",col="cyan",xlab="e-",
        cex.main=1.5,cex.lab=1.5,cex.axis=1.5)

boxplot(c(HER_ep_pos[[2]]),main="rel (Village -) \n HER model",col="cyan",xlab="e+",
        cex.main=1.5,cex.lab=1.5,cex.axis=1.5) 
boxplot(c(HER_en_pos[[2]]),main="rel (Village -) \n HER model",col="cyan",xlab="e-",
        cex.main=1.5,cex.lab=1.5,cex.axis=1.5)

boxplot(c(HER_ep_pos[[3]]),main="templecompany (Village -) \n HER model",col="cyan",xlab="e+",
        cex.main=1.5,cex.lab=1.5,cex.axis=1.5) 
boxplot(c(HER_en_pos[[3]]),main="templecompany (Village -) \n HER model",col="cyan",xlab="e-",
        cex.main=1.5,cex.lab=1.5,cex.axis=1.5)

par(mfrow=c(1,2))
boxplot(fit_BER$e_p,main="All relations (Village -) \n Butts model",col="cyan",xlab="e+",
        cex.main=1.5,cex.lab=1.5,cex.axis=1.5)
boxplot(fit_BER$e_n,main="All relations (Village -) \n Butts model",col="cyan",xlab="e-",
        cex.main=1.5,cex.lab=1.5,cex.axis=1.5)
