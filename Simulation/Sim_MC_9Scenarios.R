###----------------------------------------------------###
###  Monte Carlo Simulation Study for Nine scenarios   ###
###----------------------------------------------------###

rm(list=ls())
set.seed(1)

source("HER-function.R")
source("ER-function.R")

library(PRROC)

R <- 50                 # Monte Carlo replications
K <- 6                   # number of relations
N <- 50                  # number of actors (N=50, 100, 300)
M <- 7                   # number of knots for HER-np
draw <- 5000             # MCMC iterations
burn <- 2000             # burn-in iterations

scenarios_to_run <- 3:5  # change this to run only some scenarios (from 1 to 9 in total)

###-----------------------------###
### Basic functions             ###
###-----------------------------###

LGS <- function(x){ 1/(1+exp(-x)) }    # logistic function
Po <- function(xx){ xx*(xx>0) }

F1_fun <- function(true, pred){
  true <- as.vector(true)
  pred <- as.vector(pred)
  TP <- sum(true==1 & pred==1)
  FP <- sum(true==0 & pred==1)
  FN <- sum(true==1 & pred==0)
  if((2*TP + FP + FN)==0){ return(NA) }
  return(2*TP/(2*TP + FP + FN))
}

DE_fun <- function(true, pred){
  true <- as.vector(true)
  pred <- as.vector(pred)
  return(abs(mean(pred) - mean(true)))
}

make_spline_Z <- function(W_cont, M=7){
  # W_cont is a matrix with continuous covariates.
  # For each column v, include v^2 and (v-kappa)^2_+ terms.
  Z <- NULL
  for(ell in 1:ncol(W_cont)){
    v <- W_cont[, ell]
    knots <- quantile(v, prob=seq(0.1, 0.9, length=M), na.rm=TRUE)
    Z <- cbind(Z, v^2)
    for(j in 1:M){
      Z <- cbind(Z, Po(v - knots[j])^2)
    }
  }
  return(Z)
}

get_scenario <- function(sc){
  int <- -1.5
  gam_p <- rep(0, 7)
  gam_n <- rep(0, 7)
  gam_p[1] <- int
  gam_n[1] <- int
  node_nonlinear <- FALSE
  dyad_nonlinear <- FALSE

  if(sc==1){
    # Constant error probabilities
  }
  if(sc==2){
    # Discrete node-level effects only
    gam_p[c(2,4)] <- 2
    gam_n[c(2,4)] <- -1.5
  }
  if(sc==3){
    # Continuous node-level effects only
    gam_p[c(3,5)] <- 2
    gam_n[c(3,5)] <- -1.5
  }
  if(sc==4){
    # Linear node-level effects, both discrete and continuous
    gam_p[c(2,4)] <- 2
    gam_n[c(2,4)] <- -1.5
    gam_p[c(3,5)] <- 1
    gam_n[c(3,5)] <- -0.5
  }
  if(sc==5){
    # Linear node-level effects plus nonlinear node-level effects
    gam_p[c(2,4)] <- 1
    gam_n[c(2,4)] <- -0.5
    gam_p[c(3,5)] <- 1
    gam_n[c(3,5)] <- -0.5
    node_nonlinear <- TRUE
  }
  if(sc==6){
    # Linear dyadic-level effects only
    gam_p[c(6,7)] <- 2
    gam_n[c(6,7)] <- -1.5
  }
  if(sc==7){
    # Linear dyadic-level effects plus nonlinear dyadic-level effects
    gam_p[c(6,7)] <- 1
    gam_n[c(6,7)] <- -0.5
    dyad_nonlinear <- TRUE
  }
  if(sc==8){
    # Linear node-level and dyadic-level effects
    gam_p[c(2,4)] <- 2
    gam_n[c(2,4)] <- -1.5
    gam_p[c(3,5)] <- 1
    gam_n[c(3,5)] <- -0.5
    gam_p[c(6,7)] <- 2
    gam_n[c(6,7)] <- -1.5
  }
  if(sc==9){
    # Linear and nonlinear node-level and dyadic-level effects
    gam_p[c(2,4)] <- 1
    gam_n[c(2,4)] <- -0.5
    gam_p[c(3,5)] <- 1
    gam_n[c(3,5)] <- -0.5
    gam_p[c(6,7)] <- 0.5
    gam_n[c(6,7)] <- -0.5
    node_nonlinear <- TRUE
    dyad_nonlinear <- TRUE
  }
  return(list(gam_p=gam_p, gam_n=gam_n,
              node_nonlinear=node_nonlinear,
              dyad_nonlinear=dyad_nonlinear))
}

###-----------------------------###
### Network and ID settings     ###
###-----------------------------###

n <- N*(N-1)/2
ID_theta <- t(combn(1:N, 2))
ID <- rbind(ID_theta, ID_theta[,c(2,1)])
Phi_true <- seq(0.3, 0.8, length=K)
Meth <- c("HER-p", "HER-np", "BER", "UNO", "ITS")

###-----------------------------###
### Monte Carlo simulations     ###
###-----------------------------###

for(sc in scenarios_to_run){

  cat("\n=====================================\n")
  cat("Running Scenario", sc, "\n")
  cat("=====================================\n")

  scen <- get_scenario(sc)
  gam_p_true <- scen$gam_p
  gam_n_true <- scen$gam_n

  AB <- CP <- AL <- array(NA, c(R, K, 3))
  dimnames(AB)[[3]] <- dimnames(CP)[[3]] <- dimnames(AL)[[3]] <- Meth[1:3]

  AUC <- F1 <- DE <- array(NA, c(R, K, 5))
  dimnames(AUC)[[3]] <- dimnames(F1)[[3]] <- dimnames(DE)[[3]] <- Meth

  Runtime <- matrix(NA, R, 5)
  colnames(Runtime) <- Meth

  # ROC curves
  ROC_rep1 <- vector("list", K)
  for(k in 1:K){
    ROC_rep1[[k]] <- vector("list", length(Meth))
    names(ROC_rep1[[k]]) <- Meth
  }

  for(r in 1:R){

    ### Generate covariates
    w1 <- rbinom(N, 1, 0.75)   # binary node-level covariate
    w2 <- rnorm(N)             # continuous node-level covariate

    wi1 <- w1[ID[,1]]
    wi2 <- w2[ID[,1]]
    wj1 <- w1[ID[,2]]
    wj2 <- w2[ID[,2]]

    z1 <- rnorm(2*n)           # dyadic covariate 1
    z2 <- rnorm(2*n)           # dyadic covariate 2

    X <- cbind(1, wi1, wi2, wj1, wj2, z1, z2)

    ### Generate error probabilities
    eta_p <- c(X%*%gam_p_true)
    eta_n <- c(X%*%gam_n_true)

    if(scen$node_nonlinear){
      eta_p <- eta_p + 0.5*sin(wi2) + 0.5*wj2^2
      eta_n <- eta_n + 0.5*sin(wj2) + 0.5*wi2^2
    }
    if(scen$dyad_nonlinear){
      eta_p <- eta_p + 0.5*sin(z1) + 0.5*z2^2
      eta_n <- eta_n + 0.5*sin(z2) + 0.5*z1^2
    }

    e_p <- LGS(eta_p)
    e_n <- LGS(eta_n)

    ### Symmetric true network
    Theta_true <- matrix(NA, 2*n, K)
    for(k in 1:K){
      rn <- rbinom(n, 1, Phi_true[k])
      Theta_true[,k] <- c(rn, rn)
    }

    ### Observed reports
    Y <- matrix(NA, 2*n, K)
    for(k in 1:K){
      prob <- Theta_true[,k]*(1-e_n) + (1-Theta_true[,k])*e_p
      Y[,k] <- rbinom(2*n, 1, prob)
    }

    true <- Theta_true[1:n,]

    ### HER-p
    time0 <- Sys.time()
    fit <- HER(Y, X, ID, draw=draw, burn=burn)
    Runtime[r, "HER-p"] <- as.numeric(difftime(Sys.time(), time0, units="secs"))
    Theta_HER <- apply(fit$Theta, c(2,3), mean)

    ### HER-np
    W_cont <- cbind(wi2, wj2, z1, z2)
    Z <- make_spline_Z(W_cont, M=M)

    time0 <- Sys.time()
    fit_np <- HER(Y, X, ID, Z=Z, draw=draw, burn=burn)
    Runtime[r, "HER-np"] <- as.numeric(difftime(Sys.time(), time0, units="secs"))
    Theta_HER_np <- apply(fit_np$Theta, c(2,3), mean)

    ### BER
    time0 <- Sys.time()
    fit_BER <- ER(Y, ID, draw=draw, burn=burn)
    Runtime[r, "BER"] <- as.numeric(difftime(Sys.time(), time0, units="secs"))
    Theta_BER <- apply(fit_BER$Theta, c(2,3), mean)

    ### UNO and ITS
    time0 <- Sys.time()
    Theta_UNO <- ifelse(Y[1:n,]+Y[-(1:n),]>0, 1, 0)
    Runtime[r, "UNO"] <- as.numeric(difftime(Sys.time(), time0, units="secs"))

    time0 <- Sys.time()
    Theta_ITS <- ifelse(Y[1:n,]*Y[-(1:n),]>0, 1, 0)
    Runtime[r, "ITS"] <- as.numeric(difftime(Sys.time(), time0, units="secs"))

    ### Metrics
    Theta_HER_bin    <- ifelse(Theta_HER > 0.5, 1, 0)
    Theta_HER_np_bin <- ifelse(Theta_HER_np > 0.5, 1, 0)
    Theta_BER_bin    <- ifelse(Theta_BER > 0.5, 1, 0)
    Theta_UNO_bin    <- Theta_UNO
    Theta_ITS_bin    <- Theta_ITS

    for(sel in 1:K){
      ROC1 <- roc.curve(scores.class0=Theta_HER[true[,sel]==1,sel],
                        scores.class1=Theta_HER[true[,sel]==0,sel], curve=TRUE)
      ROC2 <- roc.curve(scores.class0=Theta_HER_np[true[,sel]==1,sel],
                        scores.class1=Theta_HER_np[true[,sel]==0,sel], curve=TRUE)
      ROC3 <- roc.curve(scores.class0=Theta_BER[true[,sel]==1,sel],
                        scores.class1=Theta_BER[true[,sel]==0,sel], curve=TRUE)
      ROC4 <- roc.curve(scores.class0=Theta_UNO[true[,sel]==1,sel],
                        scores.class1=Theta_UNO[true[,sel]==0,sel], curve=TRUE)
      ROC5 <- roc.curve(scores.class0=Theta_ITS[true[,sel]==1,sel],
                        scores.class1=Theta_ITS[true[,sel]==0,sel], curve=TRUE)

      AUC[r,sel,] <- c(ROC1$auc, ROC2$auc, ROC3$auc, ROC4$auc, ROC5$auc)

      if(r==1){
        ROC_rep1[[sel]][["HER-p"]]  <- ROC1$curve
        ROC_rep1[[sel]][["HER-np"]] <- ROC2$curve
        ROC_rep1[[sel]][["BER"]]    <- ROC3$curve
        ROC_rep1[[sel]][["UNO"]]    <- ROC4$curve
        ROC_rep1[[sel]][["ITS"]]    <- ROC5$curve
      }

      F1[r,sel,] <- c(F1_fun(true[,sel], Theta_HER_bin[,sel]),
                      F1_fun(true[,sel], Theta_HER_np_bin[,sel]),
                      F1_fun(true[,sel], Theta_BER_bin[,sel]),
                      F1_fun(true[,sel], Theta_UNO_bin[,sel]),
                      F1_fun(true[,sel], Theta_ITS_bin[,sel]))

      DE[r,sel,] <- c(DE_fun(true[,sel], Theta_HER_bin[,sel]),
                      DE_fun(true[,sel], Theta_HER_np_bin[,sel]),
                      DE_fun(true[,sel], Theta_BER_bin[,sel]),
                      DE_fun(true[,sel], Theta_UNO_bin[,sel]),
                      DE_fun(true[,sel], Theta_ITS_bin[,sel]))
    }

    ### Existing summaries for Phi
    hPhi_HER <- apply(fit$Phi, 2, mean)
    ad <- mean(Y) - mean(hPhi_HER)
    fit$Phi <- fit$Phi + ad
    hPhi_HER <- apply(fit$Phi, 2, mean)
    hPhi_HER_np <- apply(fit_np$Phi, 2, mean)
    hPhi_BER <- apply(fit_BER$Phi, 2, mean)

    AB[r,,1] <- hPhi_HER - Phi_true
    AB[r,,2] <- hPhi_HER_np - Phi_true
    AB[r,,3] <- hPhi_BER - Phi_true

    CI_HER <- apply(fit$Phi, 2, quantile, prob=c(0.025, 0.975))
    CI_HER_np <- apply(fit_np$Phi, 2, quantile, prob=c(0.025, 0.975))
    CI_BER <- apply(fit_BER$Phi, 2, quantile, prob=c(0.025, 0.975))

    CP[r,,1] <- ifelse(CI_HER[1,]<Phi_true & CI_HER[2,]>Phi_true, 1, 0)
    CP[r,,2] <- ifelse(CI_HER_np[1,]<Phi_true & CI_HER_np[2,]>Phi_true, 1, 0)
    CP[r,,3] <- ifelse(CI_BER[1,]<Phi_true & CI_BER[2,]>Phi_true, 1, 0)
    AL[r,,1] <- CI_HER[2,] - CI_HER[1,]
    AL[r,,2] <- CI_HER_np[2,] - CI_HER_np[1,]
    AL[r,,3] <- CI_BER[2,] - CI_BER[1,]

    cat("Scenario:", sc, "Replication:", r, "\n")
    print(apply(AUC[r,,], 2, mean))
    print(apply(F1[r,,], 2, mean))
    print(apply(DE[r,,], 2, mean))
    print(Runtime[r,])
  }

  if(K==6){
    save(list=ls(), file=paste0("Result-s", sc, "-N", N, ".RData"))
  }else{
    save(list=ls(), file=paste0("Result-s", sc, "-N", N, "-K", K, ".RData"))
  }
  cat("Saved Result-s", sc, ".RData\n", sep="")
}
