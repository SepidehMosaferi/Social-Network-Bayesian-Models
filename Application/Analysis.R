rm(list=ls())

source("HER-function.R")
source("ER-function.R")

num <- 15  # 1, 11, 14, 15

## read dataset
load(paste0("Dataset", num, ".RData"))
dim(Y)

## model fitting (HER)
fit_HER <- HER(Y, X, ID, Phi_est=F)   # proposed model

# model fitting (alternative methods)
fit_BER <- ER(Y, ID)      # homogeneous model 


## save 
save(list=ls(), file=paste0("Result", num, ".RData"))
