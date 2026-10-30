rm(list=ls())

num <- 15   # 1, 11, 14, 15

## read related adjacency matrices and covariates
load(file=paste0("Vil", num, ".RData"))

## network
A <- list()
A[[1]] <- Matrix_nonrel
A[[2]] <- Matrix_rel
A[[3]] <- Matrix_templecompany

for(k in 1:3){
  diag(A[[k]]) <- 0
}

## ------------------------------------------------------------
## Keep ALL village nodes, including isolates.
## The SAME full node set is used for nonrel, rel, and
## templecompany.
## ------------------------------------------------------------
Sub <- rep(TRUE, nrow(A[[1]]))

cat("Village", num, "\n")
cat("Nodes retained (including isolates):", sum(Sub), "\n")

## relation-specific covariates for all nodes
COV <- list(COV_nonrel, COV_rel, COV_templecompany)
X_ind <- vector("list", 3)

for(k in 1:3){
  X_ind[[k]] <- cbind(1, subset(COV[[k]], select=c("gender", "age")))[Sub, , drop=FALSE]
  X_ind[[k]][,2] <- ifelse(X_ind[[k]][,2] == "Male", 1, 0)
  X_ind[[k]][,3] <- scale(as.numeric(X_ind[[k]][,3]))
  X_ind[[k]] <- as.matrix(X_ind[[k]])
}

## data dimensions
K <- length(A)
ID <- vector("list", K)
n <- numeric(K)
N <- sum(Sub)

for(k in 1:K){
  n[k] <- N*(N-1)/2       # number of unordered dyads
  ID_theta <- t(combn(1:N, 2))
  ID[[k]] <- rbind(ID_theta, ID_theta[,c(2,1)])
}

## observed network
A_use <- lapply(A, function(M) M[Sub, Sub, drop=FALSE])
Y <- vector("list", K)
for(k in 1:K){
  Y[[k]] <- A_use[[k]][cbind(ID[[k]][,1], ID[[k]][,2])]
}

print(lapply(Y, sum, na.rm=TRUE))

## dyad-level design matrices
p <- ncol(X_ind[[1]])
X <- vector("list", K)
for(k in 1:K){
  X[[k]] <- matrix(NA_real_, 2*n[k], 2*p-1)
  for(i in 1:(2*n[k])){
    X[[k]][i,] <- c(X_ind[[k]][ID[[k]][i,1], 1:p],
                    X_ind[[k]][ID[[k]][i,2], -1])
  }
}

## Keep the same object names and Dataset<num>.RData filename expected
## by the downstream analysis scripts.
save(Y, X, K, ID, file=paste0("Dataset", num, ".RData"))
