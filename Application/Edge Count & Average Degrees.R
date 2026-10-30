## observed edge count or average degree
setwd("/Users/home/Library/CloudStorage/Application") 


## Village 1
load(paste0("Result", 1, "-np.RData"))
load(paste0("Result", 1, ".RData"))
relation_names <- c("nonrel", "rel", "templecompany")
N_nodes <- c(203, 203, 203)

Theta_UNO <- Theta_ITS <- list()
for(k in 1:3){
  n_k <- length(Y[[k]])/2
  Theta_UNO[[k]] <- ifelse(Y[[k]][1:n_k]+Y[[k]][-(1:n_k)]>0, 1, 0)    # Union
  Theta_ITS[[k]] <- ifelse(Y[[k]][1:n_k]*Y[[k]][-(1:n_k)]>0, 1, 0)    # Intersection
}

## For HER-p, HER-np, and BER:
## Theta is already a binary matrix:
## rows = posterior draws, columns = unordered dyads
posterior_network_summary <- function(theta, N) {
  
  # Number of edges in each posterior draw
  edges_draw <- rowSums(theta)
  
  # Posterior mean number of edges
  mean_edges <- mean(edges_draw)
  
  # Posterior mean average degree
  mean_degree <- mean(2 * edges_draw / N)
  
  c(
    Edges = mean_edges,
    Average_degree = mean_degree
  )
}

## UNO and ITS are already single binary reconstructed networks
binary_network_summary <- function(theta, N) {
  
  n_edges <- sum(theta)
  avg_degree <- 2 * n_edges / N
  
  c(
    Edges = n_edges,
    Average_degree = avg_degree
  )
}


## Collect results
results <- list()

for (k in 1:3) {
  
  results[[length(results) + 1]] <- data.frame(
    Relation = relation_names[k],
    Method = "HER-p",
    t(posterior_network_summary(
      fit_HER$Theta[[k]], N_nodes[k]
    ))
  )
  
  results[[length(results) + 1]] <- data.frame(
    Relation = relation_names[k],
    Method = "HER-np",
    t(posterior_network_summary(
      fit_HER_np$Theta[[k]], N_nodes[k]
    ))
  )
  
  results[[length(results) + 1]] <- data.frame(
    Relation = relation_names[k],
    Method = "BER",
    t(posterior_network_summary(
      fit_BER$Theta[[k]], N_nodes[k]
    ))
  )
  
  results[[length(results) + 1]] <- data.frame(
    Relation = relation_names[k],
    Method = "UNO",
    t(binary_network_summary(
      Theta_UNO[[k]], N_nodes[k]
    ))
  )
  
  results[[length(results) + 1]] <- data.frame(
    Relation = relation_names[k],
    Method = "ITS",
    t(binary_network_summary(
      Theta_ITS[[k]], N_nodes[k]
    ))
  )
}

results_v1 <- do.call(rbind, results)

results_v1$Edges <- round(results_v1$Edges, 1)
results_v1$Average_degree <- round(results_v1$Average_degree, 2)

print(results_v1)

## Village 11
load(paste0("Result", 11, "-np.RData"))
load(paste0("Result", 11, ".RData"))
relation_names <- c("nonrel", "rel", "templecompany")
N_nodes <- c(142, 142, 142)

Theta_UNO <- Theta_ITS <- list()
for(k in 1:3){
  n_k <- length(Y[[k]])/2
  Theta_UNO[[k]] <- ifelse(Y[[k]][1:n_k]+Y[[k]][-(1:n_k)]>0, 1, 0)    # Union
  Theta_ITS[[k]] <- ifelse(Y[[k]][1:n_k]*Y[[k]][-(1:n_k)]>0, 1, 0)    # Intersection
}

## For HER-p, HER-np, and BER:
## Theta is already a binary matrix:
## rows = posterior draws, columns = unordered dyads
posterior_network_summary <- function(theta, N) {
  
  # Number of edges in each posterior draw
  edges_draw <- rowSums(theta)
  
  # Posterior mean number of edges
  mean_edges <- mean(edges_draw)
  
  # Posterior mean average degree
  mean_degree <- mean(2 * edges_draw / N)
  
  c(
    Edges = mean_edges,
    Average_degree = mean_degree
  )
}

## UNO and ITS are already single binary reconstructed networks
binary_network_summary <- function(theta, N) {
  
  n_edges <- sum(theta)
  avg_degree <- 2 * n_edges / N
  
  c(
    Edges = n_edges,
    Average_degree = avg_degree
  )
}


## Collect results
results <- list()

for (k in 1:3) {
  
  results[[length(results) + 1]] <- data.frame(
    Relation = relation_names[k],
    Method = "HER-p",
    t(posterior_network_summary(
      fit_HER$Theta[[k]], N_nodes[k]
    ))
  )
  
  results[[length(results) + 1]] <- data.frame(
    Relation = relation_names[k],
    Method = "HER-np",
    t(posterior_network_summary(
      fit_HER_np$Theta[[k]], N_nodes[k]
    ))
  )
  
  results[[length(results) + 1]] <- data.frame(
    Relation = relation_names[k],
    Method = "BER",
    t(posterior_network_summary(
      fit_BER$Theta[[k]], N_nodes[k]
    ))
  )
  
  results[[length(results) + 1]] <- data.frame(
    Relation = relation_names[k],
    Method = "UNO",
    t(binary_network_summary(
      Theta_UNO[[k]], N_nodes[k]
    ))
  )
  
  results[[length(results) + 1]] <- data.frame(
    Relation = relation_names[k],
    Method = "ITS",
    t(binary_network_summary(
      Theta_ITS[[k]], N_nodes[k]
    ))
  )
}

results_v11 <- do.call(rbind, results)

results_v11$Edges <- round(results_v11$Edges, 1)
results_v11$Average_degree <- round(results_v11$Average_degree, 2)

print(results_v11)

## Village 14
load(paste0("Result", 14, "-np.RData"))
load(paste0("Result", 14, ".RData"))
relation_names <- c("nonrel", "rel", "templecompany")
N_nodes <- c(150, 150, 150)

Theta_UNO <- Theta_ITS <- list()
for(k in 1:3){
  n_k <- length(Y[[k]])/2
  Theta_UNO[[k]] <- ifelse(Y[[k]][1:n_k]+Y[[k]][-(1:n_k)]>0, 1, 0)    # Union
  Theta_ITS[[k]] <- ifelse(Y[[k]][1:n_k]*Y[[k]][-(1:n_k)]>0, 1, 0)    # Intersection
}

## For HER-p, HER-np, and BER:
## Theta is already a binary matrix:
## rows = posterior draws, columns = unordered dyads
posterior_network_summary <- function(theta, N) {
  
  # Number of edges in each posterior draw
  edges_draw <- rowSums(theta)
  
  # Posterior mean number of edges
  mean_edges <- mean(edges_draw)
  
  # Posterior mean average degree
  mean_degree <- mean(2 * edges_draw / N)
  
  c(
    Edges = mean_edges,
    Average_degree = mean_degree
  )
}

## UNO and ITS are already single binary reconstructed networks
binary_network_summary <- function(theta, N) {
  
  n_edges <- sum(theta)
  avg_degree <- 2 * n_edges / N
  
  c(
    Edges = n_edges,
    Average_degree = avg_degree
  )
}


## Collect results
results <- list()

for (k in 1:3) {
  
  results[[length(results) + 1]] <- data.frame(
    Relation = relation_names[k],
    Method = "HER-p",
    t(posterior_network_summary(
      fit_HER$Theta[[k]], N_nodes[k]
    ))
  )
  
  results[[length(results) + 1]] <- data.frame(
    Relation = relation_names[k],
    Method = "HER-np",
    t(posterior_network_summary(
      fit_HER_np$Theta[[k]], N_nodes[k]
    ))
  )
  
  results[[length(results) + 1]] <- data.frame(
    Relation = relation_names[k],
    Method = "BER",
    t(posterior_network_summary(
      fit_BER$Theta[[k]], N_nodes[k]
    ))
  )
  
  results[[length(results) + 1]] <- data.frame(
    Relation = relation_names[k],
    Method = "UNO",
    t(binary_network_summary(
      Theta_UNO[[k]], N_nodes[k]
    ))
  )
  
  results[[length(results) + 1]] <- data.frame(
    Relation = relation_names[k],
    Method = "ITS",
    t(binary_network_summary(
      Theta_ITS[[k]], N_nodes[k]
    ))
  )
}

results_v14 <- do.call(rbind, results)

results_v14$Edges <- round(results_v14$Edges, 1)
results_v14$Average_degree <- round(results_v14$Average_degree, 2)

print(results_v14)

## Village 15
load(paste0("Result", 15, "-np.RData"))
load(paste0("Result", 15, ".RData"))
relation_names <- c("nonrel", "rel", "templecompany")
N_nodes <- c(212, 212, 212)

Theta_UNO <- Theta_ITS <- list()
for(k in 1:3){
  n_k <- length(Y[[k]])/2
  Theta_UNO[[k]] <- ifelse(Y[[k]][1:n_k]+Y[[k]][-(1:n_k)]>0, 1, 0)    # Union
  Theta_ITS[[k]] <- ifelse(Y[[k]][1:n_k]*Y[[k]][-(1:n_k)]>0, 1, 0)    # Intersection
}

## For HER-p, HER-np, and BER:
## Theta is already a binary matrix:
## rows = posterior draws, columns = unordered dyads
posterior_network_summary <- function(theta, N) {
  
  # Number of edges in each posterior draw
  edges_draw <- rowSums(theta)
  
  # Posterior mean number of edges
  mean_edges <- mean(edges_draw)
  
  # Posterior mean average degree
  mean_degree <- mean(2 * edges_draw / N)
  
  c(
    Edges = mean_edges,
    Average_degree = mean_degree
  )
}

## UNO and ITS are already single binary reconstructed networks
binary_network_summary <- function(theta, N) {
  
  n_edges <- sum(theta)
  avg_degree <- 2 * n_edges / N
  
  c(
    Edges = n_edges,
    Average_degree = avg_degree
  )
}


## Collect results
results <- list()

for (k in 1:3) {
  
  results[[length(results) + 1]] <- data.frame(
    Relation = relation_names[k],
    Method = "HER-p",
    t(posterior_network_summary(
      fit_HER$Theta[[k]], N_nodes[k]
    ))
  )
  
  results[[length(results) + 1]] <- data.frame(
    Relation = relation_names[k],
    Method = "HER-np",
    t(posterior_network_summary(
      fit_HER_np$Theta[[k]], N_nodes[k]
    ))
  )
  
  results[[length(results) + 1]] <- data.frame(
    Relation = relation_names[k],
    Method = "BER",
    t(posterior_network_summary(
      fit_BER$Theta[[k]], N_nodes[k]
    ))
  )
  
  results[[length(results) + 1]] <- data.frame(
    Relation = relation_names[k],
    Method = "UNO",
    t(binary_network_summary(
      Theta_UNO[[k]], N_nodes[k]
    ))
  )
  
  results[[length(results) + 1]] <- data.frame(
    Relation = relation_names[k],
    Method = "ITS",
    t(binary_network_summary(
      Theta_ITS[[k]], N_nodes[k]
    ))
  )
}

results_v15 <- do.call(rbind, results)

results_v15$Edges <- round(results_v15$Edges, 1)
results_v15$Average_degree <- round(results_v15$Average_degree, 2)

print(results_v15)
