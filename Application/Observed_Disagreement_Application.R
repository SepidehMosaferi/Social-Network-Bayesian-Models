###----------------------------------------------------###
### Observed disagreement in self-reported networks
###----------------------------------------------------###

rm(list = ls())

relation_names <- c("nonrel", "rel", "templecompany")

disagreement_one_village <- function(Y, village){

  out <- NULL

  for(k in 1:3){

    n_k <- length(Y[[k]]) / 2

    if(n_k != floor(n_k)){
      stop(paste0("Y[[", k, "]] does not have even length."))
    }

    Y_ij <- Y[[k]][1:n_k]
    Y_ji <- Y[[k]][(n_k + 1):(2 * n_k)]

    keep <- !is.na(Y_ij) & !is.na(Y_ji)
    Y_ij <- Y_ij[keep]
    Y_ji <- Y_ji[keep]

    n_used <- length(Y_ij)

    n_00 <- sum(Y_ij == 0 & Y_ji == 0)
    n_11 <- sum(Y_ij == 1 & Y_ji == 1)
    n_10 <- sum(Y_ij == 1 & Y_ji == 0)
    n_01 <- sum(Y_ij == 0 & Y_ji == 1)

    n_disagree <- n_10 + n_01
    D_k <- n_disagree / n_used

    out <- rbind(
      out,
      data.frame(
        Village = village,
        Relation = relation_names[k],
        N_dyads = n_used,
        Both_absent_n = n_00,
        Both_present_n = n_11,
        Disagreement_n = n_disagree,
        Both_absent_prop = n_00 / n_used,
        Both_present_prop = n_11 / n_used,
        Disagreement_prop = D_k,
        Yij1_Yji0_prop = n_10 / n_used,
        Yij0_Yji1_prop = n_01 / n_used
      )
    )
  }

  return(out)
}


### All four villages
villages <- c(1, 11, 14, 15)

all_disagreement <- NULL

for(v in villages){

  cat("Processing village", v, "\n")

  E <- new.env()
  load(paste0("Result", v, ".RData"), envir = E)

  if(!exists("Y", envir = E, inherits = FALSE)){
    stop(paste0("Object Y was not found in Result", v, ".RData"))
  }

  tmp <- disagreement_one_village(
    Y = E$Y,
    village = v
  )

  all_disagreement <- rbind(all_disagreement, tmp)

  rm(E, tmp)
  gc()
}

print(all_disagreement)

write.csv(
  all_disagreement,
  "Observed_disagreement_all_villages.csv",
  row.names = FALSE
)

### Compact table in percentages
disagreement_table <- data.frame(
  Village = all_disagreement$Village,
  Relation = all_disagreement$Relation,
  Dyads = all_disagreement$N_dyads,
  Both_absent_percent =
    round(100 * all_disagreement$Both_absent_prop, 2),
  Both_present_percent =
    round(100 * all_disagreement$Both_present_prop, 2),
  Disagreement_percent =
    round(100 * all_disagreement$Disagreement_prop, 2)
)

print(disagreement_table)

write.csv(
  disagreement_table,
  "Observed_disagreement_table.csv",
  row.names = FALSE
)


### calculating the number of used nodes
K <- length(Y)

for (k in 1:K) {
  n_dyads <- length(Y[[k]]) / 2
  N_nodes <- (1 + sqrt(1 + 8 * n_dyads)) / 2
  
  cat("Relation", k,
      ": Dyads =", n_dyads,
      ", Nodes =", N_nodes, "\n")
}
