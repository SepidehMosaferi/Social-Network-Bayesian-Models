###--------------------------------------------------------------------###
###  Prior-variance sensitivity study for HER-np                       ###
###--------------------------------------------------------------------###

rm(list = ls())
set.seed(24680)

## Required packages
req_pkg <- c("pgdraw", "MCMCpack", "progress", "PRROC", "ggplot2", "reshape2")
for(pkg in req_pkg){
  if(!requireNamespace(pkg, quietly = TRUE)){
    install.packages(pkg)
  }
  library(pkg, character.only = TRUE)
}

### User settings
R <- 20                         # Monte Carlo replications for sensitivity analysis
K <- 6                          # number of relations
N <- 50                         # number of actors
scenarios <- c(3, 4)            # nonlinear scenarios
M <- 7                          # default number of knots used in the paper
prior_var_grid <- c(25, 100, 400)  # B0 = prior_var * I for gamma^+ and gamma^-
save_prefix <- "PriorVarianceSensitivity"

## MCMC settings for sensitivity analysis
draw <- 2000
burn <- 1000

### Helper functions
LGS <- function(x){ 1/(1 + exp(-x)) }
Po <- function(xx){ xx * (xx > 0) }

make_Z <- function(X, M){
  ## X has columns: intercept, x1, x2.
  ## Z includes quadratic terms and truncated quadratic spline terms for both covariates.
  x1 <- X[, 2]
  x2 <- X[, 3]
  knot1 <- as.numeric(quantile(x1, prob = seq(0.1, 0.9, length.out = M)))
  knot2 <- as.numeric(quantile(x2, prob = seq(0.1, 0.9, length.out = M)))

  Z <- cbind(x1^2, x2^2)
  colnames(Z) <- c("x1_sq", "x2_sq")

  for(j in seq_len(M)){
    Z <- cbind(Z, Po(x1 - knot1[j])^2)
    colnames(Z)[ncol(Z)] <- paste0("x1_knot", j)
  }
  for(j in seq_len(M)){
    Z <- cbind(Z, Po(x2 - knot2[j])^2)
    colnames(Z)[ncol(Z)] <- paste0("x2_knot", j)
  }
  return(Z)
}

auc_one <- function(score, truth){
  out <- PRROC::roc.curve(scores.class0 = score[truth == 1],
                          scores.class1 = score[truth == 0],
                          curve = FALSE)
  as.numeric(out$auc)
}

f1_one <- function(score, truth, cutoff = 0.5){
  pred <- as.integer(score >= cutoff)
  tp <- sum(pred == 1 & truth == 1)
  fp <- sum(pred == 1 & truth == 0)
  fn <- sum(pred == 0 & truth == 1)
  denom <- 2 * tp + fp + fn
  if(denom == 0) return(NA_real_)
  2 * tp / denom
}

get_truth_and_data <- function(sc, K, N){
  n <- N * (N - 1) / 2
  Phi_true <- seq(0.3, 0.8, length.out = K)
  int <- -1.5

  if(sc == 1){
    gam_p_true <- c(int, 2, 2)
    gam_n_true <- c(int, -1.5, -1.5)
  }
  if(sc == 2){
    gam_p_true <- c(int, 0, 0)
    gam_n_true <- c(int, 0, 0)
  }
  if(sc == 3){
    gam_p_true <- c(int, 1, 1)
    gam_n_true <- c(int, -0.5, -0.5)
  }
  if(sc == 4){
    gam_p_true <- c(int, 0, 0)
    gam_n_true <- c(int, 0, 0)
  }

  X <- cbind(1, rnorm(2 * n), rnorm(2 * n))
  colnames(X) <- c("Intercept", "x1", "x2")

  e_p <- LGS(c(X %*% gam_p_true))
  e_n <- LGS(c(X %*% gam_n_true))
  if(sc %in% c(3, 4)){
    e_p <- LGS(c(X %*% gam_p_true) + 0.5 * sin(X[, 2]) + 0.5 * X[, 3]^2)
    e_n <- LGS(c(X %*% gam_n_true) + 0.5 * sin(X[, 3]) + 0.5 * X[, 2]^2)
  }

  Theta_true <- matrix(NA, 2 * n, K)
  for(k in seq_len(K)){
    rn <- rbinom(n, 1, Phi_true[k])
    Theta_true[, k] <- c(rn, rn)
  }

  Y <- matrix(NA, 2 * n, K)
  for(k in seq_len(K)){
    prob <- Theta_true[, k] * (1 - e_n) + (1 - Theta_true[, k]) * e_p
    Y[, k] <- rbinom(2 * n, 1, prob)
  }

  list(X = X, Y = Y, Theta_true = Theta_true, Phi_true = Phi_true)
}

### HER with user-specified prior variance for gamma^+ and gamma^-
HER_prior <- function(Y, X, ID, Z = NULL, draw = 2000, burn = 1000, prior_var = 100){
  p <- dim(X)[2]
  n <- dim(X)[1] / 2
  K <- dim(Y)[2]
  N <- max(ID)
  MC <- draw - burn

  ID_theta <- t(combn(1:N, 2))
  ID_reference <- rbind(ID_theta, ID_theta[, c(2, 1)])
  if(sum(abs(ID - ID_reference)) > 0){
    stop("ID should be re-ordered")
  }

  qq <- 0
  if(!is.null(Z)){
    qq <- dim(Z)[2]
    if(dim(Z)[1] != (2 * n)){
      stop("Row dimensions of X and Z are not equal.")
    }
    lam_p <- lam_n <- 1
  }

  XZ <- cbind(X, Z)

  gam_p <- c(-2, rep(0, p - 1))
  gam_n <- c(-2, rep(0, p - 1))
  if(qq > 0){
    gam_p <- c(gam_p, rep(0, qq))
    gam_n <- c(gam_n, rep(0, qq))
  }

  Theta <- matrix(NA, n, K)
  for(k in seq_len(K)){
    Theta[, k] <- Y[1:n, k]
  }
  Phi <- apply(Y, 2, mean)

  ## Prior for gamma: gamma^+, gamma^- ~ N(0, prior_var * I).
  ## The spline-coefficient precision is updated through lam_p and lam_n below.
  A0 <- rep(0, p + qq)
  invB0 <- (1 / prior_var) * diag(p + qq)
  al <- be <- 1

  gam_p_pos <- gam_n_pos <- matrix(NA, MC, p + qq)
  Theta_pos <- array(NA, c(MC, n, K))
  Phi_pos <- matrix(NA, MC, K)

  pb <- progress::progress_bar$new(total = draw)
  for(it in seq_len(draw)){
    Theta_ex <- rbind(Theta, Theta)

    ## false negative coefficients
    Om <- matrix(NA, 2 * n, K)
    psi_n <- c(XZ %*% gam_n)
    for(k in seq_len(K)){
      Om[, k] <- Theta_ex[, k] * pgdraw::pgdraw(b = rep(1, 2 * n), c = psi_n)
    }

    invB0_n <- invB0
    if(qq > 0){ diag(invB0_n)[(p + 1):(p + qq)] <- lam_n }
    mat <- invB0_n
    for(k in seq_len(K)){
      mat <- mat + t(XZ) %*% (XZ * Om[, k])
    }
    B_n <- solve(mat)
    A_n <- B_n %*% (invB0_n %*% A0 + K * rowMeans(t(XZ) %*% ((0.5 - Y) * Theta_ex)))
    gam_n <- MCMCpack::mvrnorm(1, A_n, B_n)
    if(qq > 0){ lam_n <- rgamma(1, 1 + 0.5 * qq, 1 + 0.5 * sum(gam_n[-(1:p)]^2)) }

    ## false positive coefficients
    Om <- matrix(NA, 2 * n, K)
    psi_p <- c(XZ %*% gam_p)
    for(k in seq_len(K)){
      Om[, k] <- (1 - Theta_ex[, k]) * pgdraw::pgdraw(b = rep(1, 2 * n), c = psi_p)
    }

    invB0_p <- invB0
    if(qq > 0){ diag(invB0_p)[(p + 1):(p + qq)] <- lam_p }
    mat <- invB0_p
    for(k in seq_len(K)){
      mat <- mat + t(XZ) %*% (XZ * Om[, k])
    }
    B_p <- solve(mat)
    A_p <- B_p %*% (invB0_p %*% A0 + K * rowMeans(t(XZ) %*% ((Y - 0.5) * (1 - Theta_ex))))
    gam_p <- MCMCpack::mvrnorm(1, A_p, B_p)
    if(qq > 0){ lam_p <- rgamma(1, 1 + 0.5 * qq, 1 + 0.5 * sum(gam_p[-(1:p)]^2)) }

    ## latent network
    psi_n <- c(XZ %*% gam_n)
    psi_p <- c(XZ %*% gam_p)
    for(k in seq_len(K)){
      p1_num <- Phi[k] * exp(psi_n[1:n] * (1 - Y[1:n, k]) + psi_n[-(1:n)] * (1 - Y[-(1:n), k]))
      p1_den <- (1 + exp(psi_n[1:n])) * (1 + exp(psi_n[-(1:n)]))
      p0_num <- (1 - Phi[k]) * exp(psi_p[1:n] * Y[1:n, k] + psi_p[-(1:n)] * Y[-(1:n), k])
      p0_den <- (1 + exp(psi_p[1:n])) * (1 + exp(psi_p[-(1:n)]))
      prob <- 1 / (1 + p0_num * p1_den / (p0_den * p1_num))
      Theta[, k] <- rbinom(n, 1, prob)
    }

    ## edge probabilities
    ss <- apply(Theta, 2, sum)
    Phi <- rbeta(K, al + ss, be + (n - ss))

    if(it > burn){
      cc <- it - burn
      gam_p_pos[cc, ] <- gam_p
      gam_n_pos[cc, ] <- gam_n
      Theta_pos[cc, , ] <- Theta
      Phi_pos[cc, ] <- Phi
    }

    pb$tick()
  }

  list(gam_p = gam_p_pos, gam_n = gam_n_pos, Theta = Theta_pos, Phi = Phi_pos)
}

### ID matrix
n <- N * (N - 1) / 2
ID_theta <- t(combn(1:N, 2))
ID <- rbind(ID_theta, ID_theta[, c(2, 1)])

### Storage
AUC <- array(NA_real_, dim = c(R, length(scenarios), length(prior_var_grid), K),
             dimnames = list(rep = 1:R,
                             scenario = paste0("Scenario", scenarios),
                             prior_var = paste0("B0_", prior_var_grid),
                             relation = paste0("K", 1:K)))
F1 <- AUC
Runtime <- array(NA_real_, dim = c(R, length(scenarios), length(prior_var_grid)),
                 dimnames = list(rep = 1:R,
                                 scenario = paste0("Scenario", scenarios),
                                 prior_var = paste0("B0_", prior_var_grid)))

### Monte Carlo loop
for(s_ind in seq_along(scenarios)){
  sc <- scenarios[s_ind]
  message("Starting Scenario ", sc)

  for(r in seq_len(R)){
    dat <- get_truth_and_data(sc = sc, K = K, N = N)
    X <- dat$X
    Y <- dat$Y
    true <- dat$Theta_true[1:n, ]
    Z <- make_Z(X, M = M)

    for(v_ind in seq_along(prior_var_grid)){
      prior_var <- prior_var_grid[v_ind]
      t0 <- proc.time()["elapsed"]
      fit_np <- HER_prior(Y, X, ID, Z = Z, draw = draw, burn = burn, prior_var = prior_var)
      Runtime[r, s_ind, v_ind] <- proc.time()["elapsed"] - t0

      Theta_hat <- apply(fit_np$Theta, c(2, 3), mean)

      for(k in seq_len(K)){
        AUC[r, s_ind, v_ind, k] <- auc_one(score = Theta_hat[, k], truth = true[, k])
        F1[r, s_ind, v_ind, k] <- f1_one(score = Theta_hat[, k], truth = true[, k], cutoff = 0.5)
      }
    }

    message("Scenario ", sc, ", replication ", r, "/", R, " completed")
  }
}

### Summaries averaged over K relations within each replication
AUC_rep <- apply(AUC, c(1, 2, 3), mean, na.rm = TRUE)
F1_rep <- apply(F1, c(1, 2, 3), mean, na.rm = TRUE)

summary_list <- list()
row_id <- 1
for(s_ind in seq_along(scenarios)){
  for(v_ind in seq_along(prior_var_grid)){
    auc_vec <- AUC_rep[, s_ind, v_ind]
    f1_vec <- F1_rep[, s_ind, v_ind]
    rt_vec <- Runtime[, s_ind, v_ind]

    summary_list[[row_id]] <- data.frame(
      Scenario = scenarios[s_ind],
      PriorVariance = prior_var_grid[v_ind],
      AUC_mean = mean(auc_vec, na.rm = TRUE),
      AUC_sd = sd(auc_vec, na.rm = TRUE),
      F1_mean = mean(f1_vec, na.rm = TRUE),
      F1_sd = sd(f1_vec, na.rm = TRUE),
      Runtime_mean_sec = mean(rt_vec, na.rm = TRUE),
      Runtime_sd_sec = sd(rt_vec, na.rm = TRUE)
    )
    row_id <- row_id + 1
  }
}
summary_tab <- do.call(rbind, summary_list)

write.csv(summary_tab, paste0(save_prefix, "_summary.csv"), row.names = FALSE)
save(AUC, F1, Runtime, AUC_rep, F1_rep, summary_tab,
     file = paste0(save_prefix, "_full_results.RData"))

### Boxplot figure for AUC and F1/runtime
plot_df <- data.frame()
for(s_ind in seq_along(scenarios)){
  for(v_ind in seq_along(prior_var_grid)){
    tmp <- data.frame(
      Scenario = paste0("Scenario ", scenarios[s_ind]),
      PriorVariance = paste0("B0=", prior_var_grid[v_ind], "I"),
      AUC = AUC_rep[, s_ind, v_ind],
      F1 = F1_rep[, s_ind, v_ind],
      Runtime_sec = Runtime[, s_ind, v_ind]
    )
    plot_df <- rbind(plot_df, tmp)
  }
}
plot_long <- reshape2::melt(plot_df,
                            id.vars = c("Scenario", "PriorVariance"),
                            measure.vars = c("AUC", "F1", "Runtime_sec"),
                            variable.name = "Metric",
                            value.name = "Value")

p <- ggplot(plot_long, aes(x = PriorVariance, y = Value, fill = PriorVariance)) +
  geom_boxplot(coef = 10) +
  facet_grid(Metric ~ Scenario, scales = "free_y") +
  theme_minimal() +
  labs(x = "Prior variance", y = "Value") +
  guides(fill = "none") +
  theme(strip.text = element_text(size = 12),
        axis.text.x = element_text(size = 10, angle = 20, hjust = 1))

ggsave(paste0(save_prefix, "_boxplots.pdf"), p, width = 10, height = 8)

