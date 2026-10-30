###------------------------------------------------------------###
###  Optimized knot-number sensitivity study for HER-np        ###
###------------------------------------------------------------###

rm(list = ls())
set.seed(12345)

## Required packages
req_pkg <- c("pgdraw", "MCMCpack", "PRROC", "ggplot2", "reshape2")
for(pkg in req_pkg){
  if(!requireNamespace(pkg, quietly = TRUE)) install.packages(pkg)
  library(pkg, character.only = TRUE)
}

R <- 50
K <- 6
N <- 50
scenarios <- c(5, 7, 9)
M_grid <- c(5, 7, 10)
save_prefix <- "KnotSensitivity"

draw <- 2000
burn <- 1000

use_parallel <- FALSE
n_cores <- min(length(M_grid), max(1L, parallel::detectCores() - 1L))

### Helper functions
LGS <- function(x) plogis(x)
Po <- function(x) pmax(x, 0)

make_Z <- function(W_cont, M){
  ## Preallocate rather than repeatedly using cbind().
  q <- ncol(W_cont)
  Z <- matrix(NA_real_, nrow(W_cont), q * (M + 1L))
  znames <- character(ncol(Z))
  pos <- 1L

  probs <- seq(0.1, 0.9, length.out = M)
  for(ell in seq_len(q)){
    v <- W_cont[, ell]
    knots <- as.numeric(quantile(v, probs = probs, na.rm = TRUE, names = FALSE))

    Z[, pos] <- v^2
    znames[pos] <- paste0(colnames(W_cont)[ell], "_sq")
    pos <- pos + 1L

    for(j in seq_len(M)){
      Z[, pos] <- Po(v - knots[j])^2
      znames[pos] <- paste0(colnames(W_cont)[ell], "_knot", j)
      pos <- pos + 1L
    }
  }
  colnames(Z) <- znames
  Z
}

auc_one <- function(score, truth){
  out <- PRROC::roc.curve(
    scores.class0 = score[truth == 1],
    scores.class1 = score[truth == 0],
    curve = FALSE
  )
  as.numeric(out$auc)
}

f1_one <- function(score, truth, cutoff = 0.5){
  pred <- score >= cutoff
  tp <- sum(pred & truth == 1)
  fp <- sum(pred & truth == 0)
  fn <- sum(!pred & truth == 1)
  denom <- 2 * tp + fp + fn
  if(denom == 0) return(NA_real_)
  2 * tp / denom
}

get_scenario <- function(sc){
  gam_p <- gam_n <- numeric(7)
  gam_p[1] <- gam_n[1] <- -1.5
  node_nonlinear <- dyad_nonlinear <- FALSE

  if(sc == 5){
    gam_p[c(2,4)] <- 1
    gam_n[c(2,4)] <- -0.5
    gam_p[c(3,5)] <- 1
    gam_n[c(3,5)] <- -0.5
    node_nonlinear <- TRUE
  }
  if(sc == 7){
    gam_p[c(6,7)] <- 1
    gam_n[c(6,7)] <- -0.5
    dyad_nonlinear <- TRUE
  }
  if(sc == 9){
    gam_p[c(2,4)] <- 1
    gam_n[c(2,4)] <- -0.5
    gam_p[c(3,5)] <- 1
    gam_n[c(3,5)] <- -0.5
    gam_p[c(6,7)] <- 0.5
    gam_n[c(6,7)] <- -0.5
    node_nonlinear <- TRUE
    dyad_nonlinear <- TRUE
  }

  list(gam_p = gam_p, gam_n = gam_n,
       node_nonlinear = node_nonlinear,
       dyad_nonlinear = dyad_nonlinear)
}

get_truth_and_data <- function(sc, K, N, ID){
  n <- N * (N - 1) / 2
  Phi_true <- seq(0.3, 0.8, length.out = K)
  scen <- get_scenario(sc)

  w1 <- rbinom(N, 1, 0.75)
  w2 <- rnorm(N)

  wi1 <- w1[ID[, 1]]
  wi2 <- w2[ID[, 1]]
  wj1 <- w1[ID[, 2]]
  wj2 <- w2[ID[, 2]]
  z1 <- rnorm(2 * n)
  z2 <- rnorm(2 * n)

  X <- cbind(1, wi1, wi2, wj1, wj2, z1, z2)
  colnames(X) <- c("Intercept", "wi1", "wi2", "wj1", "wj2", "z1", "z2")

  eta_p <- drop(X %*% scen$gam_p)
  eta_n <- drop(X %*% scen$gam_n)

  if(scen$node_nonlinear){
    eta_p <- eta_p + 0.5 * sin(wi2) + 0.5 * wj2^2
    eta_n <- eta_n + 0.5 * sin(wj2) + 0.5 * wi2^2
  }
  if(scen$dyad_nonlinear){
    eta_p <- eta_p + 0.5 * sin(z1) + 0.5 * z2^2
    eta_n <- eta_n + 0.5 * sin(z2) + 0.5 * z1^2
  }

  e_p <- plogis(eta_p)
  e_n <- plogis(eta_n)

  Theta_half <- sapply(Phi_true, function(phi) rbinom(n, 1, phi))
  if(K == 1L) Theta_half <- matrix(Theta_half, ncol = 1L)
  Theta_true <- rbind(Theta_half, Theta_half)

  prob <- Theta_true * (1 - e_n) + (1 - Theta_true) * e_p
  Y <- matrix(rbinom(length(prob), 1, prob), nrow = 2 * n, ncol = K)

  W_cont <- cbind(wi2, wj2, z1, z2)
  colnames(W_cont) <- c("wi2", "wj2", "z1", "z2")

  list(X = X, W_cont = W_cont, Y = Y,
       Theta_true = Theta_true, Phi_true = Phi_true)
}


HER_np_fast <- function(Y, X, ID, Z, draw = 2000, burn = 1000){
  p <- ncol(X)
  qq <- ncol(Z)
  n <- nrow(X) / 2
  K <- ncol(Y)
  N <- max(ID)
  MC <- draw - burn

  ID_theta <- t(combn(1:N, 2))
  ID_reference <- rbind(ID_theta, ID_theta[, c(2, 1)])
  if(!identical(dim(ID), dim(ID_reference)) || any(ID != ID_reference)){
    stop("ID should be re-ordered")
  }

  XZ <- cbind(X, Z)
  d <- ncol(XZ)
  A0 <- numeric(d)

  gam_p <- c(-2, rep(0, d - 1L))
  gam_n <- c(-2, rep(0, d - 1L))
  lam_p <- lam_n <- 1

  Theta <- Y[seq_len(n), , drop = FALSE]
  Phi <- colMeans(Y)
  Theta_sum <- matrix(0, n, K)

  ## Objects reused across iterations
  Theta_ex <- matrix(0, 2 * n, K)
  Om <- matrix(0, 2 * n, K)
  ones_pg <- rep(1, 2 * n)

  for(it in seq_len(draw)){
    Theta_ex[seq_len(n), ] <- Theta
    Theta_ex[n + seq_len(n), ] <- Theta

    ## False-negative coefficients
    psi_n <- drop(XZ %*% gam_n)
    pg_n <- pgdraw::pgdraw(b = rep(ones_pg, K), c = rep(psi_n, K))
    Om[] <- matrix(pg_n, nrow = 2 * n, ncol = K) * Theta_ex

    inv_diag_n <- c(rep(0.01, p), rep(lam_n, qq))
    mat_n <- diag(inv_diag_n, d) + crossprod(XZ, XZ * rowSums(Om))
    B_n <- chol2inv(chol(mat_n))
    rhs_n <- crossprod(XZ, rowSums((0.5 - Y) * Theta_ex))
    A_n <- drop(B_n %*% rhs_n)
    gam_n <- mvrnorm(1, A_n, B_n)
    lam_n <- rgamma(1, shape = 1 + 0.5 * qq,
                    rate = 1 + 0.5 * sum(gam_n[(p + 1L):d]^2))

    ## False-positive coefficients
    psi_p <- drop(XZ %*% gam_p)
    pg_p <- pgdraw::pgdraw(b = rep(ones_pg, K), c = rep(psi_p, K))
    Om[] <- matrix(pg_p, nrow = 2 * n, ncol = K) * (1 - Theta_ex)

    inv_diag_p <- c(rep(0.01, p), rep(lam_p, qq))
    mat_p <- diag(inv_diag_p, d) + crossprod(XZ, XZ * rowSums(Om))
    B_p <- chol2inv(chol(mat_p))
    rhs_p <- crossprod(XZ, rowSums((Y - 0.5) * (1 - Theta_ex)))
    A_p <- drop(B_p %*% rhs_p)
    gam_p <- mvrnorm(1, A_p, B_p)
    lam_p <- rgamma(1, shape = 1 + 0.5 * qq,
                    rate = 1 + 0.5 * sum(gam_p[(p + 1L):d]^2))

    ## Latent symmetric networks; stable log-odds form
    psi_n <- drop(XZ %*% gam_n)
    psi_p <- drop(XZ %*% gam_p)
    idx1 <- seq_len(n)
    idx2 <- n + idx1

    for(k in seq_len(K)){
      log_lik1 <- (1 - Y[idx1, k]) * psi_n[idx1] - log1p(exp(psi_n[idx1])) +
                  (1 - Y[idx2, k]) * psi_n[idx2] - log1p(exp(psi_n[idx2]))
      log_lik0 <- Y[idx1, k] * psi_p[idx1] - log1p(exp(psi_p[idx1])) +
                  Y[idx2, k] * psi_p[idx2] - log1p(exp(psi_p[idx2]))
      log_odds <- qlogis(Phi[k]) + log_lik1 - log_lik0
      Theta[, k] <- rbinom(n, 1, plogis(log_odds))
    }

    ss <- colSums(Theta)
    Phi <- rbeta(K, 1 + ss, 1 + n - ss)

    if(it > burn) Theta_sum <- Theta_sum + Theta
  }

  list(Theta_mean = Theta_sum / MC)
}

### ID matrix
n <- N * (N - 1) / 2
ID_theta <- t(combn(1:N, 2))
ID <- rbind(ID_theta, ID_theta[, c(2, 1)])

### Storage
AUC <- array(NA_real_, c(R, length(scenarios), length(M_grid), K),
             dimnames = list(rep = seq_len(R),
                             scenario = paste0("Scenario", scenarios),
                             M = paste0("M", M_grid),
                             relation = paste0("K", seq_len(K))))
F1 <- AUC
Runtime <- array(NA_real_, c(R, length(scenarios), length(M_grid)),
                 dimnames = list(rep = seq_len(R),
                                 scenario = paste0("Scenario", scenarios),
                                 M = paste0("M", M_grid)))

fit_one_M <- function(M, X, W_cont, Y, true){
  Z <- make_Z(W_cont, M)
  t0 <- proc.time()["elapsed"]
  fit <- HER_np_fast(Y, X, ID, Z, draw = draw, burn = burn)
  runtime <- unname(proc.time()["elapsed"] - t0)
  Theta_hat <- fit$Theta_mean

  auc <- f1 <- numeric(K)
  for(k in seq_len(K)){
    auc[k] <- auc_one(Theta_hat[, k], true[, k])
    f1[k] <- f1_one(Theta_hat[, k], true[, k], cutoff = 0.5)
  }
  list(M = M, AUC = auc, F1 = f1, Runtime = runtime)
}

### Monte Carlo loop
for(s_ind in seq_along(scenarios)){
  sc <- scenarios[s_ind]
  message("Starting Scenario ", sc)

  for(r in seq_len(R)){
    dat <- get_truth_and_data(sc, K, N, ID)
    true <- dat$Theta_true[seq_len(n), , drop = FALSE]

    ## Run each knot setting.
    if(use_parallel && .Platform$OS.type != "windows"){
      fits <- parallel::mclapply(
        M_grid,
        function(m) {
          tryCatch(
            fit_one_M(m, X = dat$X, W_cont = dat$W_cont,
                      Y = dat$Y, true = true),
            error = function(e) structure(
              list(M = m, message = conditionMessage(e)),
              class = "fit_error"
            )
          )
        },
        mc.cores = n_cores,
        mc.set.seed = TRUE,
        mc.preschedule = FALSE
      )
    } else {
      fits <- lapply(
        M_grid,
        function(m) {
          tryCatch(
            fit_one_M(m, X = dat$X, W_cont = dat$W_cont,
                      Y = dat$Y, true = true),
            error = function(e) structure(
              list(M = m, message = conditionMessage(e)),
              class = "fit_error"
            )
          )
        }
      )
    }

    for(m_ind in seq_along(fits)){
      if(inherits(fits[[m_ind]], "fit_error")){
        stop(
          "HER-np fit failed for Scenario ", sc,
          ", replication ", r,
          ", M = ", fits[[m_ind]]$M,
          ": ", fits[[m_ind]]$message
        )
      }
      AUC[r, s_ind, m_ind, ] <- fits[[m_ind]][["AUC"]]
      F1[r, s_ind, m_ind, ] <- fits[[m_ind]][["F1"]]
      Runtime[r, s_ind, m_ind] <- fits[[m_ind]][["Runtime"]]
    }

    message("Scenario ", sc, ", replication ", r, "/", R, " completed")
  }
}

### Summaries averaged over K relations within each replication
AUC_rep <- apply(AUC, c(1, 2, 3), mean, na.rm = TRUE)
F1_rep <- apply(F1, c(1, 2, 3), mean, na.rm = TRUE)

summary_list <- list()
row_id <- 1L
for(s_ind in seq_along(scenarios)){
  for(m_ind in seq_along(M_grid)){
    auc_vec <- AUC_rep[, s_ind, m_ind]
    f1_vec <- F1_rep[, s_ind, m_ind]
    rt_vec <- Runtime[, s_ind, m_ind]

    summary_list[[row_id]] <- data.frame(
      Scenario = scenarios[s_ind],
      M = M_grid[m_ind],
      AUC_mean = mean(auc_vec, na.rm = TRUE),
      AUC_sd = sd(auc_vec, na.rm = TRUE),
      F1_mean = mean(f1_vec, na.rm = TRUE),
      F1_sd = sd(f1_vec, na.rm = TRUE),
      Runtime_mean_sec = mean(rt_vec, na.rm = TRUE),
      Runtime_sd_sec = sd(rt_vec, na.rm = TRUE)
    )
    row_id <- row_id + 1L
  }
}
summary_tab <- do.call(rbind, summary_list)

write.csv(summary_tab, paste0(save_prefix, "_summary.csv"), row.names = FALSE)
save(AUC, F1, Runtime, AUC_rep, F1_rep, summary_tab,
     file = paste0(save_prefix, "_full_results.RData"))
