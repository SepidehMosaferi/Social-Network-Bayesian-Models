###--------------------------------------------------------------###
### Sensitivity analysis for weakly informative hyperparameters  ###
###--------------------------------------------------------------###

rm(list = ls())
set.seed(123)

library(pgdraw)
library(MCMCpack)
library(MASS)
library(progress)
library(PRROC)

R <- 20                 
scenarios <- c(5, 7, 9) # representative nonlinear node-, dyad-, and combined scenarios
K <- 6
N <- 50
M <- 7                  # default number of knots
ndraw <- 2000
nburn <- 1000


sensitivity_groups <- c("B0", "Beta", "Spline")

### Hyperparameter settings, one factor varied at a time.
### Default: B0 = 100I, pi_k ~ Beta(1,1), spline variance lambda ~ IG(1,1).
settings_B0 <- data.frame(
  group = "B0",
  setting = c("B0=25I", "B0=100I", "B0=400I"),
  B0_var = c(25, 100, 400),
  alpha0 = 1,
  beta0 = 1,
  a_lambda = 1,
  b_lambda = 1
)

settings_Beta <- data.frame(
  group = "Beta",
  setting = c("Beta(0.5,0.5)", "Beta(1,1)", "Beta(2,2)"),
  B0_var = 100,
  alpha0 = c(0.5, 1, 2),
  beta0 = c(0.5, 1, 2),
  a_lambda = 1,
  b_lambda = 1
)

settings_Spline <- data.frame(
  group = "Spline",
  setting = c("IG(0.5,0.5)", "IG(1,1)", "IG(2,2)"),
  B0_var = 100,
  alpha0 = 1,
  beta0 = 1,
  a_lambda = c(0.5, 1, 2),
  b_lambda = c(0.5, 1, 2)
)

settings_all <- rbind(settings_B0, settings_Beta, settings_Spline)
settings_all <- settings_all[settings_all$group %in% sensitivity_groups, ]
row.names(settings_all) <- NULL

### ----------------------- Helper functions ------------------------ ###
LGS <- function(x) 1 / (1 + exp(-x))
Po <- function(x) x * (x > 0)

build_Z <- function(W_cont, M = 7) {
  ## Quadratic truncated-power spline basis for the four continuous
  ## covariates: wi2, wj2, z1, and z2.
  Z <- NULL
  for (ell in seq_len(ncol(W_cont))) {
    v <- W_cont[, ell]
    knots <- as.numeric(quantile(v,
                                 prob = seq(0.1, 0.9, length = M),
                                 na.rm = TRUE))
    Z <- cbind(Z, v^2)
    for (j in seq_len(M)) {
      Z <- cbind(Z, Po(v - knots[j])^2)
    }
  }
  return(Z)
}

compute_metrics <- function(score, truth, threshold = 0.5) {
  ## score and truth are n by K matrices. Metrics are averaged over K layers.
  K <- ncol(truth)
  auc_vec <- f1_vec <- rep(NA_real_, K)
  for (k in 1:K) {
    y_true <- truth[, k]
    sc <- score[, k]
    roc_obj <- PRROC::roc.curve(scores.class0 = sc[y_true == 1],
                                scores.class1 = sc[y_true == 0])
    auc_vec[k] <- roc_obj$auc

    pred <- ifelse(sc >= threshold, 1, 0)
    TP <- sum(pred == 1 & y_true == 1)
    FP <- sum(pred == 1 & y_true == 0)
    FN <- sum(pred == 0 & y_true == 1)
    f1_vec[k] <- ifelse((2 * TP + FP + FN) == 0, NA_real_,
                        2 * TP / (2 * TP + FP + FN))
  }
  c(AUC = mean(auc_vec, na.rm = TRUE), F1 = mean(f1_vec, na.rm = TRUE))
}

get_scenario <- function(sc) {
  int <- -1.5
  gam_p <- rep(0, 7)
  gam_n <- rep(0, 7)
  gam_p[1] <- int
  gam_n[1] <- int
  node_nonlinear <- FALSE
  dyad_nonlinear <- FALSE

  if (sc == 2) {
    ## Discrete node-level effects only
    gam_p[c(2, 4)] <- 2
    gam_n[c(2, 4)] <- -1.5
  }
  if (sc == 3) {
    ## Continuous node-level effects only
    gam_p[c(3,5)] <- 2
    gam_n[c(3,5)] <- -1.5
  }
  if (sc == 4) {
    ## Linear node-level effects: discrete and continuous
    gam_p[c(2,4)] <- 2
    gam_n[c(2,4)] <- -1.5
    gam_p[c(3,5)] <- 1
    gam_n[c(3,5)] <- -0.5
  }
  if (sc == 5) {
    ## Linear and nonlinear node-level effects
    gam_p[c(2,4)] <- 1
    gam_n[c(2,4)] <- -0.5
    gam_p[c(3,5)] <- 1
    gam_n[c(3,5)] <- -0.5
    node_nonlinear <- TRUE
  }
  if (sc == 6) {
    ## Linear dyadic-level effects only
    gam_p[c(6,7)] <- 2
    gam_n[c(6,7)] <- -1.5
  }
  if (sc == 7) {
    ## Linear and nonlinear dyadic-level effects
    gam_p[c(6,7)] <- 1
    gam_n[c(6,7)] <- -0.5
    dyad_nonlinear <- TRUE
  }
  if (sc == 8) {
    ## Linear node-level and dyadic-level effects
    gam_p[c(2,4)] <- 2
    gam_n[c(2,4)] <- -1.5
    gam_p[c(3,5)] <- 1
    gam_n[c(3,5)] <- -0.5
    gam_p[c(6,7)] <- 2
    gam_n[c(6,7)] <- -1.5
  }
  if (sc == 9) {
    ## Linear and nonlinear node-level and dyadic-level effects
    gam_p[c(2,4)] <- 1
    gam_n[c(2,4)] <- -0.5
    gam_p[c(3,5)] <- 1
    gam_n[c(3,5)] <- -0.5
    gam_p[c(6,7)] <- 0.5
    gam_n[c(6,7)] <- -0.5
    node_nonlinear <- TRUE
    dyad_nonlinear <- TRUE
  }

  return(list(gam_p = gam_p,
              gam_n = gam_n,
              node_nonlinear = node_nonlinear,
              dyad_nonlinear = dyad_nonlinear))
}

### ------------------- Configurable HER-np function ---------------- ###
HER_np_sens <- function(Y, X, ID, Z,
                        draw = 2000, burn = 1000,
                        B0_var = 100,
                        alpha0 = 1, beta0 = 1,
                        a_lambda = 1, b_lambda = 1,
                        show_progress = FALSE) {
  p <- ncol(X)
  n <- nrow(X) / 2
  K <- ncol(Y)
  N <- max(ID)
  MC <- draw - burn
  qq <- ncol(Z)

  ID_theta <- t(combn(1:N, 2))
  ID_reference <- rbind(ID_theta, ID_theta[, c(2, 1)])
  if (sum(abs(ID - ID_reference)) > 0) stop("ID should be re-ordered")
  if (nrow(Z) != 2 * n) stop("Row dimensions of X and Z are not equal.")

  XZ <- cbind(X, Z)

  gam_p <- c(-2, rep(0, p - 1), rep(0, qq))
  gam_n <- c(-2, rep(0, p - 1), rep(0, qq))

  Theta <- matrix(NA, n, K)
  for (k in 1:K) Theta[, k] <- Y[1:n, k]
  Phi <- apply(Y, 2, mean)

  ## Regression prior variance for parametric coefficients: B0_var * I.
  ## Spline coefficient prior variance is lambda, with lambda ~ IG(a_lambda,b_lambda).
  A0 <- rep(0, p + qq)
  invB0_base <- diag(0, p + qq)
  diag(invB0_base)[1:p] <- 1 / B0_var

  lambda_p <- lambda_n <- 1

  Theta_pos <- array(NA, c(MC, n, K))
  Phi_pos <- matrix(NA, MC, K)

  if (show_progress) pb <- progress_bar$new(total = draw)
  for (it in 1:draw) {
    Theta_ex <- rbind(Theta, Theta)

    ## False negative regression coefficients
    Om <- matrix(NA, 2 * n, K)
    psi_n <- c(XZ %*% gam_n)
    for (k in 1:K) {
      Om[, k] <- Theta_ex[, k] * pgdraw(b = rep(1, 2 * n), c = psi_n)
    }

    invB0 <- invB0_base
    diag(invB0)[(p + 1):(p + qq)] <- 1 / lambda_n
    mat <- invB0
    for (k in 1:K) mat <- mat + t(XZ) %*% (XZ * Om[, k])
    B_n <- solve(mat)
    A_n <- B_n %*% (invB0 %*% A0 +
                      K * rowMeans(t(XZ) %*% ((0.5 - Y) * Theta_ex)))
    gam_n <- as.numeric(MASS::mvrnorm(1, A_n, B_n))

    ## lambda_n | alpha_n under lambda_n ~ IG(a_lambda,b_lambda)
    shape_n <- a_lambda + qq / 2
    rate_n <- b_lambda + 0.5 * sum(gam_n[-(1:p)]^2)
    lambda_n <- 1 / rgamma(1, shape = shape_n, rate = rate_n)

    ## False positive regression coefficients
    Om <- matrix(NA, 2 * n, K)
    psi_p <- c(XZ %*% gam_p)
    for (k in 1:K) {
      Om[, k] <- (1 - Theta_ex[, k]) *
        pgdraw(b = rep(1, 2 * n), c = psi_p)
    }

    invB0 <- invB0_base
    diag(invB0)[(p + 1):(p + qq)] <- 1 / lambda_p
    mat <- invB0
    for (k in 1:K) mat <- mat + t(XZ) %*% (XZ * Om[, k])
    B_p <- solve(mat)
    A_p <- B_p %*% (invB0 %*% A0 +
                      K * rowMeans(t(XZ) %*% ((Y - 0.5) * (1 - Theta_ex))))
    gam_p <- as.numeric(MASS::mvrnorm(1, A_p, B_p))

    ## lambda_p | alpha_p under lambda_p ~ IG(a_lambda,b_lambda)
    shape_p <- a_lambda + qq / 2
    rate_p <- b_lambda + 0.5 * sum(gam_p[-(1:p)]^2)
    lambda_p <- 1 / rgamma(1, shape = shape_p, rate = rate_p)

    ## Latent network update
    psi_n <- c(XZ %*% gam_n)
    psi_p <- c(XZ %*% gam_p)
    for (k in 1:K) {
      p1_num <- Phi[k] * exp(psi_n[1:n] * (1 - Y[1:n, k]) +
                               psi_n[-(1:n)] * (1 - Y[-(1:n), k]))
      p1_den <- (1 + exp(psi_n[1:n])) * (1 + exp(psi_n[-(1:n)]))
      p0_num <- (1 - Phi[k]) * exp(psi_p[1:n] * Y[1:n, k] +
                                    psi_p[-(1:n)] * Y[-(1:n), k])
      p0_den <- (1 + exp(psi_p[1:n])) * (1 + exp(psi_p[-(1:n)]))
      prob <- 1 / (1 + p0_num * p1_den / (p0_den * p1_num))
      prob[!is.finite(prob)] <- 0.5
      Theta[, k] <- rbinom(n, 1, prob)
    }

    ## Edge probabilities pi_k
    ss <- apply(Theta, 2, sum)
    Phi <- rbeta(K, alpha0 + ss, beta0 + (n - ss))

    if (it > burn) {
      cc <- it - burn
      Theta_pos[cc, , ] <- Theta
      Phi_pos[cc, ] <- Phi
    }
    if (show_progress) pb$tick()
  }

  list(Theta = Theta_pos, Phi = Phi_pos)
}

### ------------------------- Simulation loop ------------------------ ###
Phi_true <- seq(0.3, 0.8, length = K)
n <- N * (N - 1) / 2
ID_theta <- t(combn(1:N, 2))
ID <- rbind(ID_theta, ID_theta[, c(2, 1)])

results <- list()
row_id <- 1

for (sc in scenarios) {
  scen <- get_scenario(sc)
  gam_p_true <- scen$gam_p
  gam_n_true <- scen$gam_n

  for (r in 1:R) {
    cat("Scenario", sc, "replication", r, "of", R, "\n")

    ## Generate node- and dyad-level covariates once per replication.
    w1 <- rbinom(N, 1, 0.75)
    w2 <- rnorm(N)

    wi1 <- w1[ID[, 1]]
    wi2 <- w2[ID[, 1]]
    wj1 <- w1[ID[, 2]]
    wj2 <- w2[ID[, 2]]

    z1 <- rnorm(2 * n)
    z2 <- rnorm(2 * n)

    X <- cbind(1, wi1, wi2, wj1, wj2, z1, z2)
    W_cont <- cbind(wi2, wj2, z1, z2)
    Z <- build_Z(W_cont, M = M)

    ## Generate error probabilities under the selected scenario.
    eta_p <- c(X %*% gam_p_true)
    eta_n <- c(X %*% gam_n_true)

    if (scen$node_nonlinear) {
      eta_p <- eta_p + 0.5 * sin(wi2) + 0.5 * wj2^2
      eta_n <- eta_n + 0.5 * sin(wj2) + 0.5 * wi2^2
    }
    if (scen$dyad_nonlinear) {
      eta_p <- eta_p + 0.5 * sin(z1) + 0.5 * z2^2
      eta_n <- eta_n + 0.5 * sin(z2) + 0.5 * z1^2
    }

    e_p <- LGS(eta_p)
    e_n <- LGS(eta_n)

    Theta_true <- matrix(NA, 2 * n, K)
    for (k in 1:K) {
      rn <- rbinom(n, 1, Phi_true[k])
      Theta_true[, k] <- c(rn, rn)
    }

    Y <- matrix(NA, 2 * n, K)
    for (k in 1:K) {
      prob <- Theta_true[, k] * (1 - e_n) +
        (1 - Theta_true[, k]) * e_p
      Y[, k] <- rbinom(2 * n, 1, prob)
    }
    truth <- Theta_true[1:n, ]

    for (s in seq_len(nrow(settings_all))) {
      st <- settings_all[s, ]
      cat("  ", st$group, st$setting, "\n")

      time_start <- Sys.time()
      fit <- HER_np_sens(Y = Y, X = X, ID = ID, Z = Z,
                         draw = ndraw, burn = nburn,
                         B0_var = st$B0_var,
                         alpha0 = st$alpha0, beta0 = st$beta0,
                         a_lambda = st$a_lambda, b_lambda = st$b_lambda,
                         show_progress = FALSE)
      runtime <- as.numeric(difftime(Sys.time(), time_start, units = "secs"))

      Theta_hat <- apply(fit$Theta, c(2, 3), mean)
      met <- compute_metrics(Theta_hat, truth)

      results[[row_id]] <- data.frame(
        Scenario = sc,
        Replication = r,
        Group = st$group,
        Setting = st$setting,
        B0_var = st$B0_var,
        alpha0 = st$alpha0,
        beta0 = st$beta0,
        a_lambda = st$a_lambda,
        b_lambda = st$b_lambda,
        AUC = met["AUC"],
        F1 = met["F1"],
        Runtime_sec = runtime
      )
      row_id <- row_id + 1
    }
  }
}

raw_results <- do.call(rbind, results)
write.csv(raw_results,
          "hyperparameter_sensitivity_raw_scenarios_5_7_9.csv",
          row.names = FALSE)

summary_results <- aggregate(cbind(AUC, F1, Runtime_sec) ~
                               Scenario + Group + Setting,
                             data = raw_results,
                             FUN = function(x) c(mean = mean(x), sd = sd(x)))

summary_table <- data.frame(
  Scenario = summary_results$Scenario,
  Group = summary_results$Group,
  Setting = summary_results$Setting,
  AUC_mean = summary_results$AUC[, "mean"],
  AUC_sd = summary_results$AUC[, "sd"],
  F1_mean = summary_results$F1[, "mean"],
  F1_sd = summary_results$F1[, "sd"],
  Runtime_mean_sec = summary_results$Runtime_sec[, "mean"],
  Runtime_sd_sec = summary_results$Runtime_sec[, "sd"]
)

write.csv(summary_table,
          "hyperparameter_sensitivity_summary_scenarios_5_7_9.csv",
          row.names = FALSE)

### Compact range table across Scenarios 5, 7, and 9 and all settings
specifications <- c(
  B0 = "$25I$, $100I$, $400I$",
  Beta = "$(0.5,0.5)$, $(1,1)$, $(2,2)$",
  Spline = "$(0.5,0.5)$, $(1,1)$, $(2,2)$"
)

group_labels <- c(
  B0 = "$B_0$",
  Beta = "$(\\alpha_0,\\beta_0)$",
  Spline = "$(a_\\lambda,b_\\lambda)$"
)

range_rows <- lapply(c("B0", "Beta", "Spline"), function(g) {
  tmp <- summary_table[summary_table$Group == g, ]
  data.frame(
    Hyperparameter = group_labels[g],
    Specifications = specifications[g],
    AUC_min = min(tmp$AUC_mean, na.rm = TRUE),
    AUC_max = max(tmp$AUC_mean, na.rm = TRUE),
    F1_min = min(tmp$F1_mean, na.rm = TRUE),
    F1_max = max(tmp$F1_mean, na.rm = TRUE),
    stringsAsFactors = FALSE
  )
})
range_table <- do.call(rbind, range_rows)

range_table$AUC_Range <- sprintf("%.3f--%.3f",
                                 range_table$AUC_min,
                                 range_table$AUC_max)
range_table$F1_Range <- sprintf("%.3f--%.3f",
                                range_table$F1_min,
                                range_table$F1_max)

range_table_output <- range_table[, c("Hyperparameter",
                                      "Specifications",
                                      "AUC_Range",
                                      "F1_Range")]

write.csv(range_table_output,
          "hyperparameter_sensitivity_ranges_scenarios_5_7_9.csv",
          row.names = FALSE)

