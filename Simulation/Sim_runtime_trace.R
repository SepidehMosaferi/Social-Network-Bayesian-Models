###----------------------------------------------------###
### Runtime and MCMC Trace Diagnostics                 ###
###----------------------------------------------------###

rm(list = ls())
set.seed(2026)

source("HER-function.R")
source("ER-function.R")

library(ggplot2)
library(reshape2)

###----------------------------------------------------###
### User settings
###----------------------------------------------------###

# Scenario 5 is used as a representative nonlinear setting,
# containing both linear and nonlinear node-level effects.
scenario_runtime <- 5
scenario_trace <- 5

M_default <- 7
draw_default <- 5000
burn_default <- 2000

###----------------------------------------------------###
### Helper functions
###----------------------------------------------------###

LGS <- function(x) 1 / (1 + exp(-x))
Po  <- function(xx) xx * (xx > 0)

make_spline_Z <- function(W_cont, M = 7){
  # W_cont contains the four continuous variables:
  # wi2, wj2, z1, and z2.
  # For each variable v, include v^2 and (v-kappa)^2_+ terms.

  Z <- NULL

  for(ell in seq_len(ncol(W_cont))){
    v <- W_cont[, ell]
    knots <- quantile(
      v,
      prob = seq(0.1, 0.9, length.out = M),
      na.rm = TRUE
    )

    Z <- cbind(Z, v^2)

    for(j in seq_len(M)){
      Z <- cbind(Z, Po(v - knots[j])^2)
    }
  }

  colnames(Z) <- paste0("Z", seq_len(ncol(Z)))
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

  if(sc == 1){
    # Constant error probabilities
  }

  if(sc == 2){
    # Discrete node-level effects only
    gam_p[c(2, 4)] <- 2
    gam_n[c(2, 4)] <- -1.5
  }

  if(sc == 3){
    # Continuous node-level effects only
    gam_p[c(3, 5)] <- 2
    gam_n[c(3, 5)] <- -1.5
  }

  if(sc == 4){
    # Linear node-level effects: discrete and continuous
    gam_p[c(2, 4)] <- 2
    gam_n[c(2, 4)] <- -1.5
    gam_p[c(3, 5)] <- 1
    gam_n[c(3, 5)] <- -0.5
  }

  if(sc == 5){
    # Linear and nonlinear node-level effects
    gam_p[c(2, 4)] <- 1
    gam_n[c(2, 4)] <- -0.5
    gam_p[c(3, 5)] <- 1
    gam_n[c(3, 5)] <- -0.5
    node_nonlinear <- TRUE
  }

  if(sc == 6){
    # Linear dyad-level effects only
    gam_p[c(6, 7)] <- 2
    gam_n[c(6, 7)] <- -1.5
  }

  if(sc == 7){
    # Linear and nonlinear dyad-level effects
    gam_p[c(6, 7)] <- 1
    gam_n[c(6, 7)] <- -0.5
    dyad_nonlinear <- TRUE
  }

  if(sc == 8){
    # Linear node-level and dyad-level effects
    gam_p[c(2, 4)] <- 2
    gam_n[c(2, 4)] <- -1.5
    gam_p[c(3, 5)] <- 1
    gam_n[c(3, 5)] <- -0.5
    gam_p[c(6, 7)] <- 2
    gam_n[c(6, 7)] <- -1.5
  }

  if(sc == 9){
    # Linear and nonlinear node-level and dyad-level effects
    gam_p[c(2, 4)] <- 1
    gam_n[c(2, 4)] <- -0.5
    gam_p[c(3, 5)] <- 1
    gam_n[c(3, 5)] <- -0.5
    gam_p[c(6, 7)] <- 0.5
    gam_n[c(6, 7)] <- -0.5
    node_nonlinear <- TRUE
    dyad_nonlinear <- TRUE
  }

  if(!sc %in% 1:9){
    stop("Scenario must be an integer from 1 to 9.")
  }

  return(list(
    gam_p = gam_p,
    gam_n = gam_n,
    node_nonlinear = node_nonlinear,
    dyad_nonlinear = dyad_nonlinear
  ))
}

generate_data <- function(N = 50, K = 6, sc = 5){
  n <- N * (N - 1) / 2

  ID_theta <- t(combn(1:N, 2))
  ID <- rbind(ID_theta, ID_theta[, c(2, 1)])

  scen <- get_scenario(sc)
  gam_p_true <- scen$gam_p
  gam_n_true <- scen$gam_n

  # Node-level covariates
  w1 <- rbinom(N, 1, 0.75)
  w2 <- rnorm(N)

  wi1 <- w1[ID[, 1]]
  wi2 <- w2[ID[, 1]]
  wj1 <- w1[ID[, 2]]
  wj2 <- w2[ID[, 2]]

  # Dyad-level covariates
  z1 <- rnorm(2 * n)
  z2 <- rnorm(2 * n)

  # Parametric design matrix used by HER-p and HER-np
  X <- cbind(1, wi1, wi2, wj1, wj2, z1, z2)
  colnames(X) <- c(
    "Intercept", "wi1", "wi2", "wj1", "wj2", "z1", "z2"
  )

  # Continuous covariates used in the spline component of HER-np
  W_cont <- cbind(wi2, wj2, z1, z2)
  colnames(W_cont) <- c("wi2", "wj2", "z1", "z2")

  eta_p <- c(X %*% gam_p_true)
  eta_n <- c(X %*% gam_n_true)

  if(scen$node_nonlinear){
    eta_p <- eta_p + 0.5 * sin(wi2) + 0.5 * wj2^2
    eta_n <- eta_n + 0.5 * sin(wj2) + 0.5 * wi2^2
  }

  if(scen$dyad_nonlinear){
    eta_p <- eta_p + 0.5 * sin(z1) + 0.5 * z2^2
    eta_n <- eta_n + 0.5 * sin(z2) + 0.5 * z1^2
  }

  e_p <- LGS(eta_p)
  e_n <- LGS(eta_n)

  Phi_true <- seq(0.3, 0.8, length.out = K)

  Theta_true <- matrix(NA, 2 * n, K)
  for(k in seq_len(K)){
    rn <- rbinom(n, 1, Phi_true[k])
    Theta_true[, k] <- c(rn, rn)
  }

  Y <- matrix(NA, 2 * n, K)
  for(k in seq_len(K)){
    prob <- Theta_true[, k] * (1 - e_n) +
      (1 - Theta_true[, k]) * e_p

    Y[, k] <- rbinom(2 * n, 1, prob)
  }

  return(list(
    Y = Y,
    X = X,
    W_cont = W_cont,
    ID = ID,
    Theta_true = Theta_true,
    Phi_true = Phi_true
  ))
}

run_fit_time <- function(
    N,
    K,
    model = c("HER-p", "HER-np"),
    sc = 5,
    M = 7,
    draw = 5000,
    burn = 2000
){
  model <- match.arg(model)

  dat <- generate_data(N = N, K = K, sc = sc)

  if(model == "HER-p"){
    start_time <- Sys.time()

    fit <- HER(
      dat$Y,
      dat$X,
      dat$ID,
      draw = draw,
      burn = burn
    )

    runtime <- as.numeric(
      difftime(Sys.time(), start_time, units = "secs")
    )
  }

  if(model == "HER-np"){
    Z <- make_spline_Z(dat$W_cont, M = M)

    start_time <- Sys.time()

    fit <- HER(
      dat$Y,
      dat$X,
      dat$ID,
      Z = Z,
      draw = draw,
      burn = burn
    )

    runtime <- as.numeric(
      difftime(Sys.time(), start_time, units = "secs")
    )
  }

  rm(fit)
  gc()

  return(runtime)
}

###----------------------------------------------------###
### Part 1: Runtime scaling with network size N
###----------------------------------------------------###

N_values <- c(50, 75, 100)
K_fixed <- 6
models <- c("HER-p", "HER-np")

runtime_N <- data.frame()

for(N in N_values){
  for(mm in models){
    cat(
      "Running Scenario", scenario_runtime,
      ": N =", N,
      ", K =", K_fixed,
      ", model =", mm, "\n"
    )

    tt <- run_fit_time(
      N = N,
      K = K_fixed,
      model = mm,
      sc = scenario_runtime,
      M = M_default,
      draw = draw_default,
      burn = burn_default
    )

    runtime_N <- rbind(
      runtime_N,
      data.frame(
        Scaling = "N",
        Scenario = scenario_runtime,
        N = N,
        K = K_fixed,
        Model = mm,
        Runtime_sec = tt
      )
    )

    write.csv(
      runtime_N,
      "runtime_scaling_N.csv",
      row.names = FALSE
    )
  }
}

###----------------------------------------------------###
### Part 2: Runtime scaling with number of relations K
###----------------------------------------------------###

K_values <- c(3, 6, 12)
N_fixed <- 50

runtime_K <- data.frame()

for(K in K_values){
  for(mm in models){
    cat(
      "Running Scenario", scenario_runtime,
      ": N =", N_fixed,
      ", K =", K,
      ", model =", mm, "\n"
    )

    tt <- run_fit_time(
      N = N_fixed,
      K = K,
      model = mm,
      sc = scenario_runtime,
      M = M_default,
      draw = draw_default,
      burn = burn_default
    )

    runtime_K <- rbind(
      runtime_K,
      data.frame(
        Scaling = "K",
        Scenario = scenario_runtime,
        N = N_fixed,
        K = K,
        Model = mm,
        Runtime_sec = tt
      )
    )

    write.csv(
      runtime_K,
      "runtime_scaling_K.csv",
      row.names = FALSE
    )
  }
}

runtime_all <- rbind(runtime_N, runtime_K)

write.csv(
  runtime_all,
  "runtime_scaling_all.csv",
  row.names = FALSE
)

save(
  runtime_all,
  runtime_N,
  runtime_K,
  file = "Runtime_results.RData"
)

### Runtime plots
pN <- ggplot(
  runtime_N,
  aes(
    x = N,
    y = Runtime_sec,
    group = Model,
    linetype = Model,
    shape = Model
  )
) +
  geom_line(linewidth = 0.7) +
  geom_point(size = 2) +
  theme_minimal() +
  labs(
    title = "Runtime scaling with network size",
    x = "Number of actors (N)",
    y = "Runtime (seconds)"
  ) +
  theme(
    plot.title = element_text(hjust = 0.5),
    legend.title = element_text(face = "bold")
  )

ggsave(
  "runtime_scaling_N.pdf",
  pN,
  width = 6,
  height = 4
)

pK <- ggplot(
  runtime_K,
  aes(
    x = K,
    y = Runtime_sec,
    group = Model,
    linetype = Model,
    shape = Model
  )
) +
  geom_line(linewidth = 0.7) +
  geom_point(size = 2) +
  theme_minimal() +
  labs(
    title = "Runtime scaling with number of relations",
    x = "Number of relations (K)",
    y = "Runtime (seconds)"
  ) +
  theme(
    plot.title = element_text(hjust = 0.5),
    legend.title = element_text(face = "bold")
  )

ggsave(
  "runtime_scaling_K.pdf",
  pK,
  width = 6,
  height = 4
)

###----------------------------------------------------###
### Part 3: Trace plots for representative parameters
###----------------------------------------------------###

cat(
  "Running one HER-np fit for trace plots under Scenario",
  scenario_trace, "...\n"
)

dat_trace <- generate_data(
  N = 50,
  K = 6,
  sc = scenario_trace
)

Z_trace <- make_spline_Z(
  dat_trace$W_cont,
  M = M_default
)

fit_trace <- HER(
  dat_trace$Y,
  dat_trace$X,
  dat_trace$ID,
  Z = Z_trace,
  draw = draw_default,
  burn = burn_default
)

save(
  fit_trace,
  scenario_trace,
  file = paste0(
    "Trace_fit_HERnp_Scenario",
    scenario_trace,
    "_N50_K6.RData"
  )
)

MC <- nrow(fit_trace$gam_p)

trace_df <- data.frame(
  Iteration = seq_len(MC),
  gamma_plus = fit_trace$gam_p[, 2],
  gamma_minus = fit_trace$gam_n[, 2],
  pi = fit_trace$Phi[, 1]
)

write.csv(
  trace_df,
  paste0(
    "trace_samples_HERnp_Scenario",
    scenario_trace,
    "_N50_K6.csv"
  ),
  row.names = FALSE
)

trace_long <- reshape2::melt(
  trace_df,
  id.vars = "Iteration",
  variable.name = "Parameter",
  value.name = "Value"
)

trace_long$Parameter <- factor(
  trace_long$Parameter,
  levels = c("gamma_plus", "gamma_minus", "pi"),
  labels = c("gamma_plus", "gamma_minus", "pi")
)

p_trace <- ggplot(
  trace_long,
  aes(x = Iteration, y = Value)
) +
  geom_line(linewidth = 0.25) +
  facet_wrap(
    ~ Parameter,
    scales = "free_y",
    ncol = 1
  ) +
  theme_minimal() +
  labs(
    title = "Representative MCMC trace plots",
    x = "Post-burn-in iteration",
    y = "Parameter value"
  ) +
  theme(
    plot.title = element_text(
      hjust = 0.5,
      size = 14,
      face = "bold"
    ),
    strip.text = element_text(
      size = 11,
      face = "bold"
    )
  )

ggsave(
  paste0(
    "trace_plots_HERnp_Scenario",
    scenario_trace,
    "_N50_K6.pdf"
  ),
  p_trace,
  width = 8,
  height = 8
)

