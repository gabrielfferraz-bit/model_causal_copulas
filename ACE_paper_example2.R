#################################################################################
# COPULA-BASED CAUSAL INFERENCE: CORRECTED IMPLEMENTATION + MONTE-CARLO
# Based on: "Modelagem causal estrutural via cópulas" (González-López, Justus & Ferraz, 2026)
# - Corrected Gaussian / FGM formulas
# - Added Monte-Carlo (compare Naive OLS, OLS+Z, Gaussian Copula, FGM Copula, BART)
#################################################################################

# ----- Setup ------------------------------------------------------------------
required_packages <- c("MASS", "copula", "ggplot2", "gridExtra", "dplyr", "tidyr",
                       "dbarts", "grid", "cowplot")
new_packages <- required_packages[!(required_packages %in% installed.packages()[,"Package"])]
if(length(new_packages)) install.packages(new_packages, repos="https://cloud.r-project.org")
invisible(lapply(required_packages, library, character.only = TRUE))

set.seed(2024)

#################################################################################
# PART 1: CORRECTED CORE FUNCTIONS (Gaussian & FGM)
#################################################################################

# Estimate K for Gaussian copula (Theorem 5.1)
# NOTE: use Pearson correlations on normal scores (qnorm of ranks) to estimate latent correlations
estimate_K_gaussian <- function(data) {
  n <- nrow(data)
  # Rank transform to uniform [0,1]
  U <- apply(data, 2, function(x) rank(x) / (n + 1))
  # Normal scores for Gaussian copula latent variables
  N <- apply(U, 2, qnorm)
  # Use Pearson correlation of normal scores to estimate latent correlations
  rho <- cor(N, use = "pairwise.complete.obs", method = "pearson")
  # partial correlation formula (latent scale)
  K <- (rho["Y","X"] - rho["Y","Z"] * rho["X","Z"]) /
    sqrt((1 - rho["Y","Z"]^2) * (1 - rho["X","Z"]^2))
  return(list(rho = rho, K = as.numeric(K)))
}

# Eq (5.17): mu(x)=Phi(K * Phi^{-1}(x))
compute_mu_gaussian <- function(x, K) {
  x <- pmin(pmax(x, 1e-6), 1 - 1e-6)
  pnorm(K * qnorm(x))
}

# Eq (5.18): ACE(x) = K * phi(K Phi^{-1}(x)) / phi(Phi^{-1}(x))
compute_ace_gaussian <- function(x, K) {
  x <- pmin(pmax(x, 1e-6), 1 - 1e-6)
  ax <- qnorm(x)
  K * dnorm(K * ax) / dnorm(ax)
}

# FGM parameter estimation (corrected)
# mapping: Kendall's tau (tau) to theta (FGM parameter) via theta = 9/2 * tau
# Spearman rho = theta / 3
estimate_fgm_params_corrected <- function(data) {
  # use Kendall's tau estimate (pairwise)
  tau_mat <- cor(data, method = "kendall")
  # convert to FGM theta parameter
  theta_XY <- pmax(pmin((9/2) * tau_mat["Y","X"], 1), -1)
  theta_XZ <- pmax(pmin((9/2) * tau_mat["X","Z"], 1), -1)
  theta_YZ <- pmax(pmin((9/2) * tau_mat["Y","Z"], 1), -1)
  # Convert to Spearman-like rho for partial-corr formula
  rho_XY <- theta_XY / 3
  rho_XZ <- theta_XZ / 3
  rho_YZ <- theta_YZ / 3
  # Partial correlation on the rho scale
  denom <- sqrt(pmax((1 - rho_YZ^2) * (1 - rho_XZ^2), 1e-6))
  alpha_rho <- (rho_XY - rho_XZ * rho_YZ) / denom
  # Convert back to FGM alpha parameter (alpha_theta = 3 * alpha_rho), clamp to [-1,1]
  alpha_theta <- pmax(pmin(3 * alpha_rho, 1), -1)
  return(list(theta_XY = theta_XY, theta_XZ = theta_XZ, theta_YZ = theta_YZ,
              rho = c(rho_XY = rho_XY, rho_XZ = rho_XZ, rho_YZ = rho_YZ),
              alpha = alpha_theta))
}

# ACE for FGM: from paper (Section 5.2): ACE(x) = alpha * (15 - eta^2) / 45  (constant in x)
# here eta is the Y-Z FGM theta parameter (theta_YZ), default zero if unknown
compute_ace_fgm_corrected <- function(alpha_theta, eta_theta = 0) {
  alpha_theta * (15 - eta_theta^2) / 45
}

#################################################################################
# PART 2: SIMULATION / DGP FUNCTIONS
#################################################################################

# helpers: build data given n, rho_ZX for latent Gaussian copula (uniform margins)
simulate_latent_uniform <- function(n, rho_ZX = 0.7) {
  Sigma <- matrix(c(1, rho_ZX, rho_ZX, 1), nrow = 2)
  latent <- MASS::mvrnorm(n, mu = c(0,0), Sigma = Sigma)
  Z_u <- pnorm(latent[,1])
  X_u <- pnorm(latent[,2])
  data.frame(X = X_u, Z = Z_u)
}

# Linear DGP: Y = beta * X + gamma * Z + noise
make_linear_outcome <- function(X, Z, beta = 2.5, gamma = 3, sigma_eps = 0.5) {
  Y <- beta * X + gamma * Z + rnorm(length(X), 0, sigma_eps)
  list(Y = Y, true_ATE = beta)
}

# Nonlinear DGP: Y = tau(X) * X + gamma * Z^2 + noise
make_nonlinear_outcome <- function(X, Z, gamma = 3, sigma_eps = 0.5) {
  tau_fun <- function(x) 2 + 1.5 * (1 - 2 * x)  # as in original script
  true_effect <- tau_fun(X)
  Y <- true_effect * X + gamma * Z^2 + rnorm(length(X), 0, sigma_eps)
  list(Y = Y, true_ATE = mean(true_effect), tau_fun = tau_fun)
}

#################################################################################
# PART 3: Estimator wrappers (produce a single estimate of average ACE)
#################################################################################

# Utility: compute average ACE over x_eval grid given ACE_x vector
avg_ace_from_grid <- function(ace_x) mean(ace_x, na.rm = TRUE)

# Estimators on one dataset (data.frame with X,Y,Z)
run_estimators_once <- function(df, 
                                x_eval = seq(0.1, 0.9, by = 0.1), 
                                bart_ndpost = 1000) {
  
  # Ensure correct names
  X <- df$X
  Y <- df$Y
  Z <- df$Z
  n <- nrow(df)
  
  ########################################
  # 1. Naive OLS
  ########################################
  ols_naive <- coef(lm(Y ~ X, data = df))["X"]
  
  ########################################
  # 2. OLS + Z
  ########################################
  ols_z <- coef(lm(Y ~ X + Z, data = df))["X"]
  
  ########################################
  # 3. Gaussian Copula (correct copula)
  ########################################
  params_gauss <- estimate_K_gaussian(df[, c("X","Y","Z")])
  ace_gauss_x <- compute_ace_gaussian(x_eval, params_gauss$K)
  est_gauss <- avg_ace_from_grid(ace_gauss_x)
  
  ########################################
  # 4. FGM Copula (misspecified copula)
  ########################################
  params_fgm <- estimate_fgm_params_corrected(df[, c("X","Y","Z")])
  est_fgm <- compute_ace_fgm_corrected(
    alpha_theta = params_fgm$alpha,
    eta_theta   = params_fgm$theta_YZ
  )
  
  ########################################
  # 5. BART (Flexible benchmark)
  ########################################
  
  bart_fit <- tryCatch({
    dbarts::bart(
      x.train  = data.frame(X = X, Z = Z),
      y.train  = Y,
      verbose  = FALSE,
      ndpost   = bart_ndpost,
      nskip    = max(200, bart_ndpost / 2),
      keeptrees = TRUE   # ← REQUIRED FIX
    )
  }, error = function(e) NULL)
  
  if (is.null(bart_fit)) {
    
    est_bart    <- NA
    ace_bart_x  <- rep(NA, length(x_eval))
    
  } else {
    
    ########################################
    # Efficient vectorized prediction
    ########################################
    
    # Build stacked prediction matrix for all x_eval at once
    newdata_all <- do.call(rbind, lapply(x_eval, function(x_val) {
      data.frame(X = rep(x_val, n), Z = Z)
    }))
    
    preds_all <- predict(bart_fit, newdata = newdata_all)
    # preds_all: ndpost × (length(x_eval)*n)
    
    # Average posterior draws first
    preds_mean <- colMeans(preds_all)
    
    # Reshape into matrix: rows = x_eval, cols = units
    preds_matrix <- matrix(preds_mean, 
                           nrow = length(x_eval), 
                           ncol = n, 
                           byrow = TRUE)
    
    # Compute μ(x) = E_Z[E[Y|X=x,Z]]
    mu_bart_x <- rowMeans(preds_matrix)
    
    ########################################
    # Numerical derivative (central difference)
    ########################################
    
    h <- 1e-3
    
    get_mu_fast <- function(x_val) {
      newdata <- data.frame(X = rep(x_val, n), Z = Z)
      preds   <- predict(bart_fit, newdata = newdata)
      mean(colMeans(preds))
    }
    
    ace_bart_x <- sapply(x_eval, function(xx) {
      xp <- min(max(xx + h, 1e-6), 1 - 1e-6)
      xm <- min(max(xx - h, 1e-6), 1 - 1e-6)
      (get_mu_fast(xp) - get_mu_fast(xm)) / (xp - xm)
    })
    
    est_bart <- avg_ace_from_grid(ace_bart_x)
  }
  
  ########################################
  # Return results
  ########################################
  
  return(list(
    ols_naive = as.numeric(ols_naive),
    ols_z     = as.numeric(ols_z),
    gauss     = as.numeric(est_gauss),
    fgm       = as.numeric(est_fgm),
    bart      = as.numeric(est_bart),
    ace_grid  = data.frame(
      x     = x_eval,
      gauss = ace_gauss_x,
      fgm   = rep(est_fgm, length(x_eval)),
      bart  = ace_bart_x,
      true_placeholder = NA
    )
  ))
}

#################################################################################
# PART 4: Monte-Carlo experiment (multiple replicates)
#################################################################################
run_montecarlo <- function(reps = 200, n = 1000, rho_ZX = 0.7, x_eval = seq(0.1, 0.9, 0.1),
                           scenario = c("linear","nonlinear"), parallel = FALSE) {
  scenario <- match.arg(scenario)
  # Storage
  store <- vector("list", reps)
  for(r in seq_len(reps)) {
    # simulate marginals
    uv <- simulate_latent_uniform(n, rho_ZX = rho_ZX)
    if(scenario == "linear") {
      outcome <- make_linear_outcome(uv$X, uv$Z, beta = 2.5, gamma = 3, sigma_eps = 0.5)
      true_ATE <- outcome$true_ATE
      df <- data.frame(X = uv$X, Y = outcome$Y, Z = uv$Z)
      true_ace_x <- rep(true_ATE, length(x_eval))   # constant in X for linear structural coefficient
    } else {
      outcome <- make_nonlinear_outcome(uv$X, uv$Z, gamma = 3, sigma_eps = 0.5)
      true_ATE <- outcome$true_ATE
      df <- data.frame(X = uv$X, Y = outcome$Y, Z = uv$Z)
      true_ace_x <- outcome$tau_fun(x_eval)  # exact true local effect function tau(x)
    }
    # run estimators
    ests <- run_estimators_once(df, x_eval = x_eval, bart_ndpost = 400)
    # attach true ACE(x) for plotting later
    ests$ace_grid$true <- true_ace_x
    store[[r]] <- list(ests = ests, true_ATE = true_ATE)
    if(r %% 10 == 0) cat(sprintf("[%s] rep %d / %d done\n", scenario, r, reps))
  }
  return(store)
}

#################################################################################
# PART 5: Run experiments (two scenarios)  -- adjust reps/n for performance
#################################################################################

# parameters - you can tune these
reps <- 200      # number of Monte-Carlo replicates (set smaller for quick tests)
n <- 1000        # sample size per replicate
rho_ZX <- 0.7
x_eval <- seq(0.1, 0.9, by = 0.1)

cat("Starting Monte-Carlo: linear scenario\n")
mc_linear <- run_montecarlo(reps = reps, n = n, rho_ZX = rho_ZX, x_eval = x_eval, scenario = "linear")
cat("Starting Monte-Carlo: nonlinear scenario\n")
mc_nonlinear <- run_montecarlo(reps = reps, n = n, rho_ZX = rho_ZX, x_eval = x_eval, scenario = "nonlinear")

#################################################################################
# PART 6: Summarize results (Bias, RMSE) across replicates
#################################################################################

summarize_mc <- function(mc_store) {
  reps <- length(mc_store)
  # methods: ols_naive, ols_z, gauss, fgm, bart
  methods <- c("ols_naive","ols_z","gauss","fgm","bart")
  res_list <- lapply(seq_len(reps), function(i) {
    s <- mc_store[[i]]
    ests <- s$ests
    est_values <- c(ests$ols_naive, ests$ols_z, ests$gauss, ests$fgm, ests$bart)
    names(est_values) <- methods
    data.frame(rep = i, method = methods, estimate = as.numeric(est_values),
               true = s$true_ATE, stringsAsFactors = FALSE)
  })
  df_all <- do.call(rbind, res_list)
  summary_stats <- df_all %>%
    group_by(method) %>%
    summarize(mean_est = mean(estimate, na.rm = TRUE),
              sd_est = sd(estimate, na.rm = TRUE),
              true_mean = mean(true),
              bias = mean(estimate - true, na.rm = TRUE),
              rmse = sqrt(mean((estimate - true)^2, na.rm = TRUE)),
              n_na = sum(is.na(estimate)),
              .groups = "drop")
  return(list(df_all = df_all, summary = summary_stats))
}

sum_lin <- summarize_mc(mc_linear)data:image/png;base64,iVBORw0KGgoAAAANSUhEUgAAABIAAAASCAYAAABWzo5XAAAAbElEQVR4Xs2RQQrAMAgEfZgf7W9LAguybljJpR3wEse5JOL3ZObDb4x1loDhHbBOFU6i2Ddnw2KNiXcdAXygJlwE8OFVBHDgKrLgSInN4WMe9iXiqIVsTMjH7z/GhNTEibOxQswcYIWYOR/zAjBJfiXh3jZ6AAAAAElFTkSuQmCC
sum_nl <- summarize_mc(mc_nonlinear)

cat("\n--- LINEAR SCENARIO SUMMARY ---\n")
print(sum_lin$summary, digits = 3)
cat("\n--- NONLINEAR SCENARIO SUMMARY ---\n")
print(sum_nl$summary, digits = 3)

#################################################################################
# PART 7: Example plots (pick the first replicate as representative)
#################################################################################

# helper to plot one replicate ACE curves for comparison
plot_one_rep <- function(store_entry, title = "ACE comparison (one replicate)") {
  df_grid <- store_entry$ests$ace_grid
  df_plot <- tidyr::pivot_longer(df_grid, cols = -x, names_to = "method", values_to = "ACE")
  # map labels
  df_plot$method <- factor(df_plot$method, levels = c("true","gauss","fgm","bart"),
                           labels = c("True ACE", "Gaussian Copula", "FGM Copula", "BART"))
  p <- ggplot(df_plot, aes(x = x, y = ACE, color = method, linetype = method)) +
    geom_line(size = 1) + geom_point(size = 1.8) +
    labs(title = title, x = "x", y = "ACE(x)") +
    theme_minimal() + theme(legend.position = "bottom")
  return(p)
}

p_lin <- plot_one_rep(mc_linear[[1]], "Linear DGP: ACE comparison (replicate 1)")
p_nl  <- plot_one_rep(mc_nonlinear[[1]], "Nonlinear DGP: ACE comparison (replicate 1)")

# display side by side
gridExtra::grid.arrange(p_lin, p_nl, ncol = 2)

#################################################################################
# PART 8: LaTeX-style tables for inclusion in draft
#################################################################################

# Prepare summary tables (linear and nonlinear)
print("LaTeX table: Linear scenario (summary)")
print(sum_lin$summary)

print("LaTeX table: Nonlinear scenario (summary)")
print(sum_nl$summary)

# Example LaTeX table output creation (linear)
cat("\n% LaTeX table: Linear Monte-Carlo summary\n")
cat("\\begin{table}[htbp]\n\\centering\n")
cat("\\caption{Monte-Carlo summary: Linear DGP}\n")
cat("\\begin{tabular}{lrrrrr}\n\\hline\n")
cat("Method & MeanEst & Bias & RMSE & SD & NA \\\\\n\\hline\n")
for(i in 1:nrow(sum_lin$summary)) {
  row <- sum_lin$summary[i,]
  cat(sprintf("%s & %.3f & %.3f & %.3f & %.3f & %d \\\\\n",
              row$method, row$mean_est, row$bias, row$rmse, row$sd_est, row$n_na))
}
cat("\\hline\n\\end{tabular}\n\\end{table}\n")

#################################################################################
# FINAL NOTE: Outputs
#################################################################################
cat("\nDone. Outputs:\n- mc_linear and mc_nonlinear (lists with per-rep results)\n- sum_lin$summary and sum_nl$summary (bias/RMSE tables)\n- p_lin and p_nl plotted above (ACE curves for a representative replicate)\n\n")
cat("Adjust reps/n or bart settings as needed. The script preserves your original structure\nand uses corrected copula formulas and BART as a flexible benchmark.\n")
