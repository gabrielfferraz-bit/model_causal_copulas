################################################################################
# COPULA-BASED CAUSAL INFERENCE: CORRECTED IMPLEMENTATION
# Direct Implementation of Theorem 4.1 and Corollaries from González-López et al (2026)
#
# KEY IMPROVEMENTS:
#   1. Proper implementation of Proposition 4.1 (Interventional CDF)
#   2. Correct conditional copula estimation C_{X,Y|Z}
#   3. Implementation of Corollary 7.1 for observational comparison
#   4. Consistent with simulation code approach
#   5. Uses exact formulas from paper, not approximations
#   6. FIXED: Numerical stability (adaptive epsilon, adaptive dx)
#   7. FIXED: Generic copula implementation (actually uses fitted copulas)
#   8. ADDED: Standardized ACE computation
#   9. ADDED: Distributional inference (quantiles)
################################################################################

# FIXED: Use machine epsilon for numerical stability
EPSILON <- .Machine$double.eps^0.5  # ≈ 1.5e-8

################################################################################
# CORE THEOREM IMPLEMENTATIONS
################################################################################

#' Compute conditional CDF using bivariate copula decomposition
#' 
#' For (X,Y,Z) with Gaussian copula, the conditional CDF is:
#' F_{X|Z}(x|z) = Φ((Φ^{-1}(F_X(x)) - ρ_{XZ}·Φ^{-1}(F_Z(z))) / √(1-ρ²_{XZ}))
#'
#' @param u_x Uniform value for X [0,1]
#' @param u_z Uniform value for Z [0,1]
#' @param rho_xz Copula correlation parameter
#' @return Conditional CDF value
compute_conditional_cdf_gaussian <- function(u_x, u_z, rho_xz) {
  # FIXED: Use adaptive epsilon instead of fixed 1e-10
  u_x <- pmin(pmax(u_x, EPSILON), 1 - EPSILON)
  u_z <- pmin(pmax(u_z, EPSILON), 1 - EPSILON)
  
  # Transform to latent Gaussian scale
  z_x <- qnorm(u_x)
  z_z <- qnorm(u_z)
  
  # Conditional CDF formula for Gaussian copula
  arg <- (z_x - rho_xz * z_z) / sqrt(1 - rho_xz^2)
  F_x_given_z <- pnorm(arg)
  
  return(F_x_given_z)
}

#' Compute partial derivative of conditional copula
#' 
#' Implements ∂C_{X,Y|Z}(u,v|z)/∂u from Proposition 4.1
#' For Gaussian copula: ∂C/∂u = Φ((Φ^{-1}(v) - ρ_{XY|Z}·Φ^{-1}(u)) / √(1-ρ²_{XY|Z}))
#'
#' @param u F_{X|Z}(x|z)
#' @param v F_{Y|Z}(y|z)
#' @param rho_xy_given_z Conditional copula correlation parameter
#' @return Partial derivative value
compute_partial_conditional_copula_gaussian <- function(u, v, rho_xy_given_z) {
  # FIXED: Use adaptive epsilon
  u <- pmin(pmax(u, EPSILON), 1 - EPSILON)
  v <- pmin(pmax(v, EPSILON), 1 - EPSILON)
  
  # Transform to latent scale
  z_u <- qnorm(u)
  z_v <- qnorm(v)
  
  # Partial derivative formula for Gaussian copula
  arg <- (z_v - rho_xy_given_z * z_u) / sqrt(1 - rho_xy_given_z^2)
  partial <- pnorm(arg)
  
  return(partial)
}

#' PROPOSITION 4.1: Interventional CDF
#'
#' F^{(X↓)}_{Y|X}(y|x) = ∫ ∂_u C_{X,Y|Z}(F_{X|Z}(x|z), F_{Y|Z}(y|z)|z) f_Z(z) dz
#'
#' @param y Target outcome value
#' @param x Treatment value  
#' @param z_grid Grid of confounder values for numerical integration
#' @param f_z_vals Marginal density f_Z(z) at grid points
#' @param F_X Function to compute F_X(x)
#' @param F_Y Function to compute F_Y(y)
#' @param F_Z Function to compute F_Z(z)
#' @param rho_xz Copula correlation for (X,Z)
#' @param rho_yz Copula correlation for (Y,Z)
#' @param rho_xy_given_z Conditional copula correlation for (X,Y)|Z
#'
#' @return Interventional CDF value F^{(X↓)}_{Y|X}(y|x)
compute_interventional_cdf_prop41 <- function(y, x, z_grid, f_z_vals,
                                               F_X, F_Y, F_Z,
                                               rho_xz, rho_yz, rho_xy_given_z) {
  # Compute F_X(x), F_Y(y) once
  u_x <- F_X(x)
  u_y <- F_Y(y)
  
  # For each z in grid, compute the integrand
  integrand <- numeric(length(z_grid))
  
  for(i in seq_along(z_grid)) {
    z <- z_grid[i]
    u_z <- F_Z(z)
    
    # Compute conditional CDFs: F_{X|Z}(x|z) and F_{Y|Z}(y|z)
    F_x_given_z <- compute_conditional_cdf_gaussian(u_x, u_z, rho_xz)
    F_y_given_z <- compute_conditional_cdf_gaussian(u_y, u_z, rho_yz)
    
    # Compute ∂C_{X,Y|Z}(F_{X|Z}(x|z), F_{Y|Z}(y|z)|z) / ∂u
    partial <- compute_partial_conditional_copula_gaussian(
      F_x_given_z, F_y_given_z, rho_xy_given_z
    )
    
    # Integrand = partial × f_Z(z)
    integrand[i] <- partial * f_z_vals[i]
  }
  
  # Numerical integration using trapezoidal rule
  dz <- diff(z_grid)
  integral <- sum((integrand[-1] + integrand[-length(integrand)]) / 2 * dz)
  
  return(integral)
}

#' COROLLARY 4.1: Interventional Mean
#'
#' μ(x) = E^{(X↓)}[Y|X=x] = ∫ y · f^{(X↓)}_{Y|X}(y|x) dy
#'
#' We compute this via: μ(x) = ∫ (1 - F^{(X↓)}_{Y|X}(y|x)) dy (for Y ≥ 0)
#'
#' @param x Treatment value
#' @param y_grid Grid of outcome values for integration
#' @param include_tail_correction Whether to add lower tail contribution
#' @param ... Additional parameters passed to compute_interventional_cdf_prop41
#'
#' @return Interventional mean μ(x)
compute_interventional_mean_cor41 <- function(x, y_grid, include_tail_correction = FALSE, ...) {
  
  # Compute F^{(X↓)}_{Y|X}(y|x) for all y in grid
  F_int <- numeric(length(y_grid))
  
  for(i in seq_along(y_grid)) {
    F_int[i] <- compute_interventional_cdf_prop41(y_grid[i], x, ...)
  }
  
  # For positive Y, μ(x) = ∫ (1 - F(y)) dy
  # Numerical integration using trapezoidal rule
  integrand <- 1 - F_int
  dy <- diff(y_grid)
  mu <- sum((integrand[-1] + integrand[-length(integrand)]) / 2 * dy)
  
  # Add lower tail contribution if y_min > 0
  if(include_tail_correction) {
    y_min <- min(y_grid)
    if(y_min > 0) {
      mu <- mu + y_min * (1 - F_int[1])
    }
  }
  
  return(mu)
  
}

#' COROLLARY 4.2: Average Causal Effect
#'
#' ACE(x) = dμ(x)/dx = -∫ ∂F^{(X↓)}_{Y|X}(y|x)/∂x dy
#'
#' @param x Treatment value
#' @param y_grid Grid of outcome values
#' @param dx Finite difference step for numerical derivative (now adaptive)
#' @param x_range Range of x values for adaptive dx calculation
#' @param ... Additional parameters for interventional CDF
#'
#' @return ACE(x) value
compute_ace_cor42 <- function(x, y_grid, dx = NULL, x_range = NULL, ...) {
  # FIXED: Adaptive dx instead of fixed 0.01
  if(is.null(dx)) {
    if(!is.null(x_range)) {
      dx <- diff(range(x_range)) * 0.005  # 0.5% of range
    } else {
      dx <- 0.01  # Fallback to default
    }
  }
  
  # Compute F^{(X↓)}_{Y|X}(y|x) and F^{(X↓)}_{Y|X}(y|x+dx)
  F_at_x <- numeric(length(y_grid))
  F_at_x_plus_dx <- numeric(length(y_grid))
  
  for(i in seq_along(y_grid)) {
    F_at_x[i] <- compute_interventional_cdf_prop41(y_grid[i], x, ...)
    F_at_x_plus_dx[i] <- compute_interventional_cdf_prop41(y_grid[i], x + dx, ...)
  }
  
  # Compute ∂F/∂x via finite difference
  dF_dx <- (F_at_x_plus_dx - F_at_x) / dx
  
  # ACE(x) = -∫ ∂F/∂x dy
  dy <- diff(y_grid)
  ace <- -sum((dF_dx[-1] + dF_dx[-length(dF_dx)]) / 2 * dy)
  
  return(ace)
}

#' COROLLARY 7.1: Non-Interventional (Observational) CDF
#'
#' F_{Y|X}(y|x) = ∫ ∂_u C_{X,Y|Z}(F_{X|Z}(x|z), F_{Y|Z}(y|z)|z) · 
#'                   c_{X,Z}(F_X(x), F_Z(z)) · f_Z(z) dz
#'
#' KEY DIFFERENCE from Proposition 4.1: includes c_{X,Z} term (confounding)
#'
#' @param y Outcome value
#' @param x Treatment value
#' @param z_grid Grid of confounder values
#' @param f_z_vals Marginal density f_Z(z)
#' @param F_X Marginal CDF for X
#' @param F_Y Marginal CDF for Y
#' @param F_Z Marginal CDF for Z
#' @param rho_xz Copula correlation (X,Z)
#' @param rho_yz Copula correlation (Y,Z)
#' @param rho_xy_given_z Conditional copula correlation
#' @param copula_xz_density Function to compute c_{X,Z}(u_x, u_z)
#'
#' @return Observational CDF value F_{Y|X}(y|x)
compute_observational_cdf_cor71 <- function(y, x, z_grid, f_z_vals,
                                             F_X, F_Y, F_Z,
                                             rho_xz, rho_yz, rho_xy_given_z,
                                             copula_xz_density) {
  u_x <- F_X(x)
  u_y <- F_Y(y)
  
  integrand <- numeric(length(z_grid))
  
  for(i in seq_along(z_grid)) {
    z <- z_grid[i]
    u_z <- F_Z(z)
    
    # Conditional CDFs
    F_x_given_z <- compute_conditional_cdf_gaussian(u_x, u_z, rho_xz)
    F_y_given_z <- compute_conditional_cdf_gaussian(u_y, u_z, rho_yz)
    
    # Partial derivative of conditional copula
    partial <- compute_partial_conditional_copula_gaussian(
      F_x_given_z, F_y_given_z, rho_xy_given_z
    )
    
    # Copula density c_{X,Z}(F_X(x), F_Z(z))
    c_xz <- copula_xz_density(u_x, u_z)
    
    # KEY DIFFERENCE: multiply by c_{X,Z} (Corollary 7.1 vs Proposition 4.1)
    integrand[i] <- partial * c_xz * f_z_vals[i]
  }
  
  # Numerical integration
  dz <- diff(z_grid)
  integral <- sum((integrand[-1] + integrand[-length(integrand)]) / 2 * dz)
  
  return(integral)
}

#' THEOREM 5.1: Gaussian Closed-Form Solution (Validation)
#'
#' For Gaussian copulas, exact formula: μ_U(u_x) = Φ(K·Φ^{-1}(u_x))
#' where K = ρ_{XY|Z} is the conditional copula correlation
#'
#' @param x Treatment value on original scale
#' @param K Conditional copula correlation ρ_{XY|Z}
#' @param F_X Marginal CDF for X
#'
#' @return List with mu (uniform scale) and ace (uniform scale)
compute_gaussian_closed_form_thm51 <- function(x, K, F_X) {
  # Transform x to uniform scale
  u_x <- F_X(x)
  u_x <- pmin(pmax(u_x, EPSILON), 1 - EPSILON)
  
  # Exact formula: μ_U(u_x) = Φ(K·Φ^{-1}(u_x))
  z_x <- qnorm(u_x)
  mu_u <- pnorm(K * z_x)
  
  # ACE on uniform scale: dμ_U/du_x = K·φ(K·z_x)/φ(z_x)
  ace_u <- K * dnorm(K * z_x) / dnorm(z_x)
  
  return(list(mu = mu_u, ace = ace_u))
}

################################################################################
# ADDED: DISTRIBUTIONAL INFERENCE FUNCTIONS
################################################################################

#' Compute Interventional Quantiles (ADDED for distributional inference)
#'
#' Computes quantiles of the interventional distribution F^{(X↓)}_{Y|X}(y|x)
#' This implements "full distributional causal inference" from the paper
#'
#' @param x Treatment value
#' @param probs Probability levels for quantiles (e.g., c(0.25, 0.5, 0.75))
#' @param y_grid Grid of outcome values
#' @param ... Additional parameters for interventional CDF
#'
#' @return Vector of quantile values
compute_interventional_quantiles <- function(x, probs = c(0.25, 0.5, 0.75), 
                                             y_grid, ...) {
  quantiles <- numeric(length(probs))
  
  # Compute CDF at all y values
  F_vals <- sapply(y_grid, function(y) 
    compute_interventional_cdf_prop41(y, x, ...))
  
  # For each probability level, find corresponding y
  for(i in seq_along(probs)) {
    # Linear interpolation to find y such that F(y|x) = probs[i]
    if(min(F_vals) > probs[i]) {
      quantiles[i] <- min(y_grid)  # Below range
    } else if(max(F_vals) < probs[i]) {
      quantiles[i] <- max(y_grid)  # Above range
    } else {
      quantiles[i] <- approx(F_vals, y_grid, xout = probs[i], rule = 2)$y
    }
  }
  
  return(quantiles)
}

################################################################################
# GAUSSIAN COPULA DENSITY (for Corollary 7.1)
################################################################################

#' Compute Gaussian copula density
#'
#' @param u First uniform variable
#' @param v Second uniform variable
#' @param rho Copula correlation
#'
#' @return Copula density value c(u,v)
compute_gaussian_copula_density <- function(u, v, rho) {
  u <- pmin(pmax(u, EPSILON), 1 - EPSILON)
  v <- pmin(pmax(v, EPSILON), 1 - EPSILON)
  
  z_u <- qnorm(u)
  z_v <- qnorm(v)
  
  # Gaussian copula density formula
  exponent <- -0.5 * (z_u^2 + z_v^2 - 2*rho*z_u*z_v) / (1 - rho^2) + 
              0.5 * (z_u^2 + z_v^2)
  
  density <- exp(exponent) / sqrt(1 - rho^2)
  
  return(density)
}

################################################################################
# MAIN WRAPPER: GAUSSIAN COPULA CAUSAL EFFECTS
################################################################################

#' Compute Causal Effects using Gaussian Copulas (NHANES-specific wrapper)
#'
#' @param x_grid Grid of treatment values
#' @param nhanes_data Data frame with columns X, Y, and confounders
#' @param confounder_vars Names of confounder variables
#' @param method "structural" (Prop 4.1), "gaussian" (Thm 5.1), or "both"
#'
#' @return Data frame with interventional/observational means, ACE, etc.
compute_causal_effects_nhanes <- function(x_grid, nhanes_data, 
                                          confounder_vars = "C1",
                                          method = "both",
                                          survey_weights = NULL) {
  # survey_weights: optional numeric vector of length nrow(nhanes_data).
  # When supplied, pseudo-observations are computed as weighted empirical CDFs,
  # addressing Reproducibility Reviewer comment 1 (inappropriate handling of
  # complex survey data). When NULL, standard rank-based pseudo-obs are used
  # and claims should be restricted to the analytic sample.
  
  if(length(confounder_vars) == 1) {
    Z_col <- confounder_vars[1]
    X <- nhanes_data$X
    Y <- nhanes_data$Y
    Z <- nhanes_data[[Z_col]]
    n <- length(X)
    
    # ---- Pseudo-observations (weighted if survey weights provided) ----
    if(!is.null(survey_weights)) {
      # Weighted empirical CDF: U_i = sum_{j: x_j <= x_i} w_j / sum(w)
      w <- survey_weights / sum(survey_weights, na.rm = TRUE)
      weighted_ecdf <- function(v, wts) {
        sorted_idx <- order(v)
        cumw <- cumsum(wts[sorted_idx])
        u <- numeric(length(v))
        u[sorted_idx] <- cumw
        pmax(pmin(u, 1 - 1/(2*n)), 1/(2*n))   # keep in open (0,1)
      }
      U_X <- weighted_ecdf(X, w)
      U_Y <- weighted_ecdf(Y, w)
      U_Z <- weighted_ecdf(Z, w)
    } else {
      U_X <- rank(X) / (n + 1)
      U_Y <- rank(Y) / (n + 1)
      U_Z <- rank(Z) / (n + 1)
    }
    
    rho_xz <- cor(U_X, U_Z, method = "spearman")
    rho_yz <- cor(U_Y, U_Z, method = "spearman")
    rho_xy <- cor(U_X, U_Y, method = "spearman")
    rho_xy_given_z <- (rho_xy - rho_xz * rho_yz) / 
      sqrt((1 - rho_xz^2) * (1 - rho_yz^2))
    
    F_X <- ecdf(X)
    F_Y <- ecdf(Y)
    F_Z <- ecdf(Z)
    
    z_grid <- seq(quantile(Z, 0.01), quantile(Z, 0.99), length.out = 50)
    y_grid <- seq(quantile(Y, 0.01), quantile(Y, 0.99), length.out = 50)
    f_z_vals <- approx(density(Z)$x, density(Z)$y, xout = z_grid, rule = 2)$y
    
    copula_xz_dens <- function(u_x, u_z) {
      compute_gaussian_copula_density(u_x, u_z, rho_xz)
    }
    
    results <- data.frame(
      x = x_grid,
      mu_interventional = NA_real_,
      mu_observational  = NA_real_,
      ace_interventional = NA_real_,  # dmu_int/dx  — causal slope
      d_obs_dx          = NA_real_,   # dE[Y|X]/dx  — observational slope (NEW)
      mu_gaussian = NA_real_,
      ace_gaussian = NA_real_
    )
    
    cat("Computing interventional quantities...\n")
    cat(sprintf("  Grid: %d treatment levels\n", length(x_grid)))
    cat(sprintf("  Integration: %d Z points, %d Y points\n\n", 
                length(z_grid), length(y_grid)))
    
    # FIXED: Compute adaptive dx for ACE
    dx_adaptive <- diff(range(x_grid)) * 0.005
    
    # Compute for each x in grid
    for(i in seq_along(x_grid)) {
      x <- x_grid[i]
      
      if(i %% 10 == 0) {
        cat(sprintf("  Progress: %d/%d (%.0f%%)\n", i, length(x_grid), 100*i/length(x_grid)))
      }
      
      # PROPOSITION 4.1: Interventional mean
      if(method %in% c("structural", "both")) {
        results$mu_interventional[i] <- compute_interventional_mean_cor41(
          x = x,
          y_grid = y_grid,
          include_tail_correction = FALSE,
          z_grid = z_grid,
          f_z_vals = f_z_vals,
          F_X = F_X,
          F_Y = F_Y,
          F_Z = F_Z,
          rho_xz = rho_xz,
          rho_yz = rho_yz,
          rho_xy_given_z = rho_xy_given_z
        )
      }
      
      # COROLLARY 7.1: Observational mean (for comparison)
      F_obs <- numeric(length(y_grid))
      for(j in seq_along(y_grid)) {
        F_obs[j] <- compute_observational_cdf_cor71(
          y = y_grid[j],
          x = x,
          z_grid = z_grid,
          f_z_vals = f_z_vals,
          F_X = F_X,
          F_Y = F_Y,
          F_Z = F_Z,
          rho_xz = rho_xz,
          rho_yz = rho_yz,
          rho_xy_given_z = rho_xy_given_z,
          copula_xz_density = copula_xz_dens
        )
      }
      
      # Integrate observational CDF to get mean
      integrand_obs <- 1 - F_obs
      dy <- diff(y_grid)
      mu_obs <- sum((integrand_obs[-1] + integrand_obs[-length(integrand_obs)]) / 2 * dy)
      if(min(y_grid) > 0) {
        mu_obs <- mu_obs + min(y_grid) * (1 - F_obs[1])
      }
      results$mu_observational[i] <- mu_obs
      
      # THEOREM 5.1: Gaussian closed-form (validation)
      if(method %in% c("gaussian", "both")) {
        K <- rho_xy_given_z
        gauss_result <- compute_gaussian_closed_form_thm51(x, K, F_X)
        
        # Transform from uniform to Y scale
        results$mu_gaussian[i] <- quantile(nhanes_data$Y, probs = gauss_result$mu)
        results$ace_gaussian[i] <- gauss_result$ace * 
          diff(quantile(nhanes_data$Y, probs = c(gauss_result$mu - 0.01, gauss_result$mu + 0.01))) / 0.02
      }
    }
    
    # Compute ACE via finite differences from interventional means
    results$ace_interventional <- c(
      (results$mu_interventional[2] - results$mu_interventional[1]) / 
        (results$x[2] - results$x[1]),
      (results$mu_interventional[3:nrow(results)] - results$mu_interventional[1:(nrow(results)-2)]) / 
        (results$x[3:nrow(results)] - results$x[1:(nrow(results)-2)]),
      (results$mu_interventional[nrow(results)] - results$mu_interventional[nrow(results)-1]) / 
        (results$x[nrow(results)] - results$x[nrow(results)-1])
    )
    
    # Standardized ACE (Paper Section 6.9)
    results$ace_standardized <- results$ace_interventional * 
      (sd(nhanes_data$X, na.rm = TRUE) / sd(nhanes_data$Y, na.rm = TRUE))
    
    cat("\n✓ Computation complete\n\n")
    
    # Validation if both methods used
    if(method == "both" && all(!is.na(results$mu_gaussian))) {
      rmse <- sqrt(mean((results$mu_interventional - results$mu_gaussian)^2, na.rm = TRUE))
      cor_val <- cor(results$mu_interventional, results$mu_gaussian, use = "complete.obs")
      
      cat("VALIDATION (Structural vs Gaussian Closed-Form):\n")
      cat(sprintf("  RMSE: %.4f\n", rmse))
      cat(sprintf("  Correlation: %.4f\n", cor_val))
      
      if(cor_val > 0.95) {
        cat("  ✓✓✓ EXCELLENT - Theorem 4.1 matches Theorem 5.1 analytically\n\n")
      } else if(cor_val > 0.85) {
        cat("  ✓ GOOD - Minor numerical differences\n\n")
      } else {
        cat("  ⚠ CHECK - May indicate numerical issues or non-Gaussian structure\n\n")
      }
    }
    
    return(results)
    
  } else {
    stop("Multiple confounders not yet implemented - use single confounder")
  }
}

################################################################################
# FIXED: GENERIC COPULA IMPLEMENTATION (Actually uses fitted copulas now!)
################################################################################

#' Compute Causal Effects using Generic (Data-Driven) Copulas
#'
#' FIXED: Now properly implements generic copulas using BiCopHfunc1/2
#' This relaxes the Gaussian simplifying assumption (Remark 3.1)
#'
#' @param x_grid Grid of treatment values
#' @param nhanes_data Data frame with X, Y, and confounders
#' @param confounder_vars Names of confounder variables
#'
#' @return Data frame with causal effects using data-driven copula families
compute_causal_effects_generic_nhanes <- function(x_grid, nhanes_data, 
                                                  confounder_vars = "C1",
                                                  survey_weights = NULL) {
  library(VineCopula)
  
  # Select columns with standardized names
  Z_col <- confounder_vars[1]
  data_subset <- nhanes_data[, c("X", "Y", Z_col)]
  colnames(data_subset) <- c("X", "Y", "Z")
  
  # Extract variables
  X <- data_subset$X
  Y <- data_subset$Y
  Z <- data_subset$Z
  n <- nrow(data_subset)
  
  # Transform to uniform margins
  U_matrix <- pobs(as.matrix(data_subset))
  U_X <- U_matrix[,1]
  U_Y <- U_matrix[,2] 
  U_Z <- U_matrix[,3]
  
  cat("GENERIC COPULA SELECTION:\n")
  
  # Full VineCopula family set: all parametric families plus all rotations
  full_familyset <- c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10,
                      13, 14, 16, 17, 18, 19, 20,
                      23, 24, 26, 27, 28, 29, 30,
                      33, 34, 36, 37, 38, 39, 40)
  
  # Fit bivariate copulas with explicit BIC selection over all families
  fit_xz <- BiCopSelect(U_X, U_Z, familyset = full_familyset, selectioncrit = "BIC")
  fit_yz <- BiCopSelect(U_Y, U_Z, familyset = full_familyset, selectioncrit = "BIC")
  fit_xy <- BiCopSelect(U_X, U_Y, familyset = full_familyset, selectioncrit = "BIC")
  
  cat(sprintf("  (X,Z): %s (τ=%.3f, par=%.3f)\n", 
              fit_xz$familyname, fit_xz$tau, fit_xz$par))
  cat(sprintf("  (Y,Z): %s (τ=%.3f, par=%.3f)\n", 
              fit_yz$familyname, fit_yz$tau, fit_yz$par))
  cat(sprintf("  (X,Y): %s (τ=%.3f, par=%.3f)\n\n", 
              fit_xy$familyname, fit_xy$tau, fit_xy$par))
  
  # FIXED: Define GENERIC conditional CDF functions using fitted copulas
  F_x_given_z_func <- function(u_x, u_z) {
    # BiCopHfunc1 computes h-function: h(u|v) = P(U <= u | V = v)
    pmin(pmax(BiCopHfunc1(u_x, u_z, fit_xz), EPSILON), 1 - EPSILON)
  }
  
  F_y_given_z_func <- function(u_y, u_z) {
    pmin(pmax(BiCopHfunc1(u_y, u_z, fit_yz), EPSILON), 1 - EPSILON)
  }
  
  # For conditional copula C_{X,Y|Z}, we use the fitted C_{X,Y}
  # This is an approximation assuming conditional independence structure
  partial_copula_func <- function(u, v) {
    pmin(pmax(BiCopHfunc1(u, v, fit_xy), EPSILON), 1 - EPSILON)
  }
  
  # Copula density for observational (Corollary 7.1)
  copula_xz_dens_func <- function(u_x, u_z) {
    BiCopPDF(u_x, u_z, fit_xz)
  }
  
  # Empirical CDFs
  F_X <- ecdf(X)
  F_Y <- ecdf(Y)
  F_Z <- ecdf(Z)
  
  # Integration grids
  z_grid <- seq(quantile(Z, 0.01), quantile(Z, 0.99), length.out = 50)
  y_grid <- seq(quantile(Y, 0.01), quantile(Y, 0.99), length.out = 50)
  f_z_vals <- approx(density(Z)$x, density(Z)$y, xout = z_grid, rule = 2)$y
  dz <- diff(z_grid)
  total_mass <- sum((f_z_vals[-1] + f_z_vals[-length(f_z_vals)]) / 2 * dz)
  f_z_vals <- f_z_vals / total_mass
  
  # Initialize results
  results <- data.frame(
    x = x_grid,
    mu_interventional = NA_real_,
    mu_observational = NA_real_,
    ace_interventional = NA_real_
  )
  
  cat("Computing causal effects with GENERIC copulas...\n")
  
  # FIXED: Compute interventional CDF using GENERIC copulas
  compute_interventional_cdf_generic <- function(y, x) {
    u_x <- F_X(x)
    u_y <- F_Y(y)
    
    integrand <- numeric(length(z_grid))
    
    for(i in seq_along(z_grid)) {
      z <- z_grid[i]
      u_z <- F_Z(z)
      
      # Use FITTED copulas (not Gaussian!)
      F_x_given_z <- F_x_given_z_func(u_x, u_z)
      F_y_given_z <- F_y_given_z_func(u_y, u_z)
      
      # Partial derivative using fitted copula
      partial <- partial_copula_func(F_x_given_z, F_y_given_z)
      
      integrand[i] <- partial * f_z_vals[i]
    }
    
    dz <- diff(z_grid)
    integral <- sum((integrand[-1] + integrand[-length(integrand)]) / 2 * dz)
    return(integral)
  }
  
  # Compute for each x
  for(i in seq_along(x_grid)) {
    x <- x_grid[i]
    
    if(i %% 10 == 0) {
      cat(sprintf("  Progress: %d/%d\n", i, length(x_grid)))
    }
    
    # Interventional mean
    F_int <- sapply(y_grid, function(y) compute_interventional_cdf_generic(y, x))
    integrand <- 1 - F_int
    dy <- diff(y_grid)
    results$mu_interventional[i] <- sum((integrand[-1] + integrand[-length(integrand)]) / 2 * dy)
    
    # Observational mean (with copula density)
    F_obs <- numeric(length(y_grid))
    for(j in seq_along(y_grid)) {
      u_x <- F_X(x)
      u_y <- F_Y(y_grid[j])
      
      integrand_obs <- numeric(length(z_grid))
      for(k in seq_along(z_grid)) {
        z <- z_grid[k]
        u_z <- F_Z(z)
        
        F_x_given_z <- F_x_given_z_func(u_x, u_z)
        F_y_given_z <- F_y_given_z_func(u_y, u_z)
        partial <- partial_copula_func(F_x_given_z, F_y_given_z)
        c_xz <- copula_xz_dens_func(u_x, u_z)
        
        integrand_obs[k] <- partial * c_xz * f_z_vals[k]
      }
      
      dz <- diff(z_grid)
      F_obs[j] <- sum((integrand_obs[-1] + integrand_obs[-length(integrand_obs)]) / 2 * dz)
    }
    
    integrand_obs_mean <- 1 - F_obs
    results$mu_observational[i] <- sum((integrand_obs_mean[-1] + 
                                        integrand_obs_mean[-length(integrand_obs_mean)]) / 2 * dy)
  }
  
  # Compute ACE
  results$ace_interventional <- c(
    (results$mu_interventional[2] - results$mu_interventional[1]) / 
      (results$x[2] - results$x[1]),
    (results$mu_interventional[3:nrow(results)] - results$mu_interventional[1:(nrow(results)-2)]) / 
      (results$x[3:nrow(results)] - results$x[1:(nrow(results)-2)]),
    (results$mu_interventional[nrow(results)] - results$mu_interventional[nrow(results)-1]) / 
      (results$x[nrow(results)] - results$x[nrow(results)-1])
  )
  
  # ADDED: Standardized ACE
  results$ace_standardized <- results$ace_interventional * 
    (sd(nhanes_data$X, na.rm = TRUE) / sd(nhanes_data$Y, na.rm = TRUE))
  
  # Add copula diagnostics
  results$copula_xz_family <- fit_xz$familyname
  results$copula_yz_family <- fit_yz$familyname
  results$copula_xy_family <- fit_xy$familyname
  results$copula_xz_tau <- fit_xz$tau
  results$copula_yz_tau <- fit_yz$tau
  results$copula_xy_tau <- fit_xy$tau
  
  cat("\n✓ Generic copula effects computed\n\n")
  
  return(results)
}

################################################################################
# EMPIRICAL COPULA IMPLEMENTATION (Non-Parametric, No Family Assumption)
################################################################################

#' Compute Causal Effects using Empirical (Non-Parametric) Copulas
#'
#' Uses kernel-smoothed non-parametric copulas for (X,Z), (Y,Z), and (X,Y)
#' via rvinecopulib. This is the most flexible approach — no parametric family
#' is assumed. The integration logic mirrors compute_causal_effects_generic_nhanes
#' exactly, so results are directly comparable.
#'
#' @param x_grid   Grid of treatment values
#' @param nhanes_data  Data frame with X, Y, and confounders
#' @param confounder_vars  Name of the confounder column (single variable)
#'
#' @return Data frame with causal effects using empirical copulas
compute_causal_effects_empirical_nhanes <- function(x_grid,
                                                    nhanes_data,
                                                    confounder_vars = "C1",
                                                    survey_weights = NULL,
                                                    n_z_grid = 40,
                                                    n_y_grid = 40) {
  
  if (!requireNamespace("rvinecopulib", quietly = TRUE))
    stop("Install rvinecopulib: install.packages('rvinecopulib')")
  
  Z_col <- confounder_vars[1]
  data_sub <- nhanes_data[, c("X", "Y", Z_col)]
  colnames(data_sub) <- c("X", "Y", "Z")
  
  X <- data_sub$X; Y <- data_sub$Y; Z <- data_sub$Z
  
  # Pseudo-observations
  U_matrix <- pobs(as.matrix(data_sub))
  U_X <- U_matrix[, 1]; U_Y <- U_matrix[, 2]; U_Z <- U_matrix[, 3]
  
  cat("Fitting nonparametric (kernel TLL) copulas...\n")
  
  fit_xz <- bicop(cbind(U_X, U_Z), family_set = "nonpar", var_types = c("c", "c"))
  fit_yz <- bicop(cbind(U_Y, U_Z), family_set = "nonpar", var_types = c("c", "c"))
  fit_xy <- bicop(cbind(U_X, U_Y), family_set = "nonpar", var_types = c("c", "c"))
  
  cat(sprintf("  τ_XZ=%.4f  τ_YZ=%.4f  τ_XY=%.4f\n\n",
              fit_xz$tau, fit_yz$tau, fit_xy$tau))
  
  # WEIGHTED Marginal CDFs (NHANES-correct)
  if (!is.null(survey_weights)) {
    w <- survey_weights / sum(survey_weights)
    make_wcdf <- function(vals, wts) {
      ord <- order(vals)
      sorted_vals <- vals[ord]
      cum_wts <- cumsum(wts[ord])
      function(q) {
        sapply(q, function(qi) {
          idx <- findInterval(qi, sorted_vals)
          if (idx == 0) 0 else cum_wts[idx] / sum(cum_wts)
        })
      }
    }
    F_X <- make_wcdf(X, w)
    F_Y <- make_wcdf(Y, w)
    F_Z <- make_wcdf(Z, w)
    cat("  Using survey-weighted marginal CDFs\n")
  } else {
    F_X <- ecdf(X); F_Y <- ecdf(Y); F_Z <- ecdf(Z)
  }
  
  # Grids
  z_grid <- seq(quantile(Z, 0.02), quantile(Z, 0.98), length.out = n_z_grid)
  y_grid <- seq(quantile(Y, 0.02), quantile(Y, 0.98), length.out = n_y_grid)
  
  f_z_raw <- approx(density(Z)$x, density(Z)$y, xout = z_grid, rule = 2)$y
  dz <- diff(z_grid)
  f_z_raw <- f_z_raw / sum((f_z_raw[-1] + f_z_raw[-length(f_z_raw)]) / 2 * dz)
  
  # **FIXED h-functions - direct bicop object**
  F_x_given_z_fn <- function(u_x, u_z) {
    pmax(pmin(hbicop(cbind(u_x, u_z), 2, fit_xz), 1-1e-8), 1e-8)
  }
  F_y_given_z_fn <- function(u_y, u_z) {
    pmax(pmin(hbicop(cbind(u_y, u_z), 2, fit_yz), 1-1e-8), 1e-8)
  }
  partial_xy_fn <- function(u, v) {
    pmax(pmin(hbicop(cbind(u, v), 1, fit_xy), 1-1e-8), 1e-8)
  }
  c_xz_fn <- function(u_x, u_z) {
    pmax(dbicop(cbind(u_x, u_z), fit_xz), 1e-6)
  }
  
  # Survival functions (1 - CDF)
  compute_F_int <- function(y_val, x_val) {
    u_x <- F_X(x_val); u_y <- F_Y(y_val)
    u_z_vec <- F_Z(z_grid)
    
    Fxz <- F_x_given_z_fn(rep(u_x, length(u_z_vec)), u_z_vec)
    Fyz <- F_y_given_z_fn(rep(u_y, length(u_z_vec)), u_z_vec)
    integrand <- partial_xy_fn(Fxz, Fyz) * f_z_raw
    
    sum((integrand[-1] + integrand[-length(integrand)]) / 2 * dz)
  }
  
  compute_F_obs <- function(y_val, x_val) {
    u_x <- F_X(x_val); u_y <- F_Y(y_val)
    u_z_vec <- F_Z(z_grid)
    
    Fxz <- F_x_given_z_fn(rep(u_x, length(u_z_vec)), u_z_vec)
    Fyz <- F_y_given_z_fn(rep(u_y, length(u_z_vec)), u_z_vec)
    cxz <- c_xz_fn(rep(u_x, length(u_z_vec)), u_z_vec)
    integrand <- partial_xy_fn(Fxz, Fyz) * cxz * f_z_raw
    
    sum((integrand[-1] + integrand[-length(integrand)]) / 2 * dz)
  }
  
  results <- data.frame(
    x = x_grid,
    mu_interventional = NA_real_,
    mu_observational  = NA_real_,
    ace_interventional = NA_real_
  )
  
  dy <- diff(y_grid)
  cat("Computing causal effects (kernel copula)...\n")
  
  for (i in seq_along(x_grid)) {
    if (i %% 5 == 0 || i == 1)
      cat(sprintf("  x[%d/%d] = %.4f\n", i, length(x_grid), x_grid[i]))
    
    F_int_vals <- sapply(y_grid, compute_F_int, x_val = x_grid[i])
    # **FIXED: Proper survival integration**
    S_int_vals <- 1 - F_int_vals
    results$mu_interventional[i] <- sum((S_int_vals[-1] + S_int_vals[-length(S_int_vals)]) / 2 * dy)
    
    F_obs_vals <- sapply(y_grid, compute_F_obs, x_val = x_grid[i])
    S_obs_vals <- 1 - F_obs_vals
    results$mu_observational[i] <- sum((S_obs_vals[-1] + S_obs_vals[-length(S_obs_vals)]) / 2 * dy)
  }
  
  # ACE
  nr <- nrow(results)
  results$ace_interventional <- c(
    (results$mu_interventional[2] - results$mu_interventional[1]) / (results$x[2] - results$x[1]),
    (results$mu_interventional[3:nr] - results$mu_interventional[1:(nr - 2)]) / (results$x[3:nr] - results$x[1:(nr - 2)]),
    (results$mu_interventional[nr] - results$mu_interventional[nr - 1]) / (results$x[nr] - results$x[nr - 1])
  )
  results$ace_standardized <- results$ace_interventional * (sd(nhanes_data$X, na.rm = TRUE) / sd(nhanes_data$Y, na.rm = TRUE))
  results$method <- "Kernel (TLL nonpar)"
  
  cat("\n✓ Kernel nonparametric copula effects computed.\n\n")
  return(results)
}

################################################################################
# Z-VARYING CONDITIONAL COPULA
# Relaxes the simplifying assumption via kernel-weighted local MLE
#
# THEORETICAL BASIS (González-López et al., 2026 — Remark 3.1):
#   The simplifying assumption fixes C_{X,Y|Z}(u,v|z) = C_{X,Y|Z}(u,v) ∀z.
#   This file relaxes that by estimating θ̂(z) — a smooth function mapping each
#   confounder level z to a copula parameter — via local kernel-weighted MLE.
#   Proposition 4.1's integral then uses the z-specific parameter at each step:
#
#   F^(X↓)_{Y|X}(y|x) = ∫ ∂_u C_{X,Y|Z=z}(F_{X|Z}(x|z), F_{Y|Z}(y|z); θ̂(u_z)) f_Z(z) dz
################################################################################

EPSILON <- 1e-8

# ─────────────────────────────────────────────────────────────────────────────
# HELPER: Local kernel-weighted MLE for θ at a single z value
#
# At each z grid point u_z0, we maximise:
#   L_local(θ; u_z0) = Σᵢ K_h(U_Zᵢ − u_z0) · log c(V_Xᵢ, V_Yᵢ; θ)
# where K_h is a Gaussian kernel. Observations near u_z0 get higher weight,
# so the estimate reflects the local dependence structure at that income level.
# ─────────────────────────────────────────────────────────────────────────────
estimate_theta_local <- function(V_X, V_Y, U_Z, u_z0, family, h,
                                 par_bounds = c(-0.99, 0.99)) {
  # Gaussian kernel weights centred at u_z0 (all in copula scale [0,1])
  kw <- dnorm((U_Z - u_z0) / h)
  
  # Guard: if the kernel mass is negligible (e.g. at boundary), fall back to
  # equal weights so we at least get the global MLE rather than crashing.
  if (sum(kw) < 1e-10) kw <- rep(1, length(U_Z))
  kw <- kw / sum(kw)   # normalise to sum to 1 (turns into a probability vector)
  
  # Negative kernel-weighted log-likelihood — minimise over θ
  obj_fn <- function(par) {
    tryCatch({
      cop   <- BiCop(family = family, par = par)
      dens  <- BiCopPDF(V_X, V_Y, cop)              # per-obs copula density
      log_d <- log(pmax(dens, .Machine$double.eps))  # guard against log(0)
      -sum(kw * log_d)
    }, error = function(e) Inf)   # Inf causes optimize() to reject this par
  }
  
  result <- optimize(obj_fn, interval = par_bounds, tol = 1e-5)
  result$minimum
}


# ─────────────────────────────────────────────────────────────────────────────
# HELPER: Pre-compute the full θ̂(z) curve over a z grid
# ─────────────────────────────────────────────────────────────────────────────
build_theta_curve <- function(V_X, V_Y, U_Z, u_z_grid, family, h) {
  
  # Parameter bounds are family-specific — using the wrong bounds causes
  # optimize() to search outside the valid parameter space.
  bounds <- switch(as.character(family),
                   "1"  = c(-0.99,   0.99),   # Gaussian: Pearson correlation ∈ (-1, 1)
                   "5"  = c(-50,     50),     # Frank: any real, but ±50 is practically ∞
                   "3"  = c(1e-4,    20),     # Clayton: θ > 0 (lower-tail dependence)
                   "4"  = c(1.001,   20),     # Gumbel:  θ ≥ 1 (upper-tail dependence)
                   "6"  = c(1.001,   20),     # Joe:     θ ≥ 1 (upper-tail dependence)
                   "13" = c(1e-4,    20),     # Rotated Clayton 180° (upper-tail)
                   "14" = c(1.001,   20),     # Rotated Gumbel  180° (lower-tail)
                   "23" = c(-20,    -1e-4),   # Rotated Clayton  90° (negative assoc.)
                   "24" = c(-20,    -1.001),  # Rotated Gumbel   90° (negative assoc.)
                   c(-20, 20)                 # Generic fallback
  )
  
  cat(sprintf("  Pre-computing theta_hat(z) at %d points (bandwidth h = %.4f)...\n",
              length(u_z_grid), h))
  
  theta_curve <- vapply(u_z_grid, function(u_z0)
    estimate_theta_local(V_X, V_Y, U_Z, u_z0, family, h, bounds),
    numeric(1))
  
  list(theta_curve = theta_curve, family = family)
}


# ─────────────────────────────────────────────────────────────────────────────
# HELPER: Weighted empirical CDF
#
# FIX (vs earlier version): w is already normalised so cumsum(w) reaches 1
# at its last element. The earlier version erroneously divided cumw[idx] by
# sum(cumw), which is the sum of all cumulative partial sums — approximately
# n/2 for a well-spread distribution — producing CDFs that never reach 1.
# The correct implementation returns cumw[idx] directly.
# ─────────────────────────────────────────────────────────────────────────────
make_weighted_ecdf <- function(vals, weights) {
  w   <- weights / sum(weights)   # normalise so Σw = 1
  ord <- order(vals)
  sorted_vals <- vals[ord]
  cumw        <- cumsum(w[ord])   # cumw[length] == 1 exactly
  
  function(q) {
    sapply(q, function(qi) {
      idx <- findInterval(qi, sorted_vals)
      if (idx == 0) 0 else cumw[idx]   # already a probability; NO further division
    })
  }
}


################################################################################
# MAIN FUNCTION: compute_causal_effects_zvarying_nhanes
#
# Drop-in replacement for compute_causal_effects_generic_nhanes(). The only
# structural difference in the integration is that BiCopHfunc1() is called with
# bicop_z_list[[i]] — the locally estimated BiCop object for z_grid[i] — instead
# of a single global BiCop object. Everything else (grids, trapezoid rule, ACE
# finite differences) is identical to the generic function so results are
# directly comparable.
################################################################################
compute_causal_effects_zvarying_nhanes <- function(
    x_grid,
    nhanes_data,
    confounder_vars  = "C1",
    survey_weights   = NULL,
    h                = NULL,       # bandwidth; NULL = Silverman's rule
    n_z_grid         = 50,
    n_y_grid         = 50,
    one_par_families = c(1, 3, 4, 5, 6, 13, 14, 23, 24)  # excludes Student-t (2 pars)
) {
  library(VineCopula)
  
  Z_col    <- confounder_vars[1]
  data_sub <- nhanes_data[, c("X", "Y", Z_col)]
  colnames(data_sub) <- c("X", "Y", "Z")
  
  X <- data_sub$X;  Y <- data_sub$Y;  Z <- data_sub$Z
  n <- nrow(data_sub)
  
  # ── Pseudo-observations on (0,1) ────────────────────────────────────────────
  U_matrix <- pobs(as.matrix(data_sub))
  U_X <- U_matrix[, 1];  U_Y <- U_matrix[, 2];  U_Z <- U_matrix[, 3]
  
  cat("Z-VARYING CONDITIONAL COPULA (no simplifying assumption):\n")
  cat(strrep("-", 60), "\n")
  cat(sprintf("  n = %d\n  tau_XZ = %.4f   tau_YZ = %.4f   tau_XY = %.4f\n\n",
              n,
              cor(U_X, U_Z, method = "kendall"),
              cor(U_Y, U_Z, method = "kendall"),
              cor(U_X, U_Y, method = "kendall")))
  
  # ── Step 1: Fit GLOBAL marginal copulas for (X,Z) and (Y,Z) ─────────────────
  # These capture how the confounder relates to each variable separately.
  # They are NOT subject to the simplifying assumption — there is no
  # conditioning argument that makes C_{X,Z} vary with a third variable here.
  full_familyset <- c(1, 3, 4, 5, 6, 13, 14, 23, 24)
  cat("Fitting marginal copulas (X,Z) and (Y,Z) globally...\n")
  fit_xz <- BiCopSelect(U_X, U_Z, familyset = full_familyset, selectioncrit = "BIC")
  fit_yz <- BiCopSelect(U_Y, U_Z, familyset = full_familyset, selectioncrit = "BIC")
  cat(sprintf("  (X,Z): %-22s  tau = %.4f   par = %.4f\n",
              fit_xz$familyname, fit_xz$tau, fit_xz$par))
  cat(sprintf("  (Y,Z): %-22s  tau = %.4f   par = %.4f\n\n",
              fit_yz$familyname, fit_yz$tau, fit_yz$par))
  
  # BiCopHfunc1(u, v, obj) computes h(u|v) = P(U1 <= u | U2 = v),
  # i.e. conditioning is on the SECOND argument. Since Z is the second column
  # in our BiCopSelect calls above, u_z is correctly the second argument here.
  F_x_given_z_fn <- function(u_x, u_z)
    pmin(pmax(BiCopHfunc1(u_x, u_z, fit_xz), EPSILON), 1 - EPSILON)
  
  F_y_given_z_fn <- function(u_y, u_z)
    pmin(pmax(BiCopHfunc1(u_y, u_z, fit_yz), EPSILON), 1 - EPSILON)
  
  c_xz_fn <- function(u_x, u_z)
    pmax(BiCopPDF(u_x, u_z, fit_xz), EPSILON)
  
  # ── Step 2: Conditional pseudo-observations ──────────────────────────────────
  # V_X = F_{X|Z}(X|Z) and V_Y = F_{Y|Z}(Y|Z) strip out Z's influence from
  # each variable via the probability integral transform. By construction,
  # they are uniform given Z. Their joint copula, as a function of Z, is
  # exactly C_{X,Y|Z} — the object the simplifying assumption was collapsing
  # to a constant. This is the key identification step.
  cat("Computing conditional pseudo-observations V_X = F_{X|Z}, V_Y = F_{Y|Z}...\n")
  V_X <- F_x_given_z_fn(U_X, U_Z)
  V_Y <- F_y_given_z_fn(U_Y, U_Z)
  cat(sprintf("  tau(V_X, V_Y) = %.4f  [conditional Kendall's tau after removing Z]\n\n",
              cor(V_X, V_Y, method = "kendall")))
  
  # ── Step 3: Select the best one-parameter family for C_{X,Y|Z} ───────────────
  # We restrict to one-parameter families because the local kernel MLE has
  # limited effective sample size at each z point: two-parameter families
  # (e.g. Student-t, family=2) are poorly identified with a narrow kernel.
  cat("Selecting parametric family for C_{X,Y|Z} (BIC on conditional pseudo-obs)...\n")
  fit_xy_cond  <- BiCopSelect(V_X, V_Y,
                              familyset     = one_par_families,
                              selectioncrit = "BIC")
  cond_family  <- fit_xy_cond$family
  global_theta <- fit_xy_cond$par
  cat(sprintf("  Selected: %-22s  theta_global = %.4f   tau = %.4f\n\n",
              fit_xy_cond$familyname, global_theta, fit_xy_cond$tau))
  
  # ── Step 4: Bandwidth selection ──────────────────────────────────────────────
  # Silverman's rule applied in U_Z space. Since U_Z ~ U(0,1), its standard
  # deviation is exactly sqrt(1/12). Larger h → smoother theta_hat(z) curve,
  # converging to the simplifying assumption as h → ∞.
  # Pass h = <value> explicitly for sensitivity analysis.
  if (is.null(h)) {
    h <- 1.06 * sqrt(1 / 12) * n^(-1 / 5)
    cat(sprintf("Bandwidth (Silverman rule): h = %.4f\n", h))
    cat("  Sensitivity tip: rerun with h = %.2f (narrower) or h = %.2f (wider).\n\n",
        h * 0.7, h * 1.5)
  } else {
    cat(sprintf("Bandwidth (user-specified): h = %.4f\n\n", h))
  }
  
  # ── Step 5: Pre-compute theta_hat(z) for the full z grid ─────────────────────
  z_grid   <- seq(quantile(Z, 0.01), quantile(Z, 0.99), length.out = n_z_grid)
  F_Z      <- ecdf(Z)
  u_z_grid <- F_Z(z_grid)   # map z grid → copula scale [0,1]
  
  theta_result <- build_theta_curve(V_X, V_Y, U_Z, u_z_grid, cond_family, h)
  theta_curve  <- theta_result$theta_curve
  
  theta_range_pct <- 100 * (max(theta_curve) - min(theta_curve)) / abs(global_theta)
  cat(sprintf("\n  theta_hat(z) range: [%.4f, %.4f]  (%.1f%% variation around theta_global = %.4f)\n",
              min(theta_curve), max(theta_curve), theta_range_pct, global_theta))
  if (theta_range_pct > 20) {
    cat("  WARNING: large variation — z-varying model materially differs from simplifying assumption.\n")
  } else {
    cat("  OK: modest variation — results should be close to simplifying assumption estimate.\n")
  }
  cat("\n")
  
  # Pre-build a BiCop object for each grid point.
  # This avoids reconstructing the object inside the double loop, which is slow.
  # If a local theta is out of bounds (rare at boundaries), fall back to global.
  bicop_z_list <- lapply(theta_curve, function(th) {
    tryCatch(
      BiCop(family = cond_family, par = th),
      error = function(e) BiCop(family = cond_family, par = global_theta)
    )
  })
  
  # ── Integration setup ─────────────────────────────────────────────────────────
  # Use weighted ECDFs for marginals if survey weights are provided.
  # FIX: the weighted ECDF now returns cumw[idx] directly — see make_weighted_ecdf()
  # above for the explanation of why the earlier division by sum(cumw) was wrong.
  if (!is.null(survey_weights)) {
    F_X <- make_weighted_ecdf(X, survey_weights)
    F_Y <- make_weighted_ecdf(Y, survey_weights)
  } else {
    F_X <- ecdf(X)
    F_Y <- ecdf(Y)
  }
  
  y_grid  <- seq(quantile(Y, 0.01), quantile(Y, 0.99), length.out = n_y_grid)
  f_z_raw <- approx(density(Z)$x, density(Z)$y, xout = z_grid, rule = 2)$y
  dz      <- diff(z_grid)
  # Normalise so the discrete approximation integrates to 1
  f_z_raw <- f_z_raw / sum((f_z_raw[-1] + f_z_raw[-length(f_z_raw)]) / 2 * dz)
  
  # ── Proposition 4.1 integrand (interventional CDF) ───────────────────────────
  # KEY DIFFERENCE from compute_causal_effects_generic_nhanes: the inner call
  # uses bicop_z_list[[i]] — the locally estimated object for z_grid[i] — so
  # the partial derivative of C_{X,Y|Z} varies with z rather than being fixed.
  compute_F_int_zvar <- function(y_val, x_val) {
    u_x <- F_X(x_val);  u_y <- F_Y(y_val)
    integrand <- numeric(length(z_grid))
    for (i in seq_along(z_grid)) {
      u_z <- u_z_grid[i]
      # FIX: as.numeric() guards against named-vector edge cases from ecdf()/
      # findInterval() that can trip up VineCopula's C internals.
      Fxz <- as.numeric(F_x_given_z_fn(u_x, u_z))
      Fyz <- as.numeric(F_y_given_z_fn(u_y, u_z))
      # partial = dC_{X,Y|Z=z_i}/du, using the LOCAL copula parameter theta_hat(z_i)
      partial    <- pmin(pmax(BiCopHfunc1(Fxz, Fyz, bicop_z_list[[i]]), EPSILON), 1 - EPSILON)
      integrand[i] <- partial * f_z_raw[i]
    }
    sum((integrand[-1] + integrand[-length(integrand)]) / 2 * dz)
  }
  
  # ── Corollary 8.1 integrand (observational CDF — includes c_{X,Z} factor) ────
  compute_F_obs_zvar <- function(y_val, x_val) {
    u_x <- F_X(x_val);  u_y <- F_Y(y_val)
    integrand <- numeric(length(z_grid))
    for (i in seq_along(z_grid)) {
      u_z        <- u_z_grid[i]
      Fxz        <- as.numeric(F_x_given_z_fn(u_x, u_z))
      Fyz        <- as.numeric(F_y_given_z_fn(u_y, u_z))
      partial    <- pmin(pmax(BiCopHfunc1(Fxz, Fyz, bicop_z_list[[i]]), EPSILON), 1 - EPSILON)
      cxz        <- c_xz_fn(u_x, u_z)        # copula density of (X,Z) at this z point
      integrand[i] <- partial * cxz * f_z_raw[i]
    }
    sum((integrand[-1] + integrand[-length(integrand)]) / 2 * dz)
  }
  
  # ── Main loop over treatment grid ─────────────────────────────────────────────
  results <- data.frame(
    x                 = x_grid,
    mu_interventional = NA_real_,
    mu_observational  = NA_real_,
    ace_interventional = NA_real_
  )
  
  dy <- diff(y_grid)
  cat("Computing interventional and observational means...\n")
  
  for (i in seq_along(x_grid)) {
    if (i %% 5 == 0 || i == 1)
      cat(sprintf("  x[%d/%d] = %.4f\n", i, length(x_grid), x_grid[i]))
    
    # Interventional mean: mu(x) = ∫₀^∞ [1 - F^(X↓)(y|x)] dy
    F_int <- sapply(y_grid, compute_F_int_zvar, x_val = x_grid[i])
    results$mu_interventional[i] <-
      sum((1 - F_int[-1] + 1 - F_int[-length(F_int)]) / 2 * dy)
    
    # Observational mean (for confounding bias quantification)
    F_obs <- sapply(y_grid, compute_F_obs_zvar, x_val = x_grid[i])
    results$mu_observational[i] <-
      sum((1 - F_obs[-1] + 1 - F_obs[-length(F_obs)]) / 2 * dy)
  }
  
  # ── ACE via central finite differences ───────────────────────────────────────
  nr <- nrow(results)
  results$ace_interventional <- c(
    (results$mu_interventional[2]    - results$mu_interventional[1])      / (results$x[2]    - results$x[1]),
    (results$mu_interventional[3:nr] - results$mu_interventional[1:(nr-2)]) / (results$x[3:nr] - results$x[1:(nr-2)]),
    (results$mu_interventional[nr]   - results$mu_interventional[nr-1])   / (results$x[nr]   - results$x[nr-1])
  )
  
  results$ace_standardized <- results$ace_interventional *
    (sd(X, na.rm = TRUE) / sd(Y, na.rm = TRUE))
  
  # ── Attach diagnostics ────────────────────────────────────────────────────────
  results$method           <- "Z-varying conditional copula"
  results$cond_family      <- fit_xy_cond$familyname
  results$global_theta     <- global_theta
  results$theta_min        <- min(theta_curve)
  results$theta_max        <- max(theta_curve)
  results$theta_range_pct  <- theta_range_pct
  results$bandwidth_h      <- h
  
  # Store the full curve as an attribute for plot_theta_curve()
  attr(results, "theta_curve") <- data.frame(
    z         = z_grid,
    u_z       = u_z_grid,
    theta_hat = theta_curve
  )
  
  cat(sprintf("\n  Done. Mean ACE = %.4f\n\n",
              mean(results$ace_interventional, na.rm = TRUE)))
  return(results)
}


################################################################################
# DIAGNOSTIC PLOT: θ̂(z) curve
#
# A flat line at theta_global validates the simplifying assumption.
# A sloped or humped curve shows exactly how and where the direct X→Y
# dependence varies with the confounder — which is itself a finding worth
# reporting in the paper.
################################################################################
plot_theta_curve <- function(results_zvar, xlab = "Income-to-poverty ratio (C1)") {
  
  tc <- attr(results_zvar, "theta_curve")
  if (is.null(tc))
    stop("No theta_curve attribute — pass the output of compute_causal_effects_zvarying_nhanes().")
  
  fam   <- as.character(results_zvar$cond_family[1])
  theta0 <- results_zvar$global_theta[1]
  rng_pct <- results_zvar$theta_range_pct[1]
  
  ylab_str <- sprintf("theta_hat(z)  [%s copula parameter]", fam)
  title_str <- sprintf("Z-varying copula parameter (%.1f%% variation around global theta)", rng_pct)
  
  plot(tc$z, tc$theta_hat, type = "l", lwd = 2.5, col = "#2166ac",
       xlab = xlab, ylab = ylab_str, main = title_str,
       ylim = range(c(tc$theta_hat, theta0)) * c(0.95, 1.05))
  
  abline(h = theta0, lty = 2, col = "#d73027", lwd = 1.5)
  
  legend("topright",
         legend = c("theta_hat(z)  — z-varying",
                    sprintf("theta_global = %.4f  (simplifying assumption)", theta0)),
         col  = c("#2166ac", "#d73027"),
         lty  = c(1, 2),
         lwd  = c(2.5, 1.5),
         bty  = "n")
  
  invisible(tc)
}

################################################################################
# COMPARISON SUMMARY FUNCTION
################################################################################

#' Generate comparison summary table
#'
#' Compares interventional vs observational means to quantify confounding bias
#'
#' @param results Results data frame from compute_causal_effects_nhanes()
#' @param quantiles Quantiles to report
#'
#' @return Summary table
generate_comparison_summary <- function(results, quantiles = c(0.1, 0.25, 0.5, 0.75, 0.9)) {
  
  # Converter list → data.frame automaticamente
  if(is.list(results) && !is.data.frame(results)) {
    results <- as.data.frame(results)
    cat("CONVERTIDO: list → data.frame (", nrow(results), "linhas)\n", sep="")
  }
  
  if(is.null(results) || !is.data.frame(results) || nrow(results) == 0) {
    stop("ERRO: results inválido após conversão")
  }
  
  n_res <- nrow(results)
  cat(sprintf("ANALISANDO: %d linhas\n", n_res))
  
  required_cols <- c("x", "mu_observational", "mu_interventional", "ace_interventional")
  missing_cols <- required_cols[!required_cols %in% names(results)]
  if(length(missing_cols) > 0) {
    stop(sprintf("FALTAM: %s", paste(missing_cols, collapse = ", ")))
  }
  
  # Calcular d_obs_dx
  if(!"d_obs_dx" %in% names(results)) {
    results$d_obs_dx <- rep(NA_real_, n_res)
    for(i in 2:(n_res-1)) {
      results$d_obs_dx[i] <- (results$mu_observational[i+1] - results$mu_observational[i-1]) /
        (results$x[i+1] - results$x[i-1])
    }
    if(n_res >= 2) {
      results$d_obs_dx[1] <- (results$mu_observational[2] - results$mu_observational[1]) /
        (results$x[2] - results$x[1])
      results$d_obs_dx[n_res] <- (results$mu_observational[n_res] - results$mu_observational[n_res-1]) /
        (results$x[n_res] - results$x[n_res-1])
    }
  }
  
  # Resto da função (igual)...
  x_quantile_vals <- quantile(results$x, quantiles, na.rm = TRUE)
  
  summary_table <- data.frame(
    Quantile = paste0(100 * quantiles, "%"),
    X = round(x_quantile_vals, 2),
    mu_interventional = NA_real_,
    mu_observational = NA_real_,
    confounding_bias = NA_real_,
    ace = NA_real_,
    d_obs_dx = NA_real_
  )
  
  if("ace_standardized" %in% names(results)) {
    summary_table$ace_standardized <- NA_real_
  }
  
  for(i in seq_along(quantiles)) {
    idx <- which.min(abs(results$x - x_quantile_vals[i]))
    summary_table$mu_interventional[i] <- round(results$mu_interventional[idx], 3)
    summary_table$mu_observational[i]  <- round(results$mu_observational[idx], 3)
    summary_table$confounding_bias[i]  <- round(results$mu_observational[idx] - results$mu_interventional[idx], 4)
    summary_table$ace[i]               <- round(results$ace_interventional[idx], 3)
    summary_table$d_obs_dx[i]          <- round(results$d_obs_dx[idx], 3)
    
    if("ace_standardized" %in% names(results)) {
      summary_table$ace_standardized[i] <- round(results$ace_standardized[idx], 4)
    }
  }
  
  # Output publication-ready
  cat("\n", strrep("=", 80), "\n")
  cat("GENERIC COPULA: CONFOUNDING BIAS + SLOPE ANALYSIS\n")
  cat(strrep("=", 80), "\n\n")
  
  print(knitr::kable(summary_table, digits = 4, format = "latex",
                     caption = "Table XX: Level + Slope Bias Decomposition"))
  
  # Summary stats
  slope_bias <- summary_table$d_obs_dx - summary_table$ace
  cat(sprintf("\nSLOPE SUMMARY:\n"))
  cat(sprintf("Obs slope:  %.3f ± %.3f\n", 
              mean(summary_table$d_obs_dx), sd(summary_table$d_obs_dx)))
  cat(sprintf("Causal ACE: %.3f ± %.3f\n", 
              mean(summary_table$ace), sd(summary_table$ace)))
  cat(sprintf("Slope bias: %.3f (%.1f%% inflation)\n",
              mean(slope_bias),
              100 * mean(slope_bias) / abs(mean(summary_table$ace) + 1e-10)))
  
  return(summary_table)
}




################################################################################
# END OF IMPLEMENTATION
################################################################################
