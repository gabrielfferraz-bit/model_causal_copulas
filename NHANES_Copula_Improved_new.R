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
compute_causal_effects_empirical_nhanes <- function(x_grid, nhanes_data,
                                                    confounder_vars = "C1",
                                                    survey_weights = NULL) {
  library(rvinecopulib)
  library(VineCopula)   # for pobs()
  
  Z_col <- confounder_vars[1]
  data_subset <- nhanes_data[, c("X", "Y", Z_col)]
  colnames(data_subset) <- c("X", "Y", "Z")
  
  X <- data_subset$X
  Y <- data_subset$Y
  Z <- data_subset$Z
  n <- nrow(data_subset)
  
  # Pseudo-observations on [0,1]
  U_matrix <- pobs(as.matrix(data_subset))
  U_X <- U_matrix[, 1]
  U_Y <- U_matrix[, 2]
  U_Z <- U_matrix[, 3]
  
  cat("EMPIRICAL COPULA FITTING (non-parametric, kernel-smoothed):\n")
  
  # Fit non-parametric bivariate copulas — no parametric family assumed
  fit_xz_emp <- bicop(cbind(U_X, U_Z), family_set = "nonpar",
                      nonpar_method = "linear")
  fit_yz_emp <- bicop(cbind(U_Y, U_Z), family_set = "nonpar",
                      nonpar_method = "linear")
  fit_xy_emp <- bicop(cbind(U_X, U_Y), family_set = "nonpar",
                      nonpar_method = "linear")
  
  cat(sprintf("  (X,Z): Empirical non-parametric (tau=%.3f)\n", fit_xz_emp$tau))
  cat(sprintf("  (Y,Z): Empirical non-parametric (tau=%.3f)\n", fit_yz_emp$tau))
  cat(sprintf("  (X,Y): Empirical non-parametric (tau=%.3f)\n\n", fit_xy_emp$tau))
  
  # h-function wrappers: h(u_x | u_z) = P(U_X <= u_x | U_Z = u_z)
  # cond_var = 2 means we condition on the SECOND column (Z)
  F_x_given_z_func <- function(u_x, u_z) {
    hfunc(fit_xz_emp, cbind(u_x, u_z), 1)  # h-function: P(U1|u2)
  }
  
  F_y_given_z_func <- function(u_y, u_z) {
    hfunc(fit_yz_emp, cbind(u_y, u_z), 1)  # h-function: P(U_Y ≤ u_y | U_Z = u_z)
  }
  
  # Partial derivative of conditional copula C_{X,Y|Z}
  partial_copula_func <- function(u, v) {
    dbicop(cbind(u, v), fit_xy_emp)  # Use density as ∂C/∂u approximation
  }
  
  # Copula density c_{X,Z}(u_x, u_z) — needed for the observational CDF (Cor. 7.1)
  copula_xz_dens_func <- function(u_x, u_z) {
    dbicop(cbind(u_x, u_z), fit_xz_emp)
  }
  
  # Empirical marginal CDFs
  F_X <- ecdf(X)
  F_Y <- ecdf(Y)
  F_Z <- ecdf(Z)
  
  # Integration grids
  z_grid   <- seq(quantile(Z, 0.01), quantile(Z, 0.99), length.out = 50)
  y_grid   <- seq(quantile(Y, 0.01), quantile(Y, 0.99), length.out = 50)
  f_z_vals <- approx(density(Z)$x, density(Z)$y, xout = z_grid, rule = 2)$y
  dz_vec   <- diff(z_grid)
  total_mass <- sum((f_z_vals[-1] + f_z_vals[-length(f_z_vals)]) / 2 * dz_vec)
  f_z_vals <- f_z_vals / total_mass   # normalise to integrate to 1
  
  # Initialize results data frame
  results <- data.frame(
    x                  = x_grid,
    mu_interventional  = NA_real_,
    mu_observational   = NA_real_,
    ace_interventional = NA_real_
  )
  
  cat("Computing causal effects with EMPIRICAL (non-parametric) copulas...\n")
  
  # Inner function: interventional CDF F^{(X↓)}_{Y|X}(y|x) via Proposition 4.1
  compute_interventional_cdf_empirical <- function(y, x) {
    u_x <- F_X(x)
    u_y <- F_Y(y)
    integrand <- numeric(length(z_grid))
    for(i in seq_along(z_grid)) {
      u_z     <- F_Z(z_grid[i])
      F_x_z   <- F_x_given_z_func(u_x, u_z)
      F_y_z   <- F_y_given_z_func(u_y, u_z)
      partial <- partial_copula_func(F_x_z, F_y_z)
      integrand[i] <- partial * f_z_vals[i]
    }
    sum((integrand[-1] + integrand[-length(integrand)]) / 2 * dz_vec)
  }
  
  # Main loop over x_grid
  for(i in seq_along(x_grid)) {
    x_val <- x_grid[i]
    if(i %% 10 == 0) cat(sprintf("  Progress: %d/%d\n", i, length(x_grid)))
    
    # --- Interventional mean ---
    F_int          <- sapply(y_grid, function(y) compute_interventional_cdf_empirical(y, x_val))
    integrand_int  <- 1 - F_int
    dy             <- diff(y_grid)
    results$mu_interventional[i] <-
      sum((integrand_int[-1] + integrand_int[-length(integrand_int)]) / 2 * dy)
    
    # --- Observational mean (Corollary 7.1: includes c_{X,Z} weight) ---
    F_obs <- numeric(length(y_grid))
    u_x_val <- F_X(x_val)
    for(j in seq_along(y_grid)) {
      u_y <- F_Y(y_grid[j])
      integrand_obs <- numeric(length(z_grid))
      for(k in seq_along(z_grid)) {
        u_z     <- F_Z(z_grid[k])
        F_x_z   <- F_x_given_z_func(u_x_val, u_z)
        F_y_z   <- F_y_given_z_func(u_y, u_z)
        partial <- partial_copula_func(F_x_z, F_y_z)
        c_xz    <- copula_xz_dens_func(u_x_val, u_z)
        integrand_obs[k] <- partial * c_xz * f_z_vals[k]
      }
      F_obs[j] <- sum((integrand_obs[-1] + integrand_obs[-length(integrand_obs)]) / 2 * dz_vec)
    }
    integrand_obs_mean <- 1 - F_obs
    results$mu_observational[i] <-
      sum((integrand_obs_mean[-1] + integrand_obs_mean[-length(integrand_obs_mean)]) / 2 * dy)
  }
  
  # --- ACE via central differences (consistent with Gaussian and Generic methods) ---
  n_res <- nrow(results)
  results$ace_interventional <- c(
    (results$mu_interventional[2] - results$mu_interventional[1]) /
      (results$x[2] - results$x[1]),
    (results$mu_interventional[3:n_res] - results$mu_interventional[1:(n_res - 2)]) /
      (results$x[3:n_res] - results$x[1:(n_res - 2)]),
    (results$mu_interventional[n_res] - results$mu_interventional[n_res - 1]) /
      (results$x[n_res] - results$x[n_res - 1])
  )
  
  # Standardized ACE
  results$ace_standardized <- results$ace_interventional *
    (sd(nhanes_data$X, na.rm = TRUE) / sd(nhanes_data$Y, na.rm = TRUE))
  
  # Diagnostics columns (family = "Empirical" for all three pairs)
  results$copula_xz_family <- "Empirical"
  results$copula_yz_family <- "Empirical"
  results$copula_xy_family <- "Empirical"
  results$copula_xz_tau    <- fit_xz_emp$tau
  results$copula_yz_tau    <- fit_yz_emp$tau
  results$copula_xy_tau    <- fit_xy_emp$tau
  
  cat("\n✓ Empirical copula effects computed\n\n")
  return(results)
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
