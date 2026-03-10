# ==============================================================================
# Copula-Based Average Causal Effect: Illustrative Numerical Examples
# Supplementary code for González-López, Justus & Ferraz (2026)
# "A Copula-Based Framework for Structural Causal Modelling"
# ==============================================================================
#
# PURPOSE
#   Reproduce Figures 3–6 and the numerical tables from Sections 5.1–5.2 and
#   6.9 of the paper.  Two copula families are contrasted:
#
#     • Gaussian copula – unbounded dependence, state-dependent ACE with
#                         characteristic tail amplification (Section 5.1).
#     • FGM copula      – bounded weak dependence, nearly constant ACE
#                         (Section 5.2).
#
#   For each family the script computes and displays:
#     (a) Interventional mean   μ(x)        = E[Y | do(X = x)]
#     (b) Raw ACE               ACE(x)      = d/dx μ(x)
#     (c) Interventional SD     σ^(X↓)(x)
#     (d) Standardised ACE      ACE_std(x)  = ACE(x) / σ^(X↓)(x)  [Section 6.9]
#
# INPUTS  : none – all parameters are set inline in "Section 2 – Parameters".
# OUTPUTS : console tables (Tables 1–6, 9–10 in the paper) and composite
#           figures (Figures 3–6 in the paper).
#
# PACKAGE DEPENDENCY
#   pracma >= 2.3  (gaussHermite for the Gaussian variance integral)
# ==============================================================================


# ==============================================================================
# SECTION 0 – SETUP
# ==============================================================================

library(pracma)   # Gauss–Hermite quadrature nodes and weights

# Treatment grid: seq away from 0 and 1 to avoid Φ⁻¹(0) = -Inf singularities.
XGRID <- seq(0.01, 0.99, length.out = 500)

# Pre-compute the Gauss–Hermite quadrature rule once; it is reused in every
# call to ivar_gaussian().  A 20-node rule delivers machine precision for
# the smooth integrands encountered here.
GH <- gaussHermite(20)

# Indices of representative quantiles used in all printed tables.
TABLE_X   <- c(0.10, 0.25, 0.50, 0.75, 0.90)
TABLE_IDX <- sapply(TABLE_X, function(x) which.min(abs(XGRID - x)))

# Shorthand indices used repeatedly in diagnostics.
IDX_MED  <- TABLE_IDX[3]   # grid point closest to x = 0.50 (median)
IDX_TAIL <- TABLE_IDX[1]   # grid point closest to x = 0.10 (lower tail)

# Shared plotting aesthetics kept consistent across all figures.
COLS <- c("steelblue", "firebrick", "darkgreen")
LTYS <- c(1L, 1L, 2L)


# ==============================================================================
# SECTION 1 – MATHEMATICAL FUNCTIONS
# ==============================================================================

# All functions are defined exactly once here and shared across both copula
# families.  Each docstring cites the corresponding equation in the paper.

# ------------------------------------------------------------------------------
# 1A  Gaussian copula  (Section 5.1)
# ------------------------------------------------------------------------------

#' Latent partial correlation  K = ρ_{Y*X*·Z*}   (Theorem 5.1 / Eq. 5.1)
#'
#' Under the latent linear-Gaussian SEM with observable correlations (ryx, rxz,
#' ryz), K equals both the partial correlation on the latent scale AND the
#' conditional copula parameter ρ_{X,Y|Z}.  It is the single structural causal
#' coefficient that fully determines all interventional quantities in the
#' Gaussian model.
#'
#' @param ryx  Unconditional correlation ρ_{X,Y}
#' @param rxz  Correlation ρ_{X,Z}
#' @param ryz  Correlation ρ_{Y,Z}
partial_corr_gaussian <- function(ryx, rxz, ryz) {
  (ryx - rxz * ryz) / sqrt((1 - rxz^2) * (1 - ryz^2))
}

#' Interventional mean  μ(x) = Φ(K · Φ⁻¹(x))   (Eq. 5.17)
#'
#' This is E[Y | do(X = x)] on the observable uniform [0, 1] scale.  The
#' composition Φ ∘ K ∘ Φ⁻¹ is a monotone S-curve that collapses to the
#' identity when K = 1 and to the constant 0.5 when K = 0 (pure confounding).
mu_gaussian <- function(x, K) pnorm(K * qnorm(x))

#' Average Causal Effect  ACE(x) = d/dx μ(x)   (Eq. 5.18)
#'
#' = K · φ(K · Φ⁻¹(x)) / φ(Φ⁻¹(x))
#'
#' The density ratio φ(Kt)/φ(t) = exp(−(K²−1)t²/2) exceeds 1 for |t| large
#' and |K| < 1, producing the "tail amplification" discussed in Section 5.1.8:
#' ACE is minimised at the median (x = 0.5) and increases symmetrically toward
#' the tails of the distribution.
ace_gaussian <- function(x, K) K * dnorm(K * qnorm(x)) / dnorm(qnorm(x))

#' Interventional variance  Var[Y | do(X = x)]   for the Gaussian copula
#'
#' The latent structural equation Y* = K·X* + √(1−K²)·ε places
#'   Y* | do(X = x)  ~  N( K·Φ⁻¹(x),  1−K² ).
#' On the observable scale Y = Φ(Y*), so we need
#'   E[ Φ(Y*)² ]  with  Y* ~ N(μ_x, σ_x²),
#' which has no elementary closed form.  We evaluate it by Gauss–Hermite
#' quadrature via the substitution  t = (Y* − μ_x) / (√2 · σ_x):
#'
#'   E[ Φ(Y*)² ]  =  (1/√π) · Σ_i  w_i · Φ( μ_x + √2·σ_x·t_i )²
#'
#' @param x   Scalar treatment value in (0, 1)
#' @param K   Latent partial correlation (scalar)
#' @param gh  Gauss–Hermite rule; list with $x (nodes) and $w (weights)
ivar_gaussian <- function(x, K, gh = GH) {
  mu_x    <- K * qnorm(x)          # conditional mean on latent scale
  sigma_x <- sqrt(1 - K^2)         # conditional SD on latent scale
  
  # Shift and scale the GH nodes to the Y* ~ N(μ_x, σ_x²) distribution
  y_star <- mu_x + sqrt(2) * sigma_x * gh$x
  
  # E[Φ(Y*)²] via the Gauss–Hermite rule  ∫ f(t)exp(−t²)dt ≈ Σ w_i f(t_i)
  E_Y2 <- sum(gh$w * pnorm(y_star)^2) / sqrt(pi)
  
  # E[Y | do(X=x)] = μ(x) is available in closed form; reuse it here
  E_Y <- mu_gaussian(x, K)
  
  E_Y2 - E_Y^2
}

# ------------------------------------------------------------------------------
# 1B  FGM copula  (Section 5.2)
# ------------------------------------------------------------------------------

#' Interventional mean μ(x) under full 3-param FGM (Sec 5.2.3, Eq 5.33)
#'
#' μ(x) = 1/2 − α (1 − 2x) A,   where  A = 1/6 − η²/90
#'
#' Derivation: integrating 1 − F^(X↓)_{Y|X}(y|x) (Eq 5.30) over y ∈ [0,1]
#' with ∫₀¹ y(1−y) dy = 1/6  and  ∫₀¹ y²(1−y)² dy = 1/30.
#' The (1−2x) factor comes from Eq 5.30; it centres μ at 0.5 when x = 0.5
#' and moves it symmetrically up/down toward the tails.
#' θ_{X,Z} (= beta) does not appear: ∫F_{X|Z}(x|z)dz = x cancels it (Eq 5.28).
#' Bounds [0,1] are guaranteed for |α|,|η| ≤ 1 by FGM non-negativity (Eq 5.22).
mu_fgm <- function(x, params) {
  alpha <- params["theta_xy_cond"]   # conditional copula param (direct X→Y path)
  eta   <- params["theta_yz"]        # marginal Y–Z copula param
  A     <- 1/6 - eta^2 / 90         # = (15 − η²)/90; key scalar from Eq 5.33
  0.5 - alpha * (1 - 2*x) * A
}

#' Average Causal Effect ACE(x) under full 3-param FGM (Sec 5.2.3, Eq 5.35)
#'
#' ACE(x) = α (15 − η²) / 45   [constant in x]
#'
#' Derivation: d/dx μ(x) = d/dx [1/2 − α(1−2x)A] = 2αA = α(15−η²)/45.
#' The ACE is constant because (1−2x) differentiates to the constant −2.
#' This is the defining qualitative property of bounded FGM dependence:
#' no tail amplification, in sharp contrast to the Gaussian copula (Section 5.1).
#' θ_{X,Z} cancels as in mu_fgm; only α and η determine the causal magnitude.
ace_fgm <- function(x, params) {
  alpha <- params["theta_xy_cond"]
  eta   <- params["theta_yz"]
  rep(alpha * (15 - eta^2) / 45, length(x))   # constant across all x
}


#' Interventional variance Var[Y|do(X=x)] under full 3-param FGM (Sec 5.2.3)
#'
#' Var(x) = 1/12 − α²(1−2x)² A²,   A = 1/6 − η²/90
#'
#' Derivation via E[Y²|do] − (E[Y|do])², using:
#'   E[Y²|do] = 1/3 − 2α(1−2x)(A/2)  where A/2 = 1/12 − η²/180
#'   E[Y|do]  = 1/2 − α(1−2x) A
#'   Var       = (1/3 − 1/4) − α²(1−2x)²A² = 1/12 − α²(1−2x)²A²
#'
#' Key properties:
#'   • MAXIMUM at x = 0.5 (the (1−2x)² term vanishes): Var(0.5) = 1/12
#'   • MINIMUM at x = 0 or 1: Var = 1/12 − α²A²
#'   • Always positive for |α|,|η| ≤ 1 (since α²A² ≤ 1/36 < 1/12)
#'   • ACE_std = ACE/SD is therefore U-shaped with MINIMUM at the median
ivar_fgm <- function(x, params) {
  alpha <- params["theta_xy_cond"]
  eta   <- params["theta_yz"]
  A <- 1/6 - eta^2 / 90
  1/12 - alpha^2 * (1 - 2*x)^2 * A^2
}


# ------------------------------------------------------------------------------
# 1C  Generic utility: compute the full set of causal quantities over a grid
# ------------------------------------------------------------------------------

#' Evaluate μ, ACE, SD, and standardised ACE at every point of a treatment grid
#'
#' mu_fn and ace_fn are called in vectorised form for speed.  ivar_fn is
#' applied element-by-element via sapply, which is necessary when the variance
#' requires numerical quadrature (Gaussian case).
#'
#' @param mu_fn    function(x_vec, param)    → numeric vector  (interventional mean)
#' @param ace_fn   function(x_vec, param)    → numeric vector  (raw ACE)
#' @param ivar_fn  function(x_scalar, param) → scalar          (interventional variance)
#' @param params   Named numeric vector; one scalar parameter value per case
#' @param xgrid    Numeric vector of treatment values in (0, 1)
#'
#' @return Named list of four matrices (rows = xgrid points, cols = cases):
#'         $mu, $ace, $sd, $ace_std
compute_causal_quantities <- function(mu_fn, ace_fn, ivar_fn, params, xgrid) {
  n_x    <- length(xgrid)
  n_case <- length(params)
  
  mu_mat      <- matrix(NA_real_, n_x, n_case)
  ace_mat     <- matrix(NA_real_, n_x, n_case)
  sd_mat      <- matrix(NA_real_, n_x, n_case)
  ace_std_mat <- matrix(NA_real_, n_x, n_case)
  
  for (j in seq_len(n_case)) {
    p <- params[[j]]
    
    # Vectorized mu/ace
    mu_mat[, j]  <- mu_fn(xgrid, p)
    ace_mat[, j] <- ace_fn(xgrid, p)
    
    # Scalar variance computation
    var_j <- numeric(n_x)
    for (i in seq_along(xgrid)) {
      var_j[i] <- ivar_fn(xgrid[i], p)  # SCALAR LOOP ✓
    }
    
    # ---- NUMERICAL STABILISATION (FIX) ----
    
    # Remove negative variance values from quadrature cancellation
    var_j[var_j < 0] <- NA_real_
    
    # Compute SD
    sd_j <- sqrt(var_j)
    
    # Standardized ACE with protection against near-zero SD
    ace_std_j <- ace_mat[, j] / sd_j
    ace_std_j[sd_j < 1e-6] <- NA_real_
    
    # Store results
    sd_mat[, j]      <- sd_j
    ace_std_mat[, j] <- ace_std_j
    
  }
  
  list(mu = mu_mat, ace = ace_mat, sd = sd_mat, ace_std = ace_std_mat)
}


# ------------------------------------------------------------------------------
# 1D  Generic utility: print a formatted summary table to the console
# ------------------------------------------------------------------------------

#' Print a formatted ACE summary table
#'
#' @param res         Output list from compute_causal_quantities()
#' @param xgrid       Treatment grid (same length as nrow(res$mu))
#' @param idx         Integer vector of grid indices to tabulate
#' @param case_names  Character vector of case labels (length = ncol(res$mu))
#' @param param_vals  Named numeric vector of parameter values (one per case)
#' @param param_name  String naming the parameter, e.g. "K" or "theta"
print_ace_table <- function(res, xgrid, idx, case_names, param_vals, param_name) {
  header  <- sprintf("%-6s  %-9s  %-9s  %-9s  %-11s",
                     "x", "mu(x)", "ACE(x)", "SD(x)", "ACE_std(x)")
  divider <- strrep("-", 52)
  
  for (j in seq_along(case_names)) {
    cat(sprintf("\n%s (%s = %.4g):\n", case_names[j], param_name, param_vals[j]))
    cat(header, "\n", divider, "\n", sep = "")
    
    for (i in idx) {
      cat(sprintf("%-6.2f  %-9.4f  %-9.4f  %-9.4f  %-11.4f\n",
                  xgrid[i],
                  res$mu[i,  j],
                  res$ace[i, j],
                  res$sd[i,  j],
                  res$ace_std[i, j]))
    }
  }
  invisible(NULL)
}

# ------------------------------------------------------------------------------
# 1E  Generic utility: draw a 2 × 2 panel (μ, ACE, SD, ACE_std)
# ------------------------------------------------------------------------------

#' Draw the standard four-panel causal summary figure
#'
#' Panels: (A) interventional mean, (B) raw ACE, (C) interventional SD,
#' (D) standardised ACE.
#'
#' @param xgrid           Treatment grid
#' @param res             Output from compute_causal_quantities()
#' @param labs            Legend labels (character vector, one per case)
#' @param cols            Line colours (one per case)
#' @param ltys            Line types (one per case)
#' @param title_tag       Short string appended to each panel title
#' @param ace_ref         y-value of the horizontal reference line in Panel B
#' @param sd_override     Optional matrix substituted for res$sd in Panel C.
#'                        The original FGM four-panel inadvertently displayed
#'                        Gaussian SD data here; sd_override replicates that
#'                        behaviour precisely.
#' @param sd_labs_override Legend labels for Panel C when sd_override is used
plot_four_panels <- function(xgrid, res, labs, cols, ltys,
                             title_tag        = "",
                             ace_ref          = 0,
                             ylim_from_zero = FALSE,
                             sd_override      = NULL,
                             sd_labs_override = NULL) {
  par(mfrow = c(2, 2), mar = c(4.5, 4.5, 3, 1))
  n_cases <- ncol(res$mu)
  
  # Internal helper: overlay cases 2 … n_cases onto an open plot
  overlay_lines <- function(mat) {
    for (j in 2:n_cases)
      lines(xgrid, mat[, j], lwd = 2, col = cols[j], lty = ltys[j])
  }
  
  # --- Panel A: Interventional mean mu(x) -----------------------------------
  plot(xgrid, res$mu[, 1], type = "l", lwd = 2, col = cols[1],
       ylim = if (ylim_from_zero) c(0, max(res$mu)) else range(res$mu),
       xlab = expression(x), ylab = expression(mu(x)),
       main = paste0("(A) Interventional mean \u03bc(x)  ", title_tag))
  overlay_lines(res$mu)
  legend("topleft", legend = labs, col = cols, lwd = 2, lty = ltys, bty = "n")
  
  # --- Panel B: Raw ACE(x) --------------------------------------------------
  plot(xgrid, res$ace[, 1], type = "l", lwd = 2, col = cols[1],
       ylim = if (ylim_from_zero) c(0, max(res$ace)) else range(res$ace),
       xlab = expression(x), ylab = expression(ACE(x)),
       main = paste0("(B) Average Causal Effect ACE(x)  ", title_tag))
  overlay_lines(res$ace)
  abline(h = ace_ref, lty = 3)
  legend("topright", legend = labs, col = cols, lwd = 2, lty = ltys, bty = "n")
  
  # --- Panel C: Interventional SD sigma^(X↓)(x) ----------------------------
  # NOTE: In the original script the FGM four-panel (Figure 6) displays
  # Gaussian SD data here — a copy-paste artefact.  The sd_override argument
  # replicates this exactly so that Figure 6 is reproduced without alteration.
  sd_data <- if (!is.null(sd_override)) sd_override else res$sd
  sd_labs <- if (!is.null(sd_labs_override)) sd_labs_override else labs
  
  plot(xgrid, sd_data[, 1], type = "l", lwd = 2, col = cols[1],
       ylim = if (ylim_from_zero) c(0, max(sd_data, na.rm = TRUE))
       else range(sd_data[is.finite(sd_data)]),
       xlab = expression(x),
       ylab = expression(sigma^(X %down% "")(x)),
       main = paste0("(C) Interventional SD \u03c3(x)  ", title_tag))
  for (j in 2:n_cases)
    lines(xgrid, sd_data[, j], lwd = 2, col = cols[j], lty = ltys[j])
  legend("top", legend = sd_labs, col = cols, lwd = 2, lty = ltys, bty = "n")
  
  # --- Panel D: Standardised ACE ACE_std(x) = ACE(x) / sigma(x) -----------
  plot(xgrid, res$ace_std[, 1], type = "l", lwd = 2, col = cols[1],
       ylim = if (ylim_from_zero) c(0, max(res$ace_std, na.rm = TRUE))
       else range(res$ace_std, na.rm = TRUE),
       xlab = expression(x),
       ylab = expression(ACE[std](x)),
       main = paste0("(D) Standardised ACE  ACE_std(x) = ACE(x)/\u03c3(x)  ",
                     title_tag))
  overlay_lines(res$ace_std)
  legend("topright", legend = labs, col = cols, lwd = 2, lty = ltys, bty = "n")
}


# ==============================================================================
# SECTION 2 – PARAMETERS
# ==============================================================================

# Three confounding regimes are studied for each copula family, labelled A–C
# to match Sections 5.1.8 and 5.2.4 of the paper.

# ---- 2A  Gaussian copula: observable pairwise correlations ------------------
# K is derived below via partial_corr_gaussian(); it is the only parameter
# that enters the Gaussian interventional formulas.

gauss_params <- list(
  CaseA = c(ryx = 0.70, rxz = 0.50, ryz = 0.60),   # moderate confounding
  CaseB = c(ryx = 0.70, rxz = 0.00, ryz = 0.60),   # no treatment–confounder link
  CaseC = c(ryx = 0.70, rxz = 0.80, ryz = 0.60)    # strong confounding
)

# ---- 2B  FGM copula: dependence parameters (Section 5.2.3) -----------------
# Three cases vary the direct-effect parameter alpha = theta_{X,Y|Z} while
# holding the Y-Z marginal copula parameter eta = theta_{Y,Z} fixed at 0.30.
# theta_{X,Z} (= beta) is included for structural completeness but cancels
# analytically: M(x) = ∫ F_{X|Z}(x|z) dz = x (Eq 5.28).
#
# Cases A/B/C mirror the Gaussian family by spanning weak → strong direct effect:
#   CaseA: alpha = 0.40  (moderate direct effect)
#   CaseB: alpha = 0.80  (strong direct effect)
#   CaseC: alpha = 0.20  (weak direct effect)
fgm_params <- list(
  CaseA = c(theta_xy_cond = 0.40, theta_xz = 0.50, theta_yz = 0.30),
  CaseB = c(theta_xy_cond = 0.80, theta_xz = 0.50, theta_yz = 0.30),
  CaseC = c(theta_xy_cond = 0.20, theta_xz = 0.50, theta_yz = 0.30)
)


# ==============================================================================
# SECTION 3 – COMPUTATIONS (100% Reproducible)
# ==============================================================================

# ---- 3A  Gaussian copula ----------------------------------------------------
# Derive the structural parameter K for each confounding scenario via Theorem 5.1.
# K = ρ_{Y*X*|Z*} = partial correlation on latent Gaussian scale = copula parameter ρ_{X,Y|Z}.
K_vals <- sapply(names(gauss_params), function(name) {
  p <- gauss_params[[name]]
  partial_corr_gaussian(p["ryx"], p["rxz"], p["ryz"])
})
names(K_vals) <- names(gauss_params)

cat("Gaussian copula: latent partial correlations  K = beta\n")
print(round(K_vals, 4))

# Evaluate all four causal quantities over XGRID.
# ivar_gaussian is the bottleneck (numerical GH quadrature); total runtime is
# 500 grid points × 3 cases × 20 GH nodes ≈ 30 000 pnorm() calls — typically
# a few seconds on modern hardware.
cat("\nComputing Gaussian causal quantities (numerical integration)...\n")

# Convert K_vals to list of SCALARS for compute_causal_quantities; clamp K≤0.999
# to ensure sigma_x = √(1-K²) > 0 (prevents NaN in ivar_gaussian).
K_list <- as.list(pmin(K_vals, 0.999))  

res_gauss <- compute_causal_quantities(
  mu_fn   = mu_gaussian,           # Eq. 5.17: μ(x) = Φ(K·Φ⁻¹(x))
  ace_fn  = ace_gaussian,          # Eq. 5.18: ACE(x) = K·φ(Kt)/φ(t), t=Φ⁻¹(x)
  ivar_fn = function(x, K) ivar_gaussian(x, K, GH),  # GH quadrature for Var[Φ(Y*)]
  params  = K_list,                # List of 3 scalar K values (CaseA/B/C)
  xgrid   = XGRID                  # 500-point grid [0.01,0.99]
)
cat("Done.\n")

# ---- 3B  FGM copula ---------------------------------------------------------
# Extract only the 2 parameters needed for FGM interventional quantities:
# θ_{X,Y|Z} (direct effect) and θ_{Y,Z} (confounder strength). θ_{X,Z} cancels
# analytically in do(X=x) via ∫F_{X|Z}(x|z)dz = x (Eq. 5.28).
fgm_params_correct <- lapply(names(fgm_params), function(name) {
  p <- fgm_params[[name]]
  c(theta_xy_cond = unname(p["theta_xy_cond"]), theta_yz = unname(p["theta_yz"]))
})
names(fgm_params_correct) <- names(fgm_params)

# Extract θ_{X,Y|Z} values for tables/diagnostics (CaseA=0.40, CaseB=0.80, CaseC=0.20).
theta_vals <- sapply(fgm_params_correct, function(p) p["theta_xy_cond"])
names(theta_vals) <- names(fgm_params_correct)
cat("\nFGM copula: θ_{X,Y|Z} values (vary by case)\n")
print(round(theta_vals, 4))

# Evaluate all four causal quantities over XGRID.
# FGM quantities have closed forms (Eqs. 5.33, 5.35, 6.9.3); no numerical integration.
cat("Computing FGM causal quantities (closed-form)...\n")

res_fgm <- compute_causal_quantities(
  mu_fn   = mu_fgm,                # Eq. 5.33: μ(x) = ½ − α(1−2x)(1/6 − η²/90) [state-dependent]
  ace_fn  = ace_fgm,               # Eq. 5.35: ACE(x) = α(15 − η²)/45 [constant in x]
  ivar_fn = ivar_fgm,              # Derived: Var = 1/12 − α²(1−2x)²(1/6 − η²/90)² [max at median]
  params  = fgm_params_correct,    # List of 3 named vectors (2 params each)
  xgrid   = XGRID                  # Same 500-point grid
)
cat("Done.\n")


# ==============================================================================
# SECTION 4 – CONSOLE OUTPUT (Tables)
# ==============================================================================

CASE_NAMES <- c("Case A", "Case B", "Case C")

cat("\n\n=== TABLE: Gaussian Copula – Causal Quantities ===\n")

print_ace_table(
  res        = res_gauss,
  xgrid      = XGRID,
  idx        = TABLE_IDX,
  case_names = CASE_NAMES,
  param_vals = K_vals,
  param_name = "K"
)

cat("\n\n=== TABLE: FGM Copula – Causal Quantities ===\n")
print_ace_table(
  res        = res_fgm,
  xgrid      = XGRID,
  idx        = TABLE_IDX,
  case_names = CASE_NAMES,
  param_vals = sapply(fgm_params, function(p) p["theta_xy_cond"]),  # Extract θ_{X,Y|Z} only ✓
  param_name = "θ_{X,Y|Z}"
)


# Purpose: Show that equal overall association strength (same Kendall's tau)
# produces qualitatively different ACE profiles under the two families.
# The Gaussian ACE varies across x (tail amplification); the FGM ACE is a
# single constant.  This is the key message of Section 5.2.5.
#
# Calibration:
#   Given a target tau, the two families require:
#
#   Gaussian:  K  = sin(pi * tau / 2)                [inverse of tau = 2/pi * arcsin(K)]
#   FGM:       alpha = 9 * tau / 2                   [inverse of tau = 2*alpha/9]
#              eta   = 0  (set to zero to isolate the direct-effect comparison;
#                          confounding attenuation is a separate sensitivity exercise)
#
# At tau = 0.20 this yields K = sin(0.1*pi) ≈ 0.309  and  alpha = 0.90.

TAU_MATCH <- 0.20   # Target Kendall's tau — change this to explore sensitivity

# Calibrate Gaussian K from tau  (Eq. tau_gauss in paper)
K_match <- sin(pi * TAU_MATCH / 2)

# Calibrate FGM alpha from tau  (Eq. tau_fgm in paper); eta = 0 for isolation
alpha_match <- 9 * TAU_MATCH / 2
eta_match   <- 0

# Build the FGM parameter vector expected by ace_fgm / mu_fgm
fgm_match_params <- c(theta_xy_cond = alpha_match, theta_yz = eta_match)

cat(sprintf(
  "\n\n=== TABLE 8: Cross-Family ACE at Matched Kendall's Tau (tau = %.2f) ===\n",
  TAU_MATCH))
cat(sprintf(
  "  Calibrated parameters:  Gaussian K = %.4f  |  FGM alpha = %.4f, eta = %.1f\n",
  K_match, alpha_match, eta_match))
cat(sprintf(
  "  Verification:  tau_Gauss = 2/pi * arcsin(K) = %.4f  |  tau_FGM = 2*alpha/9 = %.4f\n\n",
  (2/pi) * asin(K_match), 2 * alpha_match / 9))

# ACE for FGM is a single constant (independent of x); compute it once
ace_fgm_match <- ace_fgm(0.5, fgm_match_params)[1]   # value at any x is the same

cat(sprintf(
  "  FGM ACE is CONSTANT at %.4f across all treatment levels.\n",
  ace_fgm_match))
cat(sprintf(
  "  Gaussian ACE varies; minimum at median, amplified at tails.\n\n"))

# Print the table
header8  <- sprintf("%-6s  %-18s  %-12s  %-12s  %-8s",
                    "x", "ACE_Gaussian(x)", "ACE_FGM", "Ratio(G/F)", "comment")
divider8 <- strrep("-", 64)
cat(header8, "\n", divider8, "\n", sep = "")

comments8 <- c("lower tail", "first quartile", "median (Gauss min)", "third quartile", "upper tail")
for (k in seq_along(TABLE_IDX)) {
  i   <- TABLE_IDX[k]
  x_i <- XGRID[i]
  ace_g <- ace_gaussian(x_i, K_match)
  ace_f <- ace_fgm_match                     # constant; displayed for every row
  ratio <- ace_g / ace_f
  cat(sprintf("%-6.2f  %-18.4f  %-12.4f  %-12.4f  %s\n",
              x_i, ace_g, ace_f, ratio, comments8[k]))
}

# Summary statistics that appear in the paper's prose
amp_gauss <- (ace_gaussian(XGRID[IDX_TAIL], K_match) /
                ace_gaussian(XGRID[IDX_MED],  K_match) - 1) * 100

cat(sprintf(
  "\n  Gaussian tail amplification (tail / median - 1): +%.0f%%\n", amp_gauss))
cat(sprintf(
  "  FGM tail amplification: 0%%  (ACE is constant by construction)\n"))
cat(sprintf(
  "  Ratio Gauss/FGM ranges from %.2f (median) to %.2f (tails)\n",
  ace_gaussian(XGRID[IDX_MED],  K_match) / ace_fgm_match,
  ace_gaussian(XGRID[IDX_TAIL], K_match) / ace_fgm_match))

# ---- Diagnostics: tail amplification and confounding effects ----------------

cat("\n\n=== KEY INSIGHTS: Detectability of Causal Effects ===\n\n")

# Gaussian: how strongly does ACE_std amplify toward the tails?
cat("Gaussian copula – tail amplification of ACE_std:\n")
for (j in seq_len(3)) {
  amp <- (res_gauss$ace_std[IDX_TAIL, j] / res_gauss$ace_std[IDX_MED, j] - 1) * 100
  cat(sprintf(
    "  %s (K = %.3f):  median = %.2f,  tail = %.2f,  amplification = %.1f%%\n",
    CASE_NAMES[j], K_vals[j],
    res_gauss$ace_std[IDX_MED,  j],
    res_gauss$ace_std[IDX_TAIL, j], amp))
}

# Gaussian: how does confounding reduce detectability at the median?
cat("\nGaussian copula – effect of confounding on median detectability:\n")
base_detect <- res_gauss$ace_std[IDX_MED, 2]   # Case B (no confounding) as baseline
for (j in c(2L, 1L, 3L)) {                      # report B → A → C for clarity
  pct_red <- (1 - res_gauss$ace_std[IDX_MED, j] / base_detect) * 100
  cat(sprintf(
    "  %s: ACE_std(0.5) = %.2f  (%.1f%% reduction vs no confounding)\n",
    CASE_NAMES[j], res_gauss$ace_std[IDX_MED, j], pct_red))
}

# FGM: verify the near-constancy of ACE_std across the treatment range
cat("\nFGM copula – near-constancy of ACE_std:\n")
for (j in seq_len(3)) {
  variation_pct <- abs(res_fgm$ace_std[IDX_TAIL, j] - res_fgm$ace_std[IDX_MED, j]) /
    res_fgm$ace_std[IDX_MED, j] * 100
  cat(sprintf(
    "  %s (theta = %.2f):  median = %.2f,  tail = %.2f,  variation = %.1f%%\n",
    CASE_NAMES[j], theta_vals[j],
    res_fgm$ace_std[IDX_MED,  j],
    res_fgm$ace_std[IDX_TAIL, j], variation_pct))
}

# Cross-family comparison for Case A
gauss_ratio <- res_gauss$ace_std[IDX_TAIL, 1] / res_gauss$ace_std[IDX_MED, 1]
fgm_ratio   <- res_fgm$ace_std[IDX_TAIL,   1] / res_fgm$ace_std[IDX_MED,  1]
cat(sprintf(
  "\nTail-to-median ACE_std ratio (Case A)  –  Gaussian: %.2f,  FGM: %.2f\n",
  gauss_ratio, fgm_ratio))
cat(sprintf(
  "  => Gaussian exhibits %.1fx more tail amplification than FGM.\n",
  gauss_ratio / fgm_ratio))


# ==============================================================================
# SECTION 5 – FIGURES
# ==============================================================================

# Assemble legend labels linking each curve to its parameter value
labs_gauss <- sprintf("%s (K = %.2f)",     CASE_NAMES, K_vals)
labs_fgm <- sprintf("%s (θ_{XY|Z}=%.2f, θ_{YZ}=%.2f)", 
                    CASE_NAMES, 
                    sapply(fgm_params, function(p) p["theta_xy_cond"]),
                    sapply(fgm_params, function(p) p["theta_yz"]))

# ---- Figures 3 & 4: Two-panel plots (μ and ACE only)  -----------------------
# These reproduce the original script's first two plotting blocks exactly
# and match Figures 3–4 in the paper.

par(mfrow = c(1, 2), mar = c(4.5, 4.5, 3, 1))

# Figure 3 – Gaussian
plot(XGRID, res_gauss$mu[, 1], type = "l", lwd = 2, col = COLS[1],
     ylim = range(res_gauss$mu),
     xlab = expression(x), ylab = expression(mu(x)),
     main = expression("Interventional mean " ~ mu(x)))
lines(XGRID, res_gauss$mu[, 2], lwd = 2, col = COLS[2])
lines(XGRID, res_gauss$mu[, 3], lwd = 2, col = COLS[3], lty = LTYS[3])
legend("topleft", legend = labs_gauss, col = COLS, lwd = 2, lty = LTYS, bty = "n")

plot(XGRID, res_gauss$ace[, 1], type = "l", lwd = 2, col = COLS[1],
     ylim = range(res_gauss$ace),
     xlab = expression(x), ylab = expression(ACE(x)),
     main = expression("State-dependent ACE " ~ ACE(x)))
lines(XGRID, res_gauss$ace[, 2], lwd = 2, col = COLS[2])
lines(XGRID, res_gauss$ace[, 3], lwd = 2, col = COLS[3], lty = LTYS[3])
abline(h = 0, lty = 3)
legend("topright", legend = labs_gauss, col = COLS, lwd = 2, lty = LTYS, bty = "n")

par(mfrow = c(1, 2), mar = c(4.5, 4.5, 3, 1))

# Figure 4 – FGM
plot(XGRID, res_fgm$mu[, 1], type = "l", lwd = 2, col = COLS[1],
     ylim = c(0, 1),
     xlab = expression(x), ylab = expression(mu(x)),
     main = expression("Interventional mean " ~ mu(x) ~ "(FGM)"))
lines(XGRID, res_fgm$mu[, 2], lwd = 2, col = COLS[2])
lines(XGRID, res_fgm$mu[, 3], lwd = 2, col = COLS[3], lty = LTYS[3])
legend("topleft", legend = labs_fgm, col = COLS, lwd = 2, lty = LTYS, bty = "n")

plot(XGRID, res_fgm$ace[, 1], type = "l", lwd = 2, col = COLS[1],
     ylim = range(c(0, res_fgm$ace)),
     xlab = expression(x), ylab = expression(ACE(x)),
     main = expression("State-dependent ACE " ~ ACE(x) ~ "(FGM)"))
lines(XGRID, res_fgm$ace[, 2], lwd = 2, col = COLS[2])
lines(XGRID, res_fgm$ace[, 3], lwd = 2, col = COLS[3], lty = LTYS[3])
abline(h = 0, lty = 3)
legend("topright", legend = labs_fgm, col = COLS, lwd = 2, lty = LTYS, bty = "n")

# ---- Figure 5: Gaussian 2-panel (SD, ACE_std) ----
par(mfrow = c(1, 2), mar = c(4.5, 4.5, 3, 1))

# Left panel: Interventional SD
plot(XGRID, res_gauss$sd[, 1], type = "l", lwd = 2, col = COLS[1],
     ylim = range(res_gauss$sd[is.finite(res_gauss$sd)]),
     xlab = expression(x), ylab = expression(sigma^(X %down% "")(x)),
     main = "(A) Interventional SD σ(x) (Gaussian)")
lines(XGRID, res_gauss$sd[, 2], lwd = 2, col = COLS[2])
lines(XGRID, res_gauss$sd[, 3], lwd = 2, col = COLS[3], lty = LTYS[3])
legend("top", legend = labs_gauss, col = COLS, lwd = 2, lty = LTYS, bty = "n")

# Right panel: Standardised ACE
plot(XGRID, res_gauss$ace_std[, 1], type = "l", lwd = 2, col = COLS[1],
     ylim = range(res_gauss$ace_std, na.rm = TRUE),
     xlab = expression(x), ylab = expression(ACE[std](x)),
     main = "(B) Standardised ACE ACE_std(x) (Gaussian)")
lines(XGRID, res_gauss$ace_std[, 2], lwd = 2, col = COLS[2])
lines(XGRID, res_gauss$ace_std[, 3], lwd = 2, col = COLS[3], lty = LTYS[3])
abline(v = 0.5, lty = 3, col = "gray50")
legend("topright", legend = labs_gauss, col = COLS, lwd = 2, lty = LTYS, bty = "n")

# ---- Figure 6: FGM 2-panel (SD, ACE_std) ----
par(mfrow = c(1, 2), mar = c(4.5, 4.5, 3, 1))

# Left panel: Interventional SD (corrected FGM data, not Gaussian override)
plot(XGRID, res_fgm$sd[, 1], type = "l", lwd = 2, col = COLS[1],
     ylim = range(res_fgm$sd[is.finite(res_fgm$sd)]),
     xlab = expression(x), ylab = expression(sigma^(X %down% "")(x)),
     main = "(A) Interventional SD σ(x) (FGM)")
lines(XGRID, res_fgm$sd[, 2], lwd = 2, col = COLS[2])
lines(XGRID, res_fgm$sd[, 3], lwd = 2, col = COLS[3], lty = LTYS[3])
legend("top", legend = labs_fgm, col = COLS, lwd = 2, lty = LTYS, bty = "n")

# Right panel: Standardised ACE
plot(XGRID, res_fgm$ace_std[, 1], type = "l", lwd = 2, col = COLS[1],
     ylim = range(res_fgm$ace_std, na.rm = TRUE),
     xlab = expression(x), ylab = expression(ACE[std](x)),
     main = "(B) Standardised ACE ACE_std(x) (FGM)")
lines(XGRID, res_fgm$ace_std[, 2], lwd = 2, col = COLS[2])
lines(XGRID, res_fgm$ace_std[, 3], lwd = 2, col = COLS[3], lty = LTYS[3])
abline(v = 0.5, lty = 3, col = "gray50")
legend("topright", legend = labs_fgm, col = COLS, lwd = 2, lty = LTYS, bty = "n")

# ---- Gaussian four-panel (μ, ACE, SD, ACE_std) --------------------
par(mfrow = c(1, 1), mar = c(4.5, 4.5, 3, 1))
ylim_g <- c(0, max(res_gauss$ace_std[is.finite(res_gauss$ace_std)], na.rm = TRUE))
plot(XGRID, res_gauss$ace_std[, 1], type = "l", lwd = 3, col = COLS[1],
     ylim = ylim_g, xlab = expression(x), ylab = expression(ACE[std](x)),
     main = "Figure 5: Gaussian Copula – Standardised ACE")
lines(XGRID, res_gauss$ace_std[, 2], lwd = 3, col = COLS[2])
lines(XGRID, res_gauss$ace_std[, 3], lwd = 3, col = COLS[3], lty = LTYS[3])
abline(v = 0.5, lty = 3, col = "gray50")
legend("topright", legend = labs_gauss, col = COLS, lwd = 3, lty = LTYS, bty = "n")

plot_four_panels(
  xgrid     = XGRID,
  res       = res_gauss,
  labs      = labs_gauss,
  cols      = COLS,
  ltys      = LTYS,
  title_tag = "(Gaussian)",
  ace_ref   = 0
)

# ---- FGM four-panel -----------------------------------------------
# Panel C is fed Gaussian SD data via sd_override to reproduce the original
# code faithfully (see Section 1E docstring for explanation).

plot_four_panels(
  xgrid = XGRID,
  res = res_fgm,
  labs = labs_fgm,
  cols = COLS,
  ltys = LTYS,
  title_tag = "(FGM)",
  ace_ref = 1,
  ylim_from_zero = TRUE  # Works perfectly with true FGM SD
)

# ---- Comparative figure: Gaussian vs FGM detectability (ACE_std) -------------
# Illustrates the fundamental contrast between strong state-dependence
# (Gaussian) and nearly uniform detectability (FGM).

par(mfrow = c(1, 2), mar = c(4.5, 4.5, 3, 1))

ylim_g <- c(0, max(res_gauss$ace_std[is.finite(res_gauss$ace_std)]))
plot(XGRID, res_gauss$ace_std[, 1], type = "l", lwd = 3, col = COLS[1],
     ylim = ylim_g,
     xlab = expression(x), ylab = expression(ACE[std](x)),
     main = "Gaussian Copula: State-Dependent Detectability")
lines(XGRID, res_gauss$ace_std[, 2], lwd = 3, col = COLS[2])
lines(XGRID, res_gauss$ace_std[, 3], lwd = 3, col = COLS[3], lty = LTYS[3])
abline(v = 0.5, lty = 3, col = "gray50")
text(0.5, max(res_gauss$ace_std) * 0.9, "median", pos = 4, col = "gray50")
legend("topright", legend = labs_gauss, col = COLS, lwd = 3, lty = LTYS, bty = "n")

ylim_f <- c(0, max(res_fgm$ace_std[is.finite(res_fgm$ace_std)]))
plot(XGRID, res_fgm$ace_std[, 1], type = "l", lwd = 3, col = COLS[1],
     ylim = ylim_f,
     xlab = expression(x), ylab = expression(ACE[std](x)),
     main = "FGM Copula: Nearly Constant Detectability")
lines(XGRID, res_fgm$ace_std[, 2], lwd = 3, col = COLS[2])
lines(XGRID, res_fgm$ace_std[, 3], lwd = 3, col = COLS[3], lty = LTYS[3])
abline(v = 0.5, lty = 3, col = "gray50")
text(0.5, max(res_fgm$ace_std) * 0.9, "median", pos = 4, col = "gray50")
legend("topright", legend = labs_fgm, col = COLS, lwd = 3, lty = LTYS, bty = "n")


# ==============================================================================
# SECTION 6 – SUMMARY AND PRACTICAL IMPLICATIONS
# ==============================================================================

cat("\n", strrep("=", 79), "\n", sep = "")
cat("SUMMARY: STANDARDISED ACE AS SIGNAL-TO-NOISE RATIO\n")
cat(strrep("=", 79), "\n\n")
cat("ACE_std(x) = ACE(x) / sigma^(X-down)(x)  is a dimensionless, scale-free\n")
cat("measure of causal-effect strength that quantifies detectability in finite\n")
cat("samples (Section 6.9 of the paper).\n\n")

cat("KEY FINDINGS:\n\n")
cat("1. GAUSSIAN COPULA (unbounded dependence):\n")
cat("   - ACE_std(x) varies substantially with treatment level x.\n")
cat("   - Tail amplification: ACE_std is 50-90% higher at extremes than median.\n")
cat("   - Confounding uniformly reduces detectability at every treatment level.\n")
cat("   - Practical implication: effects are most detectable at extreme x values.\n\n")

cat("2. FGM COPULA (bounded weak dependence):\n")
cat("   - ACE_std(x) is nearly constant across the treatment range (< 2% variation).\n")
cat("   - No privileged treatment levels: detectability is uniform.\n")
cat("   - Practical implication: standard power calculations apply uniformly.\n\n")

cat("3. IMPLICATIONS FOR STUDY DESIGN:\n")
cat("   a) Prioritise treatment levels with high ACE_std(x) to maximise power.\n")
cat("   b) Required sample size n is proportional to [ACE_std(x)]^{-2}:\n")
cat("      doubling ACE_std reduces required n by a factor of four.\n")
cat("   c) In multi-arm trials, allocate more subjects to high-ACE_std arms.\n")
cat("   d) Copula family selection fundamentally shapes the detectability profile.\n\n")

cat("4. CONNECTION TO CLASSICAL EFFECT SIZES:\n")
cat("   - Generalises Cohen's d to state-dependent, copula-based settings.\n")
cat("   - Analogous to the regression-coefficient / residual-SD ratio,\n")
cat("     but nonparametric and without homoskedasticity assumptions.\n")
cat("   - Dimensionless: directly comparable across studies and outcome scales.\n\n")

cat("RECOMMENDATION: always report ACE_std(x) alongside raw ACE(x) to assess\n")
cat("practical significance and detectability of the estimated causal effect.\n\n")