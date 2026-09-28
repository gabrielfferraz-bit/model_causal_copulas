# ==============================================================================
# Illustrative Numerical Examples: Gaussian and FGM Copulas
# ==============================================================================
#
# This script reproduces the numerical examples in Chapter 4.
#
# The script implements two copula families:
#     - Gaussian copula: mu(x) is U-shaped with tail amplification, governed
#                         by the structural coefficient K (Eq. 5.16), which is
#                         non-monotonic in rho_XZ (Section 5.6.1).
#     - FGM copula      : mu(x) is exactly affine in x, with beta cancelling
#                          out of the mean by an antisymmetry argument
#                          (Eq. 5.23; Section 5.6.2).
#
# For each family and each confounding regime the script computes and
# tabulates, over a grid of treatment quantiles x:
#     (a) Observational mean   E[Y | X = x]           (Eqs. 5.17, 5.24)
#     (b) Interventional mean  mu(x) = E[Y | do(X=x)]  (Eqs. 5.16, 5.23)
#     (c) Confounding bias     B(x) = E[Y|X=x] - mu(x) (Eqs. 5.17, 5.25)
#
# These three quantities -- E[Y|X=x], mu(x), and the analytic confounding
# bias B(x) -- are the only causal quantities reported in this version of
# the dissertation; no ACE, interventional SD, or standardised-ACE object
# is computed.
#
# INPUTS  : none -- all parameters are set inline in "Section 2 - Parameters".
# OUTPUTS : console tables (Table 3: Gaussian; Table 4: FGM) and composite
#           figures (Figure 12: Gaussian; Figure 13: FGM), each showing
#           E[Y|X=x] and mu(x) (left panel) and B(x) (right panel).
#
# PACKAGE DEPENDENCY
#   None beyond base R (stats for pnorm/qnorm/dnorm).
# ==============================================================================


# ==============================================================================
# SECTION 0 - SETUP
# ==============================================================================

# Treatment grid: kept away from 0 and 1 to avoid Phi^{-1}(0) = -Inf
# singularities in the Gaussian family. The representative table values are
# included explicitly so that TABLE_IDX refers to the exact stated x-values.
TABLE_X <- c(0.10, 0.25, 0.50, 0.75, 0.90)

XGRID <- sort(unique(c(
  seq(0.01, 0.99, length.out = 500),
  TABLE_X
)))

# Indices of the representative quantiles used in Tables 3-4.
TABLE_IDX <- sapply(TABLE_X, function(x) which(XGRID == x)[1])

# Shared plotting aesthetics kept consistent across all figures.
COLS <- c("steelblue", "firebrick", "darkgreen")
LTYS <- c(1L, 1L, 2L)


# ==============================================================================
# SECTION 1 - MATHEMATICAL FUNCTIONS
# ==============================================================================

# All functions are defined exactly once here and shared across both copula
# families. Each docstring cites the corresponding equation in Chapter 5 of
# the dissertation.

# ------------------------------------------------------------------------------
# 1A  Gaussian copula (Section 5.6.1)
# ------------------------------------------------------------------------------

#' Structural coefficients of the Gaussian-copula interventional model
#' (dissertation Sec. 5.6.1, Step 2, Eq. 5.16).
#'
#' Regressing Y* on (X*, Z*) under the trivariate-normal representation
#' (X* = Phi^{-1}(X), etc., all marginally N(0,1) with correlation matrix
#' built from rho_XY, rho_XZ, rho_YZ -- see the "Why the asterisk" remark)
#' gives
#'   a = (rho_XY - rho_XZ*rho_YZ) / (1 - rho_XZ^2),
#'   b = (rho_YZ - rho_XY*rho_XZ) / (1 - rho_XZ^2),
#' and the interventional latent variance
#'   V = 1 - a^2 - 2*a*b*rho_XZ,
#' so that Y* | do(X* = x*) ~ N(a * Phi^{-1}(x), V) exactly. Marginalising
#' the probit link Phi over this normal (via the elementary identity Eq. 5.15)
#' gives the interventional-mean slope
#'   K = a / sqrt(1 + V),   so that   mu(x) = Phi(K * Phi^{-1}(x)).
#'
#' The no-confounding value K_0 (i.e. K evaluated at rho_XZ = 0) governs the
#' OBSERVATIONAL mean instead: E[Y|X=x] = Phi(K_0 * Phi^{-1}(x)), since
#' Y* | X* = x* ~ N(rho_XY * Phi^{-1}(x), 1 - rho_XY^2) directly, with no Z*
#' to marginalise (Step 3). K_0 = rho_XY / sqrt(2 - rho_XY^2).
gaussian_K_and_V <- function(rho_xy, rho_xz, rho_yz) {
  a <- (rho_xy - rho_xz * rho_yz) / (1 - rho_xz^2)
  b <- (rho_yz - rho_xy * rho_xz) / (1 - rho_xz^2)
  V <- 1 - a^2 - 2 * a * b * rho_xz
  K <- a / sqrt(1 + V)
  c(a = unname(a), b = unname(b), V = unname(V), K = unname(K))
}

#' No-confounding structural coefficient K_0 = K evaluated at rho_XZ = 0
#' (Eq. 5.17): K_0 = rho_XY / sqrt(2 - rho_XY^2).
gaussian_K0 <- function(rho_xy) {
  rho_xy / sqrt(2 - rho_xy^2)
}

#' Interventional mean  mu(x) = Phi(K * Phi^{-1}(x))   (Eq. 5.16)
#'
#' E[Y | do(X = x)] on the observable uniform [0, 1] scale.
mu_gaussian <- function(x, K) pnorm(K * qnorm(x))

#' Observational mean  E[Y | X = x] = Phi(K_0 * Phi^{-1}(x))   (Eq. 5.17)
#'
#' Same functional form as mu_gaussian(), but driven by K_0 rather than K.
eobs_gaussian <- function(x, K0) pnorm(K0 * qnorm(x))

#' Confounding bias  B(x) = E[Y|X=x] - mu(x)
bias_gaussian <- function(x, K, K0) {
  eobs_gaussian(x, K0) - mu_gaussian(x, K)
}


# ------------------------------------------------------------------------------
# 1B  FGM copula (Section 5.6.2)
# ------------------------------------------------------------------------------

#' Interventional mean mu(x) under the 3-parameter FGM model
#' (Eq. 5.23): mu(x) = 1/2 - alpha*(1-2x)*(1/6 - eta^2/90).
#'
#' Derivation: integrating 1 - F^{(X-down)}_{Y|X}(y|x) (Eq. 5.21) over
#' y in [0,1], using int_0^1 y(1-y) dy = 1/6 and
#' int_0^1 y^2(1-y)^2 dy = 1/30. The beta-term in the CDF contains the
#' factor y(1-y)(1-2y), whose integral vanishes by antisymmetry about
#' y = 1/2. Thus beta cancels from mu(x), although it remains in the full
#' interventional distribution through F_{X|Z}(x|z).
mu_fgm <- function(x, params) {
  alpha <- params["theta_xy_cond"]   # direct effect X -> Y
  eta   <- params["theta_yz"]        # outcome-confounder parameter
  A     <- 1/6 - eta^2 / 90
  
  unname(
    0.5 - alpha * (1 - 2*x) * A
  )
}

#' Observational mean E[Y | X = x] under the confounded FGM DAG (Eq. 5.24).
#'
#' E[Y|X=x] = 1/2 - (1-2x) * [ beta*eta/18 + alpha*(15-eta^2)/90
#'                              - alpha*beta^2*(25-3*eta^2)/225 * x*(1-x) ]
eobs_fgm <- function(x, params) {
  alpha <- params["theta_xy_cond"]
  beta  <- params["theta_xz"]
  eta   <- params["theta_yz"]
  
  bracket <- beta * eta / 18 +
    alpha * (15 - eta^2) / 90 -
    alpha * beta^2 * (25 - 3*eta^2) / 225 * x * (1 - x)
  
  unname(
    0.5 - (1 - 2*x) * bracket
  )
}

#' Confounding bias B(x) = E[Y|X=x] - mu(x) under the FGM model.
#'
#' Subtracting Eqs. 5.23 and 5.24 gives
#'
#'   B(x) = -(beta*eta/18)(1-2x)
#'          + alpha*beta^2*(25-3*eta^2)/225
#'            * x*(1-x)*(1-2x).
#'
#' The beta*eta/18 term is the first-order back-door contribution carried by
#' the Z -> X and Z -> Y legs together. The alpha*(15-eta^2)/90 contribution
#' cancels exactly against the corresponding term in mu(x). The remaining
#' alpha-dependent contribution is the higher-order alpha*beta^2 term.
#' Thus eta enters the bias both through the first-order back-door term and
#' through the higher-order alpha*beta^2 term.
bias_fgm <- function(x, params) {
  alpha <- params["theta_xy_cond"]
  beta  <- params["theta_xz"]
  eta   <- params["theta_yz"]
  
  unname(
    -beta * eta / 18 * (1 - 2*x) +
      alpha * beta^2 * (25 - 3*eta^2) / 225 *
      x * (1 - x) * (1 - 2*x)
  )
}


# ------------------------------------------------------------------------------
# 1C  Generic utility: compute E[Y|X=x], mu(x), and B(x) over a grid
# ------------------------------------------------------------------------------

#' Evaluate the observational mean, interventional mean, and confounding
#' bias at every point of a treatment grid, for a list of parameter cases.
compute_causal_quantities <- function(eobs_fn, mu_fn, bias_fn, params, xgrid) {
  n_x    <- length(xgrid)
  n_case <- length(params)
  
  eobs_mat <- matrix(NA_real_, n_x, n_case)
  mu_mat   <- matrix(NA_real_, n_x, n_case)
  bias_mat <- matrix(NA_real_, n_x, n_case)
  
  for (j in seq_len(n_case)) {
    p <- params[[j]]
    eobs_mat[, j] <- eobs_fn(xgrid, p)
    mu_mat[, j]   <- mu_fn(xgrid, p)
    bias_mat[, j] <- bias_fn(xgrid, p)
  }
  
  list(eobs = eobs_mat, mu = mu_mat, bias = bias_mat)
}


# ------------------------------------------------------------------------------
# 1D  Generic utility: print a formatted summary table to the console
# ------------------------------------------------------------------------------

#' Print a formatted E[Y|X=x] / mu(x) / B(x) summary table, reproducing the
#' layout of Tables 3-4 of the dissertation.
#'
#' @param res         Output list from compute_causal_quantities()
#' @param xgrid       Treatment grid (same length as nrow(res$mu))
#' @param idx         Integer vector of exact grid indices to tabulate
#' @param case_names  Character vector of case labels (length = ncol(res$mu))
#' @param param_label Character vector of per-case parameter annotations
print_causal_table <- function(res, xgrid, idx, case_names, param_label) {
  header  <- sprintf("%-6s  %-13s  %-9s  %-9s",
                     "x", "E[Y|X=x]", "mu(x)", "B(x)")
  divider <- strrep("-", 42)
  
  for (j in seq_along(case_names)) {
    cat(sprintf("\n%s (%s):\n", case_names[j], param_label[j]))
    cat(header, "\n", divider, "\n", sep = "")
    
    for (i in idx) {
      cat(sprintf("%-6.2f  %-13.4f  %-9.4f  %-9.4f\n",
                  xgrid[i],
                  res$eobs[i, j],
                  res$mu[i,   j],
                  res$bias[i, j]))
    }
  }
  invisible(NULL)
}


# ------------------------------------------------------------------------------
# 1E  Generic utility: draw the two-panel figure
# ------------------------------------------------------------------------------

#' Draw the standard two-panel causal summary figure used for Figures 12-13.
#'
#' Left panel:  E[Y|X=x] (dashed, per case) overlaid with mu(x) (solid, per
#'              case) -- directly visualises the confounding gap that B(x)
#'              quantifies in the right panel.
#' Right panel: B(x) = E[Y|X=x] - mu(x), the analytic confounding bias.
plot_two_panels <- function(xgrid, res, labs, cols, ltys, title_tag = "") {
  par(mfrow = c(1, 2), mar = c(4.5, 4.5, 3, 1))
  n_cases <- ncol(res$mu)
  
  # --- Left panel: E[Y|X=x] (dashed) and mu(x) (solid) -----------------------
  
  ylim_left <- range(c(res$eobs, res$mu), finite = TRUE)
  ylim_left <- c(max(0, ylim_left[1] - 0.02),
                 min(1, ylim_left[2] + 0.02))
  
  plot(xgrid, res$mu[, 1], type = "l", lwd = 2, col = cols[1],
       ylim = ylim_left,
       xlab = expression(x), ylab = "",
       main = bquote("" ~ mu(x) ~
                       "& E[Y|X=x]" ~ .(title_tag)))
  
  lines(xgrid, res$eobs[, 1], lwd = 2, col = cols[1], lty = 3)
  
  for (j in 2:n_cases) {
    lines(xgrid, res$mu[, j],   lwd = 2, col = cols[j], lty = ltys[j])
    lines(xgrid, res$eobs[, j], lwd = 2, col = cols[j], lty = 3)
  }
  
  legend("topleft", legend = labs, col = cols, lwd = 2, lty = ltys, bty = "n")
  
  legend("bottomright",
         legend = c(expression(mu(x) ~ "(interventional, solid/dashed per case)"),
                    "E[Y|X=x] (observational, dotted)"),
         lty = c(1, 3), lwd = 2, col = c("black", "black"),
         bty = "n", cex = 0.75)
  
  
  # --- Right panel: Confounding bias B(x) -----------------------------------
  
  ylim_right <- range(res$bias, finite = TRUE)
  pad <- 0.05 * diff(ylim_right)
  if (pad == 0) pad <- 0.01
  ylim_right <- ylim_right + c(-pad, pad)
  
  plot(xgrid, res$bias[, 1], type = "l", lwd = 2, col = cols[1],
       ylim = ylim_right,
       xlab = expression(x), ylab = expression(B(x)),
       main = bquote("" ~ B(x) == "E[Y|X=x] -" ~ mu(x) ~
                       .(title_tag)))
  
  for (j in 2:n_cases)
    lines(xgrid, res$bias[, j], lwd = 2, col = cols[j], lty = ltys[j])
  
  abline(h = 0, lty = 3, col = "gray50")
  
  legend("topright", legend = labs, col = cols, lwd = 2, lty = ltys, bty = "n")
}


# ==============================================================================
# SECTION 2 - PARAMETERS
# ==============================================================================

# Three confounding regimes are studied for each copula family, matching the
# regimes reported in Tables 3-4 of the dissertation (Section 5.6).

# ---- 2A  Gaussian copula: observable pairwise correlations ------------------

# All three cases share rho_XY = 0.70, rho_YZ = 0.60; only rho_XZ varies.

gauss_params <- list(
  Moderate = c(rho_xy = 0.70, rho_xz = 0.50, rho_yz = 0.60),
  None     = c(rho_xy = 0.70, rho_xz = 0.00, rho_yz = 0.60),
  Strong   = c(rho_xy = 0.70, rho_xz = 0.80, rho_yz = 0.60)
)


# ---- 2B  FGM copula: dependence parameters ---------------------------------

# Three cases vary alpha = theta_{X,Y|Z} while holding
# beta = theta_{X,Z} = 0.50 and eta = theta_{Y,Z} = 0.30 fixed.

fgm_params <- list(
  CaseA = c(theta_xy_cond = 0.40, theta_xz = 0.50, theta_yz = 0.30),
  CaseB = c(theta_xy_cond = 0.80, theta_xz = 0.50, theta_yz = 0.30),
  CaseC = c(theta_xy_cond = 0.20, theta_xz = 0.50, theta_yz = 0.30)
)


# ==============================================================================
# SECTION 3 - COMPUTATIONS
# ==============================================================================

# ---- 3A  Gaussian copula ----------------------------------------------------

gauss_coeffs <- lapply(gauss_params, function(p) {
  gaussian_K_and_V(p["rho_xy"], p["rho_xz"], p["rho_yz"])
})

K_vals  <- sapply(gauss_coeffs, function(cc) cc["K"])
K0_vals <- sapply(gauss_params, function(p) gaussian_K0(p["rho_xy"]))

names(K_vals)  <- names(gauss_params)
names(K0_vals) <- names(gauss_params)

cat("Gaussian copula: interventional-mean slope K = a / sqrt(1+V)\n")
print(round(K_vals, 4))

cat("Gaussian copula: no-confounding / observational slope K_0\n")
print(round(K0_vals, 4))

# Each case needs both K and K0 supplied to eobs_fn / mu_fn / bias_fn.
gauss_case_params <- Map(
  function(K, K0) c(K = unname(K), K0 = unname(K0)),
  K_vals,
  K0_vals
)

cat("\nComputing Gaussian causal quantities...\n")

res_gauss <- compute_causal_quantities(
  eobs_fn = function(x, p) eobs_gaussian(x, p["K0"]),
  mu_fn   = function(x, p) mu_gaussian(x, p["K"]),
  bias_fn = function(x, p) bias_gaussian(x, p["K"], p["K0"]),
  params  = gauss_case_params,
  xgrid   = XGRID
)

cat("Done.\n")

# Internal consistency check: bias_gaussian() must equal eobs - mu exactly.
stopifnot(
  max(abs(res_gauss$bias - (res_gauss$eobs - res_gauss$mu))) < 1e-10
)


# ---- 3B  FGM copula ---------------------------------------------------------

cat("\nFGM copula: alpha = theta_{X,Y|Z} values\n")

alpha_vals <- sapply(fgm_params, function(p) p["theta_xy_cond"])
names(alpha_vals) <- names(fgm_params)

print(round(alpha_vals, 4))

cat("Computing FGM causal quantities...\n")

res_fgm <- compute_causal_quantities(
  eobs_fn = eobs_fgm,
  mu_fn   = mu_fgm,
  bias_fn = bias_fgm,
  params  = fgm_params,
  xgrid   = XGRID
)

cat("Done.\n")

# Internal consistency check: bias_fgm() must equal eobs_fgm - mu_fgm exactly.
stopifnot(
  max(abs(res_fgm$bias - (res_fgm$eobs - res_fgm$mu))) < 1e-10
)

# Check the scale of the FGM bias used in the figure.
cat("\nFGM bias range over plotting grid:\n")
print(range(res_fgm$bias, finite = TRUE))

# For the chosen parameter values, the bias remains below 0.01 in absolute
# value over the plotting grid.
stopifnot(max(abs(res_fgm$bias)) < 0.01)


# ==============================================================================
# SECTION 4 - CONSOLE OUTPUT (Tables 3-4)
# ==============================================================================

CASE_NAMES_GAUSS <- c("Moderate", "None", "Strong")
CASE_NAMES_FGM   <- c("Case A", "Case B", "Case C")

cat("\n\n=== TABLE 3: Gaussian Copula -- E[Y|X=x], mu(x), B(x) ===\n")
cat("(rho_XY = 0.70, rho_YZ = 0.60 throughout; matches tab:gaussian)\n")

gauss_labels <- sprintf(
  "rho_XZ = %.2f, K = %.4f",
  sapply(gauss_params, function(p) p["rho_xz"]),
  K_vals
)

print_causal_table(
  res         = res_gauss,
  xgrid       = XGRID,
  idx         = TABLE_IDX,
  case_names  = CASE_NAMES_GAUSS,
  param_label = gauss_labels
)


cat("\n\n=== TABLE 4: FGM Copula -- E[Y|X=x], mu(x), B(x) ===\n")
cat("(beta = 0.50, eta = 0.30 throughout; matches tab:fgm)\n")

fgm_labels <- sprintf(
  "alpha = %.2f",
  alpha_vals
)

print_causal_table(
  res         = res_fgm,
  xgrid       = XGRID,
  idx         = TABLE_IDX,
  case_names  = CASE_NAMES_FGM,
  param_label = fgm_labels
)


# ==============================================================================
# SECTION 5 - FIGURES
# ==============================================================================

# Assemble legend labels linking each curve to its parameter value.
labs_gauss <- sprintf("%s (K = %.4f)", CASE_NAMES_GAUSS, K_vals)
labs_fgm   <- sprintf("%s (alpha = %.2f)", CASE_NAMES_FGM, alpha_vals)

# Create output directory if it does not already exist.
if (!dir.exists("imagens")) {
  dir.create("imagens", recursive = TRUE)
}


# ---- Figure 12: Gaussian copula ---------------------------------------------

jpeg(
  filename = "imagens/ACE_paper_example_Gauss.jpg",
  width = 1800,
  height = 900,
  res = 220,
  quality = 95
)

plot_two_panels(
  xgrid     = XGRID,
  res       = res_gauss,
  labs      = labs_gauss,
  cols      = COLS,
  ltys      = LTYS,
  title_tag = "(Gaussian)"
)

dev.off()


# ---- Figure 13: FGM copula ---------------------------------------------------

jpeg(
  filename = "imagens/ACE_paper_example_FGM.jpg",
  width = 1800,
  height = 900,
  res = 220,
  quality = 95
)

plot_two_panels(
  xgrid     = XGRID,
  res       = res_fgm,
  labs      = labs_fgm,
  cols      = COLS,
  ltys      = LTYS,
  title_tag = "(FGM)"
)

dev.off()


# ==============================================================================
# SECTION 6 - SUMMARY
# ==============================================================================

cat("\n", strrep("=", 79), "\n", sep = "")
cat("SUMMARY: E[Y|X=x], mu(x), AND THE ANALYTIC CONFOUNDING BIAS B(x)\n")
cat(strrep("=", 79), "\n\n")

cat("Both families reproduce the qualitative contrast of the numerical examples:\n\n")

cat("1. GAUSSIAN COPULA:\n")
cat("   - mu(x) and E[Y|X=x] are both S-shaped (Phi(K*Phi^{-1}(x)) family).\n")
cat("   - K is non-monotonic in |rho_XZ| through the structural coefficients.\n")
cat("   - B(x) = E[Y|X=x] - mu(x), vanishes at x = 0.5, and changes sign\n")
cat("     across the median.\n\n")

cat("2. FGM COPULA:\n")
cat("   - mu(x) is exactly affine in x (Eq. 5.23); beta cancels out of the\n")
cat("     interventional mean by antisymmetry.\n")
cat("   - beta remains in the full interventional distribution through\n")
cat("     F_{X|Z}(x|z), even though its contribution to mu(x) vanishes.\n")
cat("   - B(x) vanishes at x = 0.5 and is antisymmetric about the median.\n")
cat("   - B(x) is generally cubic in x through the term\n")
cat("     x*(1-x)*(1-2*x) in Eq. 5.25.\n\n")

cat("3. VERIFICATION:\n")
cat("   - bias_gaussian() equals E[Y|X=x] - mu(x) to machine precision.\n")
cat("   - bias_fgm() equals E[Y|X=x] - mu(x) to machine precision.\n")
cat("   - Table values are evaluated at the exact stated treatment values,\n")
cat("     including x = 0.50.\n")
cat("   - The FGM bias remains below 0.01 in absolute value over the plotting\n")
cat("     grid for the parameter values used here.\n\n")

cat("This script implements the closed-form expressions used for the numerical\n")
cat("examples and generates the corresponding tables and figures.\n\n")

cat("Figures saved as:\n")
cat("  imagens/ACE_paper_example_Gauss.jpg\n")
cat("  imagens/ACE_paper_example_FGM.jpg\n\n")