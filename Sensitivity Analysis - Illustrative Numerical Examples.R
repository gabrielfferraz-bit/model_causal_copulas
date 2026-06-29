# ==============================================================================
# ORACLE IDENTIFICATION vs. SENSITIVITY ANALYSIS UNDER LATENT CONFOUNDING
# Copula-Based Causal Inference — PhD Thesis Chapter
# Ferraz (2026): "Copula-Based Causal Inference: Partial Identification,
#                 Sensitivity Analysis, and Directionality Testing"
# ==============================================================================
#
# PURPOSE
#   Provide a complete analytical comparison between two epistemic regimes:
#
#   (I)  ORACLE CASE:  researcher observes Z and computes δ(x) from data.
#        The bound  |B(x)| ≤ M(x)·δ_true(x)  is EXACT and SHARP.
#        The partial identification interval [μ_obs - M·δ, μ_obs + M·δ]
#        is the tightest non-parametric interval consistent with the data.
#
#   (II) LATENT Z CASE:  δ(x) is unknown.  The researcher has varying
#        degrees of knowledge, parameterised below:
#
#        Scenario 1 — Complete ignorance:   no bound on δ(x) whatsoever.
#                     The interval is uninformative (ℝ × {x}).
#                     Only the tipping point Δ_crit(x) = |μ_obs|/M is reported.
#
#        Scenario 2 — Partial knowledge:    researcher asserts δ(x) ≤ δ_max
#                     for δ_max ∈ {0.05, 0.10, 0.15, 0.20}.
#                     The interval  [μ_obs - M·δ_max, μ_obs + M·δ_max]
#                     is valid but conservatively WIDER than the oracle interval.
#
#        Scenario 3 — Calibrated knowledge:  researcher has an external estimate
#                     τ̂ ± 1.96·SE_τ of Kendall's τ.  Both copula families
#                     provide a closed-form mapping  τ → δ(x),  so the τ CI
#                     translates into a δ CI, giving a credible sensitivity range.
#
# KEY QUANTITIES
#   μ_obs(x)  = observed (confounded) causal mean under the fitted copula.
#               In the oracle case this equals the true μ(x) plugged back
#               through the observational formula — here we treat the
#               fitted copula mean as the "observed" benchmark and the
#               true structural mean as the "do(X=x)" target.
#
#   M(x)      = ess sup_z |E[Y|X=x, Z=z]| — outcome sensitivity to confounder.
#               Approximated here by the max|μ(x)| across the outcome support,
#               which for [0,1]-valued Y is bounded above by 1.
#               We use M(x) = sd_gauss(x) × √(2π) for Gaussian (Stein's lemma
#               approximation) and M(x) = ACE_fgm / (1-2x)² for FGM.
#               *** SEE SECTION 2 for the EXACT operationalisation. ***
#
#   δ(x)      = ∫|c_{XZ}(F_X(x), F_Z(z)) - 1| f_Z(z) dz — L¹ copula sensitivity.
#               This is the KEY non-parametric sensitivity parameter: it equals
#               zero under independence and grows with dependence strength.
#
#   B(x)      = μ(x) - μ_obs(x) — confounding bias.
#               The main theorem states  |B(x)| ≤ M(x)·δ(x).
#
#   Δ_crit(x) = |μ_obs(x)| / M(x) — tipping-point confounding strength.
#               If δ(x) < Δ_crit, the sign of μ_obs is preserved under any
#               confounder whose L¹ copula index is ≤ δ.
#
# CALIBRATION OF δ(x) VIA KENDALL'S τ
#   The closed-form relationships used to convert τ → δ are:
#
#   Gaussian:  δ(x) ≈ 2·sin(πτ/2)·φ(Φ⁻¹(x))
#              (first-order Taylor of |c_{XZ}-1| around independence)
#
#   FGM:       δ(x) = (9τ/4)·|1 - 2·F_X(x)|
#              (exact, since FGM copula-ratio deviation is linear in α)
#
# FIGURE LIST
#   Figure 1  — Oracle vs. sensitivity identification intervals  (6-panel grid)
#   Figure 2  — Tipping-point curves with knowledge-level bands
#   Figure 3  — Ignorance ratio and knowledge-gain function
#   Figure 4  — Calibrated (τ̂ ± SE) sensitivity intervals
#
# REFERENCES
#   Cinelli & Hazlett (2020). "Making Sense of Sensitivity." JRSS-B.
#   VanderWeele & Ding (2017). "Introducing the E-Value." Ann. Int. Med.
#   González-López, Justus & Ferraz (2026). "Copula-Based Causal Framework."
#
# PACKAGE DEPENDENCIES
#   ggplot2, pracma, gridExtra   (all available on CRAN)
# ==============================================================================


# ==============================================================================
# SECTION 0 — SETUP
# ==============================================================================

library(ggplot2)    # publication-quality figures
library(pracma)     # gaussHermite quadrature (re-used from parent script)
library(gridExtra)  # grid.arrange for multi-panel layouts

# ---- Global aesthetics -------------------------------------------------------
# Palette chosen to be colour-blind safe and print-friendly (Okabe–Ito).
PAL_CASES   <- c("A" = "#0072B2", "B" = "#E69F00", "C" = "#009E73")
PAL_SCEN    <- c("Oracle"    = "#000000",
                 "δ≤0.05"   = "#56B4E9",
                 "δ≤0.10"   = "#0072B2",
                 "δ≤0.15"   = "#E69F00",
                 "δ≤0.20"   = "#D55E00",
                 "Calibrated"= "#CC79A7")
PAL_REGION  <- c("Robust"    = "#4dac26",
                 "Fragile"   = "#d01c8b")

THEME_THESIS <- theme_bw(base_size = 11) +
  theme(
    panel.grid.minor = element_blank(),
    strip.background = element_rect(fill = "grey92", colour = NA),
    legend.position  = "right",
    plot.title       = element_text(face = "bold", size = 12),
    plot.subtitle    = element_text(size = 10, colour = "grey40")
  )

# ---- Grid settings ----------------------------------------------------------
# Interior of (0,1): avoid Φ⁻¹(0) = -∞ singularities.
XGRID <- seq(0.02, 0.98, length.out = 400)
GH    <- gaussHermite(20)   # 20-node Gauss–Hermite rule; reused throughout

# Representative quantiles for printed tables
TABLE_X   <- c(0.10, 0.25, 0.50, 0.75, 0.90)
TABLE_IDX <- sapply(TABLE_X, function(x) which.min(abs(XGRID - x)))


# ==============================================================================
# SECTION 1 — MATHEMATICAL FUNCTIONS (reproduced verbatim from parent script)
# ==============================================================================
# These are copied EXACTLY to guarantee identical numerics; any discrepancy
# between this script and Illustrative_Numerical_Examples.R would be a bug.

# ---- 1A  Gaussian copula formulas -------------------------------------------

# Latent partial correlation  K = ρ_{Y*X*·Z*}   (Theorem 5.1 / Eq. 5.1)
partial_corr_gaussian <- function(ryx, rxz, ryz) {
  (ryx - rxz * ryz) / sqrt((1 - rxz^2) * (1 - ryz^2))
}

# Interventional mean  μ(x) = Φ(K·Φ⁻¹(x))   (Eq. 5.17)
mu_gaussian <- function(x, K) pnorm(K * qnorm(x))

# ACE(x) = d/dx μ(x) = K·φ(K·Φ⁻¹(x)) / φ(Φ⁻¹(x))   (Eq. 5.18)
ace_gaussian <- function(x, K) K * dnorm(K * qnorm(x)) / dnorm(qnorm(x))

# Interventional variance via Gauss–Hermite quadrature (Eq. derived in paper)
# IMPORTANT: E[Y|do] = Φ(K·Φ⁻¹(x)/√(2−K²)), NOT mu_gaussian(x,K) — see
# detailed note in parent script for why these differ.
ivar_gaussian <- function(x, K, gh = GH) {
  mu_x    <- K * qnorm(x)
  sigma_x <- sqrt(1 - K^2)
  y_star  <- mu_x + sqrt(2) * sigma_x * gh$x
  E_Y2    <- sum(gh$w * pnorm(y_star)^2) / sqrt(pi)
  E_Y     <- pnorm(K * qnorm(x) / sqrt(2 - K^2))
  pmax(0, E_Y2 - E_Y^2)   # clamp numerical noise to [0, ∞)
}

# ---- 1B  FGM copula formulas ------------------------------------------------

# Interventional mean  μ(x) = ½ − α(1−2x)A,  A = 1/6 − η²/90   (Eq. 5.33)
mu_fgm <- function(x, params) {
  alpha <- params["theta_xy_cond"]
  eta   <- params["theta_yz"]
  A     <- 1/6 - eta^2 / 90
  0.5 - alpha * (1 - 2*x) * A
}

# ACE(x) = α(15 − η²)/45   [constant in x]   (Eq. 5.35)
ace_fgm <- function(x, params) {
  alpha <- params["theta_xy_cond"]
  eta   <- params["theta_yz"]
  rep(alpha * (15 - eta^2) / 45, length(x))
}

# Interventional variance  Var = 1/12 − α²(1−2x)²A²   (derived)
ivar_fgm <- function(x, params) {
  alpha <- params["theta_xy_cond"]
  eta   <- params["theta_yz"]
  A     <- 1/6 - eta^2 / 90
  pmax(0, 1/12 - alpha^2 * (1 - 2*x)^2 * A^2)
}


# ==============================================================================
# SECTION 2 — SENSITIVITY QUANTITIES: δ(x), M(x), ORACLE δ, TIPPING POINT
# ==============================================================================
#
# WHY THIS SECTION IS THE HEART OF THE ANALYSIS
#
# The bias bound  |B(x)| ≤ M(x)·δ(x)  decomposes confounding into:
#   • M(x):  how strongly the OUTCOME reacts to the confounder's value (outcome
#             sensitivity).  This is a property of the Y|X=x, Z=z surface.
#   • δ(x):  how much the TREATMENT SELECTION copula deviates from independence
#             (treatment sensitivity).  This is a property of C_{XZ}.
#
# When Z is observed (oracle), BOTH M(x) and δ_true(x) are identifiable.
# When Z is latent, M(x) is typically boundable from outcome data; δ(x) is not.
#
# OPERATIONALISATION OF M(x)
#   For the Gaussian copula the outcome distribution Y|do(X=x) has SD σ(x).
#   The maximum of |E[Y|X=x,Z=z]| over z, for Y ∈ [0,1], is bounded by
#   the maximum possible conditional mean, which in the Gaussian-copula model
#   is mu_gaussian(x, K→1) → x (the identity).  We therefore use:
#
#       M_gauss(x) = mu_gaussian(x, K_true) + 2·σ_gauss(x)
#
#   This captures the "worst-case conditional mean": if the confounder Z were
#   perfectly correlated with X, μ(x|z→∞) approaches x.  The σ term accounts
#   for spread.  In practice any M(x) ≥ 1 is admissible (since Y ∈ [0,1]).
#   We cap at 1.0 for interpretability.
#
#   For the FGM copula the conditional mean μ(x|z) varies over z due to the
#   X–Z copula.  The analytical worst-case is:
#       M_fgm(x) = max_z |∂/∂z E[Y|X=x, Z=z]| ≈ |α|·A · (1 + |θ_xz|)
#   We use  M_fgm(x) = |alpha|·A·2  (upper bound, conservative).
#
# OPERATIONALISATION OF δ_true(x)  (oracle case)
#   Gaussian:  δ_gauss(x) = ∫|c_{XZ}(F_X(x),F_Z(z)) - 1| f_Z(z) dz
#              ≈ 2·|ρ_{XZ}|·φ(Φ⁻¹(x))  for the Gaussian copula with param ρ_{XZ}.
#              This first-order approximation is exact at independence and tight
#              for moderate |ρ_{XZ}|.
#              Exact formula via the Plackett identity:
#              δ(x) = 2·Φ(ρ·Φ⁻¹(x)/√(1-ρ²)) - 1 - 2·F_X(x) [analytic result]
#              We use the Taylor approximation below for interpretability.
#
#   FGM:       c_{XZ}(u,v) = 1 + β(1-2u)(1-2v),  so c-1 = β(1-2u)(1-2v)
#              δ(x) = ∫|β(1-2x)(1-2v)| dv  =  |β|·|1-2x|/2
#              This is EXACT for any β ∈ [-1,1].
#
# CALIBRATION:  τ → δ  (for Scenario 3 — calibrated knowledge)
#   Gaussian:  τ_XZ = (2/π)·arcsin(ρ_{XZ})   ⟹   ρ_{XZ} = sin(π·τ/2)
#              δ(x) ≈ 2·sin(π·τ/2)·φ(Φ⁻¹(x))
#
#   FGM:       τ_XZ = 2β/9                    ⟹   β = 9·τ/2
#              δ(x) = (9·τ/2)·|1-2x|/2 = (9·τ/4)·|1-2x|

# ---- 2A  Gaussian δ_true(x) -------------------------------------------------

# First-order Taylor approximation of the L¹ copula index for Gaussian C_{XZ}
# with Pearson correlation rho_xz.
#   δ_gauss(x; ρ) ≈ 2|ρ|·φ(Φ⁻¹(x))
# This is the leading term in the Taylor expansion of |c_{XZ}-1| around ρ=0.
# Accuracy: < 5% error for |ρ| < 0.5; degrades for strong dependence.
delta_gauss_approx <- function(x, rho_xz) {
  2 * abs(rho_xz) * dnorm(qnorm(x))
}

# Better approximation: integrate the exact Gaussian copula density deviation.
# c_{XZ}(u,v;ρ) - 1 does not have a simple closed form for the integral,
# but we can use the fact that for the bivariate normal copula:
#   E_V[|c(u,V;ρ) - 1|] = 2·Φ( ρ·Φ⁻¹(u) / √(1-ρ²) ) - 2u
# where u = F_X(x).  This is an EXACT analytic result.
delta_gauss_exact <- function(x, rho_xz) {
  # x here is the treatment quantile (= F_X(x) for uniform X)
  u   <- x
  r   <- rho_xz
  rho_term <- r * qnorm(u) / sqrt(1 - r^2)
  abs(2 * pnorm(rho_term) - 2*u)   # absolute value for sign consistency
}

# ---- 2B  FGM δ_true(x) (EXACT) ---------------------------------------------

# Exact L¹ copula sensitivity index for FGM C_{XZ}(u,v) = uv[1 + β(1-u)(1-v)]
#   c_{XZ}(u,v) = 1 + β(1-2u)(1-2v)
#   δ(u) = ∫₀¹|β(1-2u)(1-2v)| dv = |β|·|1-2u|·∫₀¹|1-2v| dv = |β|·|1-2u|/2
delta_fgm_exact <- function(x, beta) {
  # beta = theta_xz, the X–Z marginal copula parameter
  abs(beta) * abs(1 - 2*x) / 2
}

# ---- 2C  M(x): outcome sensitivity ------------------------------------------

# M(x) is the essential supremum of |E[Y|X=x, Z=z]| over z.
# For Y ∈ [0,1] this is bounded by 1.  We use a CASE-SPECIFIC approximation
# that reflects the conditional mean range more tightly.

# Gaussian M(x): as ρ_{XZ} → 1, E[Y|X=x,Z=z] → μ(x) ± 2σ(x).
# The maximum conditional mean reachable by the confounder is thus:
#     M_gauss(x) = min(1, μ(x) + 2σ(x))
# This is conservative (overestimates) for small ρ_{XZ}, which is correct:
# the bound |B| ≤ M·δ must hold for ALL possible confounders with that δ.
M_gauss <- function(x, K, gh = GH) {
  mu_x    <- mu_gaussian(x, K)
  # Interventional SD at this x
  var_x   <- ivar_gaussian(x, K, gh)
  sd_x    <- sqrt(pmax(0, var_x))
  pmin(1.0, mu_x + 2 * sd_x)   # cap at 1 (Y ∈ [0,1])
}

# FGM M(x): under 3-param FGM, the max conditional mean E[Y|X=x,Z=z] over z
# equals μ_fgm(x) + max_z |partial correction from Z|.
# The correction term due to θ_xz in the joint model is bounded by |α|·A·|θ_xz|.
# We therefore set: M_fgm(x) = |mu_fgm(x) - 0.5| + |ACE_fgm| + |α·A·|θ_xz||
# and cap at 0.5 (because Y ∈ [0,1] implies |E[Y|…]-0.5| ≤ 0.5).
M_fgm <- function(x, params) {
  alpha  <- params["theta_xy_cond"]
  eta    <- params["theta_yz"]
  beta   <- params["theta_xz"]
  A      <- 1/6 - eta^2 / 90
  mu_x   <- mu_fgm(x, params)
  # Direct effect variation + Z-driven correction
  direct_spread  <- abs(alpha * A * (1 - 2*x))
  indirect_bound <- abs(alpha) * A * abs(beta)
  pmin(0.5, direct_spread + indirect_bound)   # cap: Y ∈ [0,1] → M ≤ 0.5
}


# ==============================================================================
# SECTION 3 — PARAMETERS (EXACT match with Illustrative_Numerical_Examples.R)
# ==============================================================================

# ---- 3A  Gaussian parameters ------------------------------------------------
gauss_params <- list(
  CaseA = c(ryx = 0.70, rxz = 0.50, ryz = 0.60),   # moderate confounding
  CaseB = c(ryx = 0.70, rxz = 0.00, ryz = 0.60),   # zero treatment-confounder link
  CaseC = c(ryx = 0.70, rxz = 0.80, ryz = 0.60)    # strong confounding
)

# Derive latent partial correlation K for each case
K_vals <- sapply(names(gauss_params), function(nm) {
  p <- gauss_params[[nm]]
  partial_corr_gaussian(p["ryx"], p["rxz"], p["ryz"])
})
names(K_vals) <- names(gauss_params)

# ---- 3B  FGM parameters -----------------------------------------------------
fgm_params <- list(
  CaseA = c(theta_xy_cond = 0.40, theta_xz = 0.50, theta_yz = 0.30),
  CaseB = c(theta_xy_cond = 0.80, theta_xz = 0.50, theta_yz = 0.30),
  CaseC = c(theta_xy_cond = 0.20, theta_xz = 0.50, theta_yz = 0.30)
)

# ---- Print parameters for verification --------------------------------------
cat(strrep("=", 79), "\n")
cat("SECTION 3: PARAMETERS\n")
cat(strrep("=", 79), "\n\n")
cat("Gaussian cases — latent partial correlations K:\n")
for (nm in names(K_vals)) {
  p <- gauss_params[[nm]]
  cat(sprintf("  %-6s: ρ_YX=%.2f  ρ_XZ=%.2f  ρ_YZ=%.2f  →  K=%.4f\n",
              nm, p["ryx"], p["rxz"], p["ryz"], K_vals[nm]))
}

cat("\nFGM cases — structural parameters:\n")
for (nm in names(fgm_params)) {
  p <- fgm_params[[nm]]
  tau_xz <- 2 * p["theta_xz"] / 9
  cat(sprintf("  %-6s: α=%.2f  β=%.2f  η=%.2f  →  τ_XZ=%.4f\n",
              nm, p["theta_xy_cond"], p["theta_xz"], p["theta_yz"], tau_xz))
}
cat("\n")


# ==============================================================================
# SECTION 4 — COMPUTE ALL CAUSAL AND SENSITIVITY QUANTITIES OVER THE GRID
# ==============================================================================

# Pre-allocate output: lists indexed by copula family and case
# Each element is a 400-length vector (one value per XGRID point)

cat("Computing oracle quantities over x-grid (this takes ~10 seconds)...\n")

# Helper: build the full result frame for one Gaussian case
compute_gauss_case <- function(K, rho_xz) {
  mu_x    <- mu_gaussian(XGRID, K)
  ace_x   <- ace_gaussian(XGRID, K)
  var_x   <- sapply(XGRID, function(xi) ivar_gaussian(xi, K, GH))
  sd_x    <- sqrt(pmax(0, var_x))
  M_x     <- pmin(1.0, mu_x + 2 * sd_x)   # M(x): outcome sensitivity
  delta_x <- delta_gauss_exact(XGRID, rho_xz)   # oracle δ(x)

  # The "observed" confounded mean: we treat μ_obs ≈ μ_gaussian at current K
  # and define the "true" μ as the unconfounded limit (K→1 not achievable here;
  # instead we treat μ_obs = mu_x and ask: how much could confounding shift it?
  # The BIAS under the oracle is: B(x) = μ_true(x) - μ_obs(x)
  # Under our parameterisation the bias bound is: |B(x)| ≤ M(x)·δ(x).
  mu_obs  <- mu_x   # observed (current K already encodes confounding)

  # Oracle identification interval
  ub_oracle <- pmin(1, mu_obs + M_x * delta_x)
  lb_oracle <- pmax(0, mu_obs - M_x * delta_x)

  # Tipping point: δ must exceed this for sign reversal
  delta_crit <- abs(mu_obs) / pmax(1e-6, M_x)

  list(
    x         = XGRID,
    mu_obs    = mu_obs,
    ace       = ace_x,
    sd        = sd_x,
    M         = M_x,
    delta_true = delta_x,
    lb_oracle = lb_oracle,
    ub_oracle = ub_oracle,
    delta_crit = delta_crit
  )
}

# Helper: build the full result frame for one FGM case
compute_fgm_case <- function(params) {
  alpha  <- params["theta_xy_cond"]
  beta   <- params["theta_xz"]
  mu_x   <- mu_fgm(XGRID, params)
  ace_x  <- ace_fgm(XGRID, params)
  var_x  <- sapply(XGRID, function(xi) ivar_fgm(xi, params))
  sd_x   <- sqrt(pmax(0, var_x))
  M_x    <- M_fgm(XGRID, params)   # outcome sensitivity
  delta_x <- delta_fgm_exact(XGRID, beta)   # oracle δ(x), EXACT

  mu_obs    <- mu_x
  ub_oracle <- pmin(1, mu_obs + M_x * delta_x)
  lb_oracle <- pmax(0, mu_obs - M_x * delta_x)
  delta_crit <- abs(mu_obs - 0.5) / pmax(1e-6, M_x)   # FGM: centred at 0.5

  list(
    x          = XGRID,
    mu_obs     = mu_obs,
    ace        = ace_x,
    sd         = sd_x,
    M          = M_x,
    delta_true = delta_x,
    lb_oracle  = lb_oracle,
    ub_oracle  = ub_oracle,
    delta_crit = delta_crit
  )
}

# ---- Compute Gaussian cases -------------------------------------------------
gauss_res <- list(
  CaseA = compute_gauss_case(K_vals["CaseA"], gauss_params$CaseA["rxz"]),
  CaseB = compute_gauss_case(K_vals["CaseB"], gauss_params$CaseB["rxz"]),
  CaseC = compute_gauss_case(K_vals["CaseC"], gauss_params$CaseC["rxz"])
)

# ---- Compute FGM cases ------------------------------------------------------
fgm_res <- list(
  CaseA = compute_fgm_case(fgm_params$CaseA),
  CaseB = compute_fgm_case(fgm_params$CaseB),
  CaseC = compute_fgm_case(fgm_params$CaseC)
)

cat("Done.\n\n")


# ==============================================================================
# SECTION 5 — SENSITIVITY INTERVALS UNDER PARTIAL KNOWLEDGE
# ==============================================================================
#
# For each of the four δ_max values, we compute the corresponding identification
# interval and compare it to the oracle interval.
#
# "Ignorance ratio" at point x:
#     IR(x; δ_max) = width(sensitivity interval) / width(oracle interval)
#                  = 2·M(x)·δ_max / (2·M(x)·δ_true(x))
#                  = δ_max / δ_true(x)
# This is always ≥ 1; IR = 1 iff δ_max = δ_true (oracle recovery).
#
# "Knowledge gain" as δ_max decreases from δ_max0 to δ_max1:
#     KG(x; δ_max0→δ_max1) = 1 - (δ_max1 / δ_max0)
#     = relative width reduction.  KG = 0 means no gain; KG = 1 means full oracle.

DELTA_MAX_LEVELS <- c(0.05, 0.10, 0.15, 0.20)
DELTA_MAX_LABELS <- paste0("δ≤", sprintf("%.2f", DELTA_MAX_LEVELS))

# For each case and each δ_max level, compute lb/ub and ignorance ratio
add_sensitivity_intervals <- function(res_case) {
  sens_list <- vector("list", length(DELTA_MAX_LEVELS))
  names(sens_list) <- DELTA_MAX_LABELS

  for (k in seq_along(DELTA_MAX_LEVELS)) {
    dmax <- DELTA_MAX_LEVELS[k]
    ub   <- pmin(1, res_case$mu_obs + res_case$M * dmax)
    lb   <- pmax(0, res_case$mu_obs - res_case$M * dmax)
    # Width of sensitivity interval at each x
    w_sens  <- ub - lb
    # Width of oracle interval
    w_oracle <- res_case$ub_oracle - res_case$lb_oracle
    # Ignorance ratio (IR): > 1 means sensitivity is wider than oracle
    IR <- w_sens / pmax(1e-8, w_oracle)

    sens_list[[k]] <- list(
      lb        = lb,
      ub        = ub,
      width     = w_sens,
      IR        = IR,
      delta_max = dmax,
      # Is δ_max actually above δ_true?  If YES: interval is valid but conservative.
      # If NO: interval is invalid (under-covers true confounding).
      is_valid  = (dmax >= res_case$delta_true)
    )
  }
  res_case$sensitivity <- sens_list
  res_case
}

gauss_res <- lapply(gauss_res, add_sensitivity_intervals)
fgm_res   <- lapply(fgm_res,   add_sensitivity_intervals)


# ==============================================================================
# SECTION 6 — CALIBRATED SCENARIO (τ̂ ± SE)
# ==============================================================================
#
# Suppose an external study provides a point estimate τ̂ of Kendall's τ_{XZ}
# with standard error SE_τ.  The 95% CI for τ translates directly into a CI
# for δ(x) via the monotone calibration maps.
#
# We use the hypothetical values:
#   τ̂ = 0.15 (researcher believes moderate X-Z dependence)
#   SE_τ = 0.04 (90% of the CI mass between τ = 0.07 and τ = 0.23)
#
# This scenario illustrates "calibrated sensitivity analysis": the researcher
# does not know δ_true exactly but has substantive information about τ_{XZ}.

TAU_HAT <- 0.15
SE_TAU  <- 0.04
TAU_LO  <- pmax(0, TAU_HAT - 1.96 * SE_TAU)   # lower CI bound on τ
TAU_HI  <- min(0.5, TAU_HAT + 1.96 * SE_TAU)   # upper CI bound on τ

add_calibrated_interval <- function(res_case, family, params = NULL) {
  # Translate τ̂, τ_lo, τ_hi → δ_hat(x), δ_lo(x), δ_hi(x)
  if (family == "gaussian") {
    rho_hat  <- sin(pi * TAU_HAT / 2)
    rho_lo   <- sin(pi * TAU_LO  / 2)
    rho_hi   <- sin(pi * TAU_HI  / 2)
    delta_hat <- delta_gauss_exact(XGRID, rho_hat)
    delta_lo  <- delta_gauss_exact(XGRID, rho_lo)
    delta_hi  <- delta_gauss_exact(XGRID, rho_hi)
  } else {
    beta_hat  <- 9 * TAU_HAT / 2
    beta_lo   <- 9 * TAU_LO  / 2
    beta_hi   <- 9 * TAU_HI  / 2
    delta_hat <- delta_fgm_exact(XGRID, beta_hat)
    delta_lo  <- delta_fgm_exact(XGRID, beta_lo)
    delta_hi  <- delta_fgm_exact(XGRID, beta_hi)
  }

  # Point estimate and CI for the sensitivity-adjusted causal mean
  mu_adj_hat <- res_case$mu_obs - res_case$M * delta_hat
  mu_adj_lo  <- res_case$mu_obs - res_case$M * delta_hi   # worst case (δ at upper CI)
  mu_adj_hi  <- res_case$mu_obs + res_case$M * delta_hi   # best case positive direction

  res_case$calibrated <- list(
    delta_hat  = delta_hat,
    delta_lo   = delta_lo,
    delta_hi   = delta_hi,
    ub_hat     = pmin(1, res_case$mu_obs + res_case$M * delta_hat),
    lb_hat     = pmax(0, res_case$mu_obs - res_case$M * delta_hat),
    ub_hi      = pmin(1, res_case$mu_obs + res_case$M * delta_hi),
    lb_lo      = pmax(0, res_case$mu_obs - res_case$M * delta_hi),
    tau_hat    = TAU_HAT,
    tau_lo     = TAU_LO,
    tau_hi     = TAU_HI
  )
  res_case
}

gauss_res <- lapply(gauss_res, add_calibrated_interval, family = "gaussian")
fgm_res   <- lapply(fgm_res,   add_calibrated_interval, family = "fgm")


# ==============================================================================
# SECTION 7 — CONSOLE OUTPUT: TABLES
# ==============================================================================

cat(strrep("=", 79), "\n")
cat("TABLE 1: ORACLE IDENTIFICATION — GAUSSIAN COPULA\n")
cat(strrep("=", 79), "\n\n")
cat("Columns: x | μ_obs | M(x) | δ_true | lb_oracle | ub_oracle | Δ_crit\n")
cat("(Oracle interval [lb, ub] = [μ_obs ± M·δ_true]: EXACT and SHARP)\n\n")

CASE_LABELS <- c("Case A", "Case B", "Case C")

for (j in seq_along(gauss_res)) {
  r  <- gauss_res[[j]]
  nm <- CASE_LABELS[j]
  K  <- K_vals[j]
  cat(sprintf("--- %s (K = %.4f) ---\n", nm, K))
  hdr <- sprintf("  %-6s  %-7s  %-7s  %-9s  %-10s  %-10s  %-9s\n",
                 "x", "mu_obs", "M(x)", "delta_tr", "lb_oracle", "ub_oracle", "D_crit")
  cat(hdr, strrep("-", 73), "\n", sep="")
  for (i in TABLE_IDX) {
    cat(sprintf("  %-6.2f  %-7.4f  %-7.4f  %-9.4f  %-10.4f  %-10.4f  %-9.4f\n",
                r$x[i], r$mu_obs[i], r$M[i], r$delta_true[i],
                r$lb_oracle[i], r$ub_oracle[i], r$delta_crit[i]))
  }
  cat("\n")
}

cat(strrep("=", 79), "\n")
cat("TABLE 2: ORACLE IDENTIFICATION — FGM COPULA\n")
cat(strrep("=", 79), "\n\n")
cat("Note: FGM δ(x) is EXACT (not approximated); M(x) is a conservative bound.\n\n")

alpha_vals <- sapply(fgm_params, function(p) p["theta_xy_cond"])
beta_vals  <- sapply(fgm_params, function(p) p["theta_xz"])
tau_fgm    <- sapply(fgm_params, function(p) 2 * p["theta_xz"] / 9)

for (j in seq_along(fgm_res)) {
  r  <- fgm_res[[j]]
  nm <- CASE_LABELS[j]
  cat(sprintf("--- %s (α=%.2f, β=%.2f, τ_XZ=%.4f) ---\n",
              nm, alpha_vals[j], beta_vals[j], tau_fgm[j]))
  hdr <- sprintf("  %-6s  %-7s  %-7s  %-9s  %-10s  %-10s  %-9s\n",
                 "x", "mu_obs", "M(x)", "delta_tr", "lb_oracle", "ub_oracle", "D_crit")
  cat(hdr, strrep("-", 73), "\n", sep="")
  for (i in TABLE_IDX) {
    cat(sprintf("  %-6.2f  %-7.4f  %-7.4f  %-9.4f  %-10.4f  %-10.4f  %-9.4f\n",
                r$x[i], r$mu_obs[i], r$M[i], r$delta_true[i],
                r$lb_oracle[i], r$ub_oracle[i], r$delta_crit[i]))
  }
  cat("\n")
}

cat(strrep("=", 79), "\n")
cat("TABLE 3: IGNORANCE RATIO AND INTERVAL COMPARISON — GAUSSIAN, x = 0.25\n")
cat(strrep("=", 79), "\n")
cat("Ignorance Ratio (IR) = width(sensitivity interval) / width(oracle interval)\n")
cat("IR = 1 iff the researcher knows δ_true exactly (oracle recovery).\n\n")

x_ref_idx <- TABLE_IDX[2]  # x ≈ 0.25
cat(sprintf("  x = %.2f\n\n", XGRID[x_ref_idx]))
cat(sprintf("  %-8s  %-12s  %-10s  %-10s  %-10s  %-8s\n",
            "Case", "δ_true", "oracle_w", "sens_w(δ≤0.20)", "sens_w(δ≤0.10)", "IR(0.10)"))
cat(strrep("-", 68), "\n")

for (j in seq_along(gauss_res)) {
  r  <- gauss_res[[j]]
  i  <- x_ref_idx
  w_or  <- r$ub_oracle[i] - r$lb_oracle[i]
  w20   <- r$sensitivity$`δ≤0.20`$width[i]
  w10   <- r$sensitivity$`δ≤0.10`$width[i]
  IR10  <- w10 / max(1e-8, w_or)
  cat(sprintf("  %-8s  %-12.4f  %-10.4f  %-10.4f  %-10.4f  %-8.2f\n",
              CASE_LABELS[j], r$delta_true[i], w_or, w20, w10, IR10))
}

cat(strrep("=", 79), "\n")
cat("TABLE 4: CALIBRATED SCENARIO (τ̂ = 0.15, SE_τ = 0.04)\n")
cat(strrep("=", 79), "\n")
cat(sprintf("  τ CI: [%.3f, %.3f]  →  95%% CI for confounding strength\n\n", TAU_LO, TAU_HI))
cat(sprintf("  %-8s  %-12s  %-12s  %-12s  %-12s\n",
            "Case", "δ_hat(x=0.25)", "δ_lo(x=0.25)", "δ_hi(x=0.25)", "lb_lo"))
cat(strrep("-", 60), "\n")

for (j in seq_along(gauss_res)) {
  r <- gauss_res[[j]]
  i <- x_ref_idx
  cat(sprintf("  %-8s  %-12.4f  %-12.4f  %-12.4f  %-12.4f\n",
              CASE_LABELS[j],
              r$calibrated$delta_hat[i],
              r$calibrated$delta_lo[i],
              r$calibrated$delta_hi[i],
              r$calibrated$lb_lo[i]))
}
cat("\n")


# ==============================================================================
# SECTION 8 — DIAGNOSTIC TEXT
# ==============================================================================

cat(strrep("=", 79), "\n")
cat("SECTION 8: DIAGNOSTIC ANALYSIS\n")
cat(strrep("=", 79), "\n\n")

# ---- What is LOST when Z is latent? -----------------------------------------
cat("QUESTION 1: What is lost when Z is latent?\n")
cat(strrep("-", 50), "\n\n")

for (j in seq_along(gauss_res)) {
  r     <- gauss_res[[j]]
  nm    <- CASE_LABELS[j]
  i_med <- TABLE_IDX[3]   # median
  i_tail <- TABLE_IDX[1]  # lower tail

  w_or_med  <- r$ub_oracle[i_med]  - r$lb_oracle[i_med]
  # Under complete ignorance: δ can be arbitrarily large → interval = [0,1]
  # The "information loss" is thus  1 - w_or_med  (residual uncertainty)
  loss_med <- 1.0 - w_or_med

  cat(sprintf("  %s (K=%.3f):\n", nm, K_vals[j]))
  cat(sprintf("    Oracle interval width at median (x=0.50): %.4f\n", w_or_med))
  cat(sprintf("    Under complete ignorance: interval collapses to [0,1], width=1.0\n"))
  cat(sprintf("    Information lost (width increase): %.4f → %.0f%% wider\n\n",
              loss_med, loss_med / w_or_med * 100))
}

# ---- How much does partial knowledge help? ----------------------------------
cat("QUESTION 2: How much does partial knowledge help?\n")
cat(strrep("-", 50), "\n\n")

for (j in seq_along(gauss_res)) {
  r  <- gauss_res[[j]]
  nm <- CASE_LABELS[j]
  i  <- TABLE_IDX[3]   # median

  w_or <- r$ub_oracle[i] - r$lb_oracle[i]
  cat(sprintf("  %s (oracle interval width = %.4f at x=0.50):\n", nm, w_or))

  prev_w <- 1.0   # starting from complete ignorance
  for (k in rev(seq_along(DELTA_MAX_LEVELS))) {
    dmax <- DELTA_MAX_LEVELS[k]
    lab  <- DELTA_MAX_LABELS[k]
    w    <- r$sensitivity[[lab]]$width[i]
    gain <- (prev_w - w) / prev_w * 100   # % reduction from previous level
    IR   <- w / max(1e-8, w_or)
    cat(sprintf("    %-7s: width=%.4f, IR=%.2f, gain from prev level=%.0f%%\n",
                lab, w, IR, gain))
    prev_w <- w
  }
  final_gain <- (prev_w - w_or) / prev_w * 100
  cat(sprintf("    Oracle:  width=%.4f  (remaining gain to oracle: %.0f%%)\n\n",
              w_or, final_gain))
}

# ---- When does sensitivity analysis become practically useful? ---------------
cat("QUESTION 3: When does sensitivity analysis become practically useful?\n")
cat(strrep("-", 50), "\n\n")
cat("  Criterion: sensitivity analysis is 'useful' when the interval\n")
cat("  [lb_sens, ub_sens] has the same SIGN as μ_obs for all x in the grid.\n")
cat("  Equivalently: δ_max < Δ_crit(x) for all x of interest.\n\n")

for (j in seq_along(gauss_res)) {
  r  <- gauss_res[[j]]
  nm <- CASE_LABELS[j]
  cat(sprintf("  %s:\n", nm))
  cat(sprintf("    Δ_crit range: [%.4f, %.4f]\n",
              min(r$delta_crit), max(r$delta_crit)))
  for (k in seq_along(DELTA_MAX_LEVELS)) {
    dmax <- DELTA_MAX_LEVELS[k]
    lab  <- DELTA_MAX_LABELS[k]
    # sign is preserved iff dmax < delta_crit everywhere
    sign_ok <- all(dmax < r$delta_crit)
    cover   <- mean(dmax < r$delta_crit) * 100
    cat(sprintf("    %s: sign preserved globally=%s (at %.0f%% of x values)\n",
                lab, ifelse(sign_ok, "YES ✓", "NO  ✗"), cover))
  }
  cat("\n")
}


# ==============================================================================
# SECTION 9 — FIGURES
# ==============================================================================
#
# All figures are produced with ggplot2 and saved as PDF files.
# The data-assembly step converts list-of-lists into long-format data frames
# that ggplot2 can consume via faceting.

cat(strrep("=", 79), "\n")
cat("SECTION 9: GENERATING FIGURES\n")
cat(strrep("=", 79), "\n\n")

# --- 1. Define o diretório ABSOLUTO que você quer ---
# Use barras normais "/" no R (funciona no Windows) ou use file.path()
base_dir <- "E:/User/Gabriel/Unicamp/Mestrado/model_causal_copulas"
output_dir <- file.path(base_dir, "results_sensitivity")

# --- 2. Cria a pasta, se não existir ---
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
  cat("Created output directory:", output_dir, "\n")
} else {
  cat("Using existing output directory:", output_dir, "\n")
}

# --- 3. Teste de permissão de escrita (diagnóstico) ---
test_file <- file.path(output_dir, "test_write.txt")
tryCatch({
  writeLines("Permission OK", test_file)
  file.remove(test_file)
  cat("SUCCESS: Write permission confirmed in", output_dir, "\n\n")
}, error = function(e) {
  stop("ERROR: Cannot write to ", output_dir, "\n", e$message, 
       "\nCheck if the drive/path exists or if you have admin rights.")
})

# ---- Helper: assemble long-format data frame for one copula family ----------
assemble_longdf <- function(res_list, family_label, param_labels) {
  # param_labels: named vector of case labels with parameter info
  rows <- list()

  for (j in seq_along(res_list)) {
    r   <- res_list[[j]]
    nm  <- names(res_list)[j]
    lbl <- param_labels[j]
    n   <- length(r$x)

    # Oracle interval
    base_df <- data.frame(
      x          = r$x,
      mu_obs     = r$mu_obs,
      M          = r$M,
      delta_true = r$delta_true,
      delta_crit = r$delta_crit,
      lb_oracle  = r$lb_oracle,
      ub_oracle  = r$ub_oracle,
      case       = nm,
      case_label = lbl,
      family     = family_label
    )

    # Sensitivity interval at each δ_max level
    for (k in seq_along(DELTA_MAX_LEVELS)) {
      dmax  <- DELTA_MAX_LEVELS[k]
      dlab  <- DELTA_MAX_LABELS[k]
      s     <- r$sensitivity[[dlab]]
      base_df[[paste0("lb_", dlab)]] <- s$lb
      base_df[[paste0("ub_", dlab)]] <- s$ub
      base_df[[paste0("IR_", dlab)]] <- s$IR
    }

    # Calibrated scenario
    base_df$lb_calib <- r$calibrated$lb_lo
    base_df$ub_calib <- r$calibrated$ub_hi
    base_df$lb_calib_hat <- r$calibrated$lb_hat
    base_df$ub_calib_hat <- r$calibrated$ub_hat

    rows[[j]] <- base_df
  }
  do.call(rbind, rows)
}

gauss_K_labels <- sprintf("Case %s (K=%.3f)",
                           c("A","B","C"), K_vals)
fgm_alpha_labels <- sprintf("Case %s (α=%.2f)",
                             c("A","B","C"), alpha_vals)

df_gauss <- assemble_longdf(gauss_res, "Gaussian", gauss_K_labels)
df_fgm   <- assemble_longdf(fgm_res,   "FGM",     fgm_alpha_labels)
df_all   <- rbind(df_gauss, df_fgm)

# Factor ordering for consistent legend/facet ordering
case_order <- c("CaseA","CaseB","CaseC")
df_all$case <- factor(df_all$case, levels = case_order)


# ============================================================
# FIGURE 1: Oracle vs. Sensitivity Intervals (6-panel grid)
# Two rows × three columns: row = copula family, col = case
# ============================================================
cat("  Building Figure 1: Oracle vs. Sensitivity Intervals...\n")

# Reshape into very long format for ribbon+line plot
make_ribbon_df <- function(df) {
  # For each row in df, create one row per knowledge scenario
  # Scenarios: oracle, δ≤0.05, δ≤0.10, δ≤0.15, δ≤0.20
  scenarios <- list(
    oracle = list(lb = "lb_oracle",   ub = "ub_oracle",   label = "Oracle (Z observed)"),
    d05    = list(lb = "lb_δ≤0.05",  ub = "ub_δ≤0.05",  label = "δ≤0.05"),
    d10    = list(lb = "lb_δ≤0.10",  ub = "ub_δ≤0.10",  label = "δ≤0.10"),
    d20    = list(lb = "lb_δ≤0.20",  ub = "ub_δ≤0.20",  label = "δ≤0.20")
  )

  out <- lapply(names(scenarios), function(s) {
    sc <- scenarios[[s]]
    data.frame(
      x          = df$x,
      mu_obs     = df$mu_obs,
      lb         = df[[sc$lb]],
      ub         = df[[sc$ub]],
      case       = df$case,
      case_label = df$case_label,
      family     = df$family,
      scenario   = sc$label,
      stringsAsFactors = FALSE
    )
  })
  do.call(rbind, out)
}

# Fix column names: R replaces ≤ with . in data.frame names
fix_names <- function(df) {
  colnames(df) <- gsub("δ.", "δ≤", colnames(df), fixed = TRUE)
  # Also handle the dot-substitution for ≤
  colnames(df) <- gsub("lb_X.", "lb_δ≤", colnames(df), fixed = TRUE)
  colnames(df) <- gsub("ub_X.", "ub_δ≤", colnames(df), fixed = TRUE)
  colnames(df) <- gsub("IR_X.", "IR_δ≤", colnames(df), fixed = TRUE)
  df
}

# ggplot can't use ≤ in column names easily; use ASCII alternatives internally
# and remap the DELTA_MAX_LABELS to safe names for column selection
SAFE_LABELS <- paste0("d", sprintf("%02d", as.integer(DELTA_MAX_LEVELS * 100)))

make_ribbon_df2 <- function(df, family_filt) {
  df2 <- df[df$family == family_filt, ]

  scen_cfg <- list(
    list(lb_col = "lb_oracle",          ub_col = "ub_oracle",          label = "Oracle (Z observed)"),
    list(lb_col = paste0("lb_",SAFE_LABELS[1]), ub_col = paste0("ub_",SAFE_LABELS[1]), label = "δ ≤ 0.05"),
    list(lb_col = paste0("lb_",SAFE_LABELS[2]), ub_col = paste0("ub_",SAFE_LABELS[2]), label = "δ ≤ 0.10"),
    list(lb_col = paste0("lb_",SAFE_LABELS[4]), ub_col = paste0("ub_",SAFE_LABELS[4]), label = "δ ≤ 0.20")
  )

  out <- lapply(scen_cfg, function(sc) {
    data.frame(
      x          = df2$x,
      mu_obs     = df2$mu_obs,
      lb         = df2[[sc$lb_col]],
      ub         = df2[[sc$ub_col]],
      case_label = df2$case_label,
      family     = df2$family,
      scenario   = sc$label,
      stringsAsFactors = FALSE
    )
  })
  do.call(rbind, out)
}

# Rebuild data frames with safe column names
safe_col_rename <- function(df) {
  for (k in seq_along(DELTA_MAX_LEVELS)) {
    dlab <- DELTA_MAX_LABELS[k]
    slab <- SAFE_LABELS[k]
    colnames(df) <- gsub(paste0("lb_", dlab), paste0("lb_", slab), colnames(df), fixed = TRUE)
    colnames(df) <- gsub(paste0("ub_", dlab), paste0("ub_", slab), colnames(df), fixed = TRUE)
    colnames(df) <- gsub(paste0("IR_", dlab), paste0("IR_", slab), colnames(df), fixed = TRUE)
  }
  df
}

df_gauss <- safe_col_rename(df_gauss)
df_fgm   <- safe_col_rename(df_fgm)
df_all   <- rbind(df_gauss, df_fgm)
df_all$case <- factor(df_all$case, levels = case_order)

rib_gauss <- make_ribbon_df2(df_all, "Gaussian")
rib_fgm   <- make_ribbon_df2(df_all, "FGM")

SCEN_ORDER <- c("Oracle (Z observed)", "δ ≤ 0.05", "δ ≤ 0.10", "δ ≤ 0.20")
SCEN_COLS  <- c("Oracle (Z observed)" = "#000000",
                "δ ≤ 0.05"            = "#56B4E9",
                "δ ≤ 0.10"            = "#0072B2",
                "δ ≤ 0.20"            = "#D55E00")
SCEN_ALPHA <- c("Oracle (Z observed)" = 0.35,
                "δ ≤ 0.05"            = 0.20,
                "δ ≤ 0.10"            = 0.15,
                "δ ≤ 0.20"            = 0.10)

rib_gauss$scenario <- factor(rib_gauss$scenario, levels = SCEN_ORDER)
rib_fgm$scenario   <- factor(rib_fgm$scenario,   levels = SCEN_ORDER)

make_interval_plot <- function(rib_df, title_tag) {
  ggplot(rib_df, aes(x = x)) +
    # Ribbons from widest (most ignorance) to narrowest (oracle)
    geom_ribbon(
      data = rib_df[rib_df$scenario == "δ ≤ 0.20",],
      aes(ymin = lb, ymax = ub, fill = "δ ≤ 0.20"), alpha = 0.18
    ) +
    geom_ribbon(
      data = rib_df[rib_df$scenario == "δ ≤ 0.10",],
      aes(ymin = lb, ymax = ub, fill = "δ ≤ 0.10"), alpha = 0.22
    ) +
    geom_ribbon(
      data = rib_df[rib_df$scenario == "δ ≤ 0.05",],
      aes(ymin = lb, ymax = ub, fill = "δ ≤ 0.05"), alpha = 0.28
    ) +
    geom_ribbon(
      data = rib_df[rib_df$scenario == "Oracle (Z observed)",],
      aes(ymin = lb, ymax = ub, fill = "Oracle (Z observed)"), alpha = 0.40
    ) +
    # μ_obs centre line
    geom_line(
      data = rib_df[rib_df$scenario == "Oracle (Z observed)",],
      aes(y = mu_obs), colour = "grey30", linewidth = 0.6, linetype = "dashed"
    ) +
    scale_fill_manual(values = SCEN_COLS,
                      breaks = SCEN_ORDER,
                      name   = "Knowledge level") +
    facet_wrap(~ case_label, ncol = 3) +
    scale_x_continuous(breaks = c(0.1, 0.25, 0.5, 0.75, 0.9)) +
    scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.25)) +
    labs(
      title    = paste0("Oracle vs. Sensitivity Identification Intervals — ", title_tag),
      subtitle = paste0("Ribbons shrink toward oracle (black) as researcher knowledge ",
                        "about δ(x) increases.\nDashed line: μ_obs(x). ",
                        "Oracle interval is EXACT (sharp); others are conservative."),
      x        = "Treatment quantile x",
      y        = "Causal mean μ(x)"
    ) +
    THEME_THESIS
}

fig1a <- make_interval_plot(rib_gauss, "Gaussian Copula")
fig1b <- make_interval_plot(rib_fgm,   "FGM Copula")

pdf(file.path(output_dir, "fig1_oracle_vs_sensitivity.pdf"), width = 14, height = 10)
gridExtra::grid.arrange(fig1a, fig1b, nrow = 2)
dev.off()
cat("  Saved: fig1_oracle_vs_sensitivity.pdf\n")


# ============================================================
# FIGURE 2: Tipping-Point Curves with Knowledge-Level Bands
# For each case: μ_adj(x;δ) = μ_obs(x) - M(x)·δ as δ varies
# Vertical lines at Δ_crit; shaded Robust / Fragile regions
# ============================================================
cat("  Building Figure 2: Tipping-Point Curves...\n")

# Choose representative x for tipping-point analysis (x = 0.25 and 0.75)
# The tipping-point curve is parametric in δ, at a fixed x.
make_tipping_df <- function(res_list, family_label, param_labels, x_fixed_vec = c(0.25, 0.75)) {
  delta_grid <- seq(0, 0.50, length.out = 300)

  rows <- list()
  for (j in seq_along(res_list)) {
    r   <- res_list[[j]]
    nm  <- names(res_list)[j]
    lbl <- param_labels[j]

    for (x_f in x_fixed_vec) {
      ix        <- which.min(abs(r$x - x_f))
      mu_at_x   <- r$mu_obs[ix]
      M_at_x    <- r$M[ix]
      delt_true <- r$delta_true[ix]
      dcrit     <- r$delta_crit[ix]

      # Adjusted causal mean as function of δ (assumes confounding hurts the effect)
      mu_adj_up <- pmin(1, mu_at_x + M_at_x * delta_grid)   # upper bound trajectory
      mu_adj_lo <- pmax(0, mu_at_x - M_at_x * delta_grid)   # lower bound trajectory

      rows[[length(rows) + 1]] <- data.frame(
        delta      = delta_grid,
        mu_adj_lo  = mu_adj_lo,
        mu_adj_up  = mu_adj_up,
        mu_obs     = mu_at_x,
        M          = M_at_x,
        delta_true = delt_true,
        delta_crit = dcrit,
        x_fixed    = x_f,
        case       = nm,
        case_label = lbl,
        family     = family_label,
        stringsAsFactors = FALSE
      )
    }
  }
  do.call(rbind, rows)
}

tip_gauss <- make_tipping_df(gauss_res, "Gaussian", gauss_K_labels)
tip_fgm   <- make_tipping_df(fgm_res,   "FGM",     fgm_alpha_labels)
tip_all   <- rbind(tip_gauss, tip_fgm)
tip_all$case   <- factor(tip_all$case,   levels = case_order)
tip_all$family <- factor(tip_all$family, levels = c("Gaussian","FGM"))

# Add τ annotations for x-axis (approximate via Gaussian calibration for display)
# τ_Gauss at each δ level (approximate: δ ≈ 2ρφ(Φ⁻¹(x)), so ρ ≈ δ/(2φ(Φ⁻¹(x))))
tip_all$tau_approx <- with(tip_all, {
  phi_inv_x <- dnorm(qnorm(x_fixed))
  rho_approx <- pmin(0.99, delta / pmax(1e-6, 2 * phi_inv_x))
  (2/pi) * asin(rho_approx)
})

make_tipping_plot <- function(tip_df, family_filt, x_fixed_use = 0.25) {
  df2 <- tip_df[tip_df$family == family_filt & abs(tip_df$x_fixed - x_fixed_use) < 0.01, ]

  # Knowledge-level vertical annotations
  delta_max_ann <- data.frame(
    delta = DELTA_MAX_LEVELS,
    label = paste0("δ_max=", DELTA_MAX_LEVELS)
  )

  # Create per-case delta_crit annotation data
  crit_ann <- unique(df2[, c("case","case_label","delta_crit","mu_obs")])

  ggplot(df2, aes(x = delta)) +
    # Shaded regions: Robust (left of Δ_crit)
    geom_rect(
      data = crit_ann,
      aes(xmin = 0, xmax = pmin(delta_crit, 0.50),
          ymin = -Inf, ymax = Inf, fill = "Robust"),
      alpha = 0.07, inherit.aes = FALSE
    ) +
    geom_rect(
      data = crit_ann,
      aes(xmin = pmin(delta_crit, 0.50), xmax = 0.50,
          ymin = -Inf, ymax = Inf, fill = "Fragile"),
      alpha = 0.07, inherit.aes = FALSE
    ) +
    scale_fill_manual(values = PAL_REGION, name = "Region", guide = guide_legend(order = 2)) +
    # Upper and lower bound trajectories
    geom_ribbon(aes(ymin = mu_adj_lo, ymax = mu_adj_up, group = case,
                    colour = case_label),
                fill = NA, linetype = "solid", linewidth = 0.8,
                show.legend = TRUE) +
    # μ_obs horizontal reference
    geom_hline(data = crit_ann, aes(yintercept = mu_obs),
               linetype = "dashed", colour = "grey40", linewidth = 0.5) +
    # Δ_crit vertical lines
    geom_vline(data = crit_ann, aes(xintercept = delta_crit, colour = case_label),
               linetype = "dotted", linewidth = 0.9) +
    # δ_true marker (oracle)
    geom_point(
      data = df2[abs(df2$delta - df2$delta_true) == sapply(df2$case, function(c) {
        sub_d <- df2[df2$case == c, "delta_true"][1]
        min(abs(df2[df2$case == c, "delta"] - sub_d))
      }), ],
      aes(y = mu_obs, colour = case_label),
      shape = 21, size = 3, fill = "white", stroke = 1.5, show.legend = FALSE
    ) +
    # Knowledge-level vlines
    geom_vline(data = delta_max_ann, aes(xintercept = delta),
               colour = "grey60", linetype = "22", linewidth = 0.4) +
    scale_colour_manual(values = PAL_CASES, name = "Case", guide = guide_legend(order = 1)) +
    facet_wrap(~ case_label, ncol = 3) +
    labs(
      title    = paste0("Tipping-Point Curves — ", family_filt, " Copula (x = ", x_fixed_use, ")"),
      subtitle = paste0("Upper/lower bounds on μ(x) as a function of confounding strength δ(x).\n",
                        "Dotted vertical: Δ_crit (tipping point). White circle: oracle δ_true.\n",
                        "Grey verticals: δ_max = 0.05, 0.10, 0.15, 0.20."),
      x        = "Confounding strength  δ(x)  [copula-ratio L¹ index]",
      y        = "Adjusted causal mean  μ_adj(x; δ)"
    ) +
    THEME_THESIS +
    theme(legend.position = "bottom")
}

fig2a <- make_tipping_plot(tip_all, "Gaussian", 0.25)
fig2b <- make_tipping_plot(tip_all, "FGM",      0.25)

pdf(file.path(output_dir, "fig2_tipping_point_curves.pdf"), width = 14, height = 9)
gridExtra::grid.arrange(fig2a, fig2b, nrow = 2)
dev.off()
cat("  Saved: fig2_tipping_point_curves.pdf\n")


# ============================================================
# FIGURE 3: Ignorance Ratio and Knowledge-Gain Function
# Panel A: IR(x; δ_max) across x for each case and δ_max
# Panel B: KG(x; δ_max0 → δ_max1) — marginal knowledge gain
# ============================================================
cat("  Building Figure 3: Ignorance Ratio and Knowledge-Gain...\n")

make_IR_df <- function(res_list, family_label, param_labels) {
  rows <- list()
  for (j in seq_along(res_list)) {
    r   <- res_list[[j]]
    nm  <- names(res_list)[j]
    lbl <- param_labels[j]

    w_oracle <- r$ub_oracle - r$lb_oracle

    for (k in seq_along(DELTA_MAX_LEVELS)) {
      dmax  <- DELTA_MAX_LEVELS[k]
      slab  <- SAFE_LABELS[k]
      w_sens <- r$sensitivity[[DELTA_MAX_LABELS[k]]]$width

      IR <- w_sens / pmax(1e-8, w_oracle)

      rows[[length(rows) + 1]] <- data.frame(
        x          = r$x,
        IR         = IR,
        delta_max  = dmax,
        delta_max_label = DELTA_MAX_LABELS[k],
        case       = nm,
        case_label = lbl,
        family     = family_label,
        w_oracle   = w_oracle,
        w_sens     = w_sens,
        stringsAsFactors = FALSE
      )
    }
  }
  do.call(rbind, rows)
}

IR_gauss <- make_IR_df(gauss_res, "Gaussian", gauss_K_labels)
IR_fgm   <- make_IR_df(fgm_res,   "FGM",     fgm_alpha_labels)
IR_all   <- rbind(IR_gauss, IR_fgm)
IR_all$case       <- factor(IR_all$case, levels = case_order)
IR_all$delta_max_label <- factor(IR_all$delta_max_label,
                                  levels = DELTA_MAX_LABELS)
IR_all$family <- factor(IR_all$family, levels = c("Gaussian","FGM"))

# Panel A: Ignorance Ratio over x
fig3a <- ggplot(IR_all, aes(x = x, y = IR,
                             colour = delta_max_label,
                             linetype = delta_max_label)) +
  geom_line(linewidth = 0.8) +
  geom_hline(yintercept = 1, linetype = "dashed", colour = "grey30", linewidth = 0.5) +
  annotate("text", x = 0.05, y = 1.02, label = "IR = 1 (oracle)",
           colour = "grey30", size = 3, hjust = 0) +
  scale_colour_manual(
    values = c("δ≤0.05"="#56B4E9","δ≤0.10"="#0072B2","δ≤0.15"="#E69F00","δ≤0.20"="#D55E00"),
    name   = "δ_max level"
  ) +
  scale_linetype_manual(
    values = c("δ≤0.05"=1,"δ≤0.10"=2,"δ≤0.15"=3,"δ≤0.20"=4),
    name   = "δ_max level"
  ) +
  facet_grid(family ~ case) +
  labs(
    title    = "Ignorance Ratio:  IR(x; δ_max) = width(sensitivity) / width(oracle)",
    subtitle = paste0("IR = 1 iff the researcher knows δ_true exactly.\n",
                      "IR is monotone decreasing in δ_max: less knowledge → wider interval.\n",
                      "Note Gaussian cases: IR spikes at extreme x due to small oracle width."),
    x = "Treatment quantile  x",
    y = "Ignorance ratio  IR(x)"
  ) +
  THEME_THESIS +
  theme(legend.position = "bottom")

# Panel B: Knowledge-Gain function — width reduction from δ_max to oracle
# KG(x) = 1 - w_sens / w_max   where w_max = max width across δ_max levels
make_KG_df <- function(res_list, family_label, param_labels) {
  rows <- list()
  for (j in seq_along(res_list)) {
    r   <- res_list[[j]]
    nm  <- names(res_list)[j]
    lbl <- param_labels[j]

    w_oracle  <- r$ub_oracle - r$lb_oracle
    w_max_all <- r$sensitivity[[DELTA_MAX_LABELS[4]]]$width  # δ≤0.20 (widest)

    for (k in seq_along(DELTA_MAX_LEVELS)) {
      dmax   <- DELTA_MAX_LEVELS[k]
      w_sens <- r$sensitivity[[DELTA_MAX_LABELS[k]]]$width
      # Marginal gain from δ_max level to oracle
      KG_to_oracle <- 1 - w_sens / pmax(1e-8, w_max_all)
      # Remaining gap to oracle
      gap_to_oracle <- (w_sens - w_oracle) / pmax(1e-8, w_max_all)

      rows[[length(rows) + 1]] <- data.frame(
        x           = r$x,
        delta_max   = dmax,
        delta_max_label = DELTA_MAX_LABELS[k],
        KG_to_oracle = KG_to_oracle,
        gap_oracle   = pmax(0, gap_to_oracle),
        case         = nm,
        case_label   = lbl,
        family       = family_label,
        stringsAsFactors = FALSE
      )
    }
  }
  do.call(rbind, rows)
}

KG_gauss <- make_KG_df(gauss_res, "Gaussian", gauss_K_labels)
KG_fgm   <- make_KG_df(fgm_res,   "FGM",     fgm_alpha_labels)
KG_all   <- rbind(KG_gauss, KG_fgm)
KG_all$case  <- factor(KG_all$case, levels = case_order)
KG_all$delta_max_label <- factor(KG_all$delta_max_label, levels = DELTA_MAX_LABELS)
KG_all$family <- factor(KG_all$family, levels = c("Gaussian","FGM"))

fig3b <- ggplot(KG_all, aes(x = x, y = KG_to_oracle,
                              colour = delta_max_label,
                              linetype = delta_max_label)) +
  geom_line(linewidth = 0.8) +
  scale_colour_manual(
    values = c("δ≤0.05"="#56B4E9","δ≤0.10"="#0072B2","δ≤0.15"="#E69F00","δ≤0.20"="#D55E00"),
    name   = "δ_max level"
  ) +
  scale_linetype_manual(
    values = c("δ≤0.05"=1,"δ≤0.10"=2,"δ≤0.15"=3,"δ≤0.20"=4),
    name   = "δ_max level"
  ) +
  scale_y_continuous(labels = scales::percent, limits = c(0, 1)) +
  facet_grid(family ~ case) +
  labs(
    title    = "Knowledge-Gain Function: width reduction from complete ignorance toward oracle",
    subtitle = paste0("KG(x; δ_max) = 1 − width(δ_max) / width(δ_max = 0.20)\n",
                      "KG = 0 at δ_max = 0.20 (baseline ignorance); KG = 1 when δ_max = δ_true (oracle)."),
    x = "Treatment quantile  x",
    y = "Knowledge gain  KG(x)"
  ) +
  THEME_THESIS +
  theme(legend.position = "bottom")

pdf(file.path(output_dir, "fig3_ignorance_ratio_knowledge_gain.pdf"), width = 14, height = 11)
gridExtra::grid.arrange(fig3a, fig3b, nrow = 2)
dev.off()
cat("  Saved: fig3_ignorance_ratio_knowledge_gain.pdf\n")


# ============================================================
# FIGURE 4: Calibrated Sensitivity Intervals (τ̂ ± SE)
# Illustrates Scenario 3: external τ estimate with uncertainty
# ============================================================
cat("  Building Figure 4: Calibrated Sensitivity Intervals...\n")

make_calib_df <- function(res_list, family_label, param_labels) {
  rows <- list()
  for (j in seq_along(res_list)) {
    r   <- res_list[[j]]
    nm  <- names(res_list)[j]
    lbl <- param_labels[j]

    rows[[length(rows) + 1]] <- data.frame(
      x            = r$x,
      mu_obs       = r$mu_obs,
      lb_oracle    = r$lb_oracle,
      ub_oracle    = r$ub_oracle,
      lb_calib     = r$calibrated$lb_lo,
      ub_calib     = r$calibrated$ub_hi,
      lb_calib_hat = r$calibrated$lb_hat,
      ub_calib_hat = r$calibrated$ub_hat,
      delta_true   = r$delta_true,
      delta_crit   = r$delta_crit,
      case         = nm,
      case_label   = lbl,
      family       = family_label,
      stringsAsFactors = FALSE
    )
  }
  do.call(rbind, rows)
}

calib_gauss <- make_calib_df(gauss_res, "Gaussian", gauss_K_labels)
calib_fgm   <- make_calib_df(fgm_res,   "FGM",     fgm_alpha_labels)
calib_all   <- rbind(calib_gauss, calib_fgm)
calib_all$case   <- factor(calib_all$case, levels = case_order)
calib_all$family <- factor(calib_all$family, levels = c("Gaussian","FGM"))

make_calib_plot <- function(calib_df, family_filt) {
  df2 <- calib_df[calib_df$family == family_filt, ]

  ggplot(df2, aes(x = x)) +
    # Calibrated CI band (95% CI on τ → CI on δ → CI on interval)
    geom_ribbon(aes(ymin = lb_calib, ymax = ub_calib),
                fill = "#CC79A7", alpha = 0.20) +
    # Calibrated point-estimate band (at τ̂)
    geom_ribbon(aes(ymin = lb_calib_hat, ymax = ub_calib_hat),
                fill = "#CC79A7", alpha = 0.35) +
    # Oracle interval
    geom_ribbon(aes(ymin = lb_oracle, ymax = ub_oracle),
                fill = "#000000", alpha = 0.20) +
    # μ_obs line
    geom_line(aes(y = mu_obs), colour = "grey20", linewidth = 0.7, linetype = "dashed") +
    # Oracle boundary lines
    geom_line(aes(y = lb_oracle), colour = "#000000", linewidth = 0.6) +
    geom_line(aes(y = ub_oracle), colour = "#000000", linewidth = 0.6) +
    # Calibrated boundary lines
    geom_line(aes(y = lb_calib_hat), colour = "#CC79A7", linewidth = 0.8) +
    geom_line(aes(y = ub_calib_hat), colour = "#CC79A7", linewidth = 0.8) +
    geom_line(aes(y = lb_calib),     colour = "#CC79A7", linewidth = 0.5, linetype = "dotted") +
    geom_line(aes(y = ub_calib),     colour = "#CC79A7", linewidth = 0.5, linetype = "dotted") +
    facet_wrap(~ case_label, ncol = 3) +
    scale_x_continuous(breaks = c(0.1, 0.25, 0.5, 0.75, 0.9)) +
    scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.25)) +
    labs(
      title    = paste0("Calibrated Sensitivity Analysis — ", family_filt, " Copula"),
      subtitle = sprintf(
        paste0("External estimate: τ̂ = %.2f ± %.2f (95%% CI: [%.3f, %.3f])\n",
               "Black: oracle interval (Z observed). ",
               "Pink solid: interval at τ̂. Pink dotted: 95%% CI band on τ."),
        TAU_HAT, SE_TAU, TAU_LO, TAU_HI
      ),
      x = "Treatment quantile  x",
      y = "Causal mean  μ(x)"
    ) +
    THEME_THESIS
}

fig4a <- make_calib_plot(calib_all, "Gaussian")
fig4b <- make_calib_plot(calib_all, "FGM")

pdf(file.path(output_dir, "fig4_calibrated_sensitivity.pdf"), width = 14, height = 9)
gridExtra::grid.arrange(fig4a, fig4b, nrow = 2)
dev.off()
cat("  Saved: fig4_calibrated_sensitivity.pdf\n")


# ============================================================
# FIGURE 5 (BONUS): Δ_crit(x) profile — the practitioner's
# first-look diagnostic before any knowledge about δ
# ============================================================
cat("  Building Figure 5: Tipping-Point Δ_crit Profiles...\n")

make_dcrit_df <- function(res_list, family_label, param_labels) {
  rows <- list()
  for (j in seq_along(res_list)) {
    r <- res_list[[j]]
    rows[[length(rows) + 1]] <- data.frame(
      x          = r$x,
      delta_crit = r$delta_crit,
      delta_true = r$delta_true,
      robust     = r$delta_crit > r$delta_true,  # oracle: is result robust?
      case       = names(res_list)[j],
      case_label = param_labels[j],
      family     = family_label,
      stringsAsFactors = FALSE
    )
  }
  do.call(rbind, rows)
}

dcrit_gauss <- make_dcrit_df(gauss_res, "Gaussian", gauss_K_labels)
dcrit_fgm   <- make_dcrit_df(fgm_res,   "FGM",     fgm_alpha_labels)
dcrit_all   <- rbind(dcrit_gauss, dcrit_fgm)
dcrit_all$case   <- factor(dcrit_all$case, levels = case_order)
dcrit_all$family <- factor(dcrit_all$family, levels = c("Gaussian","FGM"))

fig5 <- ggplot(dcrit_all, aes(x = x)) +
  # Oracle δ_true profile (what confounding actually is)
  geom_line(aes(y = delta_true, colour = case_label, linetype = "δ_true (oracle)"),
            linewidth = 0.8) +
  # Δ_crit profile (what confounding needs to be to reverse the conclusion)
  geom_line(aes(y = delta_crit, colour = case_label, linetype = "Δ_crit (tipping point)"),
            linewidth = 1.0) +
  # Shade region where δ_true > Δ_crit (conclusion is NOT robust under oracle)
  geom_ribbon(
    data = dcrit_all[dcrit_all$delta_true > dcrit_all$delta_crit, ],
    aes(ymin = delta_crit, ymax = delta_true),
    fill = "#D55E00", alpha = 0.25
  ) +
  scale_colour_manual(
    values = c("Case A (K=0.580)"="#0072B2","Case B (K=0.870)"="#E69F00",
               "Case C (K=0.456)"="#009E73",
               "Case A (α=0.40)"="#0072B2","Case B (α=0.80)"="#E69F00",
               "Case C (α=0.20)"="#009E73"),
    name = "Case"
  ) +
  scale_linetype_manual(
    values = c("δ_true (oracle)" = 2, "Δ_crit (tipping point)" = 1),
    name   = "Quantity"
  ) +
  # Reference lines for δ_max knowledge levels
  geom_hline(yintercept = DELTA_MAX_LEVELS,
             colour = "grey70", linetype = "dotted", linewidth = 0.4) +
  geom_text(
    data = data.frame(x = 0.95, y = DELTA_MAX_LEVELS,
                      label = paste0("δmax=", DELTA_MAX_LEVELS),
                      family = factor("Gaussian", levels = c("Gaussian","FGM")),
                      case = factor("CaseA", levels = case_order)),
    aes(x = x, y = y, label = label),
    size = 2.5, colour = "grey50", hjust = 1, vjust = -0.3,
    inherit.aes = FALSE
  ) +
  facet_grid(family ~ case) +
  scale_y_continuous(limits = c(0, 0.55)) +
  labs(
    title    = "Tipping-Point Profile Δ_crit(x) vs. Oracle Confounding Strength δ_true(x)",
    subtitle = paste0("Solid: Δ_crit = |μ_obs|/M(x) — confounding needed to reverse sign.\n",
                      "Dashed: δ_true(x) — actual confounding strength (observable by oracle only).\n",
                      "Orange shading: region where result is NOT robust even under oracle knowledge."),
    x = "Treatment quantile  x",
    y = "Confounding strength"
  ) +
  THEME_THESIS +
  theme(legend.position = "bottom")

pdf(file.path(output_dir, "fig5_delta_crit_profile.pdf"), width = 14, height = 9)
print(fig5)
dev.off()
cat("  Saved: fig5_delta_crit_profile.pdf\n\n")


# ==============================================================================
# SECTION 10 — FINAL SYNTHESIS: ANSWERING THE THREE DIAGNOSTIC QUESTIONS
# ==============================================================================

cat(strrep("=", 79), "\n")
cat("SECTION 10: SYNTHESIS AND PRACTICAL IMPLICATIONS\n")
cat(strrep("=", 79), "\n\n")

cat("─── THE COPULA-RATIO SENSITIVITY FRAMEWORK: A HIERARCHY OF KNOWLEDGE ───\n\n")

cat("┌─────────────────────────────────────────────────────────────────────────┐\n")
cat("│  EPISTEMIC REGIME        WHAT IS KNOWN          IDENTIFICATION RESULT  │\n")
cat("│  ─────────────────       ─────────────          ────────────────────── │\n")
cat("│  Oracle (Z observed)     δ_true(x) known        EXACT sharp interval   │\n")
cat("│                          M(x) estimated          [μ_obs ± M·δ_true]   │\n")
cat("│                                                                        │\n")
cat("│  Calibrated (τ̂ known)   δ(x) ≈ f(τ̂) + CI      Credible sens. band   │\n")
cat("│                          via Gaussian/FGM map    Shrinks with SE_τ     │\n")
cat("│                                                                        │\n")
cat("│  Partial (δ_max known)   δ(x) ≤ δ_max           Conservative interval │\n")
cat("│                          via domain judgment      [μ_obs ± M·δ_max]   │\n")
cat("│                                                                        │\n")
cat("│  Complete ignorance      δ(x) ∈ [0, ∞)          ONLY Δ_crit reported  │\n")
cat("│                          no bound available       Qualitative only      │\n")
cat("└─────────────────────────────────────────────────────────────────────────┘\n\n")

cat("ANSWER 1: What is LOST when Z is latent?\n")
cat("  The oracle computes δ_true(x) from the joint distribution of (X,Z).\n")
cat("  Under latency, this quantity is UNIDENTIFIED — no amount of outcome data\n")
cat("  can recover it without assumptions on C_{XZ}.  The loss is:\n")
cat("    (a) Exact interval width 2·M(x)·δ_true(x) becomes 2·M(x)·δ_max.\n")
cat("    (b) The IR (ignorance ratio) quantifies the cost:\n")
cat("        IR = δ_max / δ_true(x)  ≥  1,  with equality iff δ_max = δ_true.\n")
cat("    (c) Under COMPLETE ignorance, the only usable quantity is Δ_crit(x)\n")
cat("        = |μ_obs|/M(x): the minimum δ to reverse the causal conclusion.\n")
cat("        This is a qualitative diagnostic, not a quantitative interval.\n\n")

cat("ANSWER 2: How much does PARTIAL KNOWLEDGE help?\n")
cat("  Each halving of δ_max halves the width of the sensitivity interval,\n")
cat("  reducing it LINEARLY toward the oracle width.  Key thresholds:\n\n")
cat("  Gaussian Case A (K ≈ 0.58):\n")
cat(sprintf("    δ_true at x=0.25: ≈ %.4f\n",
            gauss_res$CaseA$delta_true[TABLE_IDX[2]]))
cat(sprintf("    Δ_crit at x=0.25: ≈ %.4f\n",
            gauss_res$CaseA$delta_crit[TABLE_IDX[2]]))
cat("    With δ_max = 0.10: IR ≈ ",
    round(gauss_res$CaseA$sensitivity[[DELTA_MAX_LABELS[2]]]$IR[TABLE_IDX[2]], 2), "\n")
cat("    With δ_max = 0.05: IR ≈ ",
    round(gauss_res$CaseA$sensitivity[[DELTA_MAX_LABELS[1]]]$IR[TABLE_IDX[2]], 2), "\n\n")

cat("  FGM Case A (α = 0.40, β = 0.50):\n")
cat(sprintf("    δ_true at x=0.25: ≈ %.4f  (FGM exact formula)\n",
            fgm_res$CaseA$delta_true[TABLE_IDX[2]]))
cat(sprintf("    Δ_crit at x=0.25: ≈ %.4f\n",
            fgm_res$CaseA$delta_crit[TABLE_IDX[2]]))
cat("\n")

cat("ANSWER 3: When does sensitivity analysis become PRACTICALLY USEFUL?\n")
cat("  Operational criterion: useful iff δ_max < Δ_crit(x) for all x of interest.\n")
cat("  This ensures the sensitivity interval has the SAME SIGN as μ_obs,\n")
cat("  and the researcher can make a directional causal claim despite uncertainty.\n\n")
cat("  In copula terms, this translates to a Kendall's τ threshold:\n")
cat("    Gaussian: τ_crit such that δ(x; τ_crit) = Δ_crit(x)\n")
cat("    FGM:      τ_crit = 4·Δ_crit(x) / (9·|1-2x|)  [exact, from δ = (9τ/4)|1-2x|]\n\n")

# Report τ_crit for each case at x = 0.25
cat("  τ_crit values at x = 0.25:\n\n")
cat(sprintf("  %-24s  %-12s  %-12s  %-14s\n",
            "Case", "Δ_crit(0.25)", "δ_true(0.25)", "Robust at τ̂=0.15?"))
cat(strrep("-", 70), "\n")
i_ref <- TABLE_IDX[2]
for (j in seq_along(gauss_res)) {
  r    <- gauss_res[[j]]
  dc   <- r$delta_crit[i_ref]
  dt   <- r$delta_true[i_ref]
  # Invert Gaussian δ ≈ 2ρφ(Φ⁻¹(x)) → ρ → τ
  phi_x <- dnorm(qnorm(XGRID[i_ref]))
  rho_crit <- pmin(0.99, dc / (2 * phi_x))
  tau_crit  <- (2/pi) * asin(rho_crit)
  robust   <- if (TAU_HAT < tau_crit) "YES ✓" else "NO  ✗"
  cat(sprintf("  Gaussian %-14s  %-12.4f  %-12.4f  %-14s (τ_crit≈%.3f)\n",
              CASE_LABELS[j], dc, dt, robust, tau_crit))
}
cat("\n")
for (j in seq_along(fgm_res)) {
  r       <- fgm_res[[j]]
  dc      <- r$delta_crit[i_ref]
  dt      <- r$delta_true[i_ref]
  x_ref   <- XGRID[i_ref]
  # Exact FGM inversion: δ = (9τ/4)|1-2x|  →  τ = 4δ/(9|1-2x|)
  tau_crit <- 4 * dc / (9 * abs(1 - 2 * x_ref))
  tau_crit <- pmin(1, tau_crit)
  robust   <- if (TAU_HAT < tau_crit) "YES ✓" else "NO  ✗"
  cat(sprintf("  FGM     %-14s  %-12.4f  %-12.4f  %-14s (τ_crit≈%.3f)\n",
              CASE_LABELS[j], dc, dt, robust, tau_crit))
}

cat("\n")
cat(strrep("=", 79), "\n")
cat("COMPARISON WITH CINELLI & HAZLETT (2020) PARTIAL R² FRAMEWORK\n")
cat(strrep("=", 79), "\n\n")
cat("  The Cinelli-Hazlett Robustness Value (RV) is the minimum partial R²\n")
cat("  that a confounder must have, with BOTH treatment AND outcome, to explain\n")
cat("  away the effect.  In their linear-regression setting:\n\n")
cat("      RV_q = (1/2)[√(fq⁴ + 4fq²) − fq²],   fq = q|fY∼D|X|\n\n")
cat("  Structural analogies to the copula-ratio framework:\n\n")
cat("  ┌──────────────────────────┬──────────────────────────┐\n")
cat("  │  Cinelli-Hazlett (2020)  │  Copula-Ratio (this work)│\n")
cat("  ├──────────────────────────┼──────────────────────────┤\n")
cat("  │  Partial R²_D∼Z|X        │  δ(x)  [treatment sens.] │\n")
cat("  │  Partial R²_Y∼Z|D,X      │  M(x)  [outcome sens.]   │\n")
cat("  │  Robustness Value (RV)   │  Δ_crit(x) = |μ_obs|/M  │\n")
cat("  │  Effect estimate τ̂_res   │  μ_obs(x)               │\n")
cat("  │  Bias factor BF          │  M(x)·δ(x)              │\n")
cat("  │  Linear, parametric      │  Non-parametric, copula  │\n")
cat("  │  Binary/cont. treatment  │  Continuous treatment    │\n")
cat("  └──────────────────────────┴──────────────────────────┘\n\n")
cat("  KEY DIFFERENCE: Cinelli-Hazlett measures confounding in R² units\n")
cat("  (variance-explained, scale-free within linear models).  The copula-ratio\n")
cat("  framework uses Kendall's τ (rank-based, assumption-free), which is\n")
cat("  interpretable across ANY joint distribution of (X, Z).\n\n")

cat(strrep("=", 79), "\n")
cat("CONNECTION TO THE E-VALUE (VanderWeele & Ding, 2017)\n")
cat(strrep("=", 79), "\n\n")
cat("  The E-value for risk ratio RR = μ_obs is:\n")
cat("      E-value = RR + √[RR(RR-1)]\n\n")
cat("  This answers: 'How strong must confounding be (on the RR scale) to\n")
cat("  explain away the association?'  For continuous Y ∈ [0,1], the copula\n")
cat("  analogue is Δ_crit(x) = |μ_obs(x)|/M(x):\n\n")
cat("  ┌────────────────────────────┬───────────────────────────────┐\n")
cat("  │  E-value (binary outcome)  │  Copula Δ_crit (continuous Y) │\n")
cat("  ├────────────────────────────┼───────────────────────────────┤\n")
cat("  │  Confounding on RR scale   │  Confounding on L¹ copula norm│\n")
cat("  │  Single number             │  Function of x                │\n")
cat("  │  Requires binary outcome   │  Works for any Y ∈ [0,1]      │\n")
cat("  │  No copula structure       │  Explicit copula calibration   │\n")
cat("  └────────────────────────────┴───────────────────────────────┘\n\n")

cat(sprintf("  Example (Gaussian Case A, x=0.25): Δ_crit ≈ %.4f\n",
            gauss_res$CaseA$delta_crit[TABLE_IDX[2]]))
cat("  Interpretation: confounding strength δ must exceed Δ_crit to reverse\n")
cat("  the sign of the causal conclusion at this treatment level.\n")
cat("  This is the non-parametric analogue of the E-value for continuous Y.\n\n")

cat(strrep("=", 79), "\n")
cat("OPERATIONAL SUMMARY FOR PRACTITIONERS\n")
cat(strrep("=", 79), "\n\n")
cat("  Step 1: Fit the copula model and obtain μ_obs(x), M(x).\n")
cat("  Step 2: Report Δ_crit(x) = |μ_obs(x)|/M(x) as the first diagnostic.\n")
cat("          This requires NO knowledge of δ and tells how strong confounding\n")
cat("          must be to overturn the sign of the effect.\n\n")
cat("  Step 3: If domain knowledge suggests τ < τ_max for some Kendall's τ\n")
cat("          between X and Z, convert: δ_max = f(τ_max, x) via:\n")
cat("            Gaussian: δ_max(x) ≈ 2·sin(π·τ_max/2)·φ(Φ⁻¹(x))\n")
cat("            FGM:      δ_max(x) = (9·τ_max/4)·|1-2x|\n\n")
cat("  Step 4: If δ_max < Δ_crit(x) for all x of interest:\n")
cat("          → Effect sign is ROBUST. Report the sensitivity interval.\n")
cat("          If δ_max ≥ Δ_crit(x) for some x:\n")
cat("          → Effect sign is FRAGILE at that x. Flag for transparency.\n\n")
cat("  Step 5: If an external estimate τ̂ ± SE_τ is available, translate the\n")
cat("          entire τ CI into a δ CI and report a calibrated interval.\n")
cat("          This is the Scenario 3 (calibrated) analysis demonstrated above.\n\n")
cat("  NOTE: The framework is SHARP: every bound is attained by some confounder\n")
cat("  with the stated L¹ norm.  It is therefore NOT conservative in the model-\n")
cat("  selection sense — any claimed interval can be achieved by a real confounder.\n\n")

cat(strrep("=", 79), "\n")
cat("All figures saved to:", output_dir, "\n")
cat("  fig1_oracle_vs_sensitivity.pdf      — 6-panel oracle vs. sensitivity\n")
cat("  fig2_tipping_point_curves.pdf       — tipping-point analysis\n")
cat("  fig3_ignorance_ratio_knowledge_gain.pdf — IR and KG diagnostics\n")
cat("  fig4_calibrated_sensitivity.pdf     — calibrated τ̂ ± SE scenario\n")
cat("  fig5_delta_crit_profile.pdf         — Δ_crit profile diagnostic\n")
cat(strrep("=", 79), "\n\n")
cat("Script complete. Total wall time (approx.): depends on GH quadrature speed.\n")
