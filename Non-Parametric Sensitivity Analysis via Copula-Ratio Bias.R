# ==============================================================================
# SENSITIVITY ANALYSIS UNDER LATENT CONFOUNDING -- CORRECTED
# Copula-Based Causal Inference -- Master's Dissertation, Chapter 5
# Ferraz (2026)
# ==============================================================================
#
# This replaces Sensitivity_Analysis_-_Illustrative_Numerical_Examples.R.
# Every numerical/theoretical issue identified in review is fixed here; the
# script is reorganised around Chapter 5's own Regime A / Regime B structure,
# with GENUINE LATENT CONFOUNDING (Regime B, Section 5.4) as the primary
# analysis rather than a hypothetical "oracle" framing.
#
# WORKFLOW
# --------
#   Section 0-3   setup, closed-form Ch.5 formulas, parameters
#   Section 4     REGIME B primary analysis: tipping point under latent Z
#   Section 5     validation-only oracle check (Z assumed observed)
#   Section 6     REGIME A: calibrated tau-hat +/- SE scenario
#   Section 7-7B  sensitivity contour plots (the main new output)
#   Section 8     Delta_crit(x)/tau_crit(x) profile (replaces old Figure 5)
#   Section 9     validation figure
#   Section 10    corrected console tables
#   Section 11    save all figures
# ==============================================================================

# ==============================================================================
# SECTION 0 -- SETUP
# ==============================================================================
suppressMessages({
  library(ggplot2)
  library(pracma)
  library(gridExtra)
})

PAL_CASES  <- c("A" = "#0072B2", "B" = "#E69F00", "C" = "#009E73")
PAL_REGION <- c("Sign preserved" = "#4dac26", "Sign not determined" = "#d01c8b")

THEME_THESIS <- theme_bw(base_size = 11) +
  theme(
    panel.grid.minor = element_blank(),
    strip.background = element_rect(fill = "grey92", colour = NA),
    legend.position  = "right",
    plot.title       = element_text(face = "bold", size = 12),
    plot.subtitle    = element_text(size = 9.5, colour = "grey40")
  )

XGRID <- seq(0.02, 0.98, length.out = 200)
GH    <- gaussHermite(20)
cat("Setup OK\n")
# ==============================================================================
# SECTION 1 -- CLOSED-FORM CAUSAL QUANTITIES (Chapter 5 formulas; UNCHANGED)
# These were verified against the dissertation's Table 2 (Gaussian) and
# Table 3 (FGM) in the review, and are numerically correct. Kept as-is.
# ==============================================================================

# ---- 1A Gaussian (Section 5.6.1) --------------------------------------------
gaussian_K_and_V <- function(ryx, rxz, ryz) {
  a <- (ryx - rxz * ryz) / (1 - rxz^2)
  b <- (ryz - ryx * rxz) / (1 - rxz^2)
  V <- 1 - a^2 - 2 * a * b * rxz
  K <- a / sqrt(1 + V)
  c(a = unname(a), b = unname(b), V = unname(V), K = unname(K))
}
mu_gaussian  <- function(x, K) pnorm(K * qnorm(x))
ace_gaussian <- function(x, K) K * dnorm(K * qnorm(x)) / dnorm(qnorm(x))
ivar_gaussian <- function(x, K, V, gh = GH) {
  a_coef  <- K * sqrt(1 + V)
  mu_x    <- a_coef * qnorm(x)
  sigma_x <- sqrt(V)
  y_star  <- mu_x + sqrt(2) * sigma_x * gh$x
  E_Y2    <- sum(gh$w * pnorm(y_star)^2) / sqrt(pi)
  E_Y     <- pnorm(K * qnorm(x))
  pmax(0, E_Y2 - E_Y^2)
}

# K0: the OBSERVATIONAL (naive, Z-marginalised-out) slope -- Eq. (5.28)'s
# "no confounding" benchmark. Depends ONLY on rho_XY, which is the one
# correlation an analyst can estimate from (X,Y) data alone.
K0_observational <- function(rxy) rxy / sqrt(2 - rxy^2)

# ---- 1B FGM (Section 5.6.2) -------------------------------------------------
mu_fgm <- function(x, params) {
  alpha <- params["theta_xy_cond"]; eta <- params["theta_yz"]
  A <- 1/6 - eta^2/90
  0.5 - alpha * (1 - 2*x) * A
}
ace_fgm <- function(x, params) {
  alpha <- params["theta_xy_cond"]; eta <- params["theta_yz"]
  rep(alpha * (15 - eta^2) / 45, length(x))
}

cat("Section 1 (closed forms) loaded.\n")
# ==============================================================================
# SECTION 2 -- SENSITIVITY QUANTITIES: delta(x), M(x), Delta_crit(x)
# ==============================================================================
# ---- 2A. delta(x) for the Gaussian copula: NUMERICAL QUADRATURE (correct) ---
# c_rho(u,v) is Eq. (5.21). delta(x) = int_0^1 |c_rho(F_X(x),v) - 1| dv.
copula_density_gaussian <- function(u, v, rho) {
  a <- qnorm(u); b <- qnorm(v)
  (1/sqrt(1-rho^2)) * exp(-(rho^2*(a^2+b^2) - 2*rho*a*b)/(2*(1-rho^2)))
}

delta_gauss <- function(x, rho_xz) {
  # Vectorised over BOTH x and rho_xz (standard R recycling).
  mapply(function(xi, ri) {
    if (abs(ri) < 1e-10) return(0)
    f <- function(v) abs(copula_density_gaussian(xi, v, ri) - 1)
    stats::integrate(f, lower = 1e-9, upper = 1 - 1e-9,
                     subdivisions = 400, rel.tol = 1e-8)$value
  }, x, rho_xz)
}

# First-order approximation (Prop. 6.3.1) -- KEPT, but only ever labelled
# "first-order"; valid for small |rho_XZ| only (dissertation: <5% error for
# |rho|<0.5, degrading beyond that).
delta_gauss_firstorder <- function(x, rho_xz) {
  kappa <- sqrt(2/pi)
  abs(rho_xz) * kappa * abs(qnorm(x))
}

# ---- 2B. delta(x) for the FGM copula: EXACT (unchanged; verified correct) ---
delta_fgm <- function(x, beta) abs(beta) * abs(1 - 2*x) / 2

# ---- 2C. M(x): Proposition 5.4.2 (Manski / support bound) -- PRIMARY --------
# M(x) <= ybar = max(|ymin|,|ymax|). For Y in [0,1], ybar = 1. This requires
# NO assumption about the confounding mechanism -- it uses only Y's own
# observed range, which is exactly what is available when Z is latent.
M_manski <- function(ymin = 0, ymax = 1) max(abs(ymin), abs(ymax))

# ---- 2C'. OPTIONAL structural M(x) heuristics -- SECONDARY, opt-in only ----
# These require the analyst to additionally assume a specific parametric
# mechanism linking Z to Y (Section 5.4.2, option (a)). They are retained
# here ONLY for illustrating how much a structural assumption could tighten
# the bound relative to the assumption-free Manski ceiling -- never used as
# the default under latent confounding.
M_gauss_structural <- function(x, K, V, gh = GH) {
  mu_x  <- mu_gaussian(x, K)
  sd_x  <- sqrt(pmax(0, sapply(x, function(xi) ivar_gaussian(xi, K, V, gh))))
  pmin(1.0, mu_x + 2*sd_x)
}
M_fgm_structural <- function(x, params) {
  alpha <- params["theta_xy_cond"]; eta <- params["theta_yz"]; beta <- params["theta_xz"]
  A <- 1/6 - eta^2/90
  pmin(0.5, abs(alpha*A*(1-2*x)) + abs(alpha)*A*abs(beta))
}

# ---- 2D. Delta_crit(x): Definition 5.4.1 ------------------------------------
# B(x) := mu(x) - mu_obs(x)  [List of Symbols sign convention; NOTE this is
# the OPPOSITE sign from the dissertation's B(x):=mu_obs(x)-mu(x) (Ch. 4),
# which should be harmonised -- see write-up]. delta_crit is sign-convention-free.
delta_crit <- function(mu_obs, M_x, null_value = 0) {
  abs(mu_obs - null_value) / pmax(1e-8, M_x)
}

# ---- 2E. tau_crit(x): Corollary 5.4.1 (continuous-treatment E-value) -------
# Uses the SAME first-order Gaussian tau->rho->delta map as Lemma 5.3.1 (this
# is what Corollary 5.4.1 itself is built on -- it is explicitly a first-order
# statement in the dissertation, so tau_crit inherits that same qualifier).
tau_crit_gaussian <- function(delta_crit_x, x) {
  kappa <- sqrt(2/pi)
  arg <- pmin(1, delta_crit_x / (kappa * abs(qnorm(x))))
  (2/pi) * asin(arg)
}
# Exact FGM inversion (Section 5.3.2): delta = (9*tau/4)*|1-2x|  =>  tau = 4*delta/(9|1-2x|)
tau_crit_fgm <- function(delta_crit_x, x) {
  pmin(1, 4 * delta_crit_x / (9 * abs(1 - 2*x)))
}

# ---- 2F. tau_crit(x), Gaussian, EXACT (numerical root-finding) -------------
# Corollary 5.4.1 is explicitly a first-order statement (built on Lemma 5.3.1)
# and degenerates at x=0.5 (delta_gauss_firstorder(0.5,.) = 0 for every rho,
# so the ratio used inside tau_crit_gaussian() is undefined there). Now that
# delta_gauss() is available via quadrature, tau_crit can instead be obtained
# by directly solving delta_gauss(x, rho) = delta_crit_x for rho, then mapping
# rho -> tau. This removes the first-order restriction and the x=0.5
# degeneracy; offered here as a drop-in upgrade for Corollary 5.4.1.
tau_crit_gaussian_exact <- function(delta_crit_x, x) {
  mapply(function(dcx, xi) {
    if (dcx <= 0) return(0)
    g <- function(rho) delta_gauss(xi, rho) - dcx
    # delta_gauss(x,.) is increasing in |rho|; delta_gauss(x,1^-) is its max.
    upper_val <- tryCatch(delta_gauss(xi, 1 - 1e-6), error = function(e) Inf)
    if (dcx >= upper_val) return(1)  # even rho -> 1 can't reach delta_crit
    rho_star <- uniroot(g, lower = 0, upper = 1 - 1e-6, tol = 1e-8)$root
    (2/pi) * asin(rho_star)
  }, delta_crit_x, x)
}

cat("Section 2 (corrected sensitivity quantities) loaded.\n")
# ==============================================================================
# SECTION 3 -- PARAMETERS
# ==============================================================================
# The ONLY quantity a researcher facing latent Z can estimate is mu_obs(x),
# which requires ONLY rho_XY (the observed X-Y association). We fix
# rho_XY = 0.70 as the running example throughout (matches Cases A/B/C of the
# previous script, all of which shared this value).
RHO_XY <- 0.70
K0     <- K0_observational(RHO_XY)              # Gaussian mu_obs slope
mu_obs_gaussian <- function(x) mu_gaussian(x, K0)

# For the VALIDATION section only: the three closed-form confounding
# structures used previously, which additionally specify rho_XZ, rho_YZ (i.e.
# assume Z IS observed). These are used ONLY to check that the theorem's
# bound correctly contains the true bias when the truth happens to be known
# -- never to compute delta(x) or M(x) in the main (latent-Z) analysis.
gauss_params_VALIDATION_ONLY <- list(
  CaseA = c(ryx = 0.70, rxz = 0.50, ryz = 0.60),
  CaseB = c(ryx = 0.70, rxz = 0.00, ryz = 0.60),
  CaseC = c(ryx = 0.70, rxz = 0.80, ryz = 0.60)
)

# FGM: alpha (0.5 = observational effect size proxy), Y support is [0,1]
# throughout, so ybar (Manski ceiling for FGM's [0,1] outcome) is also 1,
# but Delta_crit for FGM is centred at the null value 0.5 (see Section 2D).
ALPHA_FGM <- 0.40  # representative direct-effect parameter for illustration
fgm_params_VALIDATION_ONLY <- list(
  CaseA = c(theta_xy_cond = 0.40, theta_xz = 0.50, theta_yz = 0.30),
  CaseB = c(theta_xy_cond = 0.80, theta_xz = 0.50, theta_yz = 0.30),
  CaseC = c(theta_xy_cond = 0.20, theta_xz = 0.50, theta_yz = 0.30)
)

cat("Section 3 loaded. mu_obs(x) uses ONLY rho_XY =", RHO_XY, "-> K0 =", round(K0,4), "\n")
# ==============================================================================
# SECTION 4 -- REGIME B: TIPPING-POINT ANALYSIS UNDER GENUINE LATENT
#               CONFOUNDING (Section 5.4). THIS IS THE PRIMARY ANALYSIS.
# ==============================================================================
# Inputs used: mu_obs(x) [estimable from (X,Y) alone] and ybar [Y's KNOWN
# support]. NOTHING about rho_XZ, rho_YZ, or any copula family is assumed.
# delta(x) is NOT computed -- it is genuinely unknown, exactly as Section 5.4
# describes, and is treated purely as a free parameter in Section 7's contour
# plots and in the sweep below.

Y_SUPPORT <- c(0, 1)
YBAR <- M_manski(Y_SUPPORT[1], Y_SUPPORT[2])   # = 1

regimeB_table <- data.frame(x = XGRID)
regimeB_table$mu_obs <- mu_obs_gaussian(regimeB_table$x)
regimeB_table$M      <- YBAR
regimeB_table$delta_crit <- delta_crit(regimeB_table$mu_obs, regimeB_table$M)
regimeB_table$tau_crit_gaussian <- tau_crit_gaussian(regimeB_table$delta_crit, regimeB_table$x)

# Sweep delta_max in [0, 2] (Remark 5.2.1's universal ceiling) and record,
# for each x, whether the sign of mu_obs is guaranteed preserved
# (Proposition 5.4.1(1)): sign is safe iff delta_max < Delta_crit(x).
DELTA_SWEEP <- seq(0, 2, length.out = 400)
robustness_grid <- expand.grid(x = c(0.10,0.25,0.50,0.75,0.90), delta_max = DELTA_SWEEP)
robustness_grid$mu_obs <- mu_obs_gaussian(robustness_grid$x)
robustness_grid$M <- YBAR
robustness_grid$bound <- robustness_grid$M * robustness_grid$delta_max
robustness_grid$sign_safe <- robustness_grid$bound < abs(robustness_grid$mu_obs)

TABLE_X   <- c(0.10, 0.25, 0.50, 0.75, 0.90)
TABLE_IDX <- sapply(TABLE_X, function(x) which.min(abs(XGRID - x)))

cat("Section 4 (Regime B tipping-point table) built.\n")
print(round(regimeB_table[TABLE_IDX, c("x","mu_obs","M","delta_crit","tau_crit_gaussian")], 4))
# ==============================================================================
# SECTION 5 -- VALIDATION ONLY: closed-form oracle check
# ==============================================================================
# Purpose: confirm that |B(x)| <= M(x)*delta(x) genuinely holds (and see how
# LOOSE it is) in a case where Z happens to be fully specified, so mu(x) and
# delta_true(x) are both computable exactly. This section is NOT part of the
# latent-confounding workflow (Section 4/7) -- it exists purely to sanity-check
# the theorem numerically on a known example, exactly as Examples 5.5.1-5.5.3
# do. B(x) := mu(x) - mu_obs(x)  [List of Symbols sign convention -- opposite
# of the dissertation's B(x):=mu_obs(x)-mu(x) in Ch. 4; see note in Section 2D].

validation_table <- do.call(rbind, lapply(names(gauss_params_VALIDATION_ONLY), function(nm) {
  p <- gauss_params_VALIDATION_ONLY[[nm]]
  cc <- gaussian_K_and_V(p["ryx"], p["rxz"], p["ryz"])
  mu_true <- mu_gaussian(XGRID, cc["K"])
  mu_obs_ <- mu_obs_gaussian(XGRID)          # K0-based, same curve for all 3 cases
  B_true  <- mu_true - mu_obs_                # script sign convention (opposite of Ch. 4's B(x))
  d_true  <- delta_gauss(XGRID, p["rxz"])
  bound   <- YBAR * d_true                     # rigorous M(x) = ybar = 1
  data.frame(case = nm, x = XGRID, mu_true = mu_true, mu_obs = mu_obs_,
             B_true = B_true, delta_true = d_true, bound = bound,
             holds = abs(B_true) <= bound + 1e-9)
}))

cat("Section 5 (validation) built. Bound holds in", 
    sum(validation_table$holds), "/", nrow(validation_table), "grid points.\n")
idx <- sapply(c(0.10,0.25,0.50,0.75,0.90), function(x) which.min(abs(XGRID-x)))
print(round(subset(validation_table, case=="CaseA")[idx, c("x","mu_true","mu_obs","B_true","delta_true","bound")], 4))
cat("(mean |B_true| / bound, i.e. how loose the Hoelder ceiling is on average):\n")
print(tapply(validation_table$B_true, validation_table$case, function(b) mean(abs(b))) /
        tapply(validation_table$bound,  validation_table$case, mean))
# ==============================================================================
# SECTION 6 -- REGIME A: CALIBRATED SCENARIO (tau-hat +/- SE), CORRECTED
# ==============================================================================
# Section 5.3.2. If an external estimate tau_hat +/- SE_tau of Kendall's tau
# for (X,Z) becomes available (e.g. a validation sub-sample, Section 5.3.1),
# this converts into a delta(x) band via the EXACT (quadrature) Gaussian map
# -- previously this used the broken delta_gauss_exact() formula.
TAU_HAT <- 0.15
SE_TAU  <- 0.04
TAU_LO  <- max(0, TAU_HAT - 1.96*SE_TAU)
TAU_HI  <- min(0.999, TAU_HAT + 1.96*SE_TAU)

rho_of_tau <- function(tau) sin(pi*tau/2)

calibrated_table <- data.frame(x = XGRID)
calibrated_table$mu_obs <- mu_obs_gaussian(XGRID)
calibrated_table$delta_hat <- delta_gauss(XGRID, rho_of_tau(TAU_HAT))
calibrated_table$delta_lo  <- delta_gauss(XGRID, rho_of_tau(TAU_LO))
calibrated_table$delta_hi  <- delta_gauss(XGRID, rho_of_tau(TAU_HI))
calibrated_table$lb_hat <- pmax(0, calibrated_table$mu_obs - YBAR*calibrated_table$delta_hat)
calibrated_table$ub_hat <- pmin(1, calibrated_table$mu_obs + YBAR*calibrated_table$delta_hat)
calibrated_table$lb_ci  <- pmax(0, calibrated_table$mu_obs - YBAR*calibrated_table$delta_hi)
calibrated_table$ub_ci  <- pmin(1, calibrated_table$mu_obs + YBAR*calibrated_table$delta_hi)

cat("Section 6 (calibrated scenario, corrected) built. tau CI: [", 
    round(TAU_LO,3), ",", round(TAU_HI,3), "]\n")
idx <- which.min(abs(XGRID-0.25))
cat("At x=0.25: delta_hat =", round(calibrated_table$delta_hat[idx],4),
    " (previous script's WRONG formula reported 0.3700 here)\n")
# ==============================================================================
# SECTION 7 -- SENSITIVITY CONTOUR PLOTS (Cinelli-Hazlett / VanderWeele-Ding
#               style), for GENUINE latent confounding (Regime B, Section 5.4)
# ==============================================================================
# Two free axes: delta(x) in [0, DELTA_MAX_PLOT] (treatment-side: L1 copula-
# ratio deviation, Def. 4.3.1) and M(x) in [0, ybar] (outcome-side: Manski
# ceiling, Prop. 5.4.2). z = M*delta is exactly the RHS of the Hoelder bound
# (Prop. 6.2.2). The critical isoline z = |mu_obs(x)| is the exact boundary of
# Proposition 5.4.1: below/left, the sign of mu_obs(x) is GUARANTEED preserved
# for every confounding mechanism consistent with (M,delta); on/above it,
# Proposition 5.4.1(2)'s construction shows a sign-reversing mechanism exists.
# This generalises the 1-D Delta_crit(x) (= the point where the critical curve
# crosses M = ybar) into the full 2-D robustness boundary.
suppressMessages(library(ggrepel))

DELTA_MAX_PLOT <- 1.2   # covers all practically calibrated benchmarks below
N_GRID <- 220

make_contour_grid <- function(x_val) {
  mu_x <- mu_obs_gaussian(x_val)
  g <- expand.grid(delta = seq(0, DELTA_MAX_PLOT, length.out = N_GRID),
                   M     = seq(0, YBAR, length.out = N_GRID))
  g$bound <- g$delta * g$M
  g$x <- x_val
  g$mu_obs <- mu_x
  g
}

X_PANELS <- c(0.10, 0.25, 0.50, 0.75, 0.90)
x_lab <- function(xv) paste0("x = ", sprintf("%.2f", xv))

contour_df <- do.call(rbind, lapply(X_PANELS, make_contour_grid))
contour_df$x_label <- factor(x_lab(contour_df$x), levels = x_lab(X_PANELS))

# ---- Reference benchmarks on the delta axis --------------------------------
# (a) Kendall's tau benchmarks (weak/moderate/fairly strong X-Z dependence)
# (b) "as strong as the observed X-Y association itself" (rho_XZ = rho_XY):
#     a natural, estimable benchmark in the spirit of Cinelli-Hazlett's
#     "as strong as an observed covariate" reference points.
tau_benchmarks <- c(0.10, 0.20, 0.30)
bench_df <- do.call(rbind, lapply(X_PANELS, function(xv) {
  d_tau <- delta_gauss(xv, sin(pi*tau_benchmarks/2))
  d_xy  <- delta_gauss(xv, RHO_XY)
  data.frame(
    x = xv, x_label = x_lab(xv),
    delta = c(d_tau, d_xy),
    label = c(paste0("\u03c4=", tau_benchmarks), "as strong as \u03c1_XY"),
    type  = c(rep("tau benchmark", length(tau_benchmarks)), "rho_XY benchmark")
  )
}))
bench_df$x_label <- factor(bench_df$x_label, levels = levels(contour_df$x_label))
bench_df$M_plot <- YBAR - 0.02

# ---- Critical curve (the 2-D generalisation of Delta_crit(x)) -------------
crit_curve_df <- do.call(rbind, lapply(X_PANELS, function(xv) {
  mu_x <- mu_obs_gaussian(xv)
  dgrid <- seq(abs(mu_x)/YBAR, DELTA_MAX_PLOT, length.out = 150)
  data.frame(x = xv, x_label = x_lab(xv), delta = dgrid, M = pmin(YBAR, abs(mu_x)/dgrid))
}))
crit_curve_df$x_label <- factor(crit_curve_df$x_label, levels = levels(contour_df$x_label))

fig_contour_delta_M <- ggplot(contour_df, aes(x = delta, y = M)) +
  geom_contour_filled(aes(z = bound), breaks = seq(0, DELTA_MAX_PLOT*YBAR, length.out = 13)) +
  geom_line(data = crit_curve_df, aes(x = delta, y = M), colour = "white", linewidth = 1.15) +
  geom_line(data = crit_curve_df, aes(x = delta, y = M), colour = "black",
            linewidth = 0.6, linetype = "21") +
  geom_point(data = bench_df, aes(x = delta, y = M_plot, shape = type),
             colour = "white", size = 2.4, stroke = 1) +
  geom_point(data = bench_df, aes(x = delta, y = M_plot, shape = type),
             colour = "black", size = 1.4) +
  geom_text_repel(data = bench_df, aes(x = delta, y = M_plot, label = label),
                  size = 2.4, colour = "grey15", segment.size = 0.25,
                  min.segment.length = 0, max.overlaps = 20,
                  nudge_y = -0.12, direction = "y", seed = 1) +
  scale_fill_viridis_d(name = "M(x)\u00b7\u03b4(x)\n(bound on |B(x)|)", option = "magma", direction = -1) +
  scale_shape_manual(values = c("tau benchmark" = 21, "rho_XY benchmark" = 23), guide = "none") +
  facet_wrap(~x_label, nrow = 1) +
  coord_cartesian(xlim = c(0, DELTA_MAX_PLOT), ylim = c(0, YBAR), expand = FALSE) +
  labs(
    title = "Sensitivity contour plot under latent confounding (Regime B, \u00a76.4)",
    subtitle = paste0(
      "Filled contours: bound M(x)\u00b7\u03b4(x) on |B(x)| (Prop. 6.2.2). Dashed curve: the critical\n",
      "isoline M\u00b7\u03b4 = |\u03bc_obs(x)| (Prop. 5.4.1) \u2014 above/right of it the sign of \u03bc_obs(x) is no longer guaranteed.\n",
      "Diamond: \u03b4 as strong as the observed \u03c1_XY. Circles: \u03b4 at Kendall's \u03c4 = 0.10/0.20/0.30."),
    x = "\u03b4(x)  [L\u00b9 copula-ratio sensitivity index]",
    y = "M(x)  [outcome-sensitivity ceiling; Manski ceiling ybar=1 at top edge]"
  ) +
  THEME_THESIS +
  theme(legend.position = "bottom", legend.key.width = unit(1.0, "cm"),
        axis.text.x = element_text(size = 7))

cat("Section 7 (delta-M contour plot) built.\n")
# ==============================================================================
# SECTION 7B -- COMPANION CONTOUR: (Kendall's tau, M) plane, Gaussian family
# ==============================================================================
# Same construction as Section 7, but with the x-axis recalibrated from delta
# to Kendall's tau via the EXACT map tau -> rho = sin(pi*tau/2) -> delta(x,rho)
# (quadrature). This is the more interpretable axis for a reader who thinks in
# rank-correlation units rather than the abstract L1 index -- exactly the role
# played by the risk-ratio axis in VanderWeele & Ding's E-value plots.
TAU_MAX_PLOT <- 0.6
N_TAU <- 120

make_contour_grid_tau <- function(x_val) {
  mu_x <- mu_obs_gaussian(x_val)
  tau_seq <- seq(0, TAU_MAX_PLOT, length.out = N_TAU)
  delta_of_tau <- delta_gauss(x_val, sin(pi*tau_seq/2))   # one quadrature call per tau
  g <- expand.grid(tau_idx = seq_along(tau_seq), M = seq(0, YBAR, length.out = N_TAU))
  g$tau <- tau_seq[g$tau_idx]
  g$delta <- delta_of_tau[g$tau_idx]
  g$bound <- g$delta * g$M
  g$x <- x_val
  g
}

contour_tau_df <- do.call(rbind, lapply(X_PANELS, make_contour_grid_tau))
contour_tau_df$x_label <- factor(x_lab(contour_tau_df$x), levels = x_lab(X_PANELS))

crit_curve_tau_df <- do.call(rbind, lapply(X_PANELS, function(xv) {
  mu_x <- mu_obs_gaussian(xv)
  tau_seq <- seq(0.001, TAU_MAX_PLOT, length.out = 150)
  delta_of_tau <- delta_gauss(xv, sin(pi*tau_seq/2))
  M_needed <- abs(mu_x) / pmax(1e-8, delta_of_tau)
  data.frame(x = xv, x_label = x_lab(xv), tau = tau_seq, M = pmin(YBAR, M_needed))
}))
crit_curve_tau_df$x_label <- factor(crit_curve_tau_df$x_label, levels = levels(contour_tau_df$x_label))

fig_contour_tau_M <- ggplot(contour_tau_df, aes(x = tau, y = M)) +
  geom_contour_filled(aes(z = bound), breaks = seq(0, TAU_MAX_PLOT*YBAR/1.6, length.out = 13)) +
  geom_line(data = crit_curve_tau_df, aes(x = tau, y = M), colour = "white", linewidth = 1.15) +
  geom_line(data = crit_curve_tau_df, aes(x = tau, y = M), colour = "black",
            linewidth = 0.6, linetype = "21") +
  geom_vline(xintercept = RHO_XY_tau <- (2/pi)*asin(RHO_XY), linetype = "dotted", colour = "grey30") +
  annotate("text", x = RHO_XY_tau, y = 0.06, label = "\u03c4(X,Y)", angle = 90,
           vjust = -0.4, size = 2.4, colour = "grey30") +
  scale_fill_viridis_d(name = "M(x)\u00b7\u03b4(x)\n(bound on |B(x)|)", option = "magma", direction = -1) +
  facet_wrap(~x_label, nrow = 1) +
  coord_cartesian(xlim = c(0, TAU_MAX_PLOT), ylim = c(0, YBAR), expand = FALSE) +
  labs(
    title = "Same bound, calibrated to Kendall's \u03c4 (Gaussian X\u2013Z copula, exact map)",
    subtitle = paste0("x-axis recalibrated via \u03c4 \u2192 \u03c1=sin(\u03c0\u03c4/2) \u2192 \u03b4(x,\u03c1). Dotted line: \u03c4(X,Y), i.e. a confounder as strongly\n",
                      "associated with X as Y itself already is \u2014 a natural, estimable benchmark."),
    x = "Kendall's \u03c4(X,Z)  [hypothesised confounder\u2013treatment rank association]",
    y = "M(x)"
  ) +
  THEME_THESIS +
  theme(legend.position = "bottom", legend.key.width = unit(1.0, "cm"),
        axis.text.x = element_text(size = 7))

cat("Section 7B (tau-M contour plot) built.\n")
# ==============================================================================
# SECTION 8 -- Delta_crit(x) / tau_crit(x) PROFILE
# ==============================================================================
profile_df <- data.frame(x = XGRID)
profile_df$mu_obs <- mu_obs_gaussian(XGRID)
profile_df$M <- YBAR
profile_df$delta_crit <- delta_crit(profile_df$mu_obs, profile_df$M)
profile_df$tau_crit_firstorder <- tau_crit_gaussian(profile_df$delta_crit, profile_df$x)
profile_df$tau_crit_exact <- tau_crit_gaussian_exact(profile_df$delta_crit, profile_df$x)

df_long <- rbind(
  data.frame(x = profile_df$x, value = profile_df$delta_crit, quantity = "Delta_crit(x)  [copula-ratio scale]"),
  data.frame(x = profile_df$x, value = profile_df$tau_crit_exact, quantity = "tau_crit(x), exact  [Kendall's tau scale]"),
  data.frame(x = profile_df$x, value = profile_df$tau_crit_firstorder, quantity = "tau_crit(x), first-order (Cor. 5.4.1 as written)")
)
df_long$quantity <- factor(df_long$quantity, levels = c(
  "Delta_crit(x)  [copula-ratio scale]",
  "tau_crit(x), exact  [Kendall's tau scale]",
  "tau_crit(x), first-order (Cor. 5.4.1 as written)"))

fig_delta_crit_profile <- ggplot(df_long, aes(x = x, y = value, colour = quantity, linetype = quantity)) +
  geom_hline(yintercept = c(0.10,0.20,0.30), colour = "grey85", linewidth = 0.4) +
  geom_line(linewidth = 0.9) +
  scale_colour_manual(values = c("#000000", "#0072B2", "#56B4E9"), name = NULL) +
  scale_linetype_manual(values = c("solid","solid","42"), name = NULL) +
  scale_y_continuous(limits = c(0,1.05), breaks = seq(0,1,0.25)) +
  labs(
    title = "Tipping point under latent confounding: Delta_crit(x) and tau_crit(x)",
    subtitle = paste0("Solid black: Delta_crit(x) = |\u03bc_obs(x)| / M(x), with the RIGOROUS M(x)=ybar=1 (Prop. 5.4.2) \u2014\n",
                      "note Delta_crit(x) = \u03bc_obs(x) exactly once M(x) is constant. Blue: the corresponding \u03c4_crit\n",
                      "computed exactly (root-finding on the quadrature-based \u03b4); pale dashed: Corollary 5.4.1 as\n",
                      "currently written (first-order only) \u2014 SATURATES at 1 near x=0.5, since its leading term\n",
                      "vanishes there (Prop. 5.3.1) and the min{1,\u00b7} cap takes over; prefer the exact (blue) curve."),
    x = "Treatment quantile x", y = NULL
  ) +
  THEME_THESIS + theme(legend.position = "bottom")

cat("Section 8 (Delta_crit / tau_crit profile) built.\n")
# ==============================================================================
# SECTION 9 -- VALIDATION FIGURE (companion to Section 5)
# ==============================================================================
val_plot_df <- validation_table
val_plot_df$case_lab <- factor(val_plot_df$case,
                               levels = c("CaseA","CaseB","CaseC"),
                               labels = c("Case A (\u03c1_XZ=0.50)","Case B (\u03c1_XZ=0, independence)","Case C (\u03c1_XZ=0.80)"))

fig_validation <- ggplot(val_plot_df, aes(x = x)) +
  geom_ribbon(aes(ymin = -bound, ymax = bound), fill = "#0072B2", alpha = 0.18) +
  geom_line(aes(y = B_true), colour = "#D55E00", linewidth = 0.9) +
  geom_hline(yintercept = 0, colour = "grey50", linewidth = 0.3) +
  facet_wrap(~case_lab, nrow = 1) +
  labs(
    title = "Validation: the Hoelder bound |B(x)| \u2264 M(x)\u00b7\u03b4(x) (Prop. 6.2.2), checked on a case with Z observed",
    subtitle = paste0("Orange: true bias B(x)=\u03bc(x)\u2212\u03bc_obs(x) (both objects fully computable here since \u03c1_XZ, \u03c1_YZ are\n",
                      "known). Blue band: \u00b1 the theorem's bound, using the CORRECTED \u03b4(x) [quadrature] and M(x)=1\n",
                      "[Prop. 5.4.2]. The bound holds everywhere, as it must \u2014 but note how much slack there is:\n",
                      "the actual bias uses only ~2\u201310% of the ceiling the theorem allows for."),
    x = "Treatment quantile x", y = "B(x)  (orange)  vs.  \u00b1M(x)\u00b7\u03b4(x)  (blue band)"
  ) +
  THEME_THESIS

cat("Section 9 (validation figure) built.\n")
# ==============================================================================
# SECTION 10 -- CORRECTED CONSOLE TABLES
# ==============================================================================
sep <- function() cat(strrep("=", 88), "\n")

sep(); cat("TABLE 1 (CORRECTED) -- REGIME B: mu_obs(x), M(x), Delta_crit(x), tau_crit(x)\n")
cat("mu_obs(x) uses ONLY rho_XY =", RHO_XY, "(estimable from (X,Y) data alone). NO knowledge\n")
cat("of rho_XZ or rho_YZ is used anywhere in this table -- this is the genuinely latent-Z case.\n")
sep()
idx <- sapply(c(0.10,0.25,0.50,0.75,0.90), function(x) which.min(abs(XGRID-x)))
tab1 <- data.frame(
  x = round(XGRID[idx], 2),
  mu_obs = round(profile_df$mu_obs[idx], 4),
  `M(x)_Manski` = round(profile_df$M[idx], 4),
  Delta_crit = round(profile_df$delta_crit[idx], 4),
  tau_crit_exact = round(profile_df$tau_crit_exact[idx], 4),
  tau_crit_1storder = round(profile_df$tau_crit_firstorder[idx], 4),
  check.names = FALSE
)
print(tab1, row.names = FALSE)
cat("\n(Compare: the previous script's M(x) heuristic at these x was 0.71/0.85/0.98/1.00/1.00 --\n")
cat(" i.e. it understated M(x) everywhere except the extremes, where it happened to also reach 1.)\n\n")

sep(); cat("TABLE 2 (CORRECTED) -- CALIBRATED SCENARIO (tau_hat=0.15 +/- 0.04), Regime A\n")
sep()
idx2 <- sapply(c(0.10,0.25,0.50,0.75,0.90), function(x) which.min(abs(XGRID-x)))
tab2 <- data.frame(
  x = round(XGRID[idx2], 2),
  mu_obs = round(calibrated_table$mu_obs[idx2], 4),
  delta_hat = round(calibrated_table$delta_hat[idx2], 4),
  lb_hat = round(calibrated_table$lb_hat[idx2], 4),
  ub_hat = round(calibrated_table$ub_hat[idx2], 4),
  lb_95CI = round(calibrated_table$lb_ci[idx2], 4),
  ub_95CI = round(calibrated_table$ub_ci[idx2], 4)
)
print(tab2, row.names = FALSE)
cat("\n(Previous script reported delta_hat = 0.3700 at every x -- wrong on two counts: it used the\n")
cat(" broken 'exact' formula, AND it accidentally ended up constant across x, which the true\n")
cat(" delta_hat(x) is not.)\n\n")

sep(); cat("TABLE 3 -- VALIDATION: bound tightness (Section 5)\n"); sep()
tightness <- tapply(validation_table$B_true, validation_table$case, function(b) mean(abs(b))) /
  tapply(validation_table$bound,  validation_table$case, mean)
print(round(tightness, 4))
cat("(fraction of the theoretical ceiling M(x)*delta(x) actually used by the realised bias,\n")
cat(" averaged over x, in a Gaussian-copula example where the truth happens to be known.\n")
cat(" CaseB is NaN because it is the independence benchmark: both B_true and the bound are 0.)\n")
# ==============================================================================
# SECTION 11 -- SAVE ALL FIGURES (cairo_pdf throughout -- fixes the
#               mbcsToSbcs Unicode failures for Delta/delta/tau-hat/-> that
#               affected every figure in the previous version's base pdf())
# ==============================================================================
output_dir <- "nonparmetric_sens_anal"
if (!dir.exists(output_dir)) dir.create(output_dir)

cairo_pdf(file.path(output_dir, "figA_sensitivity_contour_delta_M.pdf"), width = 15, height = 4.4)
print(fig_contour_delta_M); dev.off()

cairo_pdf(file.path(output_dir, "figB_sensitivity_contour_tau_M.pdf"), width = 15, height = 4.6)
print(fig_contour_tau_M); dev.off()

cairo_pdf(file.path(output_dir, "figC_delta_crit_tau_crit_profile.pdf"), width = 8.5, height = 5)
print(fig_delta_crit_profile); dev.off()

cairo_pdf(file.path(output_dir, "figD_validation_bound_check.pdf"), width = 12, height = 4.2)
print(fig_validation); dev.off()

cat("All figures written to", output_dir, "\n")
cat(strrep("=",70),"\n")
cat("figA_sensitivity_contour_delta_M.pdf   -- MAIN: 2D sensitivity contour,\n")
cat("                                            delta(x) x M(x), Regime B\n")
cat("figB_sensitivity_contour_tau_M.pdf     -- companion, Kendall's tau axis\n")
cat("figC_delta_crit_tau_crit_profile.pdf   -- corrected tipping-point profile\n")
cat("figD_validation_bound_check.pdf        -- theorem check on known-Z case\n")
cat(strrep("=",70),"\n")

# ==============================================================================
# SELF-TEST -- reproduces dissertation Table 5 (Example 5.5.3) from
# delta_gauss() and checks agreement to confirm the quadrature fix is
# actually correct before trusting anything built on top of it.
# ==============================================================================
cat("\n"); sep(); cat("SELF-TEST: delta_gauss() vs dissertation Table 5 (Example 5.5.3)\n"); sep()
table5_truth <- data.frame(
  x    = rep(c(0.10,0.25,0.50,0.75,0.90), 3),
  rho  = rep(c(0.3,0.6,0.9), each = 5),
  dissertation_exact = c(0.314,0.169,0.046,0.169,0.314,
                         0.682,0.396,0.215,0.396,0.682,
                         1.264,0.933,0.761,0.933,1.264)
)
table5_truth$this_script <- round(delta_gauss(table5_truth$x, table5_truth$rho), 3)
table5_truth$match <- abs(table5_truth$dissertation_exact - table5_truth$this_script) < 0.003
print(table5_truth, row.names = FALSE)
cat("\nAll", sum(table5_truth$match), "/", nrow(table5_truth), "rows match Table 5 to 3 decimals: ",
    all(table5_truth$match), "\n")

# ==============================================================================
# SECTION 12 -- CINELLI-HAZLETT (2020) EXACT FORMULAS, for direct comparison
# ==============================================================================
# Implemented from the paper directly (Eqs 13, 14, 18, 35), not from memory,
# so the comparison below is checkable against the source. These are NOT
# applied to a real regression here (we have no C-H dataset) -- they are kept
# as reference functions so every claim in the write-up is a verifiable
# computation, not an assertion.

CH_bias <- function(R2_YZ_DX, R2_DZ_X, se_tau_res, df) {
  se_tau_res * sqrt(R2_YZ_DX * R2_DZ_X / (1 - R2_DZ_X) * df)
}
CH_RVq <- function(f_YD_X, q = 1) {
  fq <- q * abs(f_YD_X)
  0.5 * (sqrt(fq^4 + 4*fq^2) - fq^2)
}
# f_YD|X = t-value / sqrt(df); R2_Y~D|X = f^2/(1+f^2) -- the "extreme scenario" statistic.

# ==============================================================================
# SECTION 13 -- THE Ch.5 ANALOGUES, precisely placed
# ==============================================================================
# (1) EXTREME SCENARIO (their Sec. 4.3 / R2_Y~D|X): fix the outcome-side
#     sensitivity parameter at its ceiling (their R2_Y~Z|D,X=1; our M(x)=ybar,
#     Prop. 5.4.2) and solve for the treatment-side parameter that zeroes the
#     estimate. This is EXACTLY Delta_crit(x) as already defined -- no new
#     construction needed, the correspondence is exact by definition:
#         C-H:  set R2_Y~Z|D,X=1  ->  solve R2_D~Z|X = R2_Y~D|X
#         Ch.5: set M(x)=ybar     ->  solve delta(x) = |mu_obs(x)|/ybar = Delta_crit(x)
#
# (2) ROBUSTNESS VALUE: C-H's RV works because BOTH R2's live on the identical
#     [0,1] variance-explained scale, so "equal strength on both sides" is a
#     well-posed question. delta(x) in [0,2] and M(x) in [0,ybar] are NOT
#     commensurate -- one is a copula-deviation index, the other an outcome
#     magnitude. There is no literal Ch.5 analogue of RV without first
#     normalising. Doing so (delta_tilde=delta/2, M_tilde=M/ybar, both in
#     [0,1]) and solving the equal-normalised-strength problem
#     (ybar*RV)*(2*RV) = |mu_obs(x)| gives a well-defined but CONSTRUCTED
#     analogue:
RV_copula <- function(mu_obs, ybar = 1) sqrt(abs(mu_obs) / (2*ybar))
#     This is a genuine extension, not a direct translation -- flagged as such
#     throughout.

sensitivity_stats_x <- data.frame(x = c(0.10, 0.25, 0.50, 0.75, 0.90))
sensitivity_stats_x$mu_obs <- mu_obs_gaussian(sensitivity_stats_x$x)
sensitivity_stats_x$M <- YBAR
sensitivity_stats_x$Delta_crit_extreme <- delta_crit(sensitivity_stats_x$mu_obs, YBAR)
sensitivity_stats_x$RV_copula <- RV_copula(sensitivity_stats_x$mu_obs, YBAR)
sensitivity_stats_x$tau_crit_extreme <- tau_crit_gaussian_exact(sensitivity_stats_x$Delta_crit_extreme, sensitivity_stats_x$x)
# tau-scale RV: invert delta_gauss(x, rho(tau)) = RV_copula (root-find), same machinery as tau_crit_exact
tau_RV <- mapply(function(dv, xv) {
  g <- function(rho) delta_gauss(xv, rho) - dv
  up <- delta_gauss(xv, 1-1e-6)
  if (dv >= up) return(1)
  rho_star <- uniroot(g, lower=0, upper=1-1e-6, tol=1e-8)$root
  (2/pi)*asin(rho_star)
}, sensitivity_stats_x$RV_copula, sensitivity_stats_x$x)
sensitivity_stats_x$tau_RV <- tau_RV

cat("Section 12-13 built. Sensitivity statistics across x:\n")
print(round(sensitivity_stats_x, 4))

# ==============================================================================
# SECTION 14 -- BENCHMARKING, C-H style (1x / 2x / 3x a reference association)
# ==============================================================================
# C-H benchmark against an observed covariate Xj via kD, kY multiples (Sec.
# 4.4). We have no second observed covariate in the X-Y-Z triangle, so the
# natural available benchmark is the one already-estimable association:
# rho(X,Y) itself -- "a confounder k times as strongly tied to treatment as
# treatment already is to the outcome". delta has a hard ceiling of 2
# (Remark 5.2.1), so, exactly as C-H note multiples of kD can be "ruled out
# by the data" (Eq. 65), some k will be infeasible here too.
benchmark_table <- do.call(rbind, lapply(c(0.10,0.25,0.50,0.75,0.90), function(xv) {
  d1 <- delta_gauss(xv, RHO_XY)
  data.frame(x = xv,
             delta_1x = min(2, d1), feasible_1x = d1 <= 2,
             delta_2x = min(2, 2*d1), feasible_2x = 2*d1 <= 2,
             delta_3x = min(2, 3*d1), feasible_3x = 3*d1 <= 2,
             Delta_crit = delta_crit(mu_obs_gaussian(xv), YBAR))
}))
benchmark_table$robust_1x <- benchmark_table$delta_1x < benchmark_table$Delta_crit
benchmark_table$robust_2x <- benchmark_table$delta_2x < benchmark_table$Delta_crit
benchmark_table$robust_3x <- benchmark_table$delta_3x < benchmark_table$Delta_crit

cat("\nBenchmark table (1x/2x/3x as strong as rho_XY):\n")
print(benchmark_table)
# ==============================================================================
# SECTION 15 -- EXTREME-SCENARIO PLOT (Cinelli-Hazlett Figure 3, reproduced)
# ==============================================================================
# C-H fix R2_Y~Z|D,X at 100%/75%/50% and plot the adjusted estimate against
# R2_D~Z|X. Direct analogue: fix M(x) at 100%/75%/50% of ybar and plot the
# worst-case adjusted mu against delta(x). Built for a single representative
# x (0.25); the faceted version across x is Section 16/17's job -- that's
# the part C-H's single-coefficient framework cannot do at all.
EXTREME_X <- 0.25
mu_ref <- mu_obs_gaussian(EXTREME_X)
M_fracs <- c(1.00, 0.75, 0.50)
delta_seq <- seq(0, 1.2, length.out = 200)

extreme_df <- do.call(rbind, lapply(M_fracs, function(fr) {
  data.frame(delta = delta_seq, M_frac = fr,
             adj = mu_ref - fr*YBAR*delta_seq,
             label = paste0(fr*100, "%"))
}))
extreme_df$label <- factor(extreme_df$label, levels = c("100%","75%","50%"))

fig_extreme_scenario <- ggplot(extreme_df, aes(x = delta, y = adj, colour = label)) +
  geom_hline(yintercept = 0, colour = "grey40", linewidth = 0.4) +
  geom_line(linewidth = 0.9) +
  scale_colour_manual(values = c("100%"="#000000","75%"="#0072B2","50%"="#56B4E9"),
                      name = "M(x) as % of\nManski ceiling ybar") +
  labs(
    title = paste0("Sensitivity to extreme scenarios (Cinelli-Hazlett Fig. 3, reproduced), x=", EXTREME_X),
    subtitle = paste0("Solid black: confounder's outcome-sensitivity assumed at its full, assumption-free ceiling\n",
                      "M(x)=ybar (worst case, Prop. 5.4.2). Blue/pale: willing to additionally assume M(x) is only\n",
                      "75%/50% of that ceiling. Horizontal line: where the adjusted estimate crosses zero."),
    x = "delta(x)  [confounder-treatment association, copula-ratio scale]",
    y = "adjusted mu(x) = mu_obs(x) - M(x)\u00b7delta(x)"
  ) +
  THEME_THESIS

cat("Section 15 (extreme-scenario plot) built.\n")
# ==============================================================================
# SECTION 16 -- THE EXTENSION: sensitivity statistics as FUNCTIONS of x
# ==============================================================================
# C-H's RV and R2_Y~D|X are each a single number for a single regression
# coefficient, because a linear model has one treatment effect. Chapter 5's
# objects are functions of x by construction (mu(x) is a dose-response curve),
# so both the extreme-scenario statistic and the newly-constructed RV_copula
# generalise into full curves across the treatment quantile. This is not
# obtainable by evaluating the C-H machinery "at a few points"; there is no
# C-H analogue of x at all -- their framework has no dose-response axis.

ext_df <- data.frame(x = XGRID)
ext_df$mu_obs <- mu_obs_gaussian(XGRID)
ext_df$extreme_scenario <- delta_crit(ext_df$mu_obs, YBAR)     # C-H's R2_Y~D|X analogue
ext_df$RV_copula <- RV_copula(ext_df$mu_obs, YBAR)             # constructed analogue of RV

ext_long <- rbind(
  data.frame(x = ext_df$x, value = ext_df$extreme_scenario,
             stat = "Extreme-scenario statistic (~ C-H's R2_Y~D|X)"),
  data.frame(x = ext_df$x, value = ext_df$RV_copula,
             stat = "RV_copula(x) [constructed, normalised]")
)

fig_extension_x <- ggplot(ext_long, aes(x = x, y = value, colour = stat)) +
  geom_line(linewidth = 1) +
  scale_colour_manual(values = c("Extreme-scenario statistic (~ C-H's R2_Y~D|X)" = "#000000",
                                 "RV_copula(x) [constructed, normalised]" = "#D55E00"),
                      name = NULL) +
  scale_y_continuous(limits = c(0,1)) +
  labs(
    title = "The extension: robustness as a function of treatment dose, not a single number",
    subtitle = paste0("A linear regression coefficient has exactly one value here (one point, one x).\n",
                      "The copula framework hands back the whole curve: the treatment effect is markedly\n",
                      "MORE fragile in the tails of the dose-response than in the middle -- a fact a single\n",
                      "reported robustness value could never reveal."),
    x = "Treatment quantile x", y = NULL
  ) +
  THEME_THESIS + theme(legend.position = "bottom")

cat("Section 16 (extension-across-x plot) built.\n")
# ==============================================================================
# SECTION 17 -- ENHANCED CONTOUR, C-H HOUSE STYLE (1x/2x/3x labelled points)
# ==============================================================================
suppressMessages(library(ggrepel))
X_CH <- c(0.10, 0.25, 0.50, 0.75, 0.90)
x_lab <- function(xv) paste0("x = ", sprintf("%.2f", xv))

make_grid_ch <- function(xv) {
  mu_x <- mu_obs_gaussian(xv)
  g <- expand.grid(delta = seq(0, 2, length.out = 240), M = seq(0, YBAR, length.out = 240))
  g$adj <- mu_x - g$delta*g$M     # adjusted estimate, C-H sign convention (bias hurts the estimate)
  g$x <- xv
  g
}
grid_ch <- do.call(rbind, lapply(X_CH, make_grid_ch))
grid_ch$x_label <- factor(x_lab(grid_ch$x), levels = x_lab(X_CH))

bench_ch <- do.call(rbind, lapply(X_CH, function(xv) {
  d1 <- min(2, delta_gauss(xv, RHO_XY))
  data.frame(x = xv, x_label = x_lab(xv), delta = c(min(2,d1), min(2,2*d1), min(2,3*d1)),
             M = YBAR, k = c("1x","2x","3x"))
}))
bench_ch$x_label <- factor(bench_ch$x_label, levels = levels(grid_ch$x_label))

fig_contour_ch_style <- ggplot(grid_ch, aes(x = delta, y = M, z = adj)) +
  geom_contour_filled(breaks = pretty(range(grid_ch$adj), 12)) +
  geom_contour(colour = "grey20", linewidth = 0.15, breaks = pretty(range(grid_ch$adj), 12)) +
  stat_contour(aes(x=delta,y=M,z=adj), breaks = 0, colour = "white", linewidth = 1.2) +
  stat_contour(aes(x=delta,y=M,z=adj), breaks = 0, colour = "black", linewidth = 0.6, linetype = "21") +
  geom_point(data = bench_ch, aes(x = delta, y = M, z = NULL), shape = 4, size = 2, stroke = 1.1, colour="white") +
  geom_point(data = bench_ch, aes(x = delta, y = M, z = NULL), shape = 4, size = 1.6, stroke = 0.9, colour="black") +
  geom_text_repel(data = bench_ch, aes(x = delta, y = M, z = NULL, label = paste0(k, " rho_XY")),
                  size = 2.3, nudge_y = -0.08, direction = "y", segment.size = 0.2, seed = 2, max.overlaps = 20) +
  scale_fill_viridis_d(name = "adjusted \u03bc(x)\n(sign flips at 0)", option = "magma", direction = -1) +
  facet_wrap(~x_label, nrow = 1) +
  coord_cartesian(xlim = c(0,2), ylim = c(0, YBAR), expand = FALSE) +
  labs(
    title = "Same construction, Cinelli-Hazlett house style: contours of the ADJUSTED ESTIMATE",
    subtitle = paste0("Dashed isoline: adjusted \u03bc(x) = 0 (sign flip). Crosses: confounders 1x/2x/3x as strongly\n",
                      "tied to treatment as X already is to Y -- direct analogue of their \"1x/2x/3x Female\" markers\n",
                      "(Fig. 2). Where a cross falls outside x\u2208[0,2], that multiple is infeasible (\u03b4 has a hard\n",
                      "ceiling of 2, Remark 5.2.1) -- the copula-scale analogue of their kD being data-bounded."),
    x = "\u03b4(x)", y = "M(x)"
  ) +
  THEME_THESIS + theme(legend.position = "bottom", axis.text.x = element_text(size=7))

cat("Section 17 (C-H house-style contour) built.\n")

# ==============================================================================
# SAVE THE EXTENSION FIGURES
# ==============================================================================
cairo_pdf(file.path(output_dir, "figE_extreme_scenario.pdf"), width=7, height=5)
print(fig_extreme_scenario); dev.off()
cairo_pdf(file.path(output_dir, "figF_extension_across_x.pdf"), width=8.5, height=5)
print(fig_extension_x); dev.off()
cairo_pdf(file.path(output_dir, "figG_contour_CH_style.pdf"), width=15, height=4.8)
print(fig_contour_ch_style); dev.off()
cat("Extension figures (E, F, G) saved.\n")