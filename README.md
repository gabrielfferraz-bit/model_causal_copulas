# model_causal_copulas

**Paper 2026 — González-López, Justus & Ferraz**

R code for a copula-based framework for structural causal modelling. The repository contains four scripts covering analytical examples, sensitivity analysis, and an empirical application to NHANES data.

---

## Repository structure

### `Illustrative_Numerical_Examples.R`

Reproduces the paper's numerical examples and figures for two copula families under a structural causal model with a single observed confounder.

**What it does:**
- Implements closed-form and quadrature-based formulas for the interventional mean μ(x) = E[Y | do(X = x)], the Average Causal Effect ACE(x) = dμ/dx, the interventional standard deviation, and the standardised ACE ACE_std(x) = ACE(x) / σ(x) for both the **Gaussian copula** (Sections 5.1 of the paper) and the **FGM copula** (Section 5.2).
- Contrasts the two families across three confounding regimes (Cases A–C) covering weak, moderate, and strong confounding.
- Generates Figures 3–6 and Tables 1–6, 8–10 from the paper, including a cross-family comparison at matched Kendall's τ (Table 8).
- Highlights the key qualitative difference between the families: the Gaussian ACE exhibits *tail amplification* (effect magnitude increases toward the tails of the treatment distribution), while the FGM ACE is constant across treatment levels.

**Dependencies:** `pracma` (Gauss–Hermite quadrature).

---

### `Sensitivity_Analysis_-_Illustrative_Numerical_Examples.R`

Companion to the illustrative examples, implementing a nonparametric sensitivity analysis for the case in which the confounder Z is latent.

**What it does:**
- Defines the core sensitivity parameter δ(x): the L¹ copula deviation measuring the degree of dependence between treatment and the latent confounder.
- Implements three epistemic regimes: (i) *oracle case* — Z is observed and the bias bound is exact; (ii) *partial knowledge* — the researcher asserts δ(x) ≤ δ_max for several threshold values; (iii) *calibrated knowledge* — an external estimate of Kendall's τ (plus standard error) is mapped to δ(x) via closed-form formulas for each copula family.
- Computes tipping-point curves Δ_crit(x) = |μ_obs(x)| / M(x), which give the minimum confounding strength that would overturn the sign of the estimated causal effect.
- Produces Figures 1–4 and corresponding tables comparing oracle and sensitivity intervals across copula families, confounding strengths, and knowledge levels.
- Connects to the Cinelli & Hazlett (2020) and VanderWeele & Ding (2017) sensitivity frameworks.

**Dependencies:** `ggplot2`, `pracma`, `gridExtra`.

---

### `Empirical_Illustration_NHANES_functions.R`

Function library sourced by the pipeline script. Implements the core theoretical results as reusable, numerically stable R functions.

**What it contains:**
- **Proposition 4.1 (Interventional CDF):** numerical integration of the interventional conditional CDF F^{(X↓)}_{Y|X}(y|x) over the confounder distribution.
- **Theorem 5.1 (Gaussian closed form):** analytical expressions for the interventional mean, ACE, and standardised ACE under Gaussian copulas.
- **Corollary 7.1 (Observational comparison):** computation of the observational regression slope for comparison against the causal ACE.
- **Generic copula wrapper:** BIC-based copula family selection (Gaussian, Clayton, Gumbel, Frank, Joe, Student-t) and the corresponding interventional computations.
- **Empirical copula:** nonparametric implementation requiring no parametric family assumption.
- **Z-varying conditional copula:** relaxes the simplifying assumption of a constant conditional copula, estimating C_{X,Y|Z} as a function of z.
- **Distributional inference functions:** computation of interventional quantiles and distributional causal effects.
- **Bootstrap helpers:** functions called by the pipeline for confidence interval construction.

---

### `Empirical_Illustration_NHANES_Pipeline.R`

End-to-end analysis pipeline applying the copula causal framework to NHANES 2017–2018 survey data. The causal question is the effect of **dietary quality** (Healthy Eating Index, HEI-2015) on **glycated haemoglobin** (HbA1c), with **BMI** and **income-to-poverty ratio** as potential confounders.

**Pipeline stages:**
1. **Data preparation** — downloads and merges four NHANES modules (DEMO_J, GHB_J, BMX_J, DR1TOT_J); constructs HEI-2015 scores.
2. **Causal discovery** — runs three structure-learning algorithms (PC, Hill-Climbing, RESIT) and takes their consensus to identify the confounder set.
3. **Copula estimation** — fits bivariate copulas for each pair (X, Y), (X, Z), (Y, Z) with BIC-based family selection using `rvinecopulib`.
4. **Dual causal analysis** — applies both the Gaussian closed-form (Theorem 5.1) and the generic copula estimator; compares causal ACE against the observational regression slope.
5. **Variance decomposition** — partitions variability in the outcome into direct and confounded components.
6. **Sensitivity analysis** — evaluates robustness of the ACE estimate to the inclusion of BMI as a confounder and to the simplifying copula assumption.
7. **Bootstrap inference** — constructs pointwise and uniform confidence bands for ACE(x) and μ(x) via nonparametric bootstrap.
8. **Distributional inference** — estimates interventional quantile functions and distributional causal effects with survey weights.
9. **Visualisation and tables** — produces publication-ready figures and summary tables.

**Dependencies:** `tidyverse`, `nhanesA`, `rvinecopulib`, `VineCopula`, `pcalg`, `bnlearn`, `dagitty`, `ggdag`, `boot`, `grf`, `mgcv`, `survey`, `heiscore`, and others (auto-installed if missing).

---

## Replication

Run `Illustrative_Numerical_Examples.R` and `Sensitivity_Analysis_-_Illustrative_Numerical_Examples.R` independently — they require no external data and produce all numerical and graphical outputs inline.

For the empirical application, source `Empirical_Illustration_NHANES_functions.R` first, then run `Empirical_Illustration_NHANES_Pipeline.R`. An `outputs/` directory is created automatically for saved plots and tables.
