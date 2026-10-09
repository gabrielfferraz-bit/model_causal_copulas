# ==============================================================================
# Sensitivity Analysis - Illustrative Numerical Examples
# Figures  + numerical tables
#
# Run in the R console:
#   source("Sensitivity_Analysis_Illustrative_Numerical_Examples_EXACT.R")
#
# The Figure geometry below follows the manuscript construction:
#   green region  : delta < Delta_crit
#   orange region : delta >= Delta_crit
#   blue lines    : mu_obs(0.01) +/- M(0.01) delta
#   red dotted    : Delta_crit
#   grey dashed   : FGM model-implied delta(beta)
#   black bar     : interval at the FGM model-implied delta
#   black diamond : true interventional mean mu(0.01) = 0.106
#
# ==============================================================================

# ------------------------------------------------------------------------------
# 0. Packages and output directory
# ------------------------------------------------------------------------------

if (!requireNamespace("ggplot2", quietly = TRUE)) {
  install.packages("ggplot2", repos = "https://cloud.r-project.org")
}

if (!requireNamespace("gridExtra", quietly = TRUE)) {
  install.packages("gridExtra", repos = "https://cloud.r-project.org")
}

library(ggplot2)
library(gridExtra)
library(grid)

output_dir <- file.path(getwd(), "sensitivity_example_output")
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

# ------------------------------------------------------------------------------
# 1. Exact formulas from the manuscript
# ------------------------------------------------------------------------------

mu <- function(x) {
  0.6 * x + 0.1
}

mu_obs <- function(x, beta) {
  mu(x) - beta * (1 - 2 * x) / 6
}

bias <- function(x, beta) {
  mu_obs(x, beta) - mu(x)
}

M <- function(x) {
  0.6 * x + 0.6
}

delta <- function(x, beta) {
  abs(beta * (1 - 2 * x)) / 2
}

Delta_crit <- function(x, beta) {
  abs(mu_obs(x, beta)) / M(x)
}

interval_lower <- function(x, beta, d = delta(x, beta)) {
  mu_obs(x, beta) - M(x) * d
}

interval_upper <- function(x, beta, d = delta(x, beta)) {
  mu_obs(x, beta) + M(x) * d
}

# ------------------------------------------------------------------------------
# 2. Numerical checks for Examples 7-10
# ------------------------------------------------------------------------------

x_ref <- 0.01
beta_ref <- c(0.2, 0.8)

check_table <- data.frame(
  beta = beta_ref,
  x = x_ref,
  mu = mu(x_ref),
  mu_obs = mu_obs(x_ref, beta_ref),
  B = bias(x_ref, beta_ref),
  abs_B = abs(bias(x_ref, beta_ref)),
  M = M(x_ref),
  delta = delta(x_ref, beta_ref),
  Delta_crit = Delta_crit(x_ref, beta_ref),
  lower = interval_lower(x_ref, beta_ref),
  upper = interval_upper(x_ref, beta_ref),
  width = interval_upper(x_ref, beta_ref) -
    interval_lower(x_ref, beta_ref),
  sign_guaranteed =
    delta(x_ref, beta_ref) < Delta_crit(x_ref, beta_ref),
  mu_inside =
    interval_lower(x_ref, beta_ref) <= mu(x_ref) &
    mu(x_ref) <= interval_upper(x_ref, beta_ref)
)

write.csv(
  check_table,
  file.path(output_dir, "table_example7_exact_values.csv"),
  row.names = FALSE
)

cat("\n", strrep("=", 88), "\n", sep = "")
cat("EXACT NUMERICAL CHECK: x = 0.01\n")
cat(strrep("=", 88), "\n")
print(check_table, row.names = FALSE)

# ------------------------------------------------------------------------------
# 3. Full beta sweep
# ------------------------------------------------------------------------------

beta_grid <- seq(-1, 1, by = 0.01)

beta_table <- data.frame(
  beta = beta_grid,
  x = x_ref,
  mu = mu(x_ref),
  mu_obs = mu_obs(x_ref, beta_grid),
  B = bias(x_ref, beta_grid),
  abs_B = abs(bias(x_ref, beta_grid)),
  M = M(x_ref),
  delta = delta(x_ref, beta_grid),
  Delta_crit = Delta_crit(x_ref, beta_grid),
  bound_radius = M(x_ref) * delta(x_ref, beta_grid),
  lower = interval_lower(x_ref, beta_grid),
  upper = interval_upper(x_ref, beta_grid),
  width = interval_upper(x_ref, beta_grid) -
    interval_lower(x_ref, beta_grid),
  sign_guaranteed =
    delta(x_ref, beta_grid) < Delta_crit(x_ref, beta_grid),
  mu_inside =
    interval_lower(x_ref, beta_grid) <= mu(x_ref) &
    mu(x_ref) <= interval_upper(x_ref, beta_grid)
)

write.csv(
  beta_table,
  file.path(output_dir, "table_beta_sweep_x001.csv"),
  row.names = FALSE
)

beta_selected <- c(
  -1, -0.8, -0.6, -0.4, -0.2,
   0,  0.2,  0.4,  0.6,  0.8, 1
)

selected_index <- match(
  round(beta_selected, 2),
  round(beta_table$beta, 2)
)

selected_beta_table <- beta_table[selected_index, ]

write.csv(
  selected_beta_table,
  file.path(output_dir, "table_beta_selected_x001.csv"),
  row.names = FALSE
)

cat("\n", strrep("=", 88), "\n", sep = "")
cat("SELECTED BETA VALUES: x = 0.01\n")
cat(strrep("=", 88), "\n")
print(selected_beta_table, row.names = FALSE)

# ------------------------------------------------------------------------------
# 4. Critical beta values
# ------------------------------------------------------------------------------

root_fun <- function(b) {
  delta(x_ref, b) - Delta_crit(x_ref, b)
}

root_neg <- uniroot(root_fun, c(-1, -1e-8))$root
root_pos <- uniroot(root_fun, c(1e-8, 1))$root

beta_zero_muobs <- uniroot(
  function(b) mu_obs(x_ref, b),
  c(0, 1)
)$root

critical_table <- data.frame(
  quantity = c(
    "negative beta: delta = Delta_crit",
    "positive beta: delta = Delta_crit",
    "beta: mu_obs = 0"
  ),
  value = c(
    root_neg,
    root_pos,
    beta_zero_muobs
  )
)

write.csv(
  critical_table,
  file.path(output_dir, "table_beta_critical_values_x001.csv"),
  row.names = FALSE
)

cat("\n", strrep("=", 88), "\n", sep = "")
cat("CRITICAL BETA VALUES: x = 0.01\n")
cat(strrep("=", 88), "\n")
print(critical_table, row.names = FALSE)

cat(
  sprintf(
    "\nAt x = %.2f, sign preservation is guaranteed for %.6f < beta < %.6f.\n",
    x_ref, root_neg, root_pos
  )
)

cat(
  sprintf(
    "mu_obs(0.01) crosses zero at beta = %.6f.\n",
    beta_zero_muobs
  )
)

# ------------------------------------------------------------------------------
# 5. x-profile table
# ------------------------------------------------------------------------------

x_selected <- c(0.01, 0.10, 0.25, 0.50, 0.75, 0.90, 0.99)
beta_x <- c(-0.8, -0.2, 0.2, 0.8)

x_table <- expand.grid(
  x = x_selected,
  beta = beta_x,
  KEEP.OUT.ATTRS = FALSE,
  stringsAsFactors = FALSE
)

x_table$mu <- mu(x_table$x)
x_table$mu_obs <- mu_obs(x_table$x, x_table$beta)
x_table$B <- bias(x_table$x, x_table$beta)
x_table$abs_B <- abs(x_table$B)
x_table$M <- M(x_table$x)
x_table$delta <- delta(x_table$x, x_table$beta)
x_table$Delta_crit <- Delta_crit(x_table$x, x_table$beta)
x_table$bound_radius <- x_table$M * x_table$delta
x_table$lower <- interval_lower(x_table$x, x_table$beta)
x_table$upper <- interval_upper(x_table$x, x_table$beta)
x_table$width <- x_table$upper - x_table$lower
x_table$sign_guaranteed <-
  x_table$delta < x_table$Delta_crit
x_table$mu_inside <-
  x_table$lower <= x_table$mu &
  x_table$mu <= x_table$upper

x_table <- x_table[order(x_table$beta, x_table$x), ]

write.csv(
  x_table,
  file.path(output_dir, "table_x_profiles.csv"),
  row.names = FALSE
)

# ------------------------------------------------------------------------------
# 6. Figure 2: exact manuscript geometry, implemented in ggplot2
# ------------------------------------------------------------------------------

# Colors used in the manuscript figure
COL_GREEN <- "#009E73"
COL_ORANGE <- "#E69F00"
COL_VERM <- "#D55E00"
COL_BLUE <- "#005493"

make_fig2_panel <- function(beta_value, panel_letter) {

  x <- x_ref

  obs <- mu_obs(x, beta_value)
  true_mu <- mu(x)
  m <- M(x)

  d_model <- delta(x, beta_value)
  d_crit <- Delta_crit(x, beta_value)

  lower_model <- interval_lower(
    x, beta_value, d = d_model
  )

  upper_model <- interval_upper(
    x, beta_value, d = d_model
  )

  # The manuscript uses the two straight bounds over [0, 0.5].
  d <- seq(0, 0.5, length.out = 1000)

  upper <- obs + m * d
  lower <- obs - m * d

  # --------------------------------------------------------------------------
  # Exact colored regions from the TikZ construction
  # --------------------------------------------------------------------------

  # Green region:
  # apex at delta = 0, y = mu_obs
  # and vertical segment at Delta_crit from 0 to 2*mu_obs,
  # with the ordering automatically reversed when mu_obs < 0.
  green_poly <- data.frame(
    delta = c(
      0,
      d_crit,
      d_crit
    ),
    y = c(
      obs,
      max(0, 2 * obs),
      min(0, 2 * obs)
    )
  )

  # Orange region:
  # begins at Delta_crit and extends to delta = 0.5.
  orange_poly <- data.frame(
    delta = c(
      d_crit,
      d_crit,
      0.5,
      0.5
    ),
    y = c(
      min(0, 2 * obs),
      max(0, 2 * obs),
      obs + m * 0.5,
      obs - m * 0.5
    )
  )

  band <- data.frame(
    delta = d,
    lower = lower,
    upper = upper
  )

  # --------------------------------------------------------------------------
  # Figure
  # --------------------------------------------------------------------------

  ggplot() +

    # Green = sign guaranteed
    geom_polygon(
      data = green_poly,
      aes(x = delta, y = y),
      fill = COL_GREEN,
      alpha = 0.30,
      colour = NA
    ) +

    # Orange = interval reaches zero
    geom_polygon(
      data = orange_poly,
      aes(x = delta, y = y),
      fill = COL_ORANGE,
      alpha = 0.30,
      colour = NA
    ) +

    # Zero line
    geom_hline(
      yintercept = 0,
      linetype = "dashed",
      linewidth = 0.55,
      colour = "black"
    ) +

    # Blue sensitivity bounds
    geom_line(
      data = band,
      aes(x = delta, y = upper),
      colour = COL_BLUE,
      linewidth = 1.05
    ) +

    geom_line(
      data = band,
      aes(x = delta, y = lower),
      colour = COL_BLUE,
      linewidth = 1.05
    ) +

    # mu_obs point at delta = 0
    geom_point(
      aes(x = 0, y = obs),
      shape = 16,
      size = 3.1,
      colour = COL_BLUE
    ) +

    # Critical value: red dotted line downward from zero
    geom_segment(
      aes(
        x = d_crit,
        xend = d_crit,
        y = 0,
        yend = -0.27
      ),
      linetype = "dotted",
      linewidth = 0.75,
      colour = COL_VERM
    ) +

    # Critical point
    geom_point(
      aes(x = d_crit, y = 0),
      shape = 21,
      size = 2.4,
      fill = "white",
      colour = COL_VERM,
      stroke = 0.95
    ) +

    # Critical-value annotation
    # Use mathematical Delta_{crit}(0.01) notation and anchor the label
    # at the critical-value line so the panel (b) label is never clipped.
    annotate(
      "text",
      x = d_crit,
      y = -0.285,
      label = sprintf(
        "Delta[crit](0.01) == %.4f",
        d_crit
      ),
      parse = TRUE,
      colour = COL_VERM,
      size = 3.25,
      hjust = 0,
      vjust = 1
    ) +

    # Model-implied FGM delta: grey dashed line
    geom_segment(
      aes(
        x = d_model,
        xend = d_model,
        y = -0.36,
        yend = 0.33
      ),
      linetype = "dashed",
      linewidth = 0.50,
      colour = "black"
    ) +

    # Model-implied delta annotation
    annotate(
      "text",
      x = d_model,
      y = 0.345,
      label = sprintf(
        "delta(0.01) = %.3f",
        d_model
      ),
      size = 3.5,
      hjust = 0.5,
      vjust = 0
    ) +

    # Sensitivity interval at the FGM value
    geom_errorbar(
      aes(
        x = d_model,
        ymin = lower_model,
        ymax = upper_model
      ),
      width = 0,
      linewidth = 1.25,
      colour = "black"
    ) +

    # True interventional mean
    geom_point(
      aes(
        x = d_model,
        y = true_mu
      ),
      shape = 23,
      size = 3.6,
      fill = "black",
      colour = "white",
      stroke = 0.45
    ) +

    # Exact axis range and breaks from the manuscript
    scale_x_continuous(
      limits = c(0, 0.5),
      breaks = seq(0, 0.5, by = 0.1),
      labels = function(z) sprintf("%.1f", z),
      expand = c(0, 0)
    ) +

    scale_y_continuous(
      limits = c(-0.36, 0.42),
      breaks = seq(-0.3, 0.4, by = 0.1),
      labels = function(z) sprintf("%.1f", z),
      expand = c(0, 0)
    ) +

    labs(
      title = sprintf(
        "(%s) beta = %.1f: mu_obs(0.01) = %.4f",
        panel_letter,
        beta_value,
        obs
      ),
      subtitle = if (
        d_model < d_crit
      ) {
        "delta < Delta_crit: sign guaranteed"
      } else {
        "delta > Delta_crit: sign not guaranteed"
      },
      x = "hypothesised delta(0.01)",
      y = "bound on mu(0.01)"
    ) +

    theme_bw(base_size = 10) +

    theme(
      panel.grid.minor = element_blank(),
      panel.grid.major = element_line(
        colour = "grey90",
        linewidth = 0.25
      ),
      panel.border = element_rect(
        colour = "black",
        linewidth = 0.45,
        fill = NA
      ),
      plot.title = element_text(
        size = 12.0,
        face = "plain",
        hjust = 0
      ),
      plot.subtitle = element_text(
        size = 11.0,
        hjust = 0,
        margin = margin(b = 3)
      ),
      axis.title = element_text(size = 12.5),
      axis.text = element_text(size = 11.0),
      plot.margin = margin(4, 4, 4, 4)
    )
}

fig2_a <- make_fig2_panel(
  beta_value = 0.2,
  panel_letter = "a"
)

fig2_b <- make_fig2_panel(
  beta_value = 0.8,
  panel_letter = "b"
)

# ------------------------------------------------------------------------------
# 7. Manuscript-style shared legend
#
# The previous legend used fixed text widths that were too narrow for the
# longer entries.  Here each legend entry is built as a small two-column grob
# and the six entries are arranged in a clean 3 x 2 grid, matching the layout
# of the manuscript while keeping every label readable.
# ------------------------------------------------------------------------------

legend_item <- function(icon, label, icon_width = 9, text_size = 9.0) {
  arrangeGrob(
    icon,
    textGrob(
      label,
      x = 0,
      hjust = 0,
      gp = gpar(fontsize = text_size)
    ),
    ncol = 2,
    widths = unit.c(
      unit(icon_width, "mm"),
      unit(1, "null")
    )
  )
}

icon_green <- rectGrob(
  width = unit(8, "mm"),
  height = unit(3.5, "mm"),
  gp = gpar(fill = COL_GREEN, col = NA, alpha = 0.30)
)

icon_orange <- rectGrob(
  width = unit(8, "mm"),
  height = unit(3.5, "mm"),
  gp = gpar(fill = COL_ORANGE, col = NA, alpha = 0.30)
)

icon_blue_line <- segmentsGrob(
  x0 = unit(0, "npc"),
  x1 = unit(8, "mm"),
  y0 = unit(0.5, "npc"),
  y1 = unit(0.5, "npc"),
  gp = gpar(col = COL_BLUE, lwd = 1.6 * 72 / 25.4)
)

icon_blue_point <- pointsGrob(
  x = 0.5,
  y = 0.5,
  pch = 16,
  size = unit(2.1, "mm"),
  gp = gpar(col = COL_BLUE)
)

icon_black_bar <- segmentsGrob(
  x0 = unit(0.5, "npc"),
  x1 = unit(0.5, "npc"),
  y0 = unit(0.12, "npc"),
  y1 = unit(0.88, "npc"),
  gp = gpar(col = "black", lwd = 1.8 * 72 / 25.4)
)

icon_black_diamond <- polygonGrob(
  x = c(0.5, 0.83, 0.5, 0.17),
  y = c(0.84, 0.5, 0.16, 0.5),
  gp = gpar(fill = "black", col = "white", lwd = 0.5)
)

legend_items <- list(
  legend_item(
    icon_green,
    "delta < Delta[crit]: sign guaranteed"
  ),
  legend_item(
    icon_orange,
    "delta >= Delta[crit]: interval reaches 0"
  ),
  legend_item(
    icon_blue_line,
    "mu[obs](0.01) %+-% M(0.01) * delta",
    icon_width = 12
  ),
  legend_item(
    icon_blue_point,
    "mu[obs](0.01)"
  ),
  legend_item(
    icon_black_bar,
    "interval at the FGM distribution's delta(0.01)",
    icon_width = 12
  ),
  legend_item(
    icon_black_diamond,
    "true mu(0.01) = 0.106"
  )
)

shared_legend <- arrangeGrob(
  legend_items[[1]], legend_items[[2]], legend_items[[3]],
  legend_items[[4]], legend_items[[5]], legend_items[[6]],
  ncol = 3,
  nrow = 2,
  widths = unit.c(
    unit(1, "null"),
    unit(1, "null"),
    unit(1, "null")
  ),
  heights = unit.c(
    unit(7.5, "mm"),
    unit(7.5, "mm")
  ),
  padding = unit(0, "pt")
)

# ------------------------------------------------------------------------------
# 8. Assemble Figure 2
#
# IMPORTANT:
# The previous version used a large outer canvas with a relatively short
# internal gtable, which made the figure look vertically compressed.
# Here the two panels are allowed to fill the plotting area and the legend
# receives only the space it actually needs.
# ------------------------------------------------------------------------------

fig2 <- gridExtra::arrangeGrob(
  fig2_a,
  fig2_b,
  ncol = 2,
  widths = unit.c(
    unit(1, "null"),
    unit(1, "null")
  ),
  heights = unit(1, "null")
)

fig2_complete <- gridExtra::arrangeGrob(
  fig2,
  shared_legend,
  ncol = 1,
  heights = unit.c(
    unit(4.85, "in"),
    unit(0.78, "in")
  ),
  padding = unit(0, "pt")
)

# ------------------------------------------------------------------------------
# 9. Save Figure 2
#
# The output is intentionally wide and tall enough for the labels to remain
# readable after insertion in the manuscript.
# ------------------------------------------------------------------------------

FIG2_WIDTH  <- 11.0
FIG2_HEIGHT <- 5.7

ggsave(
  filename = file.path(
    output_dir,
    "Figure2_exact_R_native_READABLE.pdf"
  ),
  plot = fig2_complete,
  width = FIG2_WIDTH,
  height = FIG2_HEIGHT,
  units = "in",
  device = cairo_pdf,
  bg = "white"
)

ggsave(
  filename = file.path(
    output_dir,
    "Figure2_exact_R_native_READABLE.png"
  ),
  plot = fig2_complete,
  width = FIG2_WIDTH,
  height = FIG2_HEIGHT,
  units = "in",
  dpi = 600,
  bg = "white"
)

# ------------------------------------------------------------------------------
# 10. Companion beta-sweep figure
# ------------------------------------------------------------------------------

beta_plot <- beta_table

p_beta_mean <- ggplot(
  beta_plot,
  aes(x = beta)
) +
  geom_rect(
    xmin = root_neg,
    xmax = root_pos,
    ymin = -Inf,
    ymax = Inf,
    fill = COL_GREEN,
    alpha = 0.10,
    inherit.aes = FALSE
  ) +
  geom_hline(
    yintercept = 0,
    linetype = "dashed",
    linewidth = 0.5
  ) +
  geom_ribbon(
    aes(
      ymin = lower,
      ymax = upper
    ),
    fill = COL_BLUE,
    alpha = 0.13
  ) +
  geom_line(
    aes(y = lower),
    colour = COL_BLUE,
    linewidth = 0.8
  ) +
  geom_line(
    aes(y = upper),
    colour = COL_BLUE,
    linewidth = 0.8
  ) +
  geom_line(
    aes(y = mu_obs),
    colour = COL_VERM,
    linewidth = 1.0
  ) +
  geom_hline(
    yintercept = mu(x_ref),
    linetype = "dotdash",
    linewidth = 0.8,
    colour = "black"
  ) +
  geom_vline(
    xintercept = c(root_neg, root_pos),
    linetype = "dotted",
    linewidth = 0.7,
    colour = COL_VERM
  ) +
  labs(
    title = "Beta sweep at x = 0.01",
    subtitle = "Sensitivity interval evaluated at the model-implied delta(beta)",
    x = "beta",
    y = "mu or sensitivity interval"
  ) +
  theme_bw(base_size = 10.5) +
  theme(
    panel.grid.minor = element_blank()
  )

p_beta_delta <- ggplot(
  beta_plot,
  aes(x = beta)
) +
  geom_rect(
    xmin = root_neg,
    xmax = root_pos,
    ymin = -Inf,
    ymax = Inf,
    fill = COL_GREEN,
    alpha = 0.10,
    inherit.aes = FALSE
  ) +
  geom_line(
    aes(y = delta),
    colour = COL_BLUE,
    linewidth = 1.0
  ) +
  geom_line(
    aes(y = Delta_crit),
    colour = COL_VERM,
    linewidth = 1.0
  ) +
  geom_vline(
    xintercept = c(root_neg, root_pos),
    linetype = "dotted",
    linewidth = 0.7,
    colour = COL_VERM
  ) +
  labs(
    title = "delta(beta) versus Delta_crit(beta)",
    subtitle = "Green band: delta(beta) < Delta_crit(beta)",
    x = "beta",
    y = "delta or Delta_crit"
  ) +
  theme_bw(base_size = 10.5) +
  theme(
    panel.grid.minor = element_blank()
  )

fig_beta_sweep <- arrangeGrob(
  p_beta_mean,
  p_beta_delta,
  ncol = 2,
  top = textGrob(
    "Companion analysis: changing beta at x = 0.01",
    gp = gpar(
      fontsize = 12,
      fontface = "bold"
    )
  )
)

ggsave(
  filename = file.path(
    output_dir,
    "Figure2_beta_sweep_x001.pdf"
  ),
  plot = fig_beta_sweep,
  width = 8.6,
  height = 4.6,
  units = "in",
  device = cairo_pdf
)

ggsave(
  filename = file.path(
    output_dir,
    "Figure2_beta_sweep_x001.png"
  ),
  plot = fig_beta_sweep,
  width = 8.6,
  height = 4.6,
  units = "in",
  dpi = 600,
  bg = "white"
)

# ------------------------------------------------------------------------------
# 11. x-profile figure
# ------------------------------------------------------------------------------

x_grid <- seq(0.01, 0.99, length.out = 500)
beta_profile <- c(-0.8, -0.2, 0.2, 0.8)

profile <- expand.grid(
  x = x_grid,
  beta = beta_profile,
  KEEP.OUT.ATTRS = FALSE,
  stringsAsFactors = FALSE
)

profile$B <- bias(profile$x, profile$beta)
profile$abs_B <- abs(profile$B)
profile$delta <- delta(profile$x, profile$beta)
profile$Delta_crit <- Delta_crit(profile$x, profile$beta)

profile$beta_label <- factor(
  profile$beta,
  levels = beta_profile,
  labels = sprintf("beta = %.1f", beta_profile)
)

p_profile_bias <- ggplot(
  profile,
  aes(
    x = x,
    y = B,
    colour = beta_label
  )
) +
  geom_hline(
    yintercept = 0,
    linetype = "dashed",
    linewidth = 0.5
  ) +
  geom_line(linewidth = 0.9) +
  geom_vline(
    xintercept = 0.5,
    linetype = "dotted",
    linewidth = 0.5
  ) +
  labs(
    title = "Confounding bias across treatment values",
    x = "x",
    y = "B(x) = mu_obs(x) - mu(x)",
    colour = NULL
  ) +
  theme_bw(base_size = 10.5) +
  theme(
    panel.grid.minor = element_blank()
  )

p_profile_delta <- ggplot(
  profile,
  aes(
    x = x,
    y = delta,
    colour = beta_label
  )
) +
  geom_line(linewidth = 0.9) +
  geom_vline(
    xintercept = 0.5,
    linetype = "dotted",
    linewidth = 0.5
  ) +
  labs(
    title = "Copula-ratio sensitivity index",
    x = "x",
    y = "delta(x) = |beta| |1 - 2x| / 2",
    colour = NULL
  ) +
  theme_bw(base_size = 10.5) +
  theme(
    panel.grid.minor = element_blank()
  )

fig_x_profiles <- arrangeGrob(
  p_profile_bias,
  p_profile_delta,
  ncol = 2,
  top = textGrob(
    "Supporting profiles for Example 7",
    gp = gpar(
      fontsize = 12,
      fontface = "bold"
    )
  )
)

ggsave(
  filename = file.path(
    output_dir,
    "Figure2_x_profiles.pdf"
  ),
  plot = fig_x_profiles,
  width = 8.3,
  height = 4.5,
  units = "in",
  device = cairo_pdf
)

ggsave(
  filename = file.path(
    output_dir,
    "Figure2_x_profiles.png"
  ),
  plot = fig_x_profiles,
  width = 8.3,
  height = 4.5,
  units = "in",
  dpi = 600,
  bg = "white"
)

# ------------------------------------------------------------------------------
# 12. Final console report
# ------------------------------------------------------------------------------

cat("\n", strrep("=", 88), "\n", sep = "")
cat("FINAL SUMMARY FOR EXAMPLES 7-10\n")
cat(strrep("=", 88), "\n")

for (b in c(0.2, 0.8)) {

  obs <- mu_obs(x_ref, b)
  d <- delta(x_ref, b)
  dc <- Delta_crit(x_ref, b)
  lo <- interval_lower(x_ref, b)
  hi <- interval_upper(x_ref, b)

  cat(
    sprintf(
      paste0(
        "beta = %4.1f | mu = %.4f | mu_obs = %.4f | ",
        "B = %.4f | delta = %.4f | Delta_crit = %.4f | ",
        "interval = [%.4f, %.4f] | sign guaranteed = %s\n"
      ),
      b,
      mu(x_ref),
      obs,
      bias(x_ref, b),
      d,
      dc,
      lo,
      hi,
      ifelse(d < dc, "YES", "NO")
    )
  )
}

cat("\nFiles written to:\n")
cat("  ", file.path(output_dir, "Figure2_exact_R_native_READABLE.pdf"), "\n", sep = "")
cat("  ", file.path(output_dir, "Figure2_exact_R_native_READABLE.png"), "\n", sep = "")
cat("  ", file.path(output_dir, "Figure2_beta_sweep_x001.pdf"), "\n", sep = "")
cat("  ", file.path(output_dir, "Figure2_beta_sweep_x001.png"), "\n", sep = "")
cat("  ", file.path(output_dir, "Figure2_x_profiles.pdf"), "\n", sep = "")
cat("  ", file.path(output_dir, "Figure2_x_profiles.png"), "\n", sep = "")
cat("  ", file.path(output_dir, "table_example7_exact_values.csv"), "\n", sep = "")
cat("  ", file.path(output_dir, "table_beta_sweep_x001.csv"), "\n", sep = "")
cat("  ", file.path(output_dir, "table_beta_selected_x001.csv"), "\n", sep = "")
cat("  ", file.path(output_dir, "table_beta_critical_values_x001.csv"), "\n", sep = "")
cat("  ", file.path(output_dir, "table_x_profiles.csv"), "\n", sep = "")
cat(strrep("=", 88), "\n", sep = "")
