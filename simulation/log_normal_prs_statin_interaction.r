# =========================================================
# LOG-NORMAL / MULTIPLICATIVE DGP VERSION
#   - baseline LDL0 is strictly positive (log-normal)
#   - follow-up LDL1 is strictly positive (log-normal)
#   - causal effects operate on the multiplicative scale
#   - All regressions up use log(Y) and log(LDL0)
#
# Core causal structure is preserved:
#   G -> log(LDL0) <- U
#   A depends on log(LDL0)
#   log(Y) depends on log(LDL0), A, G, U, and A×U (but NOT A×G)
# =========================================================

library(dplyr)
library(tidyr)
library(purrr)
library(ggplot2)
library(stringr)

set.seed(20260114)
setwd(this.path::here())
source("../figures/fps_style.r")

output_dir <- "simulation-outputs"

# ----------------------------
# DGP (log-normal LDL0 and LDL1)
# ----------------------------
simulate_1period <- function(
  n = 10000,
  aG = 0.12, aU = 0.15, sigma_logL = 0.25,         # baseline log(LDL0) drivers
  alphaA = -2.5, gammaG = 0.03, gamma_logL = 1.00, # treatment depends on log(LDL0) and G
  mu_logY = 0.0,                                   # intercept for log(Y)
  bA = log(0.70),                                  # multiplicative statin effect on LDL-C (30% reduction)
  bG = 0.00,                                       # direct effect of G on log(Y) beyond log(LDL0)
  bLogL = 1.00,                                    # carryover of baseline log(LDL0) into log(Y)
  bU = 0.20,                                       # effect of U on log(Y)
  bAU = -0.25,                                     # A×U effect modification on log scale (true effect mod by U)
  sigma_logY = 0.20                                # noise on log scale
) {
  G <- rnorm(n)
  U <- rnorm(n)

  # Baseline on log scale: collider child of G and U
  logLDL0 <- aG * G + aU * U + rnorm(n, sd = sigma_logL)
  LDL0 <- exp(logLDL0)

  # Treatment depends on observed baseline log LDL0 and G (exchangeability holds given logLDL0 and G)
  pA <- plogis(alphaA + gammaG * G + gamma_logL * logLDL0)
  A  <- rbinom(n, 1, pA)

  # Follow-up LDL-C on log scale; true effect modification is by U (A×U), not by G
  logY <- mu_logY +
    bA * A +
    bG * G +
    bLogL * logLDL0 +
    bU * U +
    bAU * A * U +
    rnorm(n, sd = sigma_logY)

  Y <- exp(logY)

  tibble(G, U, logLDL0, LDL0, A, logY, Y, pA)
}

# ----------------------------
# Correct IPTW (true propensities; no truncation)
# ----------------------------
compute_IPTW_true <- function(dat, stabilised = TRUE) {
  p_den <- dat$pA

  if (stabilised) {
    p_num <- mean(p_den)
    w <- ifelse(dat$A == 1, p_num / p_den, (1 - p_num) / (1 - p_den))
  } else {
    w <- ifelse(dat$A == 1, 1 / p_den, 1 / (1 - p_den))
  }

  w
}

# ----------------------------
# Estimators: interaction A×G (all outcome/baseline on log scale)
# ----------------------------
fit_estimators <- function(dat, stabilised = TRUE) {

  # Naive 1: omit baseline logLDL0 (confounding via baseline)
  m1 <- lm(logY ~ A * G, data = dat)

  # Naive 2: condition on baseline logLDL0 (collider bias via G -> logLDL0 <- U, and A×U in logY)
  m2 <- lm(logY ~ A * G + logLDL0, data = dat)

  # MSM/IPTW: weights use baseline logLDL0, outcome model does not condition on baseline
  w <- compute_IPTW_true(dat, stabilised = stabilised)
  m3 <- lm(logY ~ A * G, data = dat, weights = w)

  tibble(
    naive_no_ldl0   = coef(m1)[["A:G"]],
    naive_plus_ldl0 = coef(m2)[["A:G"]],
    msm_IPTW         = coef(m3)[["A:G"]]
  )
}

# ----------------------------
# Monte Carlo runner
# ----------------------------
run_mc <- function(n_sim = 500, n = 10000, stabilised = TRUE, ...) {
  map_df(1:n_sim, function(s) {
    dat <- simulate_1period(n = n, ...)
    fit_estimators(dat, stabilised = stabilised) %>% mutate(sim = s)
  })
}

summarise_mc <- function(res) {
  res_long <- res %>%
    pivot_longer(cols = c(naive_no_ldl0, naive_plus_ldl0, msm_IPTW),
                 names_to = "method", values_to = "beta_int") %>%
    mutate(method = recode(method,
                           naive_no_ldl0 = "Naive (no LDL0)",
                           naive_plus_ldl0 = "Naive (+ LDL0)",
                           msm_IPTW = "MSM/IPTW (true weights; no truncation)"))

  res_long %>%
    group_by(method) %>%
    summarise(mean = mean(beta_int),
              ci_low = quantile(beta_int, 0.025),
              ci_high = quantile(beta_int, 0.975),
              sd = sd(beta_int),
              .groups = "drop")
}

# ----------------------------
# Example run
# ----------------------------
res <- run_mc(n_sim = 1000, n = 78145, stabilised = TRUE)
sumtab <- summarise_mc(res)
print(sumtab)

# ---- Long format + MC summary (keep both, not just summary table) ----
res_long <- res %>%
  pivot_longer(cols = c(naive_no_ldl0, naive_plus_ldl0, msm_IPTW),
               names_to = "method", values_to = "beta_int") %>%
  mutate(
    method = recode(method,
                    naive_no_ldl0   = "Regression not adjusting for baseline LDL-C",
                    naive_plus_ldl0 = "Regression adjusting for baseline LDL-C",
                    msm_IPTW         = "IPTW"),
    method = factor(method, levels = c("IPTW",
                                      "Regression adjusting for baseline LDL-C",
                                      "Regression not adjusting for baseline LDL-C"))
  )

sumtab <- res_long %>%
  group_by(method) %>%
  summarise(
    mean = mean(beta_int),
    ci_low = quantile(beta_int, 0.025),
    ci_high = quantile(beta_int, 0.975),
    sd = sd(beta_int),
    .groups = "drop"
  ) %>%
  mutate(
    label = sprintf("mean = %.3f\n95%% MC CI [%.3f, %.3f]", mean, ci_low, ci_high)
  )

print(sumtab)

# ---- Colour scheme (5PS) ----
cols <- c(
  "IPTW" = fps$navy,
  "Regression adjusting for baseline LDL-C" = fps$steel,
  "Regression not adjusting for baseline LDL-C" = fps$slate,
  "IPTW (stratified)" = fps$sand
)

# ---- Choose y position for CI labels above boxes ----
y_top <- max(res_long$beta_int)
y_lab <- y_top + 0.08 * diff(range(res_long$beta_int))

sumtab <- sumtab %>%
  mutate(y = y_lab)

# =========================================================
# Plot 1: Boxplot + jitter + MC CI annotation
# =========================================================
p_box <- ggplot(res_long, aes(x = method, y = beta_int, fill = method)) +
  geom_hline(yintercept = 0, linetype = "dotted", linewidth = 0.8, colour = fps$null) +
  geom_boxplot(width = 0.55, outlier.shape = NA, alpha = 0.65) +
  geom_jitter(aes(colour = method),
              width = 0.14, height = 0,
              size = 2.2, alpha = 0.25, show.legend = FALSE) +
  geom_point(data = sumtab, aes(x = method, y = mean),
             inherit.aes = FALSE, size = 3.2, shape = 21, stroke = 0.8, fill = fps$white, colour = fps$charcoal) +
  geom_errorbar(data = sumtab,
                aes(x = method, ymin = ci_low, ymax = ci_high),
                inherit.aes = FALSE, width = 0.18, linewidth = 0.9, colour = fps$charcoal) +
  geom_text(data = sumtab,
            aes(x = method, y = y, label = label),
            inherit.aes = FALSE, size = 3.4, lineheight = 1.05) +
  scale_fill_manual(values = cols) +
  scale_colour_manual(values = cols) +
  coord_cartesian(clip = "off") +
  theme_fps(base_size = 12) +
  theme(
    legend.position = "none",
    axis.title.x = element_blank(),
    axis.text.x = element_text(angle = 15, hjust = 1),
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_blank(),
    plot.background = element_rect(fill = fps$white, colour = NA),
    plot.margin = margin(8, 12, 8, 12)
  ) +
  labs(
    y = "Estimated treatment-by-PRS interaction coefficient"
  )

# =========================================================
# Plot 2: Faceted density (one panel per method)
# =========================================================
p_hist <- ggplot(res_long, aes(x = beta_int, fill = method)) +
  geom_vline(xintercept = 0, linetype = "dotted", linewidth = 0.8, colour = fps$null) +
  geom_density(linewidth = 1.0, aes(colour = method), show.legend = FALSE) +
  facet_wrap(~ method, ncol = 1, scales = "free_y") +
  scale_fill_manual(values = cols) +
  scale_colour_manual(values = cols) +
  theme_fps(base_size = 12) +
  theme(
    legend.position = "none",
    strip.text = element_text(face = "bold"),
    panel.grid.minor = element_blank(),
    plot.background = element_rect(fill = fps$white, colour = NA),
    plot.margin = margin(8, 12, 8, 12)
  ) +
  labs(
    x = "Estimated treatment-by-PRS interaction coefficient",
    y = "Density"
  )

# ---- Print plots ----
if (interactive()) print(p_box)
if (interactive()) print(p_hist)

if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)
ggsave(file.path(output_dir, "mc_boxplot_interaction.png"), p_box,  width = 9, height = 5.5, dpi = 320, device = ragg::agg_png)
ggsave(file.path(output_dir, "mc_hist_interaction.png"),    p_hist, width = 9, height = 9.5, dpi = 320, device = ragg::agg_png)

# =========================================================
# EXTENSION: Baseline LDL-C stratification across full distribution
#
#   This stratification is on baseline log(LDL0), not LDL0.
#   Outcome model uses logY.
# =========================================================

fit_estimators_interval <- function(dat,
                                   stabilised = TRUE,
                                   ldl0_cuts = c(-Inf, Inf)) {
  stopifnot(length(ldl0_cuts) == 4)

  w_full <- compute_IPTW_true(dat, stabilised = stabilised)

  # Full cohort MSM/IPTW
  m_full <- lm(logY ~ A * G, data = dat, weights = w_full)
  out <- tibble(`IPTW` = coef(m_full)[["A:G"]])

  stratum_names <- c("IPTW, bottom LDL-C tertile",
                     "IPTW, middle LDL-C tertile",
                     "IPTW, top LDL-C tertile")

  for (k in 1:3) {
    lo <- ldl0_cuts[k]
    hi <- ldl0_cuts[k + 1]

    if (k == 1) {
      S <- (dat$logLDL0 >= lo) & (dat$logLDL0 <= hi)
    } else {
      S <- (dat$logLDL0 > lo) & (dat$logLDL0 <= hi)
    }

    if (!any(S)) stop("Interval restriction kept 0 rows for stratum ", k, "; check cut points.")

    m_k <- lm(logY ~ A * G, data = dat[S, , drop = FALSE], weights = w_full[S])
    out[[stratum_names[k]]] <- coef(m_k)[["A:G"]]
  }

  out
}

run_mc_interval <- function(n_sim = 500, n = 25000,
                            stabilised = TRUE,
                            ldl0_cuts = c(-Inf, Inf),
                            ...) {
  purrr::map_df(seq_len(n_sim), function(s) {
    dat <- simulate_1period(n = n, ...)
    fit_estimators_interval(dat,
                            stabilised = stabilised,
                            ldl0_cuts = ldl0_cuts) %>%
      mutate(sim = s)
  })
}

# Define percentile cut points for baseline log(LDL0) strata
pilot <- simulate_1period(n = 500000)
q0   <- as.numeric(quantile(pilot$logLDL0, 0.00, na.rm = TRUE))
q33  <- as.numeric(quantile(pilot$logLDL0, 0.33, na.rm = TRUE))
q66  <- as.numeric(quantile(pilot$logLDL0, 0.66, na.rm = TRUE))
q100 <- as.numeric(quantile(pilot$logLDL0, 1.00, na.rm = TRUE))
ldl0_cuts_strata <- c(q0, q33, q66, q100)

res_strata <- run_mc_interval(
  n_sim = 1000,
  n = 78145,
  stabilised = TRUE,
  ldl0_cuts = ldl0_cuts_strata
)

# x-axis labels, split across lines to match the UKB/HRN meta-analysis figure
estimator_map_strata <- c(
  "IPTW"                       = "Full\ncohort",
  "IPTW, bottom LDL-C tertile" = "Bottom\nbaseline LDL\ntertile",
  "IPTW, middle LDL-C tertile" = "Middle\nbaseline LDL\ntertile",
  "IPTW, top LDL-C tertile"    = "Top\nbaseline LDL\ntertile"
)

res_strata_long <- res_strata %>%
  pivot_longer(
    cols = -sim,
    names_to = "method",
    values_to = "beta_int"
  ) %>%
  mutate(
    method = factor(
      unname(estimator_map_strata[method]),
      levels = unname(estimator_map_strata)
    )
  )

# Two-sided MC p-value from a normal approximation to the sampling distribution
# (MC mean / MC SD), labelled as in the UKB/HRN meta-analysis figure
sumtab_strata <- res_strata_long %>%
  group_by(method) %>%
  summarise(
    mean = mean(beta_int),
    ci_low = quantile(beta_int, 0.025),
    ci_high = quantile(beta_int, 0.975),
    sd = sd(beta_int),
    .groups = "drop"
  ) %>%
  mutate(
    p = 2 * pnorm(-abs(mean / sd)),
    label = ifelse(p < 0.001, "p < 0.001", sprintf("p = %.3f", p))
  )

print(sumtab_strata)

cols_strata <- setNames(
  c(cols[["IPTW"]], rep(cols[["IPTW (stratified)"]], 3)),
  levels(res_strata_long$method)
)

alpha_map_strata <- setNames(
  c(0.65, 0.75, 0.75, 0.75),
  levels(res_strata_long$method)
)

y_top <- max(res_strata_long$beta_int)
y_lab <- y_top + 0.10 * diff(range(res_strata_long$beta_int))
sumtab_strata <- sumtab_strata %>% mutate(y = y_lab)

theme_pub <- theme_fps(base_size = 18) +
  theme(
    panel.grid.minor = element_blank(),
    panel.grid.major.y = element_blank(),
    legend.title = element_blank(),
    legend.position = "top",
    plot.title.position = "plot",
    plot.title = element_text(face = "bold"),
    plot.subtitle = element_text(colour = fps$slate),
    axis.title = element_text(face = "bold")
  )

p_box_strata <- ggplot(res_strata_long, aes(x = method, y = beta_int, fill = method)) +
  geom_hline(yintercept = 0, linetype = "dotted", linewidth = 0.8, colour = fps$null) +
  geom_boxplot(aes(alpha = method), width = 0.55, outlier.shape = NA) +
  geom_jitter(aes(colour = method, alpha = method),
              width = 0.14, height = 0, size = 2.0, show.legend = FALSE) +
  geom_point(data = sumtab_strata, aes(x = method, y = mean),
             inherit.aes = FALSE, size = 3.2, shape = 21, stroke = 0.8, fill = fps$white, colour = fps$charcoal) +
  geom_errorbar(data = sumtab_strata,
                aes(x = method, ymin = ci_low, ymax = ci_high),
                inherit.aes = FALSE, width = 0.18, linewidth = 0.9, colour = fps$charcoal) +
  geom_text(data = sumtab_strata,
            aes(x = method, y = y, label = label),
            inherit.aes = FALSE, colour = fps$charcoal, size = 5.5, lineheight = 1.05) +
  coord_cartesian(clip = "off") +
  scale_fill_manual(values = cols_strata) +
  scale_colour_manual(values = cols_strata) +
  scale_alpha_manual(values = alpha_map_strata) +
  theme_pub +
  theme(
    legend.position = "none",
    axis.title.x = element_blank(),
    axis.text.x = element_text(colour = fps$charcoal, lineheight = 1.05,
                               margin = margin(t = 6)),
    axis.text.y = element_text(colour = fps$charcoal),
    axis.title.y = element_text(size = 17, lineheight = 1.1,
                                margin = margin(r = 10)),
    plot.margin = margin(12, 16, 8, 12)
  ) +
  labs(
    y = "Change in statin effect on log LDL-C\nper SD of LDL PRS in simulation"
  )

if (interactive()) print(p_box_strata)

ggsave(file.path(output_dir, "mc_boxplot_interaction_strata_tertiles.png"), p_box_strata, width = 8, height = 6, dpi = 300, device = ragg::agg_png, bg = fps$white)

# =========================================================
# EXTENSION: Prevalent-user / change-score analyses (A = 1 only)
#
# Adds three estimators in treated (A=1) only:
#   1) Absolute change on natural scale:        (Y - LDL0) ~ G
#   2) Relative change on multiplicative scale: log(Y/LDL0) ~ G
#   3) Baseline-adjusted relative change:       log(Y/LDL0) ~ G + logLDL0
#
# Notes:
# - (2) is well-behaved because Y and LDL0 are strictly positive in the log-normal DGP.
# - (3) adjusts for baseline on the same scale that governs treatment assignment.
# - None of these estimate an A×G interaction; they are "treated-only response phenotype"
#   associations. This is deliberate for illustrating prevalent-user/change-score pitfalls.
# =========================================================

# ----------------------------
# Treated-only change-score estimators
# ----------------------------
fit_change_score_estimators <- function(dat) {
  dat1 <- dat %>% filter(A == 1)
  if (nrow(dat1) < 100) stop("Too few treated individuals for change-score analyses; check parameters.")

  # outcomes
  dat1 <- dat1 %>%
    mutate(
      abs_change = Y - LDL0,
      log_ratio  = logY - logLDL0   # = log(Y/LDL0)
    )

  m_abs   <- lm(abs_change ~ G, data = dat1)
  m_rel   <- lm(log_ratio  ~ G, data = dat1)
  m_relBL <- lm(log_ratio ~ G + logLDL0, data = dat1)

  tibble(
    `Absolute change in LDL-C`                  = coef(m_abs)[["G"]],
    `Relative change in LDL-C`                  = coef(m_rel)[["G"]],
    `Relative change in LDL-C, baseline-adjusted` = coef(m_relBL)[["G"]]
  )
}

# ----------------------------
# Monte Carlo runner for change-score designs
# ----------------------------
run_mc_change_scores <- function(n_sim = 500, n = 25000, ...) {
  map_df(seq_len(n_sim), function(s) {
    dat <- simulate_1period(n = n, ...)
    fit_change_score_estimators(dat) %>% mutate(sim = s)
  })
}

# ----------------------------
# Run MC
# ----------------------------
res_change <- run_mc_change_scores(n_sim = 1000, n = 78145)

# x-axis labels, matching the UKB/HRN change-score meta-analysis figure
estimator_map_change <- c(
  "Absolute change in LDL-C"                    = "Absolute\nchange in\nLDL-C",
  "Relative change in LDL-C"                    = "Unadjusted\nrelative change\nin LDL-C",
  "Relative change in LDL-C, baseline-adjusted" = "Adjusted\nrelative change\nin LDL-C"
)

# Long format
res_change_long <- res_change %>%
  pivot_longer(cols = -sim, names_to = "method", values_to = "beta") %>%
  mutate(
    method = factor(
      unname(estimator_map_change[method]),
      levels = unname(estimator_map_change)
    )
  )

# Summary table; two-sided MC p-value from a normal approximation to the
# sampling distribution (MC mean / MC SD), as in the strata figure
sumtab_change <- res_change_long %>%
  group_by(method) %>%
  summarise(
    mean = mean(beta),
    ci_low = quantile(beta, 0.025),
    ci_high = quantile(beta, 0.975),
    sd = sd(beta),
    .groups = "drop"
  ) %>%
  mutate(
    p = 2 * pnorm(-abs(mean / sd)),
    label = ifelse(p < 0.001, "p < 0.001", sprintf("p = %.3f", p))
  )

print(sumtab_change)

# ----------------------------
# Plot: boxplot + jitter + MC p-value annotation (same style as strata figure)
# ----------------------------
cols_change <- setNames(
  rep(cols[["IPTW (stratified)"]], 3),
  levels(res_change_long$method)
)

alpha_change <- setNames(
  rep(0.75, 3),
  levels(res_change_long$method)
)

y_top_c <- max(res_change_long$beta)
y_lab_c <- y_top_c + 0.10 * diff(range(res_change_long$beta))
sumtab_change <- sumtab_change %>% mutate(y = y_lab_c)

p_box_change <- ggplot(res_change_long, aes(x = method, y = beta, fill = method)) +
  geom_hline(yintercept = 0, linetype = "dotted", linewidth = 0.8, colour = fps$null) +
  geom_boxplot(aes(alpha = method), width = 0.55, outlier.shape = NA) +
  geom_jitter(aes(colour = method, alpha = method),
              width = 0.14, height = 0, size = 2.0, show.legend = FALSE) +
  geom_point(data = sumtab_change, aes(x = method, y = mean),
             inherit.aes = FALSE, size = 3.2, shape = 21, stroke = 0.8, fill = fps$white, colour = fps$charcoal) +
  geom_errorbar(data = sumtab_change,
                aes(x = method, ymin = ci_low, ymax = ci_high),
                inherit.aes = FALSE, width = 0.18, linewidth = 0.9, colour = fps$charcoal) +
  geom_text(data = sumtab_change,
            aes(x = method, y = y, label = label),
            inherit.aes = FALSE, colour = fps$charcoal, size = 5.5, lineheight = 1.05) +
  coord_cartesian(clip = "off") +
  scale_fill_manual(values = cols_change) +
  scale_colour_manual(values = cols_change) +
  scale_alpha_manual(values = alpha_change) +
  theme_pub +
  theme(
    legend.position = "none",
    axis.title.x = element_blank(),
    axis.text.x = element_text(colour = fps$charcoal, lineheight = 1.05,
                               margin = margin(t = 6)),
    axis.text.y = element_text(colour = fps$charcoal),
    axis.title.y = element_text(size = 17, lineheight = 1.1,
                                margin = margin(r = 10)),
    plot.margin = margin(12, 16, 8, 12)
  ) +
  labs(
    y = "Change in LDL-C after statin initiation\nper SD of LDL PRS"
  )

if (interactive()) print(p_box_change)

ggsave(file.path(output_dir, "mc_boxplot_change_score_designs.png"), p_box_change, width = 8, height = 6, dpi = 300, device = ragg::agg_png, bg = fps$white)

# =========================================================
# EXTENSION: Baseline LDL-C stratification across full distribution
# with randomised A
#
#   This stratification is on baseline log(LDL0), not LDL0.
#   Outcome model uses logY.
# =============================================

fit_estimators_interval_randomized <- function(dat,
                                   stabilised = TRUE,
                                   ldl0_cuts = c(-Inf, Inf)) {
  stopifnot(length(ldl0_cuts) == 4)

  w_full <- compute_IPTW_true(dat, stabilised = stabilised)

  # Full cohort MSM/IPTW
  m_full <- lm(logY ~ A * G, data = dat, weights = w_full)
  out <- tibble(`IPTW` = coef(m_full)[["A:G"]])

  stratum_names <- c("IPTW, bottom LDL-C tertile",
                     "IPTW, middle LDL-C tertile",
                     "IPTW, top LDL-C tertile")

  for (k in 1:3) {
    lo <- ldl0_cuts[k]
    hi <- ldl0_cuts[k + 1]

    if (k == 1) {
      S <- (dat$logLDL0 >= lo) & (dat$logLDL0 <= hi)
    } else {
      S <- (dat$logLDL0 > lo) & (dat$logLDL0 <= hi)
    }

    if (!any(S)) stop("Interval restriction kept 0 rows for stratum ", k, "; check cut points.")

    m_k <- lm(logY ~ A * G, data = dat[S, , drop = FALSE], weights = w_full[S])
    out[[stratum_names[k]]] <- coef(m_k)[["A:G"]]
  }

  out
}

run_mc_interval_randomized <- function(n_sim = 500, n = 25000,
                            stabilised = TRUE,
                            ldl0_cuts = c(-Inf, Inf),
                            ...) {
  purrr::map_df(seq_len(n_sim), function(s) {
    dat <- simulate_1period(n = n, gammaG = 0, gamma_logL = 0, ...)
    fit_estimators_interval_randomized(dat,
                            stabilised = stabilised,
                            ldl0_cuts = ldl0_cuts) %>%
      mutate(sim = s)
  })
}

# Define percentile cut points for baseline log(LDL0) strata
pilot <- simulate_1period(n = 500000)
q0   <- as.numeric(quantile(pilot$logLDL0, 0.00, na.rm = TRUE))
q33  <- as.numeric(quantile(pilot$logLDL0, 0.33, na.rm = TRUE))
q66  <- as.numeric(quantile(pilot$logLDL0, 0.66, na.rm = TRUE))
q100 <- as.numeric(quantile(pilot$logLDL0, 1.00, na.rm = TRUE))
ldl0_cuts_rand <- c(q0, q33, q66, q100)

res_rand <- run_mc_interval_randomized(
  n_sim = 1000,
  n = 78145,
  stabilised = TRUE,
  ldl0_cuts = ldl0_cuts_rand
)

res_rand_long <- res_rand %>%
  pivot_longer(
    cols = -sim,
    names_to = "method",
    values_to = "beta_int"
  ) %>%
  mutate(
    method = factor(
      method,
      levels = c(
        "IPTW",
        "IPTW, bottom LDL-C tertile",
        "IPTW, middle LDL-C tertile",
        "IPTW, top LDL-C tertile"
      )
    )
  )

sumtab_rand <- res_rand_long %>%
  group_by(method) %>%
  summarise(
    mean = mean(beta_int),
    ci_low = quantile(beta_int, 0.025),
    ci_high = quantile(beta_int, 0.975),
    sd = sd(beta_int),
    .groups = "drop"
  ) %>%
  mutate(
    label = sprintf("mean=%.3f\n95%% MC CI\n[%.3f, %.3f]", mean, ci_low, ci_high)
  )

print(sumtab_rand)

cols_rand <- c(
  "IPTW"                       = cols[["IPTW"]],
  "IPTW, bottom LDL-C tertile" = cols[["IPTW (stratified)"]],
  "IPTW, middle LDL-C tertile" = cols[["IPTW (stratified)"]],
  "IPTW, top LDL-C tertile"    = cols[["IPTW (stratified)"]]
)

alpha_map <- c(
  "IPTW"                       = 0.65,
  "IPTW, bottom LDL-C tertile" = 0.35,
  "IPTW, middle LDL-C tertile" = 0.55,
  "IPTW, top LDL-C tertile"    = 0.75
)

y_top <- max(res_rand_long$beta_int)
y_lab <- y_top + 0.08 * diff(range(res_rand_long$beta_int))
sumtab_rand <- sumtab_rand %>% mutate(y = y_lab)

p_box_rand <- ggplot(res_rand_long, aes(x = method, y = beta_int, fill = method)) +
  geom_hline(yintercept = 0, linetype = "dotted", linewidth = 0.8, colour = fps$null) +
  geom_boxplot(aes(alpha = method), width = 0.55, outlier.shape = NA) +
  geom_jitter(aes(colour = method, alpha = method),
              width = 0.14, height = 0, size = 2.0, show.legend = FALSE) +
  geom_point(data = sumtab_rand, aes(x = method, y = mean),
             inherit.aes = FALSE, size = 3.2, shape = 21, stroke = 0.8, fill = fps$white, colour = fps$charcoal) +
  geom_errorbar(data = sumtab_rand,
                aes(x = method, ymin = ci_low, ymax = ci_high),
                inherit.aes = FALSE, width = 0.18, linewidth = 0.9, colour = fps$charcoal) +
  geom_text(data = sumtab_rand,
            aes(x = method, y = y, label = label),
            inherit.aes = FALSE, size = 3.2, lineheight = 1.05) +
  coord_cartesian(clip = "off") +
  scale_fill_manual(values = cols_rand) +
  scale_colour_manual(values = cols_rand) +
  scale_alpha_manual(values = alpha_map) +
  theme_fps(base_size = 12) +
  theme(
    legend.position = "none",
    axis.title.x = element_blank(),
    axis.text.x = element_text(angle = 18, hjust = 1),
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_blank(),
    plot.background = element_rect(fill = fps$white, colour = NA),
    plot.margin = margin(8, 12, 8, 12)
  ) +
  labs(
    y = "Estimated treatment-by-PRS interaction coefficient"
  )

if (interactive()) print(p_box_rand)

ggsave(file.path(output_dir, "mc_boxplot_interaction_rand_tertiles.png"), p_box_rand, width = 9, height = 7, dpi = 320, device = ragg::agg_png)