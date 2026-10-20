# =========================================================
# 1-period simulation with corrected MSM/IPW implementation
# - True A×PRS interaction = 0 by construction
# - Naive (no LDL0) biased by confounding
# - Naive (+ LDL0) biased by collider opening
# - MSM/IPW centred ~ 0 when:
#     * no weight truncation
#     * numerator uses mean(pA) not mean(A)
# =========================================================

library(dplyr)
library(tidyr)
library(purrr)
library(ggplot2)
library(stringr)

set.seed(20260114)
setwd(this.path::here())

# ----------------------------
# DGP
# ----------------------------
simulate_1period <- function(
  n = 10000,
  aG = 1.4, aU = 1.0, sigmaL = 0.8,
  alphaA = -0.2, gammaL = 1.0,
  y0 = 0.0,
  bA = -0.6, bG = 0.2, bL = 1.6, bU = 0.2, bAU = -0.8,
  sigmaY = 1.0
) {
  G <- rnorm(n)
  U <- rnorm(n)

  # LDL0 is a collider child of G and U
  LDL0 <- aG * G + aU * U + rnorm(n, sd = sigmaL)

  # Treatment depends only on observed LDL0 (exchangeability holds given LDL0)
  pA <- plogis(alphaA + gammaL * LDL0)
  A  <- rbinom(n, 1, pA)

  # Continuous outcome; true effect modification is by U (A×U), not by G
  # Under do(A=a), E[U]=0 and U ⟂ G, so E[Y^a | G=g] has no a×g term.
  Y <- y0 + bA*A + bG*G + bL*LDL0 + bU*U + bAU*A*U + rnorm(n, sd = sigmaY)

  tibble(G, U, LDL0, A, Y, pA)
}

# ----------------------------
# Correct IPW (true propensities; no truncation)
# ----------------------------
# Key correction: stabilisation numerator uses mean(pA), not mean(A).
# This removes numerator noise. With no truncation, consistency holds.
compute_ipw_true <- function(dat, stabilised = TRUE) {
  p_den <- dat$pA

  if (stabilised) {
    p_num <- mean(p_den)   # deterministic given L0 draw; avoids mean(A) noise
    w <- ifelse(dat$A == 1, p_num / p_den, (1 - p_num) / (1 - p_den))
  } else {
    # Unstabilised IPW
    w <- ifelse(dat$A == 1, 1 / p_den, 1 / (1 - p_den))
  }

  w
}

# ----------------------------
# Estimators: interaction A×G
# ----------------------------
fit_estimators <- function(dat, stabilised = TRUE) {

  # Naive 1: omit LDL0 (confounding via LDL0)
  m1 <- lm(Y ~ A * G, data = dat)

  # Naive 2: condition on LDL0 (collider bias via G -> LDL0 <- U, and A×U in Y)
  m2 <- lm(Y ~ A * G + LDL0, data = dat)

  # MSM/IPW: weights use LDL0, outcome model does not condition on LDL0
  w <- compute_ipw_true(dat, stabilised = stabilised)
  m3 <- lm(Y ~ A * G, data = dat, weights = w)

  tibble(
    naive_no_ldl0   = coef(m1)[["A:G"]],
    naive_plus_ldl0 = coef(m2)[["A:G"]],
    msm_ipw         = coef(m3)[["A:G"]]
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
    pivot_longer(cols = c(naive_no_ldl0, naive_plus_ldl0, msm_ipw),
                 names_to = "method", values_to = "beta_int") %>%
    mutate(method = recode(method,
                           naive_no_ldl0 = "Naive (no LDL0)",
                           naive_plus_ldl0 = "Naive (+ LDL0)",
                           msm_ipw = "MSM/IPW (true weights; no truncation)"))

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
res <- run_mc(n_sim = 1000, n = 25000, stabilised = TRUE)
sumtab <- summarise_mc(res)
print(sumtab)

# ---- Long format + MC summary (keep both, not just summary table) ----
res_long <- res %>%
  pivot_longer(cols = c(naive_no_ldl0, naive_plus_ldl0, msm_ipw),
               names_to = "method", values_to = "beta_int") %>%
  mutate(
    method = recode(method,
                    naive_no_ldl0   = "Regression not adjusting for LDL",
                    naive_plus_ldl0 = "Regression adjusting for LDL",
                    msm_ipw         = "MSM via IPW"),
    method = factor(method, levels = c("MSM via IPW",
                                      "Regression adjusting for LDL",
                                      "Regression not adjusting for LDL"))
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

# ---- Colour scheme (colour-blind friendly, high contrast) ----
cols <- c(
  "MSM via IPW" = "#1B9E77",
  "Regression adjusting for LDL" = "#D95F02",
  "Regression not adjusting for LDL" = "#7570B3",
  "MSM via IPW (restrict to LDL interval)" = "#E7298A"
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
  geom_hline(yintercept = 0, linetype = "dotted", linewidth = 0.8) +
  geom_boxplot(width = 0.55, outlier.shape = NA, alpha = 0.65) +
  geom_jitter(aes(colour = method),
              width = 0.14, height = 0,
              size = 2.2, alpha = 0.25, show.legend = FALSE) +
  geom_point(data = sumtab, aes(x = method, y = mean),
             inherit.aes = FALSE, size = 3.2, shape = 21, stroke = 0.8, fill = "white") +
  geom_errorbar(data = sumtab,
                aes(x = method, ymin = ci_low, ymax = ci_high),
                inherit.aes = FALSE, width = 0.18, linewidth = 0.9) +
  geom_text(data = sumtab,
            aes(x = method, y = y, label = label),
            inherit.aes = FALSE, size = 3.4, lineheight = 1.05) +
  scale_fill_manual(values = cols) +
  scale_colour_manual(values = cols) +
  coord_cartesian(clip = "off") +
  theme_bw(base_size = 14) +
  theme(
    legend.position = "none",
    axis.title.x = element_blank(),
    axis.text.x = element_text(angle = 15, hjust = 1),
    plot.margin = margin(8, 12, 8, 12)
  ) +
  labs(
    title = "Monte Carlo distribution of estimated GxA interaction",
    subtitle = "True GxA coefficient = 0.",
    y = "Estimated interaction coefficient (GxA)"
  )

# =========================================================
# Plot 2: Faceted histogram + density (one panel per method)
# =========================================================
p_hist <- ggplot(res_long, aes(x = beta_int, fill = method)) +
  geom_vline(xintercept = 0, linetype = "dotted", linewidth = 0.8) +
  # geom_histogram(aes(y = after_stat(density)), bins = 35, alpha = 0.50) +
  geom_density(linewidth = 1.0, aes(colour = method), show.legend = FALSE) +
  facet_wrap(~ method, ncol = 1, scales = "free_y") +
  scale_fill_manual(values = cols) +
  scale_colour_manual(values = cols) +
  theme_bw(base_size = 14) +
  theme(
    legend.position = "none",
    strip.text = element_text(face = "bold"),
    plot.margin = margin(8, 12, 8, 12)
  ) +
  labs(
    title = "Sampling distribution by estimator",
    x = "Estimated interaction coefficient (GxA)",
    y = "Density"
  )

# ---- Print plots ----
p_box
p_hist

# ---- Save (optional) ----
ggsave("mc_boxplot_interaction.png", p_box, width = 9, height = 5.5, dpi = 320)
ggsave("mc_hist_interaction.png",    p_hist, width = 9, height = 9.5, dpi = 320)

# =========================================================
# EXTENSION (MODIFIED): Baseline LDL stratification across full distribution
# Continue from here without modifying anything above.
#
# Changes vs prior version:
#   - Removes the two regression baselines (adjusting / not adjusting for LDL0)
#   - Runs 3 restricted MSM/IPW analyses over LDL0 percentiles:
#       0–33rd, 34–66th, 67–100th
#   - Box plot has 4 groups total:
#       MSM via IPW (full cohort) + 3 restricted MSMs
#   - Keeps your existing colour scheme; reuses the pink for all restricted MSMs
#     and uses alpha to distinguish strata.
# =========================================================

# ----------------------------
# Estimator wrapper for interval restriction demonstrations
# Returns the A:G coefficient for full MSM and for 3 strata.
# ----------------------------
fit_estimators_interval <- function(dat,
                                   stabilised = TRUE,
                                   ldl0_cuts = c(-Inf, Inf)) {
  # ldl0_cuts is a numeric vector of cut points of length 4:
  #   c(q0, q33, q66, q100)
  stopifnot(length(ldl0_cuts) == 4)

  w_full <- compute_ipw_true(dat, stabilised = stabilised)

  # Full cohort MSM/IPW
  m_full <- lm(Y ~ A * G, data = dat, weights = w_full)
  out <- tibble(`MSM via IPW` = coef(m_full)[["A:G"]])

  # 3 strata: [q0,q33], (q33,q66], (q66,q100]
  stratum_names <- c("MSM via IPW (LDL 0-33%)",
                     "MSM via IPW (LDL 34-66%)",
                     "MSM via IPW (LDL 67-100%)")

  for (k in 1:3) {
    lo <- ldl0_cuts[k]
    hi <- ldl0_cuts[k + 1]

    # Use open/closed to avoid double-counting at boundaries
    if (k == 1) {
      S <- (dat$LDL0 >= lo) & (dat$LDL0 <= hi)
    } else {
      S <- (dat$LDL0 > lo) & (dat$LDL0 <= hi)
    }

    if (!any(S)) stop("Interval restriction kept 0 rows for stratum ", k, "; check cut points.")

    m_k <- lm(Y ~ A * G, data = dat[S, , drop = FALSE], weights = w_full[S])

    out[[stratum_names[k]]] <- coef(m_k)[["A:G"]]
  }

  out
}

# ----------------------------
# Monte Carlo runner for stratification demo
# ----------------------------
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

# ----------------------------
# Define percentile cut points for LDL0 strata: 0–33, 34–66, 67–100
# ----------------------------
pilot <- simulate_1period(n = 200000)

# Use exact percentile cuts; treat q0/q100 as min/max for the pilot draw
q0   <- as.numeric(quantile(pilot$LDL0, 0.00, na.rm = TRUE))
q33  <- as.numeric(quantile(pilot$LDL0, 0.33, na.rm = TRUE))
q66  <- as.numeric(quantile(pilot$LDL0, 0.66, na.rm = TRUE))
q100 <- as.numeric(quantile(pilot$LDL0, 1.00, na.rm = TRUE))

ldl0_cuts_strata <- c(q0, q33, q66, q100)

# ----------------------------
# Run Monte Carlo
# ----------------------------
res_strata <- run_mc_interval(
  n_sim = 1000,
  n = 25000,
  stabilised = TRUE,
  ldl0_cuts = ldl0_cuts_strata
)

# ----------------------------
# Long format
# ----------------------------
res_strata_long <- res_strata %>%
  pivot_longer(
    cols = -sim,
    names_to = "method",
    values_to = "beta_int"
  ) %>%
  mutate(
    method = factor(
      method,
      levels = c(
        "MSM via IPW",
        "MSM via IPW (LDL 0-33%)",
        "MSM via IPW (LDL 34-66%)",
        "MSM via IPW (LDL 67-100%)"
      )
    )
  )

# ----------------------------
# Summary table
# ----------------------------
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
    label = sprintf("mean=%.3f\n95%% MC CI\n[%.3f, %.3f]", mean, ci_low, ci_high)
  )

print(sumtab_strata)

# ----------------------------
# Plot: boxplot + jitter + CI annotation
# Colour scheme: keep your `cols` mapping.
# - Full cohort uses the MSM green.
# - All strata use the same pink (since cols has only one restricted-MSM colour).
#   We distinguish strata by alpha on points and box fill (minimal, but keeps palette fixed).
# ----------------------------
cols_strata <- c(
  "MSM via IPW" = cols[["MSM via IPW"]],
  "MSM via IPW (LDL 0-33%)" = cols[["MSM via IPW (restrict to LDL interval)"]],
  "MSM via IPW (LDL 34-66%)" = cols[["MSM via IPW (restrict to LDL interval)"]],
  "MSM via IPW (LDL 67-100%)" = cols[["MSM via IPW (restrict to LDL interval)"]]
)

alpha_map <- c(
  "MSM via IPW" = 0.65,
  "MSM via IPW (LDL 0-33%)" = 0.35,
  "MSM via IPW (LDL 34-66%)" = 0.55,
  "MSM via IPW (LDL 67-100%)" = 0.75
)

# label placement
y_top <- max(res_strata_long$beta_int)
y_lab <- y_top + 0.08 * diff(range(res_strata_long$beta_int))
sumtab_strata <- sumtab_strata %>% mutate(y = y_lab)

p_box_strata <- ggplot(res_strata_long, aes(x = method, y = beta_int, fill = method)) +
  geom_hline(yintercept = 0, linetype = "dotted", linewidth = 0.8) +
  geom_boxplot(aes(alpha = method), width = 0.55, outlier.shape = NA) +
  geom_jitter(aes(colour = method, alpha = method),
              width = 0.14, height = 0, size = 2.0, show.legend = FALSE) +
  geom_point(data = sumtab_strata, aes(x = method, y = mean),
             inherit.aes = FALSE, size = 3.2, shape = 21, stroke = 0.8, fill = "white") +
  geom_errorbar(data = sumtab_strata,
                aes(x = method, ymin = ci_low, ymax = ci_high),
                inherit.aes = FALSE, width = 0.18, linewidth = 0.9) +
  geom_text(data = sumtab_strata,
            aes(x = method, y = y, label = label),
            inherit.aes = FALSE, size = 3.2, lineheight = 1.05) +
  coord_cartesian(clip = "off") +
  scale_fill_manual(values = cols_strata) +
  scale_colour_manual(values = cols_strata) +
  scale_alpha_manual(values = alpha_map) +
  theme_bw(base_size = 14) +
  theme(
    legend.position = "none",
    axis.title.x = element_blank(),
    axis.text.x = element_text(angle = 18, hjust = 1),
    plot.margin = margin(8, 12, 8, 12)
  ) +
  labs(
    title = "Effect on GxA estimate of baseline LDL stratification",
    subtitle = "True A×G = 0. Baseline LDL stratified into 0–33%, 34–66%, 67–100% bins",
    y = "Estimated interaction coefficient (A×G)"
  )

p_box_strata

# Optional save
ggsave("mc_boxplot_interaction_strata_tertiles.png", p_box_strata, width = 9, height = 7, dpi = 320)
