library(ggplot2)
library(dplyr)
library(this.path)

setwd(this.path::here())
source("fps_style.r")

# ---- read cohort-level interaction estimates ----
cohort_dat <- read.csv("ukb_hrn_interaction_estimates.csv")

analysis_levels <- c(
  "Full cohort",
  "LDL tertile 1, 0-33%",
  "LDL tertile 2, 34-66%",
  "LDL tertile 3, 67-100%"
)

# x-axis labels, split across lines so they fit at a large font size
estimator_map <- c(
  "Full cohort"            = "Full\ncohort",
  "LDL tertile 1, 0-33%"   = "Bottom\nbaseline LDL\ntertile",
  "LDL tertile 2, 34-66%"  = "Middle\nbaseline LDL\ntertile",
  "LDL tertile 3, 67-100%" = "Top\nbaseline LDL\ntertile"
)

stopifnot(all(cohort_dat$analysis %in% analysis_levels))

# ---- fixed-effect inverse-variance meta-analysis per analysis ----
meta_analyse <- function(cohort_dat, analysis_levels) {
  cohort_dat %>%
    mutate(w = 1 / se^2) %>%
    group_by(analysis) %>%
    summarise(
      n_cohorts = n(),
      pooled    = sum(w * estimate) / sum(w),
      pooled_se = sqrt(1 / sum(w)),
      # Cochran's Q for heterogeneity across cohorts
      Q         = sum(w * (estimate - pooled)^2),
      .groups = "drop"
    ) %>%
    transmute(
      analysis  = factor(analysis, levels = analysis_levels),
      n_cohorts,
      estimate  = pooled,
      se        = pooled_se,
      lower_ci  = estimate - qnorm(0.975) * se,
      upper_ci  = estimate + qnorm(0.975) * se,
      p         = 2 * pnorm(-abs(estimate / se)),
      Q,
      Q_df      = n_cohorts - 1,
      Q_p       = pchisq(Q, df = Q_df, lower.tail = FALSE),
      I2        = pmax(0, (Q - Q_df) / Q)
    ) %>%
    arrange(analysis)
}

meta_dat <- meta_analyse(cohort_dat, analysis_levels)

write.csv(meta_dat, "ukb_hrn_meta_analysis_results.csv", row.names = FALSE)
print(meta_dat)

# ---- plot (same style as UKB-only stratified forest) ----
plot_dat <- meta_dat %>%
  mutate(
    estimator = factor(
      unname(estimator_map[as.character(analysis)]),
      levels = unname(estimator_map)
    ),
    label = ifelse(p < 0.001, "p < 0.001", sprintf("p = %.3f", p))
  ) %>%
  arrange(estimator)

cols_strata <- setNames(
  c(fps$navy, fps$sand, fps$sand, fps$sand),
  levels(plot_dat$estimator)
)

alpha_map <- setNames(
  c(1, 1, 1, 1),
  levels(plot_dat$estimator)
)

y_rng <- range(c(plot_dat$lower_ci, plot_dat$upper_ci), na.rm = TRUE)
y_lab <- max(plot_dat$upper_ci, na.rm = TRUE) + 0.10 * diff(y_rng)

plot_dat <- plot_dat %>%
  mutate(y = y_lab)

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

p_box_strata <- ggplot(plot_dat, aes(x = estimator, y = estimate, colour = estimator)) +
  geom_hline(yintercept = 0, linetype = "dotted", linewidth = 0.8, colour = fps$null) +
  geom_errorbar(
    aes(ymin = lower_ci, ymax = upper_ci, alpha = estimator),
    width = 0.18,
    linewidth = 1.2,
    show.legend = FALSE
  ) +
  geom_point(
    aes(fill = estimator, alpha = estimator),
    size = 4.5,
    shape = 21,
    stroke = 1,
    colour = fps$charcoal,
    show.legend = FALSE
  ) +
  geom_text(
    aes(y = y, label = label),
    colour = fps$charcoal,
    family = fps_font,
    size = 5.5,
    lineheight = 1.05,
    show.legend = FALSE
  ) +
  coord_cartesian(clip = "off") +
  scale_fill_manual(values = cols_strata) +
  scale_colour_manual(values = cols_strata) +
  scale_alpha_manual(values = alpha_map) +
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
    y = "Change in statin effect on log LDL-C\nper SD of LDL PRS"
  )

ggsave(
  plot = p_box_strata,
  width = 8,
  height = 6,
  dpi = 300,
  device = ragg::agg_png,
  filename = "statin_ldl_stratified_forest_ukb_hrn_meta.png"
)

# =========================================================
# Change-score (initiator-only) designs
# =========================================================

# ---- read cohort-level change-score estimates ----
change_dat <- read.csv("ukb_hrn_change_score_estimates.csv")

change_levels <- c(
  "Absolute change in LDL-C",
  "Unadjusted relative change in LDL-C",
  "Adjusted relative change in LDL-C"
)

# x-axis labels, split across lines so they fit at a large font size
change_map <- c(
  "Absolute change in LDL-C"            = "Absolute\nchange in\nLDL-C",
  "Unadjusted relative change in LDL-C" = "Unadjusted\nrelative change\nin LDL-C",
  "Adjusted relative change in LDL-C"   = "Adjusted\nrelative change\nin LDL-C"
)

stopifnot(all(change_dat$analysis %in% change_levels))

meta_change <- meta_analyse(change_dat, change_levels) %>%
  left_join(
    change_dat %>%
      group_by(analysis) %>%
      summarise(n_total = sum(n), .groups = "drop") %>%
      mutate(analysis = factor(analysis, levels = change_levels)),
    by = "analysis"
  )

write.csv(meta_change, "ukb_hrn_change_score_meta_analysis_results.csv",
          row.names = FALSE)
print(meta_change)

# ---- plot (same style as stratified forest) ----
plot_dat_change <- meta_change %>%
  mutate(
    estimator = factor(
      unname(change_map[as.character(analysis)]),
      levels = unname(change_map)
    ),
    label = ifelse(p < 0.001, "p < 0.001", sprintf("p = %.3f", p))
  ) %>%
  arrange(estimator)

cols_change <- setNames(
  rep(fps$sand, 3),
  levels(plot_dat_change$estimator)
)

alpha_change <- setNames(
  rep(1, 3),
  levels(plot_dat_change$estimator)
)

y_rng_c <- range(c(plot_dat_change$lower_ci, plot_dat_change$upper_ci),
                 na.rm = TRUE)
y_lab_c <- max(plot_dat_change$upper_ci, na.rm = TRUE) + 0.10 * diff(y_rng_c)

plot_dat_change <- plot_dat_change %>%
  mutate(y = y_lab_c)

p_box_change <- ggplot(plot_dat_change,
                       aes(x = estimator, y = estimate, colour = estimator)) +
  geom_hline(yintercept = 0, linetype = "dotted", linewidth = 0.8, colour = fps$null) +
  geom_errorbar(
    aes(ymin = lower_ci, ymax = upper_ci, alpha = estimator),
    width = 0.18,
    linewidth = 1.2,
    show.legend = FALSE
  ) +
  geom_point(
    aes(fill = estimator, alpha = estimator),
    size = 4.5,
    shape = 21,
    stroke = 1,
    colour = fps$charcoal,
    show.legend = FALSE
  ) +
  geom_text(
    aes(y = y, label = label),
    colour = fps$charcoal,
    family = fps_font,
    size = 5.5,
    lineheight = 1.05,
    show.legend = FALSE
  ) +
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

ggsave(
  plot = p_box_change,
  width = 8,
  height = 6,
  dpi = 300,
  device = ragg::agg_png,
  filename = "statin_ldl_change_score_forest_ukb_hrn_meta.png"
)
