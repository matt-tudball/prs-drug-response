# ============================================================
# Trial-specific one-way sensitivity analyses
# One figure per trial, three panels per figure
#
# Requires that the main simulation script has already run and
# created:
#   - fourier_spec
#   - odyssey_spec
#   - woscops_spec
#   - pars_best
#   - simulate_once()
# ============================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(broom)
  library(survival)
  library(purrr)
})

setwd(this.path::here())
source("clinical_trial_sim.r")

# ------------------------------------------------------------
# 1. Update the shared "best-guess" parameters
# ------------------------------------------------------------
# If you settled on slightly different updated values for beta_G_B
# and beta_U_B, change them here only.

pars_best$beta_G_B <- 0.35
pars_best$beta_U_B <- 0.40

# ------------------------------------------------------------
# 2. Sensitivity-analysis controls
# ------------------------------------------------------------

n_rep_sens <- 25
sens_seed <- 20260501

# You said you already increased the simulated sample size.
# If sample_multiplier is already large in the base simulation
# script, you can keep trial specs unchanged here.
trial_specs_sens <- list(
  FOURIER = fourier_spec,
  ODYSSEY = odyssey_spec,
  WOSCOPS = woscops_spec
)

# ------------------------------------------------------------
# 3. Real-trial target contrasts
#    We plot the HR ratio:
#      HR_high / HR_comparator
#    < 1 means stronger apparent treatment benefit in high PRS.
# ------------------------------------------------------------

approx_ratio_ci <- function(hr_high, lo_high, hi_high,
                            hr_ref, lo_ref, hi_ref) {
  se_high <- (log(hi_high) - log(lo_high)) / (2 * 1.96)
  se_ref  <- (log(hi_ref)  - log(lo_ref))  / (2 * 1.96)

  log_ratio <- log(hr_high) - log(hr_ref)
  se_ratio  <- sqrt(se_high^2 + se_ref^2)

  tibble(
    observed_ratio = exp(log_ratio),
    observed_low   = exp(log_ratio - 1.96 * se_ratio),
    observed_high  = exp(log_ratio + 1.96 * se_ratio)
  )
}

real_trial_targets <- bind_rows(
  approx_ratio_ci(
    hr_high = 0.75, lo_high = 0.60, hi_high = 0.94,
    hr_ref  = 1.02, lo_ref  = 0.80, hi_ref  = 1.30
  ) %>% mutate(trial_name = "FOURIER"),
  approx_ratio_ci(
    hr_high = 0.63, lo_high = 0.46, hi_high = 0.86,
    hr_ref  = 0.87, lo_ref  = 0.78, hi_ref  = 0.98
  ) %>% mutate(trial_name = "ODYSSEY"),
  approx_ratio_ci(
    hr_high = 0.56, lo_high = 0.40, hi_high = 0.78,
    hr_ref  = 0.76, lo_ref  = 0.63, hi_ref  = 0.92
  ) %>% mutate(trial_name = "WOSCOPS")
)

# ------------------------------------------------------------
# 4. Trial-specific parameter grids
# ------------------------------------------------------------
# FOURIER and ODYSSEY:
#   vary beta_G_B, beta_U_B, beta_AU
#
# WOSCOPS:
#   vary beta_G_L, beta_U_L, beta_AU
#
# beta_AU is the parameter name in the simulation code
# (not beta_A_U).

grid_fourier_odyssey <- bind_rows(
  tibble(
    parameter = "beta_G_B",
    parameter_label = "CAD PRS effect on cardiovascular history",
    parameter_value = seq(0.05, 0.55, by = 0.05),
    reference_value = pars_best$beta_G_B
  ),
  tibble(
    parameter = "beta_U_B",
    parameter_label = "Comorbidity effect on cardiovascular history",
    parameter_value = seq(0.05, 0.60, by = 0.05),
    reference_value = pars_best$beta_U_B
  ),
  tibble(
    parameter = "beta_AU",
    parameter_label = "Treatment-response heterogeneity",
    parameter_value = seq(0.00, 0.40, by = 0.05),
    reference_value = pars_best$beta_AU
  )
) %>%
  mutate(is_reference = dplyr::near(parameter_value, reference_value))

grid_woscops <- bind_rows(
  tibble(
    parameter = "beta_G_L",
    parameter_label = "CAD PRS effect on LDL-C",
    parameter_value = seq(0.00, 0.06, by = 0.01),
    reference_value = pars_best$beta_G_L
  ),
  tibble(
    parameter = "beta_U_L",
    parameter_label = "Comorbidity effect on LDL-C",
    parameter_value = seq(0.00, 0.20, by = 0.02),
    reference_value = pars_best$beta_U_L
  ),
  tibble(
    parameter = "beta_AU",
    parameter_label = "Treatment-response heterogeneity",
    parameter_value = seq(0.00, 0.40, by = 0.05),
    reference_value = pars_best$beta_AU
  )
) %>%
  mutate(is_reference = dplyr::near(parameter_value, reference_value))

trial_grids <- list(
  FOURIER = grid_fourier_odyssey,
  ODYSSEY = grid_fourier_odyssey,
  WOSCOPS = grid_woscops
)

# ------------------------------------------------------------
# 5. Helpers
# ------------------------------------------------------------

update_pars <- function(pars, parameter, value) {
  out <- pars
  if (!parameter %in% names(out)) {
    stop("Parameter not found in pars: ", parameter)
  }
  out[[parameter]] <- value
  out
}

safe_glm_or <- function(dat, response, exposure, n_max = 200000) {
  if (!response %in% names(dat)) stop("Missing response: ", response)
  if (!exposure %in% names(dat)) stop("Missing exposure: ", exposure)

  if (length(unique(dat[[response]])) < 2) {
    return(NA_real_)
  }

  dat_sub <- dat %>%
    dplyr::select(all_of(c(response, exposure)))

  n_to_sample <- min(nrow(dat_sub), n_max)

  dat_sub <- dat_sub %>%
    slice_sample(n = n_to_sample)

  fit <- tryCatch(
    glm(
      reformulate(exposure, response = response),
      data = dat_sub,
      family = binomial()
    ),
    error = function(e) NULL
  )

  if (is.null(fit)) return(NA_real_)

  est <- coef(fit)[[exposure]]
  if (is.na(est) || !is.finite(est)) return(NA_real_)

  exp(est)
}

extract_sensitivity_metrics <- function(sim) {
  trial_spec <- sim$trial_spec
  cox_tab <- sim$analysis$cox_by_grs

  if (trial_spec$grs_scheme == "fourier_3") {
    reference_group <- "Low"
    high_group <- "High"
  } else if (trial_spec$grs_scheme == "top20_vs_rest") {
    reference_group <- "All others"
    high_group <- "High"
  } else if (trial_spec$grs_scheme == "top10_vs_rest") {
    reference_group <- "Lower"
    high_group <- "High"
  } else {
    stop("Unknown grs_scheme: ", trial_spec$grs_scheme)
  }

  hr_ref <- cox_tab %>%
    filter(as.character(GRS_group) == reference_group) %>%
    pull(HR)

  hr_high <- cox_tab %>%
    filter(as.character(GRS_group) == high_group) %>%
    pull(HR)

  if (length(hr_ref) != 1 || length(hr_high) != 1) {
    stop("Could not identify PRS comparison groups for ", trial_spec$name)
  }

  source <- sim$source %>%
    mutate(
      lipid_eligible_binary = as.integer(LDL >= trial_spec$ldl_cut)
    )

  tibble(
    trial_name = trial_spec$name,
    reference_group = reference_group,
    high_group = high_group,
    HR_reference = hr_ref,
    HR_high = hr_high,
    HR_ratio_high_vs_reference = hr_high / hr_ref,

    cor_G_U_trial = sim$analysis$summary$cor_G_U_trial,
    slope_U_on_G_trial = sim$analysis$summary$slope_U_on_G_trial,

    OR_D_per_SD_G_source = safe_glm_or(source, "D", "G"),
    OR_D_per_SD_U_source = safe_glm_or(source, "D", "U"),
    OR_LDLcut_per_SD_G_source = safe_glm_or(source, "lipid_eligible_binary", "G"),
    OR_LDLcut_per_SD_U_source = safe_glm_or(source, "lipid_eligible_binary", "U")
  )
}

trial_colours <- c(
  FOURIER = "#4C78A8",
  ODYSSEY = "#E45756",
  WOSCOPS = "#54A24B"
)

# ------------------------------------------------------------
# 6. Run one trial-specific sensitivity analysis
# ------------------------------------------------------------

run_trial_sensitivity <- function(trial_name,
                                  trial_spec,
                                  sensitivity_grid,
                                  n_rep = 25,
                                  seed = 1) {
  jobs <- tidyr::crossing(
    grid_id = seq_len(nrow(sensitivity_grid)),
    rep_id = seq_len(n_rep)
  ) %>%
    mutate(job_id = row_number())

  out <- vector("list", nrow(jobs))

  for (i in seq_len(nrow(jobs))) {
    g <- sensitivity_grid[jobs$grid_id[i], ]

    pars_i <- update_pars(
      pars = pars_best,
      parameter = g$parameter,
      value = g$parameter_value
    )

    set.seed(seed + i)

    message(
      "[", trial_name, "] ",
      i, "/", nrow(jobs), " | ",
      g$parameter, " = ", signif(g$parameter_value, 4),
      " | rep ", jobs$rep_id[i]
    )

    sim_i <- simulate_once(
      trial_spec = trial_spec,
      pars = pars_i,
      keep_source = TRUE
    )

    out[[i]] <- extract_sensitivity_metrics(sim_i) %>%
      mutate(
        parameter = g$parameter,
        parameter_label = g$parameter_label,
        parameter_value = g$parameter_value,
        reference_value = g$reference_value,
        is_reference = g$is_reference,
        rep_id = jobs$rep_id[i]
      )
  }

  bind_rows(out)
}

# ------------------------------------------------------------
# 7. Convert raw parameter values to interpretable x-axes
# ------------------------------------------------------------

add_interpretable_x <- function(dat) {
  dat %>%
    mutate(
      x_value = case_when(
        parameter == "beta_G_B" ~ OR_D_per_SD_G_source,
        parameter == "beta_U_B" ~ OR_D_per_SD_U_source,
        parameter == "beta_G_L" ~ 100 * (exp(parameter_value) - 1),
        parameter == "beta_U_L" ~ 100 * (exp(parameter_value) - 1),
        parameter == "beta_AU"  ~ exp(parameter_value),
        TRUE ~ parameter_value
      ),
      x_axis_label = case_when(
        parameter == "beta_G_B" ~ "OR of qualifying cardiovascular history\nper 1-SD higher CAD PRS",
        parameter == "beta_U_B" ~ "OR of qualifying cardiovascular history\nper 1-SD higher comorbidity",
        parameter == "beta_G_L" ~ "% higher LDL-C\nper 1-SD higher CAD PRS",
        parameter == "beta_U_L" ~ "% higher LDL-C\nper 1-SD higher comorbidity",
        parameter == "beta_AU"  ~ "Treatment HR multiplier\nper 1-SD higher comorbidity",
        TRUE ~ "Parameter value"
      )
    )
}

summarise_sensitivity <- function(raw_dat) {
  raw_dat %>%
    group_by(
      trial_name,
      parameter,
      parameter_label,
      parameter_value,
      is_reference
    ) %>%
    summarise(
      x_median = median(x_value, na.rm = TRUE),
      x_low = quantile(x_value, 0.025, na.rm = TRUE),
      x_high = quantile(x_value, 0.975, na.rm = TRUE),
      y_median = median(HR_ratio_high_vs_reference, na.rm = TRUE),
      y_low = quantile(HR_ratio_high_vs_reference, 0.025, na.rm = TRUE),
      y_high = quantile(HR_ratio_high_vs_reference, 0.975, na.rm = TRUE),
      x_axis_label = first(x_axis_label),
      .groups = "drop"
    ) %>%
    arrange(parameter, parameter_value)
}

# ------------------------------------------------------------
# 8. Plotting
# ------------------------------------------------------------

make_target_band_df <- function(trial_name, summary_dat) {
  target <- real_trial_targets %>%
    filter(trial_name == !!trial_name)

  tibble(
    parameter_label = unique(summary_dat$parameter_label),
    observed_ratio = target$observed_ratio,
    observed_low = target$observed_low,
    observed_high = target$observed_high
  )
}

plot_trial_sensitivity <- function(raw_dat, summary_dat, trial_name) {
  trial_colour <- trial_colours[[trial_name]]
  target_df <- make_target_band_df(trial_name, summary_dat)

  ggplot() +
    geom_rect(
      data = target_df,
      aes(
        xmin = -Inf, xmax = Inf,
        ymin = observed_low, ymax = observed_high
      ),
      inherit.aes = FALSE,
      fill = trial_colour,
      alpha = 0.08
    ) +
    geom_hline(
      yintercept = 1,
      linetype = "dashed",
      linewidth = 0.4,
      colour = "grey45"
    ) +
    geom_hline(
      data = target_df,
      aes(yintercept = observed_ratio),
      linetype = "dashed",
      linewidth = 0.7,
      colour = trial_colour
    ) +
    geom_point(
      data = raw_dat,
      aes(x = x_value, y = HR_ratio_high_vs_reference),
      colour = trial_colour,
      alpha = 0.18,
      size = 1.8
    ) +
    geom_ribbon(
      data = summary_dat,
      aes(
        x = x_median,
        ymin = y_low,
        ymax = y_high,
        group = 1
      ),
      fill = trial_colour,
      alpha = 0.12
    ) +
    geom_line(
      data = summary_dat,
      aes(x = x_median, y = y_median, group = 1),
      colour = trial_colour,
      linewidth = 1.0
    ) +
    geom_point(
      data = summary_dat,
      aes(x = x_median, y = y_median),
      colour = trial_colour,
      fill = "white",
      shape = 21,
      stroke = 0.9,
      size = 2.8
    ) +
    geom_point(
      data = summary_dat %>% filter(is_reference),
      aes(x = x_median, y = y_median),
      colour = trial_colour,
      fill = "white",
      shape = 21,
      stroke = 1.4,
      size = 4.0
    ) +
    facet_wrap(
      ~ parameter_label,
      scales = "free_x",
      nrow = 1
    ) +
    labs(
      title = paste0(trial_name, ": sensitivity of apparent high-PRS treatment benefit"),
      subtitle = paste0(
        "Outcome is the treatment HR ratio (high PRS / comparator). ",
        "Values below 1 indicate greater apparent benefit in the high-PRS group. ",
        "The shaded band and dashed coloured line show the published trial contrast."
      ),
      x = NULL,
      y = "Apparent treatment HR ratio\n(high PRS / comparator)"
    ) +
    theme_classic(base_size = 12) +
    theme(
      strip.background = element_blank(),
      strip.text = element_text(face = "bold"),
      plot.title = element_text(face = "bold"),
      plot.subtitle = element_text(size = 11),
      panel.grid.major.y = element_line(colour = "grey90", linewidth = 0.35),
      axis.title.y = element_text(margin = margin(r = 8))
    )
}

# ------------------------------------------------------------
# 9. Run all three trial-specific analyses
# ------------------------------------------------------------

raw_fourier <- run_trial_sensitivity(
  trial_name = "FOURIER",
  trial_spec = trial_specs_sens$FOURIER,
  sensitivity_grid = trial_grids$FOURIER,
  n_rep = n_rep_sens,
  seed = sens_seed + 1000
) %>%
  add_interpretable_x()

raw_odyssey <- run_trial_sensitivity(
  trial_name = "ODYSSEY",
  trial_spec = trial_specs_sens$ODYSSEY,
  sensitivity_grid = trial_grids$ODYSSEY,
  n_rep = n_rep_sens,
  seed = sens_seed + 2000
) %>%
  add_interpretable_x()

raw_woscops <- run_trial_sensitivity(
  trial_name = "WOSCOPS",
  trial_spec = trial_specs_sens$WOSCOPS,
  sensitivity_grid = trial_grids$WOSCOPS,
  n_rep = n_rep_sens,
  seed = sens_seed + 3000
) %>%
  add_interpretable_x()

sum_fourier <- summarise_sensitivity(raw_fourier)
sum_odyssey <- summarise_sensitivity(raw_odyssey)
sum_woscops <- summarise_sensitivity(raw_woscops)

# ------------------------------------------------------------
# 10. Save raw and summary results
# ------------------------------------------------------------

raw_all <- bind_rows(raw_fourier, raw_odyssey, raw_woscops)
sum_all <- bind_rows(sum_fourier, sum_odyssey, sum_woscops)

# write.csv(raw_all, "trial_specific_sensitivity_raw.csv", row.names = FALSE)
# write.csv(sum_all, "trial_specific_sensitivity_summary.csv", row.names = FALSE)

# ------------------------------------------------------------
# 11. Make figures
# ------------------------------------------------------------

p_fourier <- plot_trial_sensitivity(raw_fourier, sum_fourier, "FOURIER")
p_odyssey <- plot_trial_sensitivity(raw_odyssey, sum_odyssey, "ODYSSEY")
p_woscops <- plot_trial_sensitivity(raw_woscops, sum_woscops, "WOSCOPS")

print(p_fourier)
print(p_odyssey)
print(p_woscops)

ggsave(
  filename = "fourier_sensitivity_hr_ratio.png",
  plot = p_fourier,
  width = 11.5,
  height = 4.2,
  dpi = 320
)

ggsave(
  filename = "odyssey_sensitivity_hr_ratio.png",
  plot = p_odyssey,
  width = 11.5,
  height = 4.2,
  dpi = 320
)

ggsave(
  filename = "woscops_sensitivity_hr_ratio.png",
  plot = p_woscops,
  width = 11.5,
  height = 4.2,
  dpi = 320
)

# ------------------------------------------------------------
# 12. Small manuscript-friendly summary table
# ------------------------------------------------------------

reference_summary <- sum_all %>%
  filter(is_reference) %>%
  select(
    trial_name,
    parameter_label,
    x_median,
    y_median,
    y_low,
    y_high
  ) %>%
  arrange(trial_name, parameter_label)

print(reference_summary)
# write.csv(reference_summary, "trial_specific_sensitivity_reference_summary.csv", row.names = FALSE)