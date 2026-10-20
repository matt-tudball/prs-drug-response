# ============================================================
# Unified FOURIER / WOSCOPS / ODYSSEY parsimonious simulation
#
# Shared causal structure:
#   G -> LDL
#   U -> LDL
#
#   G -> latent atherosclerotic burden B
#   U -> latent atherosclerotic burden B
#   LDL -> latent atherosclerotic burden B
#
#   D = 1(B exceeds threshold), where D is established ASCVD / prior event.
#
#   Eligibility conditions on:
#     trial-specific LDL-C threshold
#     trial-specific disease-history rule:
#       FOURIER:  D == 1
#       WOSCOPS:  D == 0
#       ODYSSEY:  D == 1
#
#   Treatment A is randomised after eligibility.
#
#   Outcome contains A x U treatment-response heterogeneity,
#   but NO A x G term.
#
# Therefore:
#   Any estimated A x G effect is artefactual.
#
# Note on event-rate parameterisation:
#   We collapse alpha_Y and followup_y into baseline_cumhaz.
#   Since the simulation only models cumulative event occurrence over
#   a fixed horizon, exp(alpha_Y) * followup_y is the identified quantity.
# ============================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(purrr)
  library(ggplot2)
  library(broom)
  library(survival)
})

set.seed(20260425)
setwd(this.path::here())

# ----------------------------
# Helpers
# ----------------------------

std <- function(x) {
  as.numeric(scale(x))
}

safe_exp <- function(x) {
  exp(pmin(pmax(x, -20), 20))
}

safe_cor <- function(x, y) {
  if (length(x) < 2 || length(y) < 2) return(NA_real_)
  cor(x, y, use = "complete.obs")
}

rr <- function(r_treat, r_control) {
  ifelse(r_control == 0, NA_real_, r_treat / r_control)
}

fmt_p <- function(p) {
  ifelse(
    is.na(p), "NA",
    ifelse(p < 0.001, "<0.001", sprintf("%.3f", p))
  )
}

clean_group_name <- function(x) {
  make.names(as.character(x))
}

# ----------------------------
# Shared GRS grouping function
# ----------------------------

assign_grs_group <- function(g, trial_spec) {
  if (trial_spec$grs_scheme == "fourier_3") {
    qs <- quantile(g, probs = c(0.20, 0.80), na.rm = TRUE)

    out <- case_when(
      g <= qs[1] ~ "Low",
      g >  qs[2] ~ "High",
      TRUE       ~ "Intermediate"
    )

    factor(out, levels = c("Low", "Intermediate", "High"))

  } else if (trial_spec$grs_scheme == "top20_vs_rest") {
    q80 <- quantile(g, probs = 0.80, na.rm = TRUE)

    out <- ifelse(g > q80, "High", "All others")

    factor(out, levels = c("All others", "High"))

  } else if (trial_spec$grs_scheme == "top10_vs_rest") {
    q90 <- quantile(g, probs = 0.90, na.rm = TRUE)

    out <- ifelse(g > q90, "High", "Lower")

    factor(out, levels = c("Lower", "High"))

  } else {
    stop("Unknown grs_scheme: ", trial_spec$grs_scheme)
  }
}

# ----------------------------
# Trial specifications
# ----------------------------
#
# The same source-population DGP is used for all trials.
# Only trial_spec changes.
#
# baseline_cumhaz:
#   Trial-specific baseline cumulative event hazard over the analysis horizon.
#   This replaces alpha_Y + followup_y.
#
# For quick testing, set sample_multiplier <- 1.
# For smoother figures, increase to 5 or 10.

sample_multiplier <- 10

fourier_spec <- list(
  name = "FOURIER",

  n_source  = sample_multiplier * 300000,
  n_trial   = sample_multiplier * 14298,
  n_treat   = sample_multiplier * 7163,
  n_control = sample_multiplier * 7135,

  treatment_label = "Evolocumab",
  control_label   = "Placebo",
  outcome_label   = "major vascular events",

  # Previously equivalent to alpha_Y = log(0.032), followup_y = 2.3.
  baseline_cumhaz = 0.032 * 2.3,

  # FOURIER-like eligibility:
  # established ASCVD plus LDL-C >= 70 mg/dL.
  ldl_cut = 140,
  disease_rule = "D_eq_1",

  # Generate D by thresholding latent burden B so source prevalence is 25%.
  disease_threshold_basis = "source",
  p_disease_source = 0.25,

  # FOURIER genetic grouping:
  # low = quintile 1, intermediate = quintiles 2-4, high = quintile 5.
  grs_scheme = "fourier_3",

  plot_filename = "fourier_like_sim_figure2.png"
)

woscops_spec <- list(
  name = "WOSCOPS",

  n_source  = sample_multiplier * 300000,
  n_trial   = sample_multiplier * 4892,
  n_treat   = sample_multiplier * 2452,
  n_control = sample_multiplier * 2440,

  treatment_label = "Pravastatin",
  control_label   = "Placebo",
  outcome_label   = "incident CHD events",

  # Working cumulative event scale. Adjust this one value if the
  # WOSCOPS-like event rates are too high/low.
  baseline_cumhaz = 0.032 * 2.3,

  # Working WOSCOPS-like hypercholesterolaemia threshold.
  # 174 mg/dL is approximately 4.5 mmol/L.
  ldl_cut = 174,

  # WOSCOPS excluded prior MI / prior clinical coronary disease.
  disease_rule = "D_eq_0",

  # For WOSCOPS, calibrate prior-event prevalence among high-LDL people,
  # not in the whole source population.
  disease_threshold_basis = "lipid_eligible",
  p_disease_given_lipid = 0.1,

  # WOSCOPS headline genetic grouping:
  # high genetic risk = top quintile; all others = lower 80%.
  grs_scheme = "top20_vs_rest",

  plot_filename = "woscops_like_sim_figure1.png"
)

odyssey_spec <- list(
  name = "ODYSSEY",

  n_source  = sample_multiplier * 1000000,
  n_trial   = sample_multiplier * 11953,
  n_treat   = sample_multiplier * 5967,
  n_control = sample_multiplier * 5986,

  treatment_label = "Alirocumab",
  control_label   = "Placebo",
  outcome_label   = "MACE",

  # First-pass event scale.
  # Roughly equivalent to alpha_Y = log(0.038), followup_y = 2.8.
  # Adjust this one value if ODYSSEY-like marginal MACE rates are too high/low.
  baseline_cumhaz = 0.038 * 2.8,

  # ODYSSEY required LDL-C >=70 OR ApoB >=80 OR non-HDL-C >=100.
  # This LDL-only simulation uses LDL-C >=70 as the parsimonious proxy.
  ldl_cut = 140,

  # ODYSSEY was a post-ACS trial:
  # recent ACS 1-12 months before randomisation.
  disease_rule = "D_eq_1",

  # Recent ACS is sharper / more selective than broad stable ASCVD.
  disease_threshold_basis = "source",
  p_disease_source = 0.05,

  # ODYSSEY high PRS = top decile.
  grs_scheme = "top10_vs_rest",

  plot_filename = "odyssey_like_sim_figure3.png"
)

# ----------------------------
# Shared best-guess plausible parameters
# ----------------------------
#
# These parameters define the common source-population DGP.
# They are not fitted separately to FOURIER, WOSCOPS, or ODYSSEY.
#
# G = CAD genetic risk score
# U = unobserved vascular / treatment-response modifier
# LDL = baseline LDL-C
# B = latent atherosclerotic burden
# D = established ASCVD / prior event, defined by thresholding B
# A = randomised treatment
#
# Important:
#   There is no A:G term anywhere.
#
# beta_AU > 0 means high-U individuals receive less relative benefit
# from lipid lowering. If eligibility induces lower U among higher-G people,
# high-G individuals can appear to benefit more even when true A:G = 0.

pars_best <- list(
  # Source LDL-C distribution
  log_ldl_mu    = log(95),
  log_ldl_sigma = 0.32,

  # CAD GRS has weak direct effect on LDL-C.
  # U has stronger effect on LDL-C.
  beta_G_L = 0.02,
  beta_U_L = 0.12,

  # Latent atherosclerotic burden.
  beta_G_B = 0.70,
  beta_U_B = 0.80,
  beta_L_B = 0.25,
  sigma_B  = 0.80,

  # Event-risk effects after randomisation.
  # The baseline event scale is trial-specific via trial_spec$baseline_cumhaz.
  beta_G_Y = log(1.20),
  beta_U_Y = log(1.05),
  beta_L_Y = log(1.05),
  beta_B_Y = log(1.20),

  # Average treatment effect and latent response heterogeneity.
  # U is centred within the selected randomised trial before entering A:U.
  beta_A   = log(0.75),
  beta_AU  = 0.30
)

# ----------------------------
# Simulate source population
# ----------------------------

simulate_source <- function(n, pars, trial_spec) {
  G <- rnorm(n)
  U <- rnorm(n)

  log_LDL <- pars$log_ldl_mu +
    pars$beta_G_L * G +
    pars$beta_U_L * U +
    rnorm(n, sd = pars$log_ldl_sigma)

  LDL <- exp(log_LDL)
  log_LDL_z <- std(log(LDL))

  B_raw <- pars$beta_G_B * G +
    pars$beta_U_B * U +
    pars$beta_L_B * log_LDL_z +
    rnorm(n, sd = pars$sigma_B)

  B_z <- std(B_raw)

  if (trial_spec$disease_threshold_basis == "source") {
    B_threshold <- quantile(
      B_z,
      probs = 1 - trial_spec$p_disease_source,
      na.rm = TRUE
    )

  } else if (trial_spec$disease_threshold_basis == "lipid_eligible") {
    lipid_eligible <- LDL >= trial_spec$ldl_cut

    if (sum(lipid_eligible) < 100) {
      stop("Too few lipid-eligible individuals to define disease threshold.")
    }

    B_threshold <- quantile(
      B_z[lipid_eligible],
      probs = 1 - trial_spec$p_disease_given_lipid,
      na.rm = TRUE
    )

  } else {
    stop("Unknown disease_threshold_basis: ", trial_spec$disease_threshold_basis)
  }

  D <- as.integer(B_z >= B_threshold)

  tibble(
    G = G,
    U = U,
    LDL = LDL,
    log_LDL_z = log_LDL_z,
    B_z = B_z,
    B_threshold = as.numeric(B_threshold),
    D = D
  )
}

# ----------------------------
# Apply trial-specific eligibility
# ----------------------------

apply_trial_eligibility <- function(dat, trial_spec) {
  dat <- dat %>%
    mutate(
      lipid_eligible = LDL >= trial_spec$ldl_cut,
      disease_eligible = case_when(
        trial_spec$disease_rule == "D_eq_1" ~ D == 1,
        trial_spec$disease_rule == "D_eq_0" ~ D == 0,
        trial_spec$disease_rule == "none"   ~ TRUE,
        TRUE ~ NA
      )
    )

  if (any(is.na(dat$disease_eligible))) {
    stop("Unknown disease_rule: ", trial_spec$disease_rule)
  }

  dat %>%
    mutate(
      eligible = lipid_eligible & disease_eligible
    )
}

# ----------------------------
# Draw selected trial and randomise treatment
# ----------------------------

draw_trial <- function(source, trial_spec) {
  eligible <- source %>% filter(eligible)

  if (nrow(eligible) < trial_spec$n_trial) {
    stop(
      "Too few eligible individuals for ", trial_spec$name, ": ",
      nrow(eligible), " available; ", trial_spec$n_trial, " required. ",
      "Increase n_source or loosen trial_spec/DGP parameters."
    )
  }

  trial <- eligible %>%
    slice_sample(n = trial_spec$n_trial, replace = FALSE) %>%
    mutate(
      GRS_group = assign_grs_group(G, trial_spec)
    )

  A <- c(
    rep(1, trial_spec$n_treat),
    rep(0, trial_spec$n_control)
  )

  if (length(A) != trial_spec$n_trial) {
    stop("n_treat + n_control must equal n_trial for ", trial_spec$name)
  }

  A <- sample(A, length(A), replace = FALSE)

  trial %>%
    mutate(
      A = A,
      treatment = ifelse(A == 1, trial_spec$treatment_label, trial_spec$control_label),
      treatment = factor(
        treatment,
        levels = c(trial_spec$control_label, trial_spec$treatment_label)
      )
    )
}

# ----------------------------
# Simulate event outcome
# ----------------------------

simulate_events <- function(trial, pars, trial_spec) {
  # No A:G term.
  #
  # We use a unit analysis horizon. trial_spec$baseline_cumhaz already
  # combines the previous alpha_Y and followup_y terms.
  #
  # Event probability:
  #   1 - exp(-exp(lp))
  #
  # where lp includes log(baseline_cumhaz), covariates, and treatment.

  trial <- trial %>%
    mutate(
      U_response = U - mean(U)
    )

  lp <- with(
    trial,
    log(trial_spec$baseline_cumhaz) +
      pars$beta_G_Y * G +
      pars$beta_U_Y * U +
      pars$beta_L_Y * log_LDL_z +
      pars$beta_B_Y * B_z +
      A * (pars$beta_A + pars$beta_AU * U_response)
  )

  cumhaz <- safe_exp(lp)

  # Unit-horizon exponential survival.
  event_time <- rexp(nrow(trial), rate = cumhaz)
  event <- as.integer(event_time <= 1)
  time <- pmin(event_time, 1)

  trial %>%
    mutate(
      cumhaz = cumhaz,
      time = time,
      event = event
    )
}

# ----------------------------
# Analysis functions
# ----------------------------

event_table <- function(trial) {
  trial %>%
    group_by(GRS_group) %>%
    summarise(
      n_treat = sum(A == 1),
      n_control = sum(A == 0),
      events_treat = sum(event[A == 1]),
      events_control = sum(event[A == 0]),
      risk_treat = mean(event[A == 1]),
      risk_control = mean(event[A == 0]),
      ARR = risk_control - risk_treat,
      RR = rr(risk_treat, risk_control),
      .groups = "drop"
    )
}

cox_by_grs <- function(trial) {
  trial %>%
    group_split(GRS_group) %>%
    map_df(function(d) {
      fit <- coxph(Surv(time, event) ~ A, data = d)

      broom::tidy(fit, exponentiate = TRUE, conf.int = TRUE) %>%
        filter(term == "A") %>%
        transmute(
          GRS_group = unique(d$GRS_group),
          HR = estimate,
          HR_low = conf.low,
          HR_high = conf.high,
          p = p.value
        )
    })
}

interaction_models <- function(trial) {
  fit_cont <- coxph(Surv(time, event) ~ A * G, data = trial)
  fit_cat  <- coxph(Surv(time, event) ~ A * GRS_group, data = trial)

  list(
    continuous = fit_cont,
    categorical = fit_cat
  )
}

diagnostics_by_grs <- function(trial) {
  trial %>%
    group_by(GRS_group) %>%
    summarise(
      n = n(),
      mean_G = mean(G),
      mean_U = mean(U),
      mean_U_response = mean(U_response),
      mean_LDL = mean(LDL),
      mean_B_z = mean(B_z),
      event_rate = mean(event),
      control_event_rate = mean(event[A == 0]),
      treatment_event_rate = mean(event[A == 1]),
      .groups = "drop"
    )
}

trial_pattern <- function(event_tab, cox_tab, trial_spec) {
  x <- event_tab %>%
    left_join(cox_tab %>% select(GRS_group, HR), by = "GRS_group")

  if (trial_spec$grs_scheme == "fourier_3") {
    low <- x %>% filter(GRS_group == "Low")
    int <- x %>% filter(GRS_group == "Intermediate")
    high <- x %>% filter(GRS_group == "High")

    tibble(
      risk_gradient =
        high$risk_control > int$risk_control &&
        int$risk_control > low$risk_control,
      arr_gradient =
        high$ARR > int$ARR &&
        int$ARR > low$ARR,
      hr_gradient =
        high$HR < int$HR &&
        int$HR < low$HR,
      headline_like =
        risk_gradient && (arr_gradient || hr_gradient)
    )

  } else if (trial_spec$grs_scheme %in% c("top20_vs_rest", "top10_vs_rest")) {
    lower_label <- ifelse(
      trial_spec$grs_scheme == "top20_vs_rest",
      "All others",
      "Lower"
    )

    rest <- x %>% filter(GRS_group == lower_label)
    high <- x %>% filter(GRS_group == "High")

    tibble(
      risk_gradient = high$risk_control > rest$risk_control,
      arr_gradient = high$ARR > rest$ARR,
      hr_gradient = high$HR < rest$HR,
      headline_like = risk_gradient && (arr_gradient || hr_gradient)
    )

  } else {
    stop("Unknown grs_scheme: ", trial_spec$grs_scheme)
  }
}

summarise_trial <- function(trial, source = NULL, pars = NULL, trial_spec) {
  tab <- event_table(trial)
  cox_tab <- cox_by_grs(trial)
  fits <- interaction_models(trial)

  cont_int <- broom::tidy(
    fits$continuous,
    exponentiate = FALSE,
    conf.int = TRUE
  ) %>%
    filter(term == "A:G") %>%
    transmute(
      beta_AxG = estimate,
      beta_AxG_low = conf.low,
      beta_AxG_high = conf.high,
      p_AxG = p.value
    )

  tab_wide <- tab %>%
    left_join(cox_tab %>% select(GRS_group, HR), by = "GRS_group") %>%
    mutate(GRS_group_id = clean_group_name(GRS_group)) %>%
    select(GRS_group_id, risk_control, risk_treat, ARR, RR, HR) %>%
    pivot_wider(
      names_from = GRS_group_id,
      values_from = c(risk_control, risk_treat, ARR, RR, HR),
      names_sep = "_"
    )

  lm_u_g <- lm(U_response ~ G, data = trial)
  slope_u_on_g <- coef(lm_u_g)[["G"]]

  out <- tibble(
    trial_name = trial_spec$name,
    n = nrow(trial),
    ldl_cut = trial_spec$ldl_cut,
    disease_rule = trial_spec$disease_rule,
    baseline_cumhaz = trial_spec$baseline_cumhaz,
    cor_G_U_trial = cor(trial$G, trial$U),
    slope_U_on_G_trial = slope_u_on_g,
    approx_induced_AxG = if (!is.null(pars)) pars$beta_AU * slope_u_on_g else NA_real_,
    mean_LDL = mean(trial$LDL),
    median_LDL = median(trial$LDL),
    prop_LDL_above_cut = mean(trial$LDL >= trial_spec$ldl_cut),
    prop_D = mean(trial$D),
    mean_B_z = mean(trial$B_z),
    event_rate_control = mean(trial$event[trial$A == 0]),
    event_rate_treat = mean(trial$event[trial$A == 1])
  ) %>%
    bind_cols(cont_int) %>%
    bind_cols(tab_wide)

  if (!is.null(source)) {
    out <- out %>%
      mutate(
        cor_G_U_source = safe_cor(source$G, source$U),
        cor_G_U_D1 = safe_cor(source$G[source$D == 1], source$U[source$D == 1]),
        cor_G_U_D0 = safe_cor(source$G[source$D == 0], source$U[source$D == 0]),
        cor_G_U_eligible = safe_cor(source$G[source$eligible], source$U[source$eligible]),
        mean_G_source = mean(source$G),
        mean_G_D1 = mean(source$G[source$D == 1]),
        mean_G_D0 = mean(source$G[source$D == 0]),
        mean_G_eligible = mean(source$G[source$eligible]),
        eligible_fraction = mean(source$eligible),
        disease_fraction_source = mean(source$D == 1),
        disease_fraction_lipid_eligible = mean(source$D[source$LDL >= trial_spec$ldl_cut] == 1)
      )
  }

  pattern <- trial_pattern(tab, cox_tab, trial_spec)

  list(
    summary = out,
    event_table = tab,
    cox_by_grs = cox_tab,
    continuous_interaction = fits$continuous,
    categorical_interaction = fits$categorical,
    pattern = pattern
  )
}

# ----------------------------
# One full simulation
# ----------------------------

simulate_once <- function(
  trial_spec,
  pars = pars_best,
  keep_source = TRUE
) {
  source <- simulate_source(trial_spec$n_source, pars, trial_spec) %>%
    apply_trial_eligibility(trial_spec)

  trial <- source %>%
    draw_trial(trial_spec) %>%
    simulate_events(pars, trial_spec)

  analysed <- summarise_trial(
    trial = trial,
    source = if (keep_source) source else NULL,
    pars = pars,
    trial_spec = trial_spec
  )

  list(
    source = if (keep_source) source else NULL,
    trial = trial,
    analysis = analysed,
    diagnostics_by_grs = diagnostics_by_grs(trial),
    pattern = analysed$pattern,
    pars = pars,
    trial_spec = trial_spec
  )
}

# ----------------------------
# Printing helper
# ----------------------------

print_simulation_results <- function(sim) {
  cat("\n============================================================\n")
  cat(sim$trial_spec$name, "simulation summary\n")
  cat("============================================================\n")

  cat("\nSummary:\n")
  print(sim$analysis$summary)

  cat("\nDiagnostics by GRS category:\n")
  print(sim$diagnostics_by_grs)

  cat("\nEvent table by GRS category:\n")
  print(sim$analysis$event_table)

  cat("\nCox HRs by GRS category:\n")
  print(sim$analysis$cox_by_grs)

  cat("\nContinuous A x G interaction model:\n")
  print(summary(sim$analysis$continuous_interaction))

  cat("\nCategorical A x GRS-group interaction model:\n")
  print(summary(sim$analysis$categorical_interaction))

  cat("\nHeadline-like pattern flags:\n")
  print(sim$pattern)
}

# ============================================================
# Generic Figure-style plot
# Relative and absolute risk reduction by genetic-risk group
# ============================================================

plot_risk_reduction_by_grs <- function(sim) {
  trial <- sim$trial %>%
    mutate(
      GRS_ord = as.numeric(GRS_group)
    )

  trial_spec <- sim$trial_spec

  event_tab <- sim$analysis$event_table %>%
    left_join(sim$analysis$cox_by_grs, by = "GRS_group") %>%
    mutate(
      RRR = 1 - HR,
      label = sprintf(
        "HR: %.2f (%.2f-%.2f)\nRRR: %.1f%%\nARR: %.1f%%",
        HR, HR_low, HR_high, 100 * RRR, 100 * ARR
      )
    )

  plot_dat <- event_tab %>%
    select(GRS_group, risk_control, risk_treat) %>%
    pivot_longer(
      cols = c(risk_control, risk_treat),
      names_to = "arm",
      values_to = "risk"
    ) %>%
    mutate(
      arm = recode(
        arm,
        "risk_control" = trial_spec$control_label,
        "risk_treat" = trial_spec$treatment_label
      ),
      arm = factor(
        arm,
        levels = c(trial_spec$control_label, trial_spec$treatment_label)
      ),
      risk_pct = 100 * risk
    )

  ann_dat <- event_tab %>%
    mutate(
      y = 100 * pmax(risk_control, risk_treat) + 3.0
    )

  # Trend/difference test for relative effect.
  # For 3-level FOURIER grouping this is an ordered trend;
  # for 2-level WOSCOPS/ODYSSEY grouping this is equivalent to a two-group interaction.
  cox_trend <- coxph(
    Surv(time, event) ~ A * GRS_ord,
    data = trial
  )

  p_trend_hr <- broom::tidy(cox_trend) %>%
    filter(term == "A:GRS_ord") %>%
    pull(p.value)

  # Trend/difference test for absolute risk reduction.
  arr_trend <- lm(
    event ~ A * GRS_ord,
    data = trial
  )

  p_trend_arr <- broom::tidy(arr_trend) %>%
    filter(term == "A:GRS_ord") %>%
    pull(p.value)

  trend_label <- sprintf(
    "P-trend/difference (HR) = %s\nP-trend/difference (ARR) = %s",
    fmt_p(p_trend_hr),
    fmt_p(p_trend_arr)
  )

  y_max <- max(c(plot_dat$risk_pct, ann_dat$y), na.rm = TRUE) + 3

  ggplot(plot_dat, aes(x = GRS_group, y = risk_pct, fill = arm)) +
    geom_col(
      position = position_dodge(width = 0.72),
      width = 0.62,
      colour = "grey30",
      linewidth = 0.25
    ) +
    geom_text(
      data = ann_dat,
      aes(x = GRS_group, y = y, label = label),
      inherit.aes = FALSE,
      size = 3.4,
      lineheight = 0.92
    ) +
    annotate(
      "text",
      x = 0.55,
      y = y_max,
      label = trend_label,
      hjust = 0,
      vjust = 1,
      size = 3.5
    ) +
    scale_fill_manual(
      values = setNames(
        c("#4C78A8", "#D62728"),
        c(trial_spec$control_label, trial_spec$treatment_label)
      )
    ) +
    coord_cartesian(ylim = c(0, y_max), clip = "off") +
    theme_bw(base_size = 13) +
    theme(
      legend.position = "top",
      legend.title = element_blank(),
      panel.grid.major.x = element_blank(),
      panel.grid.minor = element_blank(),
      axis.title.x = element_blank(),
      plot.title = element_text(face = "bold"),
      plot.margin = margin(8, 16, 8, 8)
    ) +
    labs(
      title = paste0(
        "Relative and absolute risk reduction of ",
        trial_spec$outcome_label
      ),
      subtitle = paste0(
        trial_spec$name,
        "-like simulation; true treatment-by-PRS interaction = 0"
      ),
      y = paste0(trial_spec$outcome_label, " (%)")
    )
}

# ============================================================
# Run all trial specifications
# ============================================================

sim_fourier <- simulate_once(
  trial_spec = fourier_spec,
  pars = pars_best,
  keep_source = TRUE
)

sim_woscops <- simulate_once(
  trial_spec = woscops_spec,
  pars = pars_best,
  keep_source = TRUE
)

sim_odyssey <- simulate_once(
  trial_spec = odyssey_spec,
  pars = pars_best,
  keep_source = TRUE
)

print_simulation_results(sim_fourier)
print_simulation_results(sim_woscops)
print_simulation_results(sim_odyssey)

# ============================================================
# Make plots
# ============================================================

p_fourier <- plot_risk_reduction_by_grs(sim_fourier)
p_woscops <- plot_risk_reduction_by_grs(sim_woscops)
p_odyssey <- plot_risk_reduction_by_grs(sim_odyssey)

print(p_fourier)
print(p_woscops)
print(p_odyssey)

ggsave(
  filename = fourier_spec$plot_filename,
  plot = p_fourier,
  width = 8.5,
  height = 5.5,
  dpi = 320
)

ggsave(
  filename = woscops_spec$plot_filename,
  plot = p_woscops,
  width = 8.5,
  height = 5.5,
  dpi = 320
)

ggsave(
  filename = odyssey_spec$plot_filename,
  plot = p_odyssey,
  width = 8.5,
  height = 5.5,
  dpi = 320
)