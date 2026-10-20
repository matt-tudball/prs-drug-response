# ============================================================
# Parsimonious simulation of trial-eligibility-induced
# treatment-by-PRS interactions in lipid-lowering trials
#
# Variant: includes trial-analysis-style Cox adjustment for
# baseline LDL-C and a noisy observed proxy for atherosclerotic burden
#
# Trials represented:
#   FOURIER-like: established ASCVD + elevated LDL-C on statin
#   ODYSSEY-like: recent ACS + elevated LDL-C on intensive statin
#   WOSCOPS-like: high untreated LDL-C + no prior MI/clinical CHD
#
# Shared causal structure:
#   G -> LDL
#   U -> LDL
#   G -> B
#   U -> B
#   LDL -> B
#   eligibility = f(B, LDL), trial-specific
#   A randomized after eligibility
#   Y depends on A, G, U, LDL, B, and A:U
#   Y does NOT contain A:G
#
# Any apparent treatment-by-PRS interaction is therefore artefactual.
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
# User controls
# ----------------------------

# Number of independent Monte Carlo replicates to run at the actual trial sample size.
# Use a small value such as 5 for quick checks and a larger value such as 100-500 for manuscript figures.
n_reps <- 250L

write_plots <- TRUE
run_beta_AU_sensitivity <- TRUE
output_dir <- "simulation-outputs"

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

safe_slope <- function(y, x) {
  if (length(x) < 2 || sd(x, na.rm = TRUE) == 0) return(NA_real_)
  unname(coef(lm(y ~ x))[2])
}

rr <- function(r_treat, r_control) {
  ifelse(r_control == 0, NA_real_, r_treat / r_control)
}

risk_to_cumhaz <- function(p) {
  p <- pmin(pmax(p, 0), 1 - 1e-12)
  -log1p(-p)
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
# Eligibility builders
# ----------------------------

make_high_burden_eligibility <- function(p_high_burden, ldl_cut) {
  force(p_high_burden)
  force(ldl_cut)

  list(
    type = "high_burden_plus_ldl",
    p_high_burden = p_high_burden,
    p_excluded_burden = NA_real_,
    ldl_cut = ldl_cut,
    fn = function(dat) {
      burden_cut <- quantile(dat$B_z, probs = 1 - p_high_burden, na.rm = TRUE)
      dat$B_z >= burden_cut & dat$LDL >= ldl_cut
    },
    burden_cut_fn = function(dat) {
      as.numeric(quantile(dat$B_z, probs = 1 - p_high_burden, na.rm = TRUE))
    }
  )
}

make_no_prior_event_eligibility <- function(p_excluded_burden, ldl_cut) {
  force(p_excluded_burden)
  force(ldl_cut)

  list(
    type = "ldl_plus_no_prior_event",
    p_high_burden = NA_real_,
    p_excluded_burden = p_excluded_burden,
    ldl_cut = ldl_cut,
    fn = function(dat) {
      lipid_eligible <- dat$LDL >= ldl_cut
      if (sum(lipid_eligible) < 100) {
        stop("Too few lipid-eligible individuals to define the WOSCOPS burden exclusion.")
      }
      burden_cut <- quantile(
        dat$B_z[lipid_eligible],
        probs = 1 - p_excluded_burden,
        na.rm = TRUE
      )
      lipid_eligible & dat$B_z < burden_cut
    },
    burden_cut_fn = function(dat) {
      lipid_eligible <- dat$LDL >= ldl_cut
      as.numeric(quantile(
        dat$B_z[lipid_eligible],
        probs = 1 - p_excluded_burden,
        na.rm = TRUE
      ))
    }
  )
}

# ----------------------------
# Shared GRS grouping function
# ----------------------------

assign_grs_group <- function(g, grs_scheme) {
  if (grs_scheme == "fourier_3") {
    qs <- quantile(g, probs = c(0.20, 0.80), na.rm = TRUE)

    out <- case_when(
      g <= qs[1] ~ "Low",
      g >  qs[2] ~ "High",
      TRUE       ~ "Intermediate"
    )

    factor(out, levels = c("Low", "Intermediate", "High"))

  } else if (grs_scheme == "top20_vs_rest") {
    q80 <- quantile(g, probs = 0.80, na.rm = TRUE)
    out <- ifelse(g > q80, "High", "All others")
    factor(out, levels = c("All others", "High"))

  } else if (grs_scheme == "top10_vs_rest") {
    q90 <- quantile(g, probs = 0.90, na.rm = TRUE)
    out <- ifelse(g > q90, "High", "Lower")
    factor(out, levels = c("Lower", "High"))

  } else {
    stop("Unknown grs_scheme: ", grs_scheme)
  }
}

# ----------------------------
# Trial specifications
# ----------------------------
#
# LDL-C is interpreted as an untreated/off-statin-equivalent LDL-C liability.
# For FOURIER/ODYSSEY, the LDL thresholds are therefore back-calculated from
# the on-statin eligibility thresholds. WOSCOPS uses the observed untreated
# protocol threshold directly.

fourier_spec <- list(
  name = "FOURIER",

  n_source  = 600000L,
  n_trial   = 14298L,
  n_treat   = 7163L,
  n_control = 7135L,

  treatment_label = "Evolocumab",
  control_label   = "Placebo",
  outcome_label   = "major vascular events",

  target_control_event_rate = 0.115,
  baseline_cumhaz = NULL,

  eligibility = make_high_burden_eligibility(
    p_high_burden = 0.05,
    ldl_cut = 130
  ),

  grs_scheme = "fourier_3",
  plot_filename = "fourier_like_simulation.png"
)

odyssey_spec <- list(
  name = "ODYSSEY",

  n_source  = 1500000L,
  n_trial   = 11953L,
  n_treat   = 5967L,
  n_control = 5986L,

  treatment_label = "Alirocumab",
  control_label   = "Placebo",
  outcome_label   = "MACE",

  target_control_event_rate = 0.120,
  baseline_cumhaz = NULL,

  eligibility = make_high_burden_eligibility(
    p_high_burden = 0.02,
    ldl_cut = 140
  ),

  grs_scheme = "top10_vs_rest",
  plot_filename = "odyssey_like_simulation.png"
)

woscops_spec <- list(
  name = "WOSCOPS",

  n_source  = 300000L,
  n_trial   = 4892L,
  n_treat   = 2452L,
  n_control = 2440L,

  treatment_label = "Pravastatin",
  control_label   = "Placebo",
  outcome_label   = "incident CHD events",

  target_control_event_rate = 0.142,
  baseline_cumhaz = NULL,

  eligibility = make_no_prior_event_eligibility(
    p_excluded_burden = 0.05,
    ldl_cut = 174
  ),

  grs_scheme = "top20_vs_rest",
  plot_filename = "woscops_like_simulation.png"
)

trial_specs <- list(
  FOURIER = fourier_spec,
  ODYSSEY = odyssey_spec,
  WOSCOPS = woscops_spec
)

# ----------------------------
# Shared source-population DGP parameters
# ----------------------------
#
# These are deliberately parsimonious and are held constant across trials.
# Only eligibility, target placebo event rate, sample size, and GRS grouping
# differ by trial.

pars_anchor <- list(
  # LDL-C distribution. The intercept is a log-scale median, not a natural-scale mean.
  ldl_median_mgdl = 137.5,
  log_ldl_error_sd = 0.25,

  # Weak CAD PRS effect on clinical LDL-C; stronger non-genetic/metabolic component.
  beta_G_L = 0.03,
  beta_U_L = 0.12,

  # Latent atherosclerotic burden.
  beta_G_B = 0.55,
  beta_U_B = 0.70,
  beta_L_B = 0.20,
  sigma_B  = 0.80,

  # Prognostic effects after randomization.
  beta_G_Y = log(1.20),
  beta_U_Y = log(1.10),
  beta_L_Y = log(1.10),
  beta_B_Y = log(1.30),

  # Observed proxy for latent atherosclerotic burden. This is deliberately
  # noisy: trial analyses can adjust for clinical risk factors or disease
  # history, but not for true latent burden B itself.
  sigma_B_proxy = 0.70,

  # Treatment effect and response heterogeneity.
  # U_response is centered within each randomized trial, so beta_A is the
  # trial-average log-HR and beta_AU is heterogeneity around that average.
  # There is no beta_A_G term.
  beta_A   = log(0.80),
  beta_AU  = 0.30
)

# ----------------------------
# Source population simulation
# ----------------------------

simulate_source <- function(n, pars) {
  G <- rnorm(n)
  U <- rnorm(n)

  log_LDL <- log(pars$ldl_median_mgdl) +
    pars$beta_G_L * G +
    pars$beta_U_L * U +
    rnorm(n, sd = pars$log_ldl_error_sd)

  LDL <- exp(log_LDL)
  log_LDL_z <- std(log_LDL)

  B_raw <- pars$beta_G_B * G +
    pars$beta_U_B * U +
    pars$beta_L_B * log_LDL_z +
    rnorm(n, sd = pars$sigma_B)

  B_z <- std(B_raw)

  B_proxy_raw <- B_z + rnorm(n, sd = pars$sigma_B_proxy)
  B_proxy_z <- std(B_proxy_raw)

  tibble(
    G = G,
    U = U,
    LDL = LDL,
    log_LDL = log_LDL,
    log_LDL_z = log_LDL_z,
    B_z = B_z,
    B_proxy_z = B_proxy_z
  )
}

apply_trial_eligibility <- function(source, trial_spec) {
  source %>%
    mutate(
      eligible = trial_spec$eligibility$fn(source),
      lipid_eligible = LDL >= trial_spec$eligibility$ldl_cut
    )
}

# ----------------------------
# Draw selected trial and randomize treatment
# ----------------------------

draw_trial <- function(source, trial_spec) {
  eligible <- source %>% filter(eligible)

  if (nrow(eligible) < trial_spec$n_trial) {
    stop(
      "Too few eligible individuals for ", trial_spec$name, ": ",
      nrow(eligible), " available; ", trial_spec$n_trial, " required. ",
      "Increase n_source or loosen the eligibility rule."
    )
  }

  trial <- eligible %>%
    slice_sample(n = trial_spec$n_trial, replace = FALSE) %>%
    mutate(
      GRS_group = assign_grs_group(G, trial_spec$grs_scheme)
    )

  A <- c(
    rep(1L, trial_spec$n_treat),
    rep(0L, trial_spec$n_control)
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
# Outcome model
# ----------------------------

linear_predictor_without_treatment <- function(trial, pars) {
  with(
    trial,
    pars$beta_G_Y * G +
      pars$beta_U_Y * U +
      pars$beta_L_Y * log_LDL_z +
      pars$beta_B_Y * B_z
  )
}

calibrate_baseline_cumhaz <- function(trial, pars, target_control_event_rate) {
  lp0 <- linear_predictor_without_treatment(trial, pars)
  lp0_control <- lp0[trial$A == 0]

  mean_event_rate <- function(ch) {
    mean(1 - exp(-ch * safe_exp(lp0_control)))
  }

  f <- function(ch) mean_event_rate(ch) - target_control_event_rate

  if (f(1e-8) > 0) return(1e-8)
  upper <- 1
  while (f(upper) < 0 && upper < 100) {
    upper <- upper * 2
  }

  if (upper >= 100 && f(upper) < 0) {
    stop("Could not calibrate baseline cumulative hazard.")
  }

  uniroot(f, lower = 1e-8, upper = upper)$root
}

simulate_events <- function(trial, pars, trial_spec) {
  trial <- trial %>%
    mutate(
      U_response = U - mean(U)
    )

  baseline_cumhaz <- trial_spec$baseline_cumhaz
  if (is.null(baseline_cumhaz)) {
    baseline_cumhaz <- calibrate_baseline_cumhaz(
      trial = trial,
      pars = pars,
      target_control_event_rate = trial_spec$target_control_event_rate
    )
  }

  lp_no_treatment <- linear_predictor_without_treatment(trial, pars)

  lp <- log(baseline_cumhaz) +
    lp_no_treatment +
    trial$A * (pars$beta_A + pars$beta_AU * trial$U_response)

  cumhaz <- safe_exp(lp)

  event_time <- rexp(nrow(trial), rate = cumhaz)
  event <- as.integer(event_time <= 1)
  time <- pmin(event_time, 1)

  trial %>%
    mutate(
      baseline_cumhaz_used = baseline_cumhaz,
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
      n = n(),
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
  list(
    continuous = coxph(Surv(time, event) ~ A * G, data = trial),
    categorical = coxph(Surv(time, event) ~ A * GRS_group, data = trial)
  )
}

# ----------------------------
# Trial-analysis-style adjusted Cox models
# ----------------------------
#
# These are intentionally not part of the data-generating model. They are
# analysis-layer sensitivity models that mimic the fact that trial reports
# often estimate PRS subgroup treatment HRs after adjustment for baseline
# biomarkers and clinical risk factors. We include only two parsimonious
# measured covariates:
#   1. log_LDL_z: baseline LDL-C biomarker
#   2. B_proxy_z: noisy observed proxy for latent atherosclerotic burden
#
# We never adjust for the true latent B_z in the main adjusted models, because
# true burden is not observed in real trial analyses.

adjustment_sets <- list(
  unadjusted = character(0),
  ldl_adjusted = c("log_LDL_z"),
  burden_proxy_adjusted = c("B_proxy_z"),
  ldl_burden_proxy_adjusted = c("log_LDL_z", "B_proxy_z")
)

make_survival_formula <- function(base_rhs, adjust_vars = character(0)) {
  rhs <- paste(c(base_rhs, adjust_vars), collapse = " + ")
  as.formula(paste("Surv(time, event) ~", rhs))
}

cox_by_grs_adjusted_one <- function(trial, adjust_vars, model_set) {
  trial %>%
    group_split(GRS_group) %>%
    map_df(function(d) {
      fit <- coxph(
        make_survival_formula("A", adjust_vars),
        data = d
      )

      broom::tidy(fit, exponentiate = TRUE, conf.int = TRUE) %>%
        filter(term == "A") %>%
        transmute(
          model_set = model_set,
          adjustment = ifelse(length(adjust_vars) == 0, "None", paste(adjust_vars, collapse = " + ")),
          GRS_group = unique(d$GRS_group),
          HR = estimate,
          HR_low = conf.low,
          HR_high = conf.high,
          p = p.value
        )
    })
}

cox_by_grs_adjusted_all <- function(trial) {
  imap_dfr(
    adjustment_sets,
    ~ cox_by_grs_adjusted_one(trial, adjust_vars = .x, model_set = .y)
  )
}

continuous_interaction_adjusted_one <- function(trial, adjust_vars, model_set) {
  fit <- coxph(
    make_survival_formula("A * G", adjust_vars),
    data = trial
  )

  broom::tidy(fit, exponentiate = FALSE, conf.int = TRUE) %>%
    filter(term == "A:G") %>%
    transmute(
      model_set = model_set,
      adjustment = ifelse(length(adjust_vars) == 0, "None", paste(adjust_vars, collapse = " + ")),
      beta_AxG = estimate,
      beta_AxG_low = conf.low,
      beta_AxG_high = conf.high,
      p_AxG = p.value
    )
}

continuous_interactions_adjusted_all <- function(trial) {
  imap_dfr(
    adjustment_sets,
    ~ continuous_interaction_adjusted_one(trial, adjust_vars = .x, model_set = .y)
  )
}

categorical_interaction_adjusted_one <- function(trial, adjust_vars, model_set) {
  fit <- coxph(
    make_survival_formula("A * GRS_group", adjust_vars),
    data = trial
  )

  broom::tidy(fit, exponentiate = FALSE, conf.int = TRUE) %>%
    filter(grepl("^A:GRS_group", term)) %>%
    transmute(
      model_set = model_set,
      adjustment = ifelse(length(adjust_vars) == 0, "None", paste(adjust_vars, collapse = " + ")),
      interaction_term = term,
      beta_interaction = estimate,
      beta_interaction_low = conf.low,
      beta_interaction_high = conf.high,
      p_interaction = p.value
    )
}

categorical_interactions_adjusted_all <- function(trial) {
  imap_dfr(
    adjustment_sets,
    ~ categorical_interaction_adjusted_one(trial, adjust_vars = .x, model_set = .y)
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
      mean_B_proxy_z = mean(B_proxy_z),
      control_event_rate = mean(event[A == 0]),
      treatment_event_rate = mean(event[A == 1]),
      .groups = "drop"
    )
}

eligibility_diagnostics <- function(source, trial_spec) {
  elig <- source$eligible
  lip <- source$lipid_eligible

  tibble(
    eligibility_type = trial_spec$eligibility$type,
    ldl_cut = trial_spec$eligibility$ldl_cut,
    p_high_burden = trial_spec$eligibility$p_high_burden,
    p_excluded_burden = trial_spec$eligibility$p_excluded_burden,
    burden_cut_used = trial_spec$eligibility$burden_cut_fn(source),
    source_n = nrow(source),
    eligible_n = sum(elig),
    eligible_fraction = mean(elig),
    lipid_eligible_fraction = mean(lip),
    mean_LDL_source = mean(source$LDL),
    median_LDL_source = median(source$LDL),
    mean_LDL_eligible = mean(source$LDL[elig]),
    cor_G_U_source = safe_cor(source$G, source$U),
    cor_G_U_eligible = safe_cor(source$G[elig], source$U[elig]),
    cor_G_logLDL_source = safe_cor(source$G, source$log_LDL),
    cor_G_logLDL_eligible = safe_cor(source$G[elig], source$log_LDL[elig]),
    cor_G_Bproxy_source = safe_cor(source$G, source$B_proxy_z),
    cor_G_Bproxy_eligible = safe_cor(source$G[elig], source$B_proxy_z[elig]),
    slope_U_on_G_eligible = safe_slope(source$U[elig], source$G[elig]),
    slope_logLDL_on_G_eligible = safe_slope(source$log_LDL[elig], source$G[elig]),
    slope_Bproxy_on_G_eligible = safe_slope(source$B_proxy_z[elig], source$G[elig]),
    mean_G_eligible = mean(source$G[elig]),
    mean_U_eligible = mean(source$U[elig]),
    mean_B_z_eligible = mean(source$B_z[elig])
  )
}

trial_pattern <- function(event_tab, cox_tab, grs_scheme) {
  x <- event_tab %>%
    left_join(cox_tab %>% select(GRS_group, HR), by = "GRS_group")

  if (grs_scheme == "fourier_3") {
    low <- x %>% filter(GRS_group == "Low")
    int <- x %>% filter(GRS_group == "Intermediate")
    high <- x %>% filter(GRS_group == "High")

    tibble(
      risk_gradient = high$risk_control > int$risk_control && int$risk_control > low$risk_control,
      arr_gradient = high$ARR > int$ARR && int$ARR > low$ARR,
      hr_gradient = high$HR < int$HR && int$HR < low$HR,
      headline_like = risk_gradient && (arr_gradient || hr_gradient)
    )

  } else if (grs_scheme %in% c("top20_vs_rest", "top10_vs_rest")) {
    lower_label <- ifelse(grs_scheme == "top20_vs_rest", "All others", "Lower")

    rest <- x %>% filter(GRS_group == lower_label)
    high <- x %>% filter(GRS_group == "High")

    tibble(
      risk_gradient = high$risk_control > rest$risk_control,
      arr_gradient = high$ARR > rest$ARR,
      hr_gradient = high$HR < rest$HR,
      headline_like = risk_gradient && (arr_gradient || hr_gradient)
    )

  } else {
    stop("Unknown grs_scheme: ", grs_scheme)
  }
}


# ----------------------------
# Analytic risks implied by the data-generating model
# ----------------------------

analytic_group_risks <- function(trial, pars) {
  lp0 <- linear_predictor_without_treatment(trial, pars)
  H0 <- unique(trial$baseline_cumhaz_used)

  if (length(H0) != 1L) {
    stop("Expected one baseline cumulative hazard value per simulated trial.")
  }

  trial %>%
    mutate(
      p_placebo_analytic =
        1 - exp(-H0 * safe_exp(lp0)),
      p_treatment_analytic =
        1 - exp(-H0 * safe_exp(lp0 + pars$beta_A + pars$beta_AU * U_response)),
      individual_treatment_HR =
        exp(pars$beta_A + pars$beta_AU * U_response)
    ) %>%
    group_by(GRS_group) %>%
    summarise(
      n = n(),
      analytic_placebo_risk = mean(p_placebo_analytic),
      analytic_treatment_risk = mean(p_treatment_analytic),
      analytic_ARR = analytic_placebo_risk - analytic_treatment_risk,

      # This is the group-level cumulative-hazard ratio over the fixed follow-up
      # horizon, not a fitted Cox HR. It is the clean analytic analogue after
      # averaging heterogeneous individual risks within a PRS group.
      analytic_treatment_HR_cumhaz =
        risk_to_cumhaz(analytic_treatment_risk) /
        risk_to_cumhaz(analytic_placebo_risk),

      # This is the arithmetic mean of the individual-level treatment HRs from
      # the DGP. It is useful as a diagnostic but is not equal to a marginal Cox HR.
      analytic_mean_individual_treatment_HR = mean(individual_treatment_HR),

      analytic_mean_U_response = mean(U_response),
      analytic_mean_G = mean(G),
      .groups = "drop"
    )
}

lower_group_label <- function(grs_scheme) {
  if (grs_scheme == "fourier_3") {
    "Low"
  } else if (grs_scheme == "top20_vs_rest") {
    "All others"
  } else if (grs_scheme == "top10_vs_rest") {
    "Lower"
  } else {
    stop("Unknown grs_scheme: ", grs_scheme)
  }
}

analytic_high_low_contrasts <- function(analytic_tab, grs_scheme) {
  low_label <- lower_group_label(grs_scheme)

  low <- analytic_tab %>% filter(GRS_group == low_label)
  high <- analytic_tab %>% filter(GRS_group == "High")

  if (nrow(low) != 1L || nrow(high) != 1L) {
    stop("Expected exactly one lower-risk group and one high-risk group.")
  }

  tibble(
    lower_group = low_label,
    high_group = "High",

    analytic_placebo_risk_lower = low$analytic_placebo_risk,
    analytic_placebo_risk_high = high$analytic_placebo_risk,
    analytic_treatment_risk_lower = low$analytic_treatment_risk,
    analytic_treatment_risk_high = high$analytic_treatment_risk,

    # High-vs-lower PRS risk gradients within each randomized arm.
    analytic_placebo_high_vs_lower_cumhaz_HR =
      risk_to_cumhaz(high$analytic_placebo_risk) /
      risk_to_cumhaz(low$analytic_placebo_risk),
    analytic_treatment_high_vs_lower_cumhaz_HR =
      risk_to_cumhaz(high$analytic_treatment_risk) /
      risk_to_cumhaz(low$analytic_treatment_risk),

    # Treatment HRs within lower and high PRS groups.
    analytic_treatment_HR_lower = low$analytic_treatment_HR_cumhaz,
    analytic_treatment_HR_high = high$analytic_treatment_HR_cumhaz,
    analytic_treatment_HR_ratio_high_vs_lower =
      high$analytic_treatment_HR_cumhaz /
      low$analytic_treatment_HR_cumhaz,

    analytic_ARR_lower = low$analytic_ARR,
    analytic_ARR_high = high$analytic_ARR,
    analytic_ARR_difference_high_minus_lower =
      high$analytic_ARR - low$analytic_ARR
  )
}

summarise_trial <- function(trial, source, pars, trial_spec) {
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

  slope_u_on_g <- safe_slope(trial$U_response, trial$G)

  summary <- tibble(
    trial_name = trial_spec$name,
    n = nrow(trial),
    baseline_cumhaz_used = unique(trial$baseline_cumhaz_used),
    target_control_event_rate = trial_spec$target_control_event_rate,
    event_rate_control = mean(trial$event[trial$A == 0]),
    event_rate_treat = mean(trial$event[trial$A == 1]),
    cor_G_U_trial = safe_cor(trial$G, trial$U),
    cor_G_logLDL_trial = safe_cor(trial$G, trial$log_LDL),
    cor_G_Bproxy_trial = safe_cor(trial$G, trial$B_proxy_z),
    slope_U_on_G_trial = slope_u_on_g,
    slope_logLDL_on_G_trial = safe_slope(trial$log_LDL, trial$G),
    slope_Bproxy_on_G_trial = safe_slope(trial$B_proxy_z, trial$G),
    approx_induced_AxG = pars$beta_AU * slope_u_on_g,
    mean_LDL_trial = mean(trial$LDL),
    median_LDL_trial = median(trial$LDL),
    mean_B_z_trial = mean(trial$B_z),
    mean_B_proxy_z_trial = mean(trial$B_proxy_z)
  ) %>%
    bind_cols(cont_int) %>%
    bind_cols(tab_wide) %>%
    bind_cols(eligibility_diagnostics(source, trial_spec))

  pattern <- trial_pattern(tab, cox_tab, trial_spec$grs_scheme)
  analytic_tab <- analytic_group_risks(trial, pars)
  analytic_contrasts <- analytic_high_low_contrasts(analytic_tab, trial_spec$grs_scheme)

  list(
    summary = summary,
    event_table = tab,
    cox_by_grs = cox_tab,
    diagnostics_by_grs = diagnostics_by_grs(trial),
    analytic_group_risks = analytic_tab,
    analytic_high_low_contrasts = analytic_contrasts,
    continuous_interaction = fits$continuous,
    categorical_interaction = fits$categorical,
    cox_by_grs_adjusted = cox_by_grs_adjusted_all(trial),
    continuous_interactions_adjusted = continuous_interactions_adjusted_all(trial),
    categorical_interactions_adjusted = categorical_interactions_adjusted_all(trial),
    pattern = pattern
  )
}

# ----------------------------
# One full simulation
# ----------------------------

simulate_once <- function(trial_spec, pars = pars_anchor, seed = NULL, keep_source = TRUE) {
  if (!is.null(seed)) set.seed(seed)

  source <- simulate_source(trial_spec$n_source, pars) %>%
    apply_trial_eligibility(trial_spec)

  trial <- source %>%
    draw_trial(trial_spec) %>%
    simulate_events(pars, trial_spec)

  analysis <- summarise_trial(
    trial = trial,
    source = source,
    pars = pars,
    trial_spec = trial_spec
  )

  list(
    source = if (keep_source) source else NULL,
    trial = trial,
    analysis = analysis,
    pars = pars,
    trial_spec = trial_spec
  )
}


# ----------------------------
# Monte Carlo replicate helpers
# ----------------------------

mc_interval <- function(x, probs = c(0.025, 0.975)) {
  if (all(is.na(x))) return(c(NA_real_, NA_real_))
  as.numeric(quantile(x, probs = probs, na.rm = TRUE, names = FALSE))
}

mc_mean <- function(x) {
  mean(x, na.rm = TRUE)
}

run_trial_replicates <- function(
  trial_spec,
  pars = pars_anchor,
  n_reps = n_reps,
  base_seed = 20260425
) {
  rep_ids <- seq_len(n_reps)
  sim_list <- vector("list", length(rep_ids))

  for (r in rep_ids) {
    sim_list[[r]] <- simulate_once(
      trial_spec = trial_spec,
      pars = pars,
      seed = base_seed + r,
      keep_source = FALSE
    )
  }

  list(
    trial_spec = trial_spec,
    pars = pars,
    n_reps = n_reps,
    summaries = map2_dfr(sim_list, rep_ids, ~ mutate(.x$analysis$summary, rep = .y, .before = 1)),
    event_tables = map2_dfr(sim_list, rep_ids, ~ mutate(.x$analysis$event_table, rep = .y, .before = 1)),
    cox_tables = map2_dfr(sim_list, rep_ids, ~ mutate(.x$analysis$cox_by_grs, rep = .y, .before = 1)),
    diagnostics_by_grs = map2_dfr(sim_list, rep_ids, ~ mutate(.x$analysis$diagnostics_by_grs, rep = .y, .before = 1)),
    analytic_tables = map2_dfr(sim_list, rep_ids, ~ mutate(.x$analysis$analytic_group_risks, rep = .y, .before = 1)),
    analytic_contrasts = map2_dfr(sim_list, rep_ids, ~ mutate(.x$analysis$analytic_high_low_contrasts, rep = .y, .before = 1)),
    cox_tables_adjusted = map2_dfr(sim_list, rep_ids, ~ mutate(.x$analysis$cox_by_grs_adjusted, rep = .y, .before = 1)),
    continuous_interactions_adjusted = map2_dfr(sim_list, rep_ids, ~ mutate(.x$analysis$continuous_interactions_adjusted, rep = .y, .before = 1)),
    categorical_interactions_adjusted = map2_dfr(sim_list, rep_ids, ~ mutate(.x$analysis$categorical_interactions_adjusted, rep = .y, .before = 1)),
    patterns = map2_dfr(sim_list, rep_ids, ~ mutate(.x$analysis$pattern, rep = .y, .before = 1))
  )
}

summarise_analytic_group_risks_replicates <- function(rep_result) {
  rep_result$analytic_tables %>%
    group_by(GRS_group) %>%
    summarise(
      analytic_placebo_risk_low = mc_interval(analytic_placebo_risk)[1],
      analytic_placebo_risk_high = mc_interval(analytic_placebo_risk)[2],
      analytic_placebo_risk = mc_mean(analytic_placebo_risk),

      analytic_treatment_risk_low = mc_interval(analytic_treatment_risk)[1],
      analytic_treatment_risk_high = mc_interval(analytic_treatment_risk)[2],
      analytic_treatment_risk = mc_mean(analytic_treatment_risk),

      analytic_ARR_low = mc_interval(analytic_ARR)[1],
      analytic_ARR_high = mc_interval(analytic_ARR)[2],
      analytic_ARR = mc_mean(analytic_ARR),

      analytic_treatment_HR_cumhaz_low = mc_interval(analytic_treatment_HR_cumhaz)[1],
      analytic_treatment_HR_cumhaz_high = mc_interval(analytic_treatment_HR_cumhaz)[2],
      analytic_treatment_HR_cumhaz = mc_mean(analytic_treatment_HR_cumhaz),

      analytic_mean_individual_treatment_HR = mc_mean(analytic_mean_individual_treatment_HR),
      analytic_mean_U_response = mc_mean(analytic_mean_U_response),
      analytic_mean_G = mc_mean(analytic_mean_G),
      .groups = "drop"
    )
}

summarise_risk_reduction_replicates <- function(rep_result) {
  event_rep <- rep_result$event_tables
  cox_rep <- rep_result$cox_tables

  risk_summary <- event_rep %>%
    group_by(GRS_group) %>%
    summarise(
      n_mean = mc_mean(n),
      risk_control_low = mc_interval(risk_control)[1],
      risk_control_high = mc_interval(risk_control)[2],
      risk_control = mc_mean(risk_control),
      risk_treat_low = mc_interval(risk_treat)[1],
      risk_treat_high = mc_interval(risk_treat)[2],
      risk_treat = mc_mean(risk_treat),
      ARR_low = mc_interval(ARR)[1],
      ARR_high = mc_interval(ARR)[2],
      ARR = mc_mean(ARR),
      .groups = "drop"
    )

  hr_summary <- cox_rep %>%
    group_by(GRS_group) %>%
    summarise(
      HR_low = mc_interval(HR)[1],
      HR_high = mc_interval(HR)[2],
      HR = mc_mean(HR),
      .groups = "drop"
    )

  analytic_summary <- summarise_analytic_group_risks_replicates(rep_result)

  risk_summary %>%
    left_join(hr_summary, by = "GRS_group") %>%
    left_join(analytic_summary, by = "GRS_group") %>%
    mutate(
      RRR = 1 - HR,
      label = sprintf(
        "HR %.2f (%.2f-%.2f)\nRisk difference %.1f%%",
        HR, HR_low, HR_high, 100 * ARR
      )
    )
}


summarise_anchor_replicates <- function(rep_result) {
  s <- rep_result$summaries
  rep_ids <- sort(unique(rep_result$event_tables$rep))

  contrast_rep <- map_dfr(rep_ids, function(r) {
    event_df <- rep_result$event_tables %>% filter(rep == r)
    cox_df <- rep_result$cox_tables %>% filter(rep == r)

    fake_sim <- list(
      analysis = list(event_table = event_df, cox_by_grs = cox_df),
      trial_spec = rep_result$trial_spec
    )

    extract_contrast(fake_sim) %>%
      mutate(rep = r, .before = 1)
  })

  analytic_contrast_rep <- rep_result$analytic_contrasts

  s %>%
    select(
      rep,
      trial_name,
      beta_AxG,
      beta_AxG_low,
      beta_AxG_high,
      p_AxG,
      cor_G_U_trial,
      cor_G_logLDL_trial,
      slope_U_on_G_trial,
      approx_induced_AxG,
      event_rate_control,
      event_rate_treat,
      baseline_cumhaz_used,
      target_control_event_rate,
      eligible_fraction
    ) %>%
    left_join(contrast_rep, by = "rep") %>%
    left_join(analytic_contrast_rep, by = "rep", suffix = c("", "_analytic")) %>%
    group_by(trial_name) %>%
    summarise(
      n_reps = n(),
      beta_AxG = mc_mean(beta_AxG),
      beta_AxG_mc_low = mc_interval(beta_AxG)[1],
      beta_AxG_mc_high = mc_interval(beta_AxG)[2],
      approx_induced_AxG = mc_mean(approx_induced_AxG),
      approx_induced_AxG_mc_low = mc_interval(approx_induced_AxG)[1],
      approx_induced_AxG_mc_high = mc_interval(approx_induced_AxG)[2],
      cor_G_U_trial = mc_mean(cor_G_U_trial),
      cor_G_logLDL_trial = mc_mean(cor_G_logLDL_trial),
      slope_U_on_G_trial = mc_mean(slope_U_on_G_trial),
      event_rate_control = mc_mean(event_rate_control),
      event_rate_treat = mc_mean(event_rate_treat),
      HR_lower = mc_mean(HR_lower),
      HR_high = mc_mean(HR_high),
      HR_ratio_high_vs_lower = mc_mean(HR_ratio_high_vs_lower),
      ARR_lower = mc_mean(ARR_lower),
      ARR_high = mc_mean(ARR_high),
      ARR_difference_high_minus_lower = mc_mean(ARR_difference_high_minus_lower),

      analytic_placebo_risk_lower = mc_mean(analytic_placebo_risk_lower),
      analytic_placebo_risk_high = mc_mean(analytic_placebo_risk_high),
      analytic_treatment_risk_lower = mc_mean(analytic_treatment_risk_lower),
      analytic_treatment_risk_high = mc_mean(analytic_treatment_risk_high),
      analytic_placebo_high_vs_lower_cumhaz_HR =
        mc_mean(analytic_placebo_high_vs_lower_cumhaz_HR),
      analytic_treatment_high_vs_lower_cumhaz_HR =
        mc_mean(analytic_treatment_high_vs_lower_cumhaz_HR),
      analytic_treatment_HR_lower = mc_mean(analytic_treatment_HR_lower),
      analytic_treatment_HR_high = mc_mean(analytic_treatment_HR_high),
      analytic_treatment_HR_ratio_high_vs_lower =
        mc_mean(analytic_treatment_HR_ratio_high_vs_lower),
      analytic_ARR_lower = mc_mean(analytic_ARR_lower),
      analytic_ARR_high = mc_mean(analytic_ARR_high),
      analytic_ARR_difference_high_minus_lower =
        mc_mean(analytic_ARR_difference_high_minus_lower),

      baseline_cumhaz_used = mc_mean(baseline_cumhaz_used),
      target_control_event_rate = mc_mean(target_control_event_rate),
      eligible_fraction = mc_mean(eligible_fraction),
      .groups = "drop"
    )
}

summarise_anchor_analytic_tables <- function(rep_result) {
  summarise_analytic_group_risks_replicates(rep_result) %>%
    mutate(trial_name = rep_result$trial_spec$name, .before = 1)
}

summarise_anchor_analytic_contrasts <- function(rep_result) {
  rep_result$analytic_contrasts %>%
    group_by(lower_group, high_group) %>%
    summarise(
      n_reps = n(),
      analytic_placebo_risk_lower = mc_mean(analytic_placebo_risk_lower),
      analytic_placebo_risk_high = mc_mean(analytic_placebo_risk_high),
      analytic_treatment_risk_lower = mc_mean(analytic_treatment_risk_lower),
      analytic_treatment_risk_high = mc_mean(analytic_treatment_risk_high),
      analytic_placebo_high_vs_lower_cumhaz_HR =
        mc_mean(analytic_placebo_high_vs_lower_cumhaz_HR),
      analytic_placebo_high_vs_lower_cumhaz_HR_low =
        mc_interval(analytic_placebo_high_vs_lower_cumhaz_HR)[1],
      analytic_placebo_high_vs_lower_cumhaz_HR_high =
        mc_interval(analytic_placebo_high_vs_lower_cumhaz_HR)[2],
      analytic_treatment_high_vs_lower_cumhaz_HR =
        mc_mean(analytic_treatment_high_vs_lower_cumhaz_HR),
      analytic_treatment_high_vs_lower_cumhaz_HR_low =
        mc_interval(analytic_treatment_high_vs_lower_cumhaz_HR)[1],
      analytic_treatment_high_vs_lower_cumhaz_HR_high =
        mc_interval(analytic_treatment_high_vs_lower_cumhaz_HR)[2],
      analytic_treatment_HR_lower = mc_mean(analytic_treatment_HR_lower),
      analytic_treatment_HR_high = mc_mean(analytic_treatment_HR_high),
      analytic_treatment_HR_ratio_high_vs_lower =
        mc_mean(analytic_treatment_HR_ratio_high_vs_lower),
      analytic_ARR_lower = mc_mean(analytic_ARR_lower),
      analytic_ARR_high = mc_mean(analytic_ARR_high),
      analytic_ARR_difference_high_minus_lower =
        mc_mean(analytic_ARR_difference_high_minus_lower),
      .groups = "drop"
    ) %>%
    mutate(trial_name = rep_result$trial_spec$name, .before = 1)
}


summarise_anchor_event_tables <- function(rep_result) {
  summarise_risk_reduction_replicates(rep_result) %>%
    mutate(trial_name = rep_result$trial_spec$name, .before = 1)
}

# ----------------------------
# Plotting
# ----------------------------

plot_risk_reduction_by_grs <- function(rep_result) {
  trial_spec <- rep_result$trial_spec
  event_tab <- summarise_risk_reduction_replicates(rep_result)

  arm_levels <- c(trial_spec$control_label, trial_spec$treatment_label)
  arm_colours <- setNames(c("#4C78A8", "#D62728"), arm_levels)

  plot_dat <- event_tab %>%
    select(
      GRS_group,
      risk_control,
      risk_control_low,
      risk_control_high,
      risk_treat,
      risk_treat_low,
      risk_treat_high
    ) %>%
    pivot_longer(
      cols = -GRS_group,
      names_to = "quantity",
      values_to = "value"
    ) %>%
    mutate(
      arm = case_when(
        grepl("control", quantity) ~ trial_spec$control_label,
        grepl("treat", quantity) ~ trial_spec$treatment_label,
        TRUE ~ NA_character_
      ),
      statistic = case_when(
        grepl("_low$", quantity) ~ "low",
        grepl("_high$", quantity) ~ "high",
        TRUE ~ "mean"
      )
    ) %>%
    select(GRS_group, arm, statistic, value) %>%
    pivot_wider(names_from = statistic, values_from = value) %>%
    mutate(
      arm = factor(arm, levels = arm_levels),
      risk_pct = 100 * mean,
      risk_low_pct = 100 * low,
      risk_high_pct = 100 * high
    )

  ann_dat <- event_tab %>%
    mutate(y = 100 * pmax(risk_control_high, risk_treat_high) + 3.0)

  # Analytic event risks implied by the individual-level data-generating model.
  # These values average the model-implied risks within each PRS group before
  # binomial event sampling, and therefore show the expected pattern without
  # Monte Carlo outcome noise.
  analytic_pred_both <- bind_rows(
    event_tab %>%
      transmute(
        GRS_group,
        risk_pred_pct = 100 * analytic_placebo_risk,
        arm = factor(trial_spec$control_label, levels = arm_levels)
      ),
    event_tab %>%
      transmute(
        GRS_group,
        risk_pred_pct = 100 * analytic_treatment_risk,
        arm = factor(trial_spec$treatment_label, levels = arm_levels)
      )
  )

  y_max <- max(
    c(plot_dat$risk_high_pct, ann_dat$y, analytic_pred_both$risk_pred_pct),
    na.rm = TRUE
  ) + 5

  ggplot(plot_dat, aes(x = GRS_group, y = risk_pct, fill = arm)) +
    geom_col(
      position = position_dodge(width = 0.72),
      width = 0.62,
      colour = "grey30",
      linewidth = 0.25
    ) +
    geom_errorbar(
      aes(ymin = risk_low_pct, ymax = risk_high_pct, colour = arm),
      position = position_dodge(width = 0.72),
      width = 0.18,
      linewidth = 0.55,
      show.legend = FALSE
    ) +
    geom_text(
      data = ann_dat,
      aes(x = GRS_group, y = y, label = label),
      inherit.aes = FALSE,
      size = 3.1,
      lineheight = 0.92
    ) +
    geom_point(
      data = analytic_pred_both,
      aes(x = GRS_group, y = risk_pred_pct, group = arm),
      position = position_dodge(width = 0.72),
      shape = 23,
      size = 3.8,
      fill = "#B04A1F",
      colour = "#B04A1F",
      inherit.aes = FALSE
    ) +
    scale_fill_manual(values = arm_colours) +
    scale_colour_manual(values = arm_colours) +
    coord_cartesian(ylim = c(0, y_max), clip = "off") +
    theme_bw(base_size = 13) +
    theme(
      legend.position = "top",
      legend.title = element_blank(),
      panel.grid.major.x = element_blank(),
      panel.grid.minor = element_blank(),
      axis.title.x = element_blank(),
      plot.title = element_text(face = "bold"),
      plot.subtitle = element_text(size = 10.5),
      plot.caption = element_text(hjust = 0, size = 9),
      plot.margin = margin(8, 80, 8, 8)
    ) +
    labs(
      title = paste0(trial_spec$name, "-like simulation: event risk by CAD PRS group"),
      subtitle = paste0(
        "Mean event risk across ", rep_result$n_reps,
        " Monte Carlo replicates at the actual trial sample size; error bars show simulation-to-simulation uncertainty."
      ),
      y = paste0(trial_spec$outcome_label, " risk (%)"),
      caption = "Orange diamonds show model-implied event risks before outcome sampling. They are calculated directly from the individual-level cumulative-hazard model, with zero true PRS-by-treatment interaction."
    )
}


plot_adjusted_hr_comparison <- function(sim) {
  trial_spec <- sim$trial_spec

  plot_dat <- sim$analysis$cox_by_grs_adjusted %>%
    mutate(
      model_set = factor(
        model_set,
        levels = names(adjustment_sets),
        labels = c(
          "Unadjusted",
          "LDL adjusted",
          "Burden-proxy adjusted",
          "LDL + burden-proxy adjusted"
        )
      )
    )

  ggplot(plot_dat, aes(x = GRS_group, y = HR, ymin = HR_low, ymax = HR_high)) +
    geom_hline(yintercept = 1, linetype = "dashed", colour = "grey40") +
    geom_pointrange(position = position_dodge(width = 0.55)) +
    facet_wrap(~ model_set, ncol = 2) +
    coord_cartesian(ylim = c(0.35, 1.25)) +
    theme_bw(base_size = 12) +
    theme(
      panel.grid.major.x = element_blank(),
      panel.grid.minor = element_blank(),
      axis.title.x = element_blank(),
      plot.title = element_text(face = "bold")
    ) +
    labs(
      title = paste0(trial_spec$name, "-like subgroup HRs under Cox adjustment"),
      subtitle = "Treatment HRs by PRS group; true treatment-by-PRS interaction = 0",
      y = "Treatment hazard ratio"
    )
}

# ----------------------------
# Compact printing
# ----------------------------

print_simulation_results <- function(sim) {
  cat("\n============================================================\n")
  cat(sim$trial_spec$name, "simulation summary\n")
  cat("============================================================\n")

  cat("\nHeadline summary:\n")
  print(sim$analysis$summary)

  cat("\nDiagnostics by GRS group:\n")
  print(sim$analysis$diagnostics_by_grs)

  cat("\nEvent table by GRS group:\n")
  print(sim$analysis$event_table)

  cat("\nCox HRs by GRS group, unadjusted:\n")
  print(sim$analysis$cox_by_grs)

  cat("\nCox HRs by GRS group, adjusted model sets:\n")
  print(sim$analysis$cox_by_grs_adjusted)

  cat("\nContinuous A x G interaction, adjusted model sets:\n")
  print(sim$analysis$continuous_interactions_adjusted)

  cat("\nCategorical A x GRS interaction, adjusted model sets:\n")
  print(sim$analysis$categorical_interactions_adjusted)

  cat("\nHeadline-like pattern flags based on unadjusted subgroup HRs:\n")
  print(sim$analysis$pattern)
}

# ----------------------------
# Sensitivity to treatment-by-U heterogeneity
# ----------------------------

extract_contrast <- function(sim) {
  x <- sim$analysis$event_table %>%
    left_join(sim$analysis$cox_by_grs, by = "GRS_group")

  scheme <- sim$trial_spec$grs_scheme

  if (scheme == "fourier_3") {
    low <- x %>% filter(GRS_group == "Low")
    high <- x %>% filter(GRS_group == "High")
    tibble(
      lower_group = "Low",
      high_group = "High",
      HR_lower = low$HR,
      HR_high = high$HR,
      HR_ratio_high_vs_lower = high$HR / low$HR,
      ARR_lower = low$ARR,
      ARR_high = high$ARR,
      ARR_difference_high_minus_lower = high$ARR - low$ARR
    )

  } else if (scheme == "top20_vs_rest") {
    rest <- x %>% filter(GRS_group == "All others")
    high <- x %>% filter(GRS_group == "High")
    tibble(
      lower_group = "All others",
      high_group = "High",
      HR_lower = rest$HR,
      HR_high = high$HR,
      HR_ratio_high_vs_lower = high$HR / rest$HR,
      ARR_lower = rest$ARR,
      ARR_high = high$ARR,
      ARR_difference_high_minus_lower = high$ARR - rest$ARR
    )

  } else if (scheme == "top10_vs_rest") {
    rest <- x %>% filter(GRS_group == "Lower")
    high <- x %>% filter(GRS_group == "High")
    tibble(
      lower_group = "Lower",
      high_group = "High",
      HR_lower = rest$HR,
      HR_high = high$HR,
      HR_ratio_high_vs_lower = high$HR / rest$HR,
      ARR_lower = rest$ARR,
      ARR_high = high$ARR,
      ARR_difference_high_minus_lower = high$ARR - rest$ARR
    )

  } else {
    stop("Unknown grs_scheme: ", scheme)
  }
}

run_beta_AU_grid <- function(
  beta_AU_grid = c(0, 0.10, 0.20, 0.30, 0.40),
  specs = trial_specs,
  pars = pars_anchor,
  n_reps = n_reps,
  base_seed = 41000
) {
  raw <- list()
  k <- 1L

  for (j in seq_along(beta_AU_grid)) {
    pars_j <- pars
    pars_j$beta_AU <- beta_AU_grid[j]

    for (i in seq_along(specs)) {
      spec_i <- specs[[i]]

      for (r in seq_len(n_reps)) {
        sim_i <- simulate_once(
          trial_spec = spec_i,
          pars = pars_j,
          seed = base_seed + 100000L * j + 1000L * i + r,
          keep_source = FALSE
        )

        raw[[k]] <- sim_i$analysis$summary %>%
          select(
            trial_name,
            beta_AxG,
            beta_AxG_low_model = beta_AxG_low,
            beta_AxG_high_model = beta_AxG_high,
            p_AxG,
            cor_G_U_trial,
            cor_G_logLDL_trial,
            slope_U_on_G_trial,
            approx_induced_AxG,
            event_rate_control,
            event_rate_treat
          ) %>%
          bind_cols(extract_contrast(sim_i)) %>%
          mutate(
            beta_AU = beta_AU_grid[j],
            rep = r,
            .before = 1
          )

        k <- k + 1L
      }
    }
  }

  bind_rows(raw) %>%
    group_by(beta_AU, trial_name) %>%
    summarise(
      n_reps = n(),
      beta_AxG = mc_mean(beta_AxG),
      beta_AxG_low = mc_interval(beta_AxG)[1],
      beta_AxG_high = mc_interval(beta_AxG)[2],
      prob_apparent_benefit_higher_PRS = mean(beta_AxG < 0, na.rm = TRUE),
      approx_induced_AxG = mc_mean(approx_induced_AxG),
      approx_induced_AxG_low = mc_interval(approx_induced_AxG)[1],
      approx_induced_AxG_high = mc_interval(approx_induced_AxG)[2],
      cor_G_U_trial = mc_mean(cor_G_U_trial),
      cor_G_logLDL_trial = mc_mean(cor_G_logLDL_trial),
      slope_U_on_G_trial = mc_mean(slope_U_on_G_trial),
      event_rate_control = mc_mean(event_rate_control),
      event_rate_treat = mc_mean(event_rate_treat),
      HR_lower = mc_mean(HR_lower),
      HR_high = mc_mean(HR_high),
      HR_ratio_high_vs_lower = mc_mean(HR_ratio_high_vs_lower),
      ARR_lower = mc_mean(ARR_lower),
      ARR_high = mc_mean(ARR_high),
      ARR_difference_high_minus_lower = mc_mean(ARR_difference_high_minus_lower),
      .groups = "drop"
    )
}

run_beta_AU_grid_adjusted <- function(
  beta_AU_grid = c(0, 0.10, 0.20, 0.30, 0.40),
  specs = trial_specs,
  pars = pars_anchor,
  n_reps = n_reps,
  base_seed = 51000
) {
  raw <- list()
  k <- 1L

  for (j in seq_along(beta_AU_grid)) {
    pars_j <- pars
    pars_j$beta_AU <- beta_AU_grid[j]

    for (i in seq_along(specs)) {
      spec_i <- specs[[i]]

      for (r in seq_len(n_reps)) {
        sim_i <- simulate_once(
          trial_spec = spec_i,
          pars = pars_j,
          seed = base_seed + 100000L * j + 1000L * i + r,
          keep_source = FALSE
        )

        raw[[k]] <- sim_i$analysis$continuous_interactions_adjusted %>%
          mutate(
            beta_AU = beta_AU_grid[j],
            trial_name = spec_i$name,
            rep = r,
            .before = 1
          )

        k <- k + 1L
      }
    }
  }

  bind_rows(raw) %>%
    group_by(beta_AU, trial_name, model_set, adjustment) %>%
    summarise(
      n_reps = n(),
      beta_AxG = mc_mean(beta_AxG),
      beta_AxG_low = mc_interval(beta_AxG)[1],
      beta_AxG_high = mc_interval(beta_AxG)[2],
      prob_apparent_benefit_higher_PRS = mean(beta_AxG < 0, na.rm = TRUE),
      median_model_p = median(p_AxG, na.rm = TRUE),
      .groups = "drop"
    )
}


# ----------------------------
# Manuscript-grade sensitivity figures
# ----------------------------

label_model_set <- function(x) {
  dplyr::recode(
    x,
    unadjusted = "Unadjusted",
    ldl_adjusted = "LDL adjusted",
    burden_proxy_adjusted = "Burden proxy adjusted",
    ldl_burden_proxy_adjusted = "LDL + burden proxy adjusted",
    .default = x
  )
}

plot_beta_AU_sensitivity <- function(beta_AU_sensitivity) {
  observed_dat <- beta_AU_sensitivity %>%
    transmute(
      trial_name,
      hidden_heterogeneity = beta_AU,
      estimate = beta_AxG,
      conf.low = beta_AxG_low,
      conf.high = beta_AxG_high,
      series = "Full simulated trial analysis"
    )

  approx_dat <- beta_AU_sensitivity %>%
    transmute(
      trial_name,
      hidden_heterogeneity = beta_AU,
      estimate = approx_induced_AxG,
      conf.low = approx_induced_AxG_low,
      conf.high = approx_induced_AxG_high,
      series = "Mechanism-only prediction"
    )

  plot_dat <- bind_rows(observed_dat, approx_dat) %>%
    mutate(
      trial_name = factor(trial_name, levels = c("FOURIER", "ODYSSEY", "WOSCOPS")),
      series = factor(
        series,
        levels = c("Full simulated trial analysis", "Mechanism-only prediction")
      )
    )

  ggplot(plot_dat, aes(x = hidden_heterogeneity, y = estimate, group = series, colour = series, linetype = series)) +
    geom_hline(yintercept = 0, linewidth = 0.4, colour = "grey55") +
    geom_ribbon(
      data = plot_dat %>% filter(series == "Mechanism-only prediction"),
      aes(ymin = conf.low, ymax = conf.high, fill = series),
      alpha = 0.10,
      colour = NA,
      show.legend = FALSE
    ) +
    geom_errorbar(
      data = plot_dat %>% filter(series == "Full simulated trial analysis"),
      aes(ymin = conf.low, ymax = conf.high),
      width = 0.015,
      linewidth = 0.45,
      alpha = 0.9
    ) +
    geom_line(linewidth = 0.85) +
    geom_point(size = 2.2) +
    scale_x_continuous(breaks = sort(unique(plot_dat$hidden_heterogeneity))) +
    scale_colour_manual(values = c(
      "Full simulated trial analysis" = "#1B3A57",
      "Mechanism-only prediction" = "#B04A1F"
    )) +
    scale_fill_manual(values = c("Mechanism-only prediction" = "#B04A1F")) +
    scale_linetype_manual(values = c(
      "Full simulated trial analysis" = "solid",
      "Mechanism-only prediction" = "22"
    )) +
    facet_wrap(~ trial_name, nrow = 1) +
    labs(
      title = "Hidden treatment-response differences can create an apparent PRS treatment effect",
      subtitle = "Negative values mean that higher CAD PRS appears to receive greater relative treatment benefit, even though the true PRS-by-treatment effect is zero.",
      x = "Strength of hidden treatment-response heterogeneity",
      y = "Apparent PRS-related treatment benefit\n(log hazard-ratio change per PRS SD)",
      colour = NULL,
      linetype = NULL,
      caption = "Points show Monte Carlo means; bars and bands show simulation-to-simulation uncertainty. The mechanism-only line uses the eligibility-induced PRS-comorbidity imbalance."
    ) +
    theme_minimal(base_size = 12) +
    theme(
      plot.title = element_text(face = "bold"),
      plot.subtitle = element_text(size = 10.5),
      strip.text = element_text(face = "bold"),
      panel.grid.minor = element_blank(),
      legend.position = "bottom",
      legend.box = "vertical",
      axis.title = element_text(face = "bold"),
      plot.caption = element_text(hjust = 0, size = 9)
    )
}


plot_beta_AU_sensitivity_adjusted <- function(beta_AU_sensitivity_adjusted) {
  plot_dat <- beta_AU_sensitivity_adjusted %>%
    mutate(
      trial_name = factor(trial_name, levels = c("FOURIER", "ODYSSEY", "WOSCOPS")),
      model_set = factor(
        label_model_set(model_set),
        levels = c(
          "Unadjusted",
          "LDL adjusted",
          "Burden proxy adjusted",
          "LDL + burden proxy adjusted"
        )
      )
    )

  ggplot(plot_dat, aes(x = beta_AU, y = beta_AxG, colour = model_set, group = model_set)) +
    geom_hline(yintercept = 0, linewidth = 0.4, colour = "grey55") +
    geom_errorbar(aes(ymin = beta_AxG_low, ymax = beta_AxG_high), width = 0.015, linewidth = 0.42, alpha = 0.9) +
    geom_line(linewidth = 0.8) +
    geom_point(size = 1.9) +
    scale_x_continuous(breaks = sort(unique(plot_dat$beta_AU))) +
    scale_colour_manual(values = c(
      "Unadjusted" = "#1B3A57",
      "LDL adjusted" = "#2A7F62",
      "Burden proxy adjusted" = "#B04A1F",
      "LDL + burden proxy adjusted" = "#6B5B95"
    )) +
    facet_wrap(~ trial_name, nrow = 1) +
    labs(
      title = "Adjustment does not necessarily remove the apparent PRS treatment effect",
      subtitle = "Each line repeats the simulated trial analysis after adjusting for baseline LDL-C, a measured burden proxy, or both.",
      x = "Strength of hidden treatment-response heterogeneity",
      y = "Apparent PRS-related treatment benefit\n(log hazard-ratio change per PRS SD)",
      colour = NULL,
      caption = "Points show Monte Carlo means; bars show simulation-to-simulation uncertainty. Negative values imply apparently stronger relative treatment benefit at higher CAD PRS."
    ) +
    theme_minimal(base_size = 12) +
    theme(
      plot.title = element_text(face = "bold"),
      plot.subtitle = element_text(size = 10.5),
      strip.text = element_text(face = "bold"),
      panel.grid.minor = element_blank(),
      legend.position = "bottom",
      legend.box = "vertical",
      axis.title = element_text(face = "bold"),
      plot.caption = element_text(hjust = 0, size = 9)
    )
}


# ============================================================
# Run anchor Monte Carlo simulations
# ============================================================

replicate_results <- imap(
  trial_specs,
  ~ run_trial_replicates(
    trial_spec = .x,
    pars = pars_anchor,
    n_reps = n_reps,
    base_seed = 20260425 + 10000L * match(.y, names(trial_specs))
  )
)

anchor_summary <- map_dfr(replicate_results, summarise_anchor_replicates)
anchor_event_tables <- map_dfr(replicate_results, summarise_anchor_event_tables)
anchor_diagnostics_by_grs <- imap_dfr(replicate_results, ~ mutate(.x$diagnostics_by_grs, trial_name = .y, .before = 1))
anchor_analytic_tables <- map_dfr(replicate_results, summarise_anchor_analytic_tables)
anchor_analytic_contrasts <- map_dfr(replicate_results, summarise_anchor_analytic_contrasts)
anchor_cox_tables_adjusted <- imap_dfr(replicate_results, ~ mutate(.x$cox_tables_adjusted, trial_name = .y, .before = 1))
anchor_continuous_interactions_adjusted <- imap_dfr(replicate_results, ~ mutate(.x$continuous_interactions_adjusted, trial_name = .y, .before = 1))
anchor_categorical_interactions_adjusted <- imap_dfr(replicate_results, ~ mutate(.x$categorical_interactions_adjusted, trial_name = .y, .before = 1))

print(anchor_summary)
print(anchor_event_tables)
print(anchor_diagnostics_by_grs)
print(anchor_analytic_tables)
print(anchor_analytic_contrasts)
print(anchor_cox_tables_adjusted)
print(anchor_continuous_interactions_adjusted)
print(anchor_categorical_interactions_adjusted)

if (write_plots) {
  if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

  plots <- map(replicate_results, plot_risk_reduction_by_grs)

  walk2(plots, trial_specs, function(p, spec) {
    ggsave(
      filename = file.path(output_dir, spec$plot_filename),
      plot = p,
      width = 9.2,
      height = 5.8,
      dpi = 320
    )
  })

  write.csv(anchor_summary, file.path(output_dir, "anchor_summary.csv"), row.names = FALSE)
  write.csv(anchor_event_tables, file.path(output_dir, "anchor_event_tables.csv"), row.names = FALSE)
  write.csv(anchor_diagnostics_by_grs, file.path(output_dir, "anchor_diagnostics_by_grs.csv"), row.names = FALSE)
  write.csv(anchor_analytic_tables, file.path(output_dir, "anchor_analytic_group_risks.csv"), row.names = FALSE)
  write.csv(anchor_analytic_contrasts, file.path(output_dir, "anchor_analytic_high_low_contrasts.csv"), row.names = FALSE)
  write.csv(anchor_cox_tables_adjusted, file.path(output_dir, "anchor_cox_tables_adjusted_raw.csv"), row.names = FALSE)
  write.csv(anchor_continuous_interactions_adjusted, file.path(output_dir, "anchor_continuous_interactions_adjusted_raw.csv"), row.names = FALSE)
  write.csv(anchor_categorical_interactions_adjusted, file.path(output_dir, "anchor_categorical_interactions_adjusted_raw.csv"), row.names = FALSE)
}

if (run_beta_AU_sensitivity) {
  beta_AU_sensitivity <- run_beta_AU_grid(n_reps = n_reps)
  beta_AU_sensitivity_adjusted <- run_beta_AU_grid_adjusted(n_reps = n_reps)

  print(beta_AU_sensitivity)
  print(beta_AU_sensitivity_adjusted)

  if (write_plots) {
    if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

    sensitivity_plot <- plot_beta_AU_sensitivity(beta_AU_sensitivity)
    sensitivity_plot_adjusted <- plot_beta_AU_sensitivity_adjusted(beta_AU_sensitivity_adjusted)

    ggsave(
      filename = file.path(output_dir, "beta_AU_sensitivity_overview.png"),
      plot = sensitivity_plot,
      width = 10.8,
      height = 5.1,
      dpi = 320
    )

    ggsave(
      filename = file.path(output_dir, "beta_AU_sensitivity_adjusted_overview.png"),
      plot = sensitivity_plot_adjusted,
      width = 10.8,
      height = 5.1,
      dpi = 320
    )

    write.csv(beta_AU_sensitivity, file.path(output_dir, "beta_AU_sensitivity.csv"), row.names = FALSE)
    write.csv(beta_AU_sensitivity_adjusted, file.path(output_dir, "beta_AU_sensitivity_adjusted_models.csv"), row.names = FALSE)
  }
}
