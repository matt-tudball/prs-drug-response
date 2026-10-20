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

setwd(this.path::here())
source("../figures/fps_style.r")

set.seed(20260425)

# ----------------------------
# User controls
# ----------------------------

# Number of independent Monte Carlo replicates to run at the actual trial sample size.
# Use a small value such as 5 for quick checks and a larger value such as 100-500 for manuscript figures.
n_reps <- 1000L

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
    p_excluded_burden = 0.10,
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

# ----------------------------
# Non-genetic risk by PRS decile: general population vs trial-selected sample
# ----------------------------
#
# Eligibility depends on B (which depends on both G and U), so among those
# selected for the trial G and U are no longer independent, even though they
# are independent in the general population. This table makes that selection
# effect visible directly: mean U (the non-genetic risk factor) by decile of
# G (the CAD PRS), in the general source population versus the trial sample.
# Deciles are defined from the general population's distribution of G.

prs_decile_table <- function(trial, source) {
  decile_breaks <- quantile(source$G, probs = seq(0, 1, by = 0.1), na.rm = TRUE)
  decile_breaks[1] <- -Inf
  decile_breaks[length(decile_breaks)] <- Inf

  assign_decile <- function(dat) {
    dat %>%
      mutate(
        G_decile = as.integer(cut(G, breaks = decile_breaks, labels = FALSE, include.lowest = TRUE))
      )
  }

  source_tab <- source %>%
    assign_decile() %>%
    group_by(G_decile) %>%
    summarise(n = n(), mean_U = mean(U), .groups = "drop") %>%
    mutate(population = "General population", .before = 1)

  trial_tab <- trial %>%
    assign_decile() %>%
    group_by(G_decile) %>%
    summarise(n = n(), mean_U = mean(U), .groups = "drop") %>%
    mutate(population = "Selected for trial", .before = 1)

  bind_rows(source_tab, trial_tab)
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
    prs_decile_table = prs_decile_table(trial = trial, source = source),
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
    prs_decile_tables = map2_dfr(sim_list, rep_ids, ~ mutate(.x$analysis$prs_decile_table, rep = .y, .before = 1)),
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


summarise_prs_decile_replicates <- function(rep_result) {
  rep_result$prs_decile_tables %>%
    group_by(population, G_decile) %>%
    summarise(
      n_mean = mc_mean(n),
      mean_U_low = mc_interval(mean_U)[1],
      mean_U_high = mc_interval(mean_U)[2],
      mean_U = mc_mean(mean_U),
      .groups = "drop"
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
      beta_AxG_mc_low = mc_interval(beta_AxG)[1],
      beta_AxG_mc_high = mc_interval(beta_AxG)[2],
      beta_AxG = mc_mean(beta_AxG),
      approx_induced_AxG_mc_low = mc_interval(approx_induced_AxG)[1],
      approx_induced_AxG_mc_high = mc_interval(approx_induced_AxG)[2],
      approx_induced_AxG = mc_mean(approx_induced_AxG),
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
      analytic_placebo_high_vs_lower_cumhaz_HR_low =
        mc_interval(analytic_placebo_high_vs_lower_cumhaz_HR)[1],
      analytic_placebo_high_vs_lower_cumhaz_HR_high =
        mc_interval(analytic_placebo_high_vs_lower_cumhaz_HR)[2],
      analytic_placebo_high_vs_lower_cumhaz_HR =
        mc_mean(analytic_placebo_high_vs_lower_cumhaz_HR),
      analytic_treatment_high_vs_lower_cumhaz_HR_low =
        mc_interval(analytic_treatment_high_vs_lower_cumhaz_HR)[1],
      analytic_treatment_high_vs_lower_cumhaz_HR_high =
        mc_interval(analytic_treatment_high_vs_lower_cumhaz_HR)[2],
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

# Publication theme shared with the UKB/HRN meta-analysis figure
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

theme_pub_axes <- theme(
  axis.text.x = element_text(colour = fps$charcoal, lineheight = 1.05,
                             margin = margin(t = 6)),
  axis.text.y = element_text(colour = fps$charcoal),
  axis.title.x = element_text(size = 17, margin = margin(t = 10)),
  axis.title.y = element_text(size = 17, lineheight = 1.1,
                              margin = margin(r = 10)),
  strip.text = element_text(face = "bold", size = 18, colour = fps$charcoal,
                            margin = margin(b = 8)),
  panel.spacing.x = unit(1.4, "lines"),
  plot.margin = margin(12, 16, 8, 12)
)

plot_risk_reduction_by_grs <- function(rep_results_list,
                                        published_contrasts = published_trial_contrasts) {
  lower_labels <- setNames(
    published_contrasts$lower_group,
    published_contrasts$trial_name
  )

  sim_rows <- imap_dfr(rep_results_list, function(rep_result, tname) {
    low_label <- lower_labels[[tname]]
    summarise_risk_reduction_replicates(rep_result) %>%
      filter(GRS_group %in% c(low_label, "High")) %>%
      transmute(
        trial     = tname,
        GRS_group = ifelse(as.character(GRS_group) == "High",
                           "Higher PRS", "Lower PRS"),
        type      = "Simulated",
        HR        = HR,
        HR_low    = HR_low,
        HR_high   = HR_high
      )
  })

  pub_rows <- bind_rows(
    published_contrasts %>%
      transmute(
        trial     = trial_name,
        GRS_group = "Lower PRS",
        type      = "Observed",
        HR        = HR_lower_obs,
        HR_low    = HR_lower_obs_ci_lower,
        HR_high   = HR_lower_obs_ci_upper
      ),
    published_contrasts %>%
      transmute(
        trial     = trial_name,
        GRS_group = "Higher PRS",
        type      = "Observed",
        HR        = HR_high_obs,
        HR_low    = HR_high_obs_ci_lower,
        HR_high   = HR_high_obs_ci_upper
      )
  )

  n_sim <- rep_results_list[[1]]$n_reps

  all_dat <- bind_rows(sim_rows, pub_rows) %>%
    mutate(
      trial     = factor(trial,     levels = names(rep_results_list)),
      type      = factor(type,      levels = c("Observed", "Simulated")),
      GRS_group = factor(GRS_group, levels = c("Lower PRS", "Higher PRS"))
    )

  ggplot(all_dat, aes(x = GRS_group, y = HR, colour = type)) +
    geom_hline(yintercept = 1, linetype = "dotted", colour = fps$null, linewidth = 0.8) +
    geom_errorbar(
      aes(ymin = HR_low, ymax = HR_high),
      position    = position_dodge(width = 0.5),
      width       = 0.18,
      linewidth   = 1.2,
      show.legend = FALSE
    ) +
    geom_point(
      aes(fill = type),
      position = position_dodge(width = 0.5),
      shape    = 21,
      size     = 4.5,
      stroke   = 1,
      colour   = fps$charcoal
    ) +
    facet_wrap(~ trial, ncol = 3L) +
    scale_colour_manual(
      values = c("Observed" = fps$sand, "Simulated" = fps$navy),
      labels = c(
        "Observed"  = "Observed trial results",
        "Simulated" = paste0("Simulated trial results (", n_sim, " simulations)")
      ),
      aesthetics = c("colour", "fill")
    ) +
    scale_y_continuous(
      breaks = seq(0.4, 1.4, by = 0.2),
      labels = function(x) sprintf("%.1f", x)
    ) +
    coord_cartesian(ylim = c(0.30, 1.45)) +
    theme_pub +
    theme_pub_axes +
    theme(
      axis.title.x = element_blank(),
      legend.text  = element_text(size = 16, colour = fps$charcoal)
    ) +
    labs(
      y = "Treatment hazard ratio"
    )
}


plot_nongenetic_risk_by_prs_decile <- function(replicate_results) {
  plot_dat <- imap_dfr(replicate_results, function(rep_result, tname) {
    summarise_prs_decile_replicates(rep_result) %>%
      mutate(trial = tname, .before = 1)
  }) %>%
    mutate(
      trial      = factor(trial, levels = names(replicate_results)),
      population = factor(population, levels = c("General population", "Selected for trial"))
    )

  n_sim <- replicate_results[[1]]$n_reps

  ggplot(plot_dat, aes(x = G_decile, y = mean_U, colour = population, fill = population)) +
    geom_hline(yintercept = 0, linetype = "dashed", colour = fps$null, linewidth = 0.35) +
    geom_ribbon(
      aes(ymin = mean_U_low, ymax = mean_U_high),
      colour = NA,
      alpha  = 0.15
    ) +
    geom_line(linewidth = 0.85) +
    geom_point(size = 1.9) +
    scale_x_continuous(breaks = 1:10) +
    scale_colour_manual(
      values = c("General population" = fps$navy, "Selected for trial" = fps$sand),
      name   = NULL
    ) +
    scale_fill_manual(
      values = c("General population" = fps$navy, "Selected for trial" = fps$sand),
      name   = NULL
    ) +
    facet_wrap(~ trial, nrow = 1) +
    theme_fps(base_size = 12) +
    theme(
      legend.position    = "top",
      strip.text         = element_text(face = "bold"),
      panel.grid.minor   = element_blank(),
      axis.title         = element_text(face = "bold", size = 10),
      plot.caption       = element_text(hjust = 0, size = 9),
      plot.background    = element_rect(fill = fps$white, colour = NA)
    ) +
    labs(
      # title = "Selection into the trial induces a PRS-non-genetic-risk association\nthat is absent in the general population",
      # subtitle = "Mean non-genetic risk factor U by decile of the CAD PRS, in the general\npopulation and among those selected for the simulated trial.",
      x = "CAD PRS decile (general population)",
      y = "Average other risk factors"
      # caption = paste0(
      #   "Deciles of the CAD PRS are defined within each simulation's general population. ",
      #   "Points show the mean non-genetic risk factor within each decile, averaged across ",
      #   n_sim, " simulations; shaded bands show the 2.5-97.5th percentile interval across simulations."
      # )
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
    geom_hline(yintercept = 1, linetype = "dashed", colour = fps$null) +
    geom_pointrange(position = position_dodge(width = 0.55), colour = fps$charcoal) +
    facet_wrap(~ model_set, ncol = 2) +
    coord_cartesian(ylim = c(0.35, 1.25)) +
    theme_fps(base_size = 12) +
    theme(
      panel.grid.major.x = element_blank(),
      panel.grid.minor = element_blank(),
      axis.title.x = element_blank(),
      plot.title = element_text(face = "bold"),
      plot.subtitle = element_text(size = 10.5),
      plot.caption = element_text(hjust = 0, size = 9),
      plot.background = element_rect(fill = fps$white, colour = NA)
    ) +
    labs(
      # title = paste0(trial_spec$name, "-like subgroup HRs under Cox adjustment"),
      # subtitle = "Treatment HRs by PRS group; true treatment-by-PRS interaction = 0",
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
  base_seed = 41000,
  published_contrasts = published_trial_contrasts
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
    left_join(
      published_contrasts %>%
        select(
          trial_name,
          HR_lower_obs,
          HR_high_obs,
          observed_HR_ratio_high_vs_lower,
          observed_log_HR_ratio_high_vs_lower,
          observed_HR_ratio_high_vs_lower_ci_lower,
          observed_HR_ratio_high_vs_lower_ci_upper,
          contrast_note
        ),
      by = "trial_name"
    ) %>%
    mutate(
      # A replicate is compatible when its simulated high-vs-lower PRS treatment
      # HR ratio falls inside the published trial's own 95% CI for that
      # contrast, rather than merely being at least as extreme as the point
      # estimate. This reflects the statistical precision of the published
      # contrast: a wide (e.g. non-significant) published CI is easy to fall
      # within, exactly as it should be.
      compatible_HR_contrast =
        HR_ratio_high_vs_lower >= observed_HR_ratio_high_vs_lower_ci_lower &
        HR_ratio_high_vs_lower <= observed_HR_ratio_high_vs_lower_ci_upper,
      directional_HR_contrast = HR_ratio_high_vs_lower < 1
    ) %>%
    group_by(beta_AU, trial_name, lower_group, high_group) %>%
    summarise(
      n_reps = n(),

      observed_HR_lower = first(HR_lower_obs),
      observed_HR_high = first(HR_high_obs),
      observed_HR_ratio_high_vs_lower = first(observed_HR_ratio_high_vs_lower),
      observed_log_HR_ratio_high_vs_lower = first(observed_log_HR_ratio_high_vs_lower),
      observed_HR_ratio_high_vs_lower_ci_lower = first(observed_HR_ratio_high_vs_lower_ci_lower),
      observed_HR_ratio_high_vs_lower_ci_upper = first(observed_HR_ratio_high_vs_lower_ci_upper),

      n_compatible_HR = sum(compatible_HR_contrast, na.rm = TRUE),
      compatibility_probability_HR = mean(compatible_HR_contrast, na.rm = TRUE),
      compatibility_probability_HR_low =
        prop_ci_exact(sum(compatible_HR_contrast, na.rm = TRUE), n())[1],
      compatibility_probability_HR_high =
        prop_ci_exact(sum(compatible_HR_contrast, na.rm = TRUE), n())[2],
      probability_directional_HR = mean(directional_HR_contrast, na.rm = TRUE),

      prob_apparent_benefit_higher_PRS = mean(beta_AxG < 0, na.rm = TRUE),
      beta_AxG_low = mc_interval(beta_AxG)[1],
      beta_AxG_high = mc_interval(beta_AxG)[2],
      beta_AxG = mc_mean(beta_AxG),
      approx_induced_AxG_low = mc_interval(approx_induced_AxG)[1],
      approx_induced_AxG_high = mc_interval(approx_induced_AxG)[2],
      approx_induced_AxG = mc_mean(approx_induced_AxG),

      cor_G_U_trial = mc_mean(cor_G_U_trial),
      cor_G_logLDL_trial = mc_mean(cor_G_logLDL_trial),
      slope_U_on_G_trial = mc_mean(slope_U_on_G_trial),
      event_rate_control = mc_mean(event_rate_control),
      event_rate_treat = mc_mean(event_rate_treat),

      HR_lower = mc_mean(HR_lower),
      HR_high = mc_mean(HR_high),
      HR_ratio_high_vs_lower_low = mc_interval(HR_ratio_high_vs_lower)[1],
      HR_ratio_high_vs_lower_high = mc_interval(HR_ratio_high_vs_lower)[2],
      HR_ratio_high_vs_lower = mc_mean(HR_ratio_high_vs_lower),

      ARR_lower = mc_mean(ARR_lower),
      ARR_high = mc_mean(ARR_high),
      ARR_difference_high_minus_lower = mc_mean(ARR_difference_high_minus_lower),
      contrast_note = first(contrast_note),
      .groups = "drop"
    )
}

extract_adjusted_contrasts <- function(sim) {
  sim$analysis$cox_by_grs_adjusted %>%
    group_by(model_set, adjustment) %>%
    group_modify(function(d, key) {
      fake_sim <- list(
        analysis = list(
          event_table = sim$analysis$event_table,
          cox_by_grs = d %>% select(GRS_group, HR, HR_low, HR_high, p)
        ),
        trial_spec = sim$trial_spec
      )
      extract_contrast(fake_sim)
    }) %>%
    ungroup()
}

run_beta_AU_grid_adjusted <- function(
  beta_AU_grid = c(0, 0.10, 0.20, 0.30, 0.40),
  specs = trial_specs,
  pars = pars_anchor,
  n_reps = n_reps,
  base_seed = 51000,
  published_contrasts = published_trial_contrasts
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

        adjusted_contrasts <- extract_adjusted_contrasts(sim_i)

        raw[[k]] <- sim_i$analysis$continuous_interactions_adjusted %>%
          left_join(
            adjusted_contrasts,
            by = c("model_set", "adjustment")
          ) %>%
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
    left_join(
      published_contrasts %>%
        select(
          trial_name,
          HR_lower_obs,
          HR_high_obs,
          observed_HR_ratio_high_vs_lower,
          observed_log_HR_ratio_high_vs_lower,
          observed_HR_ratio_high_vs_lower_ci_lower,
          observed_HR_ratio_high_vs_lower_ci_upper,
          contrast_note
        ),
      by = "trial_name"
    ) %>%
    mutate(
      # See run_beta_AU_grid() for why compatibility is defined against the
      # published trial's own 95% CI rather than its point estimate.
      compatible_HR_contrast =
        HR_ratio_high_vs_lower >= observed_HR_ratio_high_vs_lower_ci_lower &
        HR_ratio_high_vs_lower <= observed_HR_ratio_high_vs_lower_ci_upper,
      directional_HR_contrast = HR_ratio_high_vs_lower < 1
    ) %>%
    group_by(beta_AU, trial_name, model_set, adjustment, lower_group, high_group) %>%
    summarise(
      n_reps = n(),

      observed_HR_lower = first(HR_lower_obs),
      observed_HR_high = first(HR_high_obs),
      observed_HR_ratio_high_vs_lower = first(observed_HR_ratio_high_vs_lower),
      observed_log_HR_ratio_high_vs_lower = first(observed_log_HR_ratio_high_vs_lower),
      observed_HR_ratio_high_vs_lower_ci_lower = first(observed_HR_ratio_high_vs_lower_ci_lower),
      observed_HR_ratio_high_vs_lower_ci_upper = first(observed_HR_ratio_high_vs_lower_ci_upper),

      n_compatible_HR = sum(compatible_HR_contrast, na.rm = TRUE),
      compatibility_probability_HR = mean(compatible_HR_contrast, na.rm = TRUE),
      compatibility_probability_HR_low =
        prop_ci_exact(sum(compatible_HR_contrast, na.rm = TRUE), n())[1],
      compatibility_probability_HR_high =
        prop_ci_exact(sum(compatible_HR_contrast, na.rm = TRUE), n())[2],
      probability_directional_HR = mean(directional_HR_contrast, na.rm = TRUE),

      HR_lower = mc_mean(HR_lower),
      HR_high = mc_mean(HR_high),
      HR_ratio_high_vs_lower_low = mc_interval(HR_ratio_high_vs_lower)[1],
      HR_ratio_high_vs_lower_high = mc_interval(HR_ratio_high_vs_lower)[2],
      HR_ratio_high_vs_lower = mc_mean(HR_ratio_high_vs_lower),

      prob_apparent_benefit_higher_PRS = mean(beta_AxG < 0, na.rm = TRUE),
      beta_AxG_low = mc_interval(beta_AxG)[1],
      beta_AxG_high = mc_interval(beta_AxG)[2],
      beta_AxG = mc_mean(beta_AxG),
      median_model_p = median(p_AxG, na.rm = TRUE),
      contrast_note = first(contrast_note),
      .groups = "drop"
    )
}

# ----------------------------
# Manuscript-grade sensitivity figures
# ----------------------------

label_model_set <- function(x) {
  case_match(
    x,
    "unadjusted"                ~ "Main analysis",
    "ldl_adjusted"              ~ "Adjusted for LDL cholesterol",
    "burden_proxy_adjusted"     ~ "Adjusted for measured risk proxy",
    "ldl_burden_proxy_adjusted" ~ "Adjusted for LDL and risk proxy",
    .default = x
  )
}

make_beta_AU_compatibility_table_pretty <- function(beta_AU_sensitivity) {
  beta_AU_sensitivity %>%
    transmute(
      Trial = trial_name,
      `Degree of treatment response modification by other risk factors` = beta_AU,
      `Observed published HR contrast` = sprintf(
        "HR %.2f vs %.2f (ratio %.2f, 95%% CI %.2f-%.2f)",
        observed_HR_high,
        observed_HR_lower,
        observed_HR_ratio_high_vs_lower,
        observed_HR_ratio_high_vs_lower_ci_lower,
        observed_HR_ratio_high_vs_lower_ci_upper
      ),
      `Mean simulated HR contrast` = sprintf(
        "ratio %.2f (%.2f-%.2f)",
        HR_ratio_high_vs_lower,
        HR_ratio_high_vs_lower_low,
        HR_ratio_high_vs_lower_high
      ),
      `Probability within observed 95% CI` = sprintf(
        "%s (%s-%s)",
        pct_label(compatibility_probability_HR),
        pct_label(compatibility_probability_HR_low),
        pct_label(compatibility_probability_HR_high)
      ),
      `Probability high PRS appears to benefit more` = pct_label(probability_directional_HR),
      Note = contrast_note
    )
}

make_beta_AU_adjusted_compatibility_table_pretty <- function(beta_AU_sensitivity_adjusted) {
  beta_AU_sensitivity_adjusted %>%
    mutate(model_set_pretty = label_model_set(model_set)) %>%
    transmute(
      Trial = trial_name,
      Analysis = model_set_pretty,
      `Degree of treatment response modification by other risk factors` = beta_AU,
      `Observed published HR contrast` = sprintf(
        "HR %.2f vs %.2f (ratio %.2f, 95%% CI %.2f-%.2f)",
        observed_HR_high,
        observed_HR_lower,
        observed_HR_ratio_high_vs_lower,
        observed_HR_ratio_high_vs_lower_ci_lower,
        observed_HR_ratio_high_vs_lower_ci_upper
      ),
      `Mean simulated HR contrast` = sprintf(
        "ratio %.2f (%.2f-%.2f)",
        HR_ratio_high_vs_lower,
        HR_ratio_high_vs_lower_low,
        HR_ratio_high_vs_lower_high
      ),
      `Probability within observed 95% CI` = sprintf(
        "%s (%s-%s)",
        pct_label(compatibility_probability_HR),
        pct_label(compatibility_probability_HR_low),
        pct_label(compatibility_probability_HR_high)
      ),
      `Probability high PRS appears to benefit more` = pct_label(probability_directional_HR),
      Note = contrast_note
    )
}

plot_beta_AU_sensitivity <- function(beta_AU_sensitivity) {
  plot_dat <- beta_AU_sensitivity %>%
    mutate(
      trial_name = factor(trial_name, levels = c("FOURIER", "ODYSSEY", "WOSCOPS")),
      hidden_response = beta_AU,
      prob_pct = 100 * compatibility_probability_HR,
      prob_low_pct = 100 * compatibility_probability_HR_low,
      prob_high_pct = 100 * compatibility_probability_HR_high,
      label = pct_label(compatibility_probability_HR)
    )

  ggplot(plot_dat, aes(x = hidden_response, y = prob_pct, group = trial_name)) +
    geom_hline(yintercept = 5, linewidth = 0.4, colour = fps$null) +
    geom_hline(yintercept = 10, linewidth = 0.4, colour = fps$null, linetype = "22") +
    geom_ribbon(
      aes(ymin = prob_low_pct, ymax = prob_high_pct),
      fill = fps$navy,
      alpha = 0.12,
      colour = NA
    ) +
    geom_line(colour = fps$navy, linewidth = 1.2) +
    geom_point(colour = fps$navy, fill = fps$white, shape = 21, size = 4.5, stroke = 1) +
    geom_text(
      data = plot_dat %>% filter(hidden_response == max(hidden_response, na.rm = TRUE)),
      aes(label = label),
      nudge_y = 7,
      size = 5.5,
      colour = fps$navy
    ) +
    scale_x_continuous(breaks = sort(unique(plot_dat$hidden_response))) +
    scale_y_continuous(labels = function(x) paste0(x, "%")) +
    coord_cartesian(ylim = c(0, 100), clip = "off") +
    facet_wrap(~ trial_name, nrow = 1) +
    labs(
      # title = "How much treatment response modification by the non-genetic risk\n is needed to fall outside the published subgroup contrasts' 95% CI?",
      # subtitle = "Each point shows the percentage of simulated trials whose high-PRS treatment\ncontrast falls within the published trial's own 95% CI for that contrast.",
      x = "Degree of treatment response modification by other risk factors",
      y = "Probability of falling within\nthe observed 95% CI",
      # caption = "A simulation is counted as compatible when the high-PRS versus lower-PRS\ntreatment contrast falls within the published trial's 95% CI. Shaded\nbands show Monte Carlo uncertainty; the true PRS-by-treatment effect is zero."
    ) +
    theme_pub +
    theme_pub_axes
}

plot_beta_AU_sensitivity_adjusted <- function(beta_AU_sensitivity_adjusted) {
  plot_dat <- beta_AU_sensitivity_adjusted %>%
    mutate(
      trial_name = factor(trial_name, levels = c("FOURIER", "ODYSSEY", "WOSCOPS")),
      hidden_response = beta_AU,
      model_set = factor(
        label_model_set(model_set),
        levels = c(
          "Main analysis",
          "Adjusted for LDL cholesterol",
          "Adjusted for measured risk proxy",
          "Adjusted for LDL and risk proxy"
        )
      ),
      prob_pct = 100 * compatibility_probability_HR,
      prob_low_pct = 100 * compatibility_probability_HR_low,
      prob_high_pct = 100 * compatibility_probability_HR_high
    )

  ggplot(plot_dat, aes(x = hidden_response, y = prob_pct, colour = model_set, group = model_set)) +
    geom_hline(yintercept = 5, linewidth = 0.3, colour = fps$null) +
    geom_hline(yintercept = 10, linewidth = 0.3, colour = fps$null, linetype = "22") +
    geom_errorbar(
      aes(ymin = prob_low_pct, ymax = prob_high_pct),
      width = 0.012,
      linewidth = 0.45,
      alpha = 0.75
    ) +
    geom_line(linewidth = 0.82) +
    geom_point(size = 2.0) +
    scale_x_continuous(breaks = sort(unique(plot_dat$hidden_response))) +
    scale_y_continuous(labels = function(x) paste0(x, "%")) +
    coord_cartesian(ylim = c(0, 100)) +
    scale_colour_manual(values = c(
      "Main analysis" = fps$navy,
      "Adjusted for LDL cholesterol" = fps$steel,
      "Adjusted for measured risk proxy" = fps$sand,
      "Adjusted for LDL and risk proxy" = fps$slate
    )) +
    facet_wrap(~ trial_name, nrow = 1) +
    labs(
      # title = "Adjustment has limited impact on compatibility with the published subgroup contrasts",
      # subtitle = "Lines show how often simulated trials' high-PRS treatment contrast falls within\nthe published trial's 95% CI, after different covariate adjustments.",
      x = "Degree of treatment response modification by other risk factors",
      y = "Probability of falling within the observed 95% CI",
      colour = NULL,
      # caption = "A simulation is counted as compatible when the adjusted high-PRS versus lower-PRS\ntreatment contrast falls within the published trial's 95% CI. The true\nPRS-by-treatment effect is zero in every simulation."
    ) +
    theme_fps(base_size = 12) +
    theme(
      plot.title = element_text(face = "bold"),
      plot.subtitle = element_text(size = 10.5),
      strip.text = element_text(face = "bold"),
      panel.grid.minor = element_blank(),
      legend.position = "bottom",
      legend.box = "vertical",
      axis.title = element_text(face = "bold", size = 10),
      plot.caption = element_text(hjust = 0, size = 9),
      plot.background = element_rect(fill = fps$white, colour = NA)
    )
}

plot_beta_AU_HR_ratio_sensitivity <- function(beta_AU_sensitivity) {
  plot_dat <- beta_AU_sensitivity %>%
    mutate(
      trial_name = factor(trial_name, levels = c("FOURIER", "ODYSSEY", "WOSCOPS")),
      hidden_response = beta_AU,
      observed_label = sprintf("Published contrast: %.2f", observed_HR_ratio_high_vs_lower)
    )

  obs_dat <- plot_dat %>%
    distinct(
      trial_name,
      observed_HR_ratio_high_vs_lower,
      observed_HR_ratio_high_vs_lower_ci_lower,
      observed_HR_ratio_high_vs_lower_ci_upper,
      observed_label
    )

  ggplot(plot_dat, aes(x = hidden_response, y = HR_ratio_high_vs_lower)) +
    geom_hline(yintercept = 1, linewidth = 0.35, colour = fps$null) +
    geom_rect(
      data = obs_dat,
      aes(
        ymin = observed_HR_ratio_high_vs_lower_ci_lower,
        ymax = observed_HR_ratio_high_vs_lower_ci_upper
      ),
      xmin = -Inf,
      xmax = Inf,
      inherit.aes = FALSE,
      fill = fps$sand,
      alpha = 0.08
    ) +
    geom_hline(
      data = obs_dat,
      aes(yintercept = observed_HR_ratio_high_vs_lower),
      colour = fps$sand,
      linetype = "22",
      linewidth = 0.8
    ) +
    geom_ribbon(
      aes(ymin = HR_ratio_high_vs_lower_low, ymax = HR_ratio_high_vs_lower_high),
      fill = fps$navy,
      alpha = 0.12,
      colour = NA
    ) +
    geom_line(colour = fps$navy, linewidth = 0.9) +
    geom_point(colour = fps$navy, fill = fps$white, shape = 21, size = 2.4, stroke = 0.8) +
    scale_x_continuous(breaks = sort(unique(plot_dat$hidden_response))) +
    facet_wrap(~ trial_name, nrow = 1) +
    labs(
      # title = "Simulated high-PRS treatment contrasts become stronger as treatment\nresponse modification by non-genetic risk increases",
      # subtitle = "Lower values mean the high-PRS subgroup appears to benefit more than the lower-PRS\nsubgroup. Dashed lines show the published point estimate; the orange band shows its 95% CI.",
      x = "Degree of treatment response modification by other risk factors",
      y = "Ratio of high-PRS subgroup HR to low-PRS subgroup HR",
      # caption = "Values below 1 mean stronger apparent treatment benefit in the high-PRS subgroup. Blue\nbands show simulation-to-simulation uncertainty; the orange band is the published 95% CI."
    ) +
    theme_fps(base_size = 12) +
    theme(
      plot.title = element_text(face = "bold"),
      plot.subtitle = element_text(size = 10.5),
      strip.text = element_text(face = "bold"),
      panel.grid.minor = element_blank(),
      axis.title = element_text(face = "bold", size = 10),
      plot.caption = element_text(hjust = 0, size = 9),
      plot.background = element_rect(fill = fps$white, colour = NA)
    )
}

# ----------------------------
# Compatibility with published PRS subgroup contrasts
# ----------------------------
#
# These functions ask a finite-sample question that is deliberately different
# from comparing the mean simulated contrast with the published point estimate:
#
#   Across Monte Carlo replicates at the actual trial sample size, what
#   proportion of the simulated high-vs-lower PRS treatment HR contrasts fall
#   within the published trial-reanalysis contrast's own 95% CI?
#
# The HR contrast is evaluated on the high-vs-lower treatment HR ratio:
#   HR_ratio = HR_high_PRS / HR_lower_PRS.
# Values below 1 imply apparently greater relative treatment benefit in the
# high-PRS group. A simulation is called compatible if its HR ratio falls
# within [observed_HR_ratio_high_vs_lower_ci_lower, observed_HR_ratio_high_vs_lower_ci_upper],
# the published CI for that ratio (derived from the two subgroups' published
# HR CIs via the delta method; see ci_hr_ratio() below). This is a more
# honest standard than comparing against the point estimate alone: when the
# published high-vs-lower contrast is itself imprecise or not statistically
# significant (e.g. WOSCOPS), a wide published CI is easy for a null
# simulation to fall inside, and the compatibility probability should be
# correspondingly high.
#
# ARR compatibility still uses the older "at least as large as the published
# point estimate" rule, because the published papers do not report CIs for
# the group-level ARRs (or their difference), so there is no published ARR
# CI to fall within. ARR compatibility is optional because not all observed
# ARR contrasts are always reported in the same way across papers. Leave ARR
# values as NA when you do not want to evaluate ARR compatibility for that
# trial.

se_log_hr_from_ci <- function(hr_ci_lower, hr_ci_upper, conf_level = 0.95) {
  z <- qnorm(1 - (1 - conf_level) / 2)
  (log(hr_ci_upper) - log(hr_ci_lower)) / (2 * z)
}

ci_log_hr_ratio <- function(
  hr_lower,
  hr_lower_ci_lower,
  hr_lower_ci_upper,
  hr_high,
  hr_high_ci_lower,
  hr_high_ci_upper,
  conf_level = 0.95
) {
  z <- qnorm(1 - (1 - conf_level) / 2)

  log_hr_ratio <- log(hr_high) - log(hr_lower)

  se_lower <- se_log_hr_from_ci(
    hr_ci_lower = hr_lower_ci_lower,
    hr_ci_upper = hr_lower_ci_upper,
    conf_level = conf_level
  )

  se_high <- se_log_hr_from_ci(
    hr_ci_lower = hr_high_ci_lower,
    hr_ci_upper = hr_high_ci_upper,
    conf_level = conf_level
  )

  se_log_hr_ratio <- sqrt(se_lower^2 + se_high^2)

  tibble::tibble(
    observed_log_HR_ratio_high_vs_lower_ci_lower =
      log_hr_ratio - z * se_log_hr_ratio,
    observed_log_HR_ratio_high_vs_lower_ci_upper =
      log_hr_ratio + z * se_log_hr_ratio
  )
}

ci_hr_ratio <- function(
  hr_lower,
  hr_lower_ci_lower,
  hr_lower_ci_upper,
  hr_high,
  hr_high_ci_lower,
  hr_high_ci_upper,
  conf_level = 0.95
) {
  log_ci <- ci_log_hr_ratio(
    hr_lower = hr_lower,
    hr_lower_ci_lower = hr_lower_ci_lower,
    hr_lower_ci_upper = hr_lower_ci_upper,
    hr_high = hr_high,
    hr_high_ci_lower = hr_high_ci_lower,
    hr_high_ci_upper = hr_high_ci_upper,
    conf_level = conf_level
  )

  tibble::tibble(
    observed_HR_ratio_high_vs_lower_ci_lower =
      exp(log_ci$observed_log_HR_ratio_high_vs_lower_ci_lower),
    observed_HR_ratio_high_vs_lower_ci_upper =
      exp(log_ci$observed_log_HR_ratio_high_vs_lower_ci_upper)
  )
}

published_trial_contrasts <- tibble::tribble(
  ~trial_name, ~lower_group, ~high_group,
  ~HR_lower_obs, ~HR_lower_obs_ci_lower, ~HR_lower_obs_ci_upper,
  ~HR_high_obs,  ~HR_high_obs_ci_lower,  ~HR_high_obs_ci_upper,
  ~ARR_lower_obs, ~ARR_high_obs, ~contrast_note,
  "FOURIER", "Low",        "High",
      1.02, 0.80, 1.30,
      0.75, 0.60, 0.94,
      NA_real_, NA_real_,
      "FOURIER genetic-risk reanalysis: low vs high GRS subgroup HRs.",
  "ODYSSEY", "Lower",      "High",
      0.87, 0.78, 0.98,
      0.63, 0.46, 0.86,
      0.015, 0.060,
      "ODYSSEY OUTCOMES pharmacogenomic analysis: lower 90% vs top decile PRS.",
  "WOSCOPS", "All others", "High",
      0.76, 0.63, 0.92,
      0.56, 0.40, 0.78,
      0.028, 0.079,
      "WOSCOPS genetic-risk reanalysis: all others vs top quintile PRS; ARR values use long-term event risks."
) %>%
  dplyr::mutate(
    observed_HR_ratio_high_vs_lower = HR_high_obs / HR_lower_obs,
    observed_log_HR_ratio_high_vs_lower = log(observed_HR_ratio_high_vs_lower),
    observed_ARR_difference_high_minus_lower = ARR_high_obs - ARR_lower_obs
  ) %>%
  dplyr::bind_cols(
    ci_log_hr_ratio(
      hr_lower = .$HR_lower_obs,
      hr_lower_ci_lower = .$HR_lower_obs_ci_lower,
      hr_lower_ci_upper = .$HR_lower_obs_ci_upper,
      hr_high = .$HR_high_obs,
      hr_high_ci_lower = .$HR_high_obs_ci_lower,
      hr_high_ci_upper = .$HR_high_obs_ci_upper
    ),
    ci_hr_ratio(
      hr_lower = .$HR_lower_obs,
      hr_lower_ci_lower = .$HR_lower_obs_ci_lower,
      hr_lower_ci_upper = .$HR_lower_obs_ci_upper,
      hr_high = .$HR_high_obs,
      hr_high_ci_lower = .$HR_high_obs_ci_lower,
      hr_high_ci_upper = .$HR_high_obs_ci_upper
    )
  )

prop_ci_exact <- function(k, n, conf.level = 0.95) {
  if (is.na(k) || is.na(n) || n == 0) return(c(NA_real_, NA_real_))

  alpha <- 1 - conf.level

  lower <- if (k == 0) {
    0
  } else {
    qbeta(alpha / 2, k, n - k + 1)
  }

  upper <- if (k == n) {
    1
  } else {
    qbeta(1 - alpha / 2, k + 1, n - k)
  }

  c(lower, upper)
}

pct_label <- function(x, digits = 0) {
  ifelse(is.na(x), "NA", paste0(sprintf(paste0("%.", digits, "f"), 100 * x), "%"))
}

make_compatibility_replicates_one <- function(
  rep_result,
  published_contrasts = published_trial_contrasts
) {
  trial_name_i <- rep_result$trial_spec$name
  obs_i <- published_contrasts %>% filter(trial_name == trial_name_i)

  if (nrow(obs_i) != 1L) {
    stop("Expected exactly one published contrast row for ", trial_name_i)
  }

  # Do not join the observed lower/high group labels directly, because
  # simulated_contrasts already contains lower_group and high_group. If those
  # columns are present on both sides of the join, dplyr creates lower_group.x /
  # lower_group.y, which later breaks group_by(lower_group, high_group).
  obs_i_for_join <- obs_i %>%
    select(-lower_group, -high_group)

  rep_ids <- sort(unique(rep_result$event_tables$rep))

  simulated_contrasts <- map_dfr(rep_ids, function(r) {
    event_df <- rep_result$event_tables %>% filter(rep == r)
    cox_df <- rep_result$cox_tables %>% filter(rep == r)

    fake_sim <- list(
      analysis = list(event_table = event_df, cox_by_grs = cox_df),
      trial_spec = rep_result$trial_spec
    )

    extract_contrast(fake_sim) %>%
      mutate(rep = r, .before = 1)
  })

  analytic_contrasts <- rep_result$analytic_contrasts %>%
    select(
      rep,
      analytic_treatment_HR_ratio_high_vs_lower,
      analytic_ARR_difference_high_minus_lower
    )

  simulated_contrasts %>%
    left_join(analytic_contrasts, by = "rep") %>%
    mutate(trial_name = trial_name_i, .before = 1) %>%
    left_join(obs_i_for_join, by = "trial_name") %>%
    mutate(
      lower_group = as.character(lower_group),
      high_group = as.character(high_group),
      simulated_log_HR_ratio_high_vs_lower = log(HR_ratio_high_vs_lower),
      analytic_log_HR_ratio_high_vs_lower = log(analytic_treatment_HR_ratio_high_vs_lower),

      # A replicate is compatible when its simulated high-vs-lower PRS
      # treatment HR ratio falls inside the published trial's own 95% CI for
      # that contrast, rather than merely being at least as extreme as the
      # point estimate. This reflects the statistical precision of the
      # published contrast: a wide (e.g. non-significant) published CI is
      # easy to fall within, exactly as it should be.
      compatible_HR_contrast =
        HR_ratio_high_vs_lower >= observed_HR_ratio_high_vs_lower_ci_lower &
        HR_ratio_high_vs_lower <= observed_HR_ratio_high_vs_lower_ci_upper,
      compatible_analytic_HR_contrast =
        analytic_treatment_HR_ratio_high_vs_lower >= observed_HR_ratio_high_vs_lower_ci_lower &
        analytic_treatment_HR_ratio_high_vs_lower <= observed_HR_ratio_high_vs_lower_ci_upper,

      compatible_ARR_contrast = ifelse(
        is.na(observed_ARR_difference_high_minus_lower),
        NA,
        ARR_difference_high_minus_lower >= observed_ARR_difference_high_minus_lower
      ),
      compatible_analytic_ARR_contrast = ifelse(
        is.na(observed_ARR_difference_high_minus_lower),
        NA,
        analytic_ARR_difference_high_minus_lower >= observed_ARR_difference_high_minus_lower
      ),

      directional_HR_contrast = HR_ratio_high_vs_lower < 1,
      directional_ARR_contrast = ARR_difference_high_minus_lower > 0
    )
}

make_compatibility_replicates <- function(
  replicate_results,
  published_contrasts = published_trial_contrasts
) {
  map_dfr(
    replicate_results,
    make_compatibility_replicates_one,
    published_contrasts = published_contrasts
  )
}

summarise_compatibility_probabilities <- function(compatibility_replicates) {
  # Backwards-compatible repair for compatibility_replicates objects created
  # before make_compatibility_replicates_one() stopped joining duplicate
  # lower_group/high_group columns.
  dat <- compatibility_replicates

  if (!"lower_group" %in% names(dat) && "lower_group.x" %in% names(dat)) {
    dat <- dat %>% rename(lower_group = lower_group.x)
  }
  if (!"high_group" %in% names(dat) && "high_group.x" %in% names(dat)) {
    dat <- dat %>% rename(high_group = high_group.x)
  }

  dat <- dat %>%
    select(-any_of(c("lower_group.y", "high_group.y")))

  dat %>%
    group_by(trial_name, lower_group, high_group) %>%
    summarise(
      n_reps = n(),

      observed_HR_lower = first(HR_lower_obs),
      observed_HR_high = first(HR_high_obs),
      observed_HR_ratio_high_vs_lower = first(observed_HR_ratio_high_vs_lower),
      observed_log_HR_ratio_high_vs_lower = first(observed_log_HR_ratio_high_vs_lower),
      observed_HR_ratio_high_vs_lower_ci_lower = first(observed_HR_ratio_high_vs_lower_ci_lower),
      observed_HR_ratio_high_vs_lower_ci_upper = first(observed_HR_ratio_high_vs_lower_ci_upper),

      simulated_HR_lower = mc_mean(HR_lower),
      simulated_HR_high = mc_mean(HR_high),
      simulated_HR_ratio_high_vs_lower = mc_mean(HR_ratio_high_vs_lower),
      simulated_HR_ratio_low = mc_interval(HR_ratio_high_vs_lower)[1],
      simulated_HR_ratio_high = mc_interval(HR_ratio_high_vs_lower)[2],

      analytic_HR_ratio_high_vs_lower = mc_mean(analytic_treatment_HR_ratio_high_vs_lower),
      analytic_HR_ratio_low = mc_interval(analytic_treatment_HR_ratio_high_vs_lower)[1],
      analytic_HR_ratio_high = mc_interval(analytic_treatment_HR_ratio_high_vs_lower)[2],

      n_compatible_HR = sum(compatible_HR_contrast, na.rm = TRUE),
      compatibility_probability_HR = mean(compatible_HR_contrast, na.rm = TRUE),
      compatibility_probability_HR_low = prop_ci_exact(n_compatible_HR, n_reps)[1],
      compatibility_probability_HR_high = prop_ci_exact(n_compatible_HR, n_reps)[2],

      n_compatible_analytic_HR = sum(compatible_analytic_HR_contrast, na.rm = TRUE),
      compatibility_probability_analytic_HR = mean(compatible_analytic_HR_contrast, na.rm = TRUE),
      compatibility_probability_analytic_HR_low = prop_ci_exact(n_compatible_analytic_HR, n_reps)[1],
      compatibility_probability_analytic_HR_high = prop_ci_exact(n_compatible_analytic_HR, n_reps)[2],

      probability_directional_HR = mean(directional_HR_contrast, na.rm = TRUE),

      observed_ARR_lower = first(ARR_lower_obs),
      observed_ARR_high = first(ARR_high_obs),
      observed_ARR_difference_high_minus_lower =
        first(observed_ARR_difference_high_minus_lower),

      simulated_ARR_lower = mc_mean(ARR_lower),
      simulated_ARR_high = mc_mean(ARR_high),
      simulated_ARR_difference_high_minus_lower =
        mc_mean(ARR_difference_high_minus_lower),
      simulated_ARR_difference_low = mc_interval(ARR_difference_high_minus_lower)[1],
      simulated_ARR_difference_high = mc_interval(ARR_difference_high_minus_lower)[2],

      analytic_ARR_difference_high_minus_lower =
        mc_mean(analytic_ARR_difference_high_minus_lower),
      analytic_ARR_difference_low = mc_interval(analytic_ARR_difference_high_minus_lower)[1],
      analytic_ARR_difference_high = mc_interval(analytic_ARR_difference_high_minus_lower)[2],

      n_compatible_ARR = sum(compatible_ARR_contrast, na.rm = TRUE),
      n_evaluable_ARR = sum(!is.na(compatible_ARR_contrast)),
      compatibility_probability_ARR = mean(compatible_ARR_contrast, na.rm = TRUE),
      compatibility_probability_ARR_low = prop_ci_exact(n_compatible_ARR, n_evaluable_ARR)[1],
      compatibility_probability_ARR_high = prop_ci_exact(n_compatible_ARR, n_evaluable_ARR)[2],

      probability_directional_ARR = mean(directional_ARR_contrast, na.rm = TRUE),
      contrast_note = first(contrast_note),
      .groups = "drop"
    )
}

make_compatibility_probability_table_pretty <- function(compatibility_probabilities) {
  compatibility_probabilities %>%
    transmute(
      Trial = trial_name,
      `PRS contrast` = paste0(high_group, " vs ", lower_group),
      `Observed HR contrast` = sprintf(
        "%.2f vs %.2f (ratio %.2f, 95%% CI %.2f-%.2f)",
        observed_HR_high,
        observed_HR_lower,
        observed_HR_ratio_high_vs_lower,
        observed_HR_ratio_high_vs_lower_ci_lower,
        observed_HR_ratio_high_vs_lower_ci_upper
      ),
      `Simulated HR ratio, mean (95% range)` = sprintf(
        "%.2f (%.2f-%.2f)",
        simulated_HR_ratio_high_vs_lower,
        simulated_HR_ratio_low,
        simulated_HR_ratio_high
      ),
      `Probability within observed 95% CI (HR contrast)` = sprintf(
        "%s (%s-%s)",
        pct_label(compatibility_probability_HR),
        pct_label(compatibility_probability_HR_low),
        pct_label(compatibility_probability_HR_high)
      ),
      `Observed ARR contrast` = ifelse(
        is.na(observed_ARR_difference_high_minus_lower),
        "Not evaluated",
        sprintf(
          "%.1f%% vs %.1f%% (difference %.1f%%)",
          100 * observed_ARR_high,
          100 * observed_ARR_lower,
          100 * observed_ARR_difference_high_minus_lower
        )
      ),
      `Simulated ARR difference, mean (95% range)` = sprintf(
        "%.1f%% (%.1f%%-%.1f%%)",
        100 * simulated_ARR_difference_high_minus_lower,
        100 * simulated_ARR_difference_low,
        100 * simulated_ARR_difference_high
      ),
      `Compatibility probability for ARR contrast` = ifelse(
        is.na(compatibility_probability_ARR),
        "Not evaluated",
        sprintf(
          "%s (%s-%s)",
          pct_label(compatibility_probability_ARR),
          pct_label(compatibility_probability_ARR_low),
          pct_label(compatibility_probability_ARR_high)
        )
      ),
      Note = contrast_note
    )
}

plot_hr_compatibility_probabilities <- function(compatibility_probabilities) {
  plot_dat <- compatibility_probabilities %>%
    mutate(
      trial_name = factor(trial_name, levels = c("FOURIER", "ODYSSEY", "WOSCOPS")),
      prob_pct = 100 * compatibility_probability_HR,
      prob_low_pct = 100 * compatibility_probability_HR_low,
      prob_high_pct = 100 * compatibility_probability_HR_high,
      label = pct_label(compatibility_probability_HR)
    )

  ggplot(plot_dat, aes(x = trial_name, y = prob_pct)) +
    geom_col(width = 0.62, fill = fps$navy, colour = fps$charcoal, linewidth = 0.25) +
    geom_errorbar(
      aes(ymin = prob_low_pct, ymax = prob_high_pct),
      width = 0.16,
      colour = fps$navy,
      linewidth = 0.65
    ) +
    geom_text(aes(x = trial_name, y = prob_high_pct, label = label), inherit.aes = FALSE, vjust = -0.55, size = 3.7) +
    coord_cartesian(ylim = c(0, max(100, max(plot_dat$prob_high_pct, na.rm = TRUE) + 8))) +
    theme_fps(base_size = 12) +
    theme(
      plot.title = element_text(face = "bold"),
      plot.subtitle = element_text(size = 10.5),
      panel.grid.major.x = element_blank(),
      panel.grid.minor = element_blank(),
      axis.title.x = element_blank(),
      axis.title.y = element_text(face = "bold"),
      plot.caption = element_text(hjust = 0, size = 9),
      plot.background = element_rect(fill = fps$white, colour = NA)
    ) +
    labs(
      # title = "Proportion of simulations whose HR contrast falls within the published trial's 95% CI",
      # subtitle = "Bars show the percentage of Monte Carlo replicates with a high-vs-lower PRS treatment\nHR contrast inside the published trial's own 95% CI for that contrast.",
      y = "Probability within observed 95% CI (%)",
      # caption = "A replicate is counted when the simulated high-PRS/lower-PRS treatment HR ratio falls\nwithin the published trial's 95% CI for that ratio. Error bars are exact binomial 95% intervals."
    )
}

plot_hr_contrast_distribution <- function(compatibility_replicates) {
  obs_dat <- compatibility_replicates %>%
    distinct(
      trial_name,
      observed_HR_ratio_high_vs_lower,
      observed_HR_ratio_high_vs_lower_ci_lower,
      observed_HR_ratio_high_vs_lower_ci_upper
    ) %>%
    mutate(trial_name = factor(trial_name, levels = c("FOURIER", "ODYSSEY", "WOSCOPS")))

  plot_dat <- compatibility_replicates %>%
    mutate(trial_name = factor(trial_name, levels = c("FOURIER", "ODYSSEY", "WOSCOPS")))

  ggplot(plot_dat, aes(x = HR_ratio_high_vs_lower)) +
    geom_rect(
      data = obs_dat,
      aes(
        xmin = observed_HR_ratio_high_vs_lower_ci_lower,
        xmax = observed_HR_ratio_high_vs_lower_ci_upper,
        ymin = -Inf,
        ymax = Inf
      ),
      inherit.aes = FALSE,
      fill = fps$sand,
      alpha = 0.08
    ) +
    geom_histogram(
      bins = 34,
      fill = fps$navy,
      colour = fps$white,
      alpha = 0.88
    ) +
    geom_vline(
      data = obs_dat,
      aes(xintercept = observed_HR_ratio_high_vs_lower),
      inherit.aes = FALSE,
      colour = fps$sand,
      linetype = "22",
      linewidth = 0.9
    ) +
    facet_wrap(~ trial_name, scales = "free_y", nrow = 1) +
    theme_fps(base_size = 12) +
    theme(
      plot.title = element_text(face = "bold"),
      plot.subtitle = element_text(size = 10.5),
      strip.text = element_text(face = "bold"),
      panel.grid.minor = element_blank(),
      axis.title = element_text(face = "bold"),
      plot.caption = element_text(hjust = 0, size = 9),
      plot.background = element_rect(fill = fps$white, colour = NA)
    ) +
    labs(
      # title = "Distribution of high-vs-lower PRS HR ratios across simulations",
      # subtitle = "Dashed line shows the published point estimate; shaded band shows the published 95% CI.",
      x = "High-PRS treatment HR divided by lower-PRS treatment HR",
      y = "Number of simulations",
      # caption = "The shaded band shows the published trial's 95% CI for the high-vs-lower PRS HR ratio."
    )
}

plot_arr_compatibility_probabilities <- function(compatibility_probabilities) {
  plot_dat <- compatibility_probabilities %>%
    filter(!is.na(compatibility_probability_ARR)) %>%
    mutate(
      trial_name = factor(trial_name, levels = c("FOURIER", "ODYSSEY", "WOSCOPS")),
      prob_pct = 100 * compatibility_probability_ARR,
      prob_low_pct = 100 * compatibility_probability_ARR_low,
      prob_high_pct = 100 * compatibility_probability_ARR_high,
      label = pct_label(compatibility_probability_ARR)
    )

  if (nrow(plot_dat) == 0L) {
    return(NULL)
  }

  ggplot(plot_dat, aes(x = trial_name, y = prob_pct)) +
    geom_col(width = 0.62, fill = fps$steel, colour = fps$charcoal, linewidth = 0.25) +
    geom_errorbar(
      aes(ymin = prob_low_pct, ymax = prob_high_pct),
      width = 0.16,
      colour = fps$steel,
      linewidth = 0.65
    ) +
    geom_text(aes(x = trial_name, y = prob_high_pct, label = label), inherit.aes = FALSE, vjust = -0.55, size = 3.7) +
    coord_cartesian(ylim = c(0, max(100, max(plot_dat$prob_high_pct, na.rm = TRUE) + 8))) +
    theme_fps(base_size = 12) +
    theme(
      plot.title = element_text(face = "bold"),
      plot.subtitle = element_text(size = 10.5),
      panel.grid.major.x = element_blank(),
      panel.grid.minor = element_blank(),
      axis.title.x = element_blank(),
      axis.title.y = element_text(face = "bold"),
      plot.caption = element_text(hjust = 0, size = 9),
      plot.background = element_rect(fill = fps$white, colour = NA)
    ) +
    labs(
      # title = "Proportion of simulations matching or exceeding the absolute HR differences in the published trials",
      # subtitle = "Bars show the percentage of Monte Carlo replicates with a high-minus-lower PRS\nabsolute risk reduction at least as large as the published contrast.",
      y = "Compatibility probability (%)",
      # caption = "ARR compatibility is evaluated only for trials with observed ARR contrasts in\npublished_trial_contrasts. Error bars are exact binomial 95% intervals."
    )
}

plot_arr_contrast_distribution <- function(compatibility_replicates) {
  plot_dat <- compatibility_replicates %>%
    filter(!is.na(observed_ARR_difference_high_minus_lower)) %>%
    mutate(trial_name = factor(trial_name, levels = c("FOURIER", "ODYSSEY", "WOSCOPS")))

  if (nrow(plot_dat) == 0L) {
    return(NULL)
  }

  obs_dat <- plot_dat %>%
    distinct(trial_name, observed_ARR_difference_high_minus_lower)

  ggplot(plot_dat, aes(x = 100 * ARR_difference_high_minus_lower)) +
    geom_rect(
      data = obs_dat,
      aes(
        xmin = 100 * observed_ARR_difference_high_minus_lower,
        xmax = Inf,
        ymin = -Inf,
        ymax = Inf
      ),
      inherit.aes = FALSE,
      fill = fps$sand,
      alpha = 0.08
    ) +
    geom_histogram(
      bins = 34,
      fill = fps$navy,
      colour = fps$white,
      alpha = 0.88
    ) +
    geom_vline(
      data = obs_dat,
      aes(xintercept = 100 * observed_ARR_difference_high_minus_lower),
      inherit.aes = FALSE,
      colour = fps$sand,
      linetype = "22",
      linewidth = 0.9
    ) +
    facet_wrap(~ trial_name, scales = "free_y", nrow = 1) +
    theme_fps(base_size = 12) +
    theme(
      plot.title = element_text(face = "bold"),
      plot.subtitle = element_text(size = 10.5),
      strip.text = element_text(face = "bold"),
      panel.grid.minor = element_blank(),
      axis.title = element_text(face = "bold"),
      plot.caption = element_text(hjust = 0, size = 9),
      plot.background = element_rect(fill = fps$white, colour = NA)
    ) +
    labs(
      # title = "Distribution of high-vs-lower PRS absolute HR differences across simulations",
      # subtitle = "Dashed line shows the published high-minus-lower PRS ARR contrast; shaded area\nshows simulated contrasts at least as large.",
      x = "High-minus-lower PRS absolute risk reduction (percentage points)",
      y = "Number of simulations",
      # caption = "Positive values imply larger absolute benefit in the high-PRS subgroup."
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
anchor_prs_decile_tables <- imap_dfr(replicate_results, ~ mutate(summarise_prs_decile_replicates(.x), trial_name = .y, .before = 1))
anchor_analytic_tables <- map_dfr(replicate_results, summarise_anchor_analytic_tables)
anchor_analytic_contrasts <- map_dfr(replicate_results, summarise_anchor_analytic_contrasts)
anchor_cox_tables_adjusted <- imap_dfr(replicate_results, ~ mutate(.x$cox_tables_adjusted, trial_name = .y, .before = 1))
anchor_continuous_interactions_adjusted <- imap_dfr(replicate_results, ~ mutate(.x$continuous_interactions_adjusted, trial_name = .y, .before = 1))
anchor_categorical_interactions_adjusted <- imap_dfr(replicate_results, ~ mutate(.x$categorical_interactions_adjusted, trial_name = .y, .before = 1))

compatibility_replicates <- make_compatibility_replicates(
  replicate_results,
  published_contrasts = published_trial_contrasts
)
compatibility_probabilities <- summarise_compatibility_probabilities(compatibility_replicates)
compatibility_probability_table <- make_compatibility_probability_table_pretty(compatibility_probabilities)

print(anchor_summary)
print(anchor_event_tables)
print(anchor_diagnostics_by_grs)
print(anchor_prs_decile_tables)
print(anchor_analytic_tables)
print(anchor_analytic_contrasts)
print(anchor_cox_tables_adjusted)
print(anchor_continuous_interactions_adjusted)
print(anchor_categorical_interactions_adjusted)
print(compatibility_probabilities)
print(compatibility_probability_table)

if (write_plots) {
  if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

  risk_reduction_plot <- plot_risk_reduction_by_grs(replicate_results, published_trial_contrasts)

  ggsave(
    filename = file.path(output_dir, "risk_reduction_by_grs.png"),
    plot = risk_reduction_plot,
    width = 12,
    height = 7.5,
    dpi = 300,
    device = ragg::agg_png,
    bg = fps$white
  )

  prs_decile_plot <- plot_nongenetic_risk_by_prs_decile(replicate_results)

  ggsave(
    filename = file.path(output_dir, "nongenetic_risk_by_prs_decile.png"),
    plot = prs_decile_plot,
    width = 10.8,
    height = 4.6,
    dpi = 320,
    device = ragg::agg_png
  )

  compatibility_plot_HR <- plot_hr_compatibility_probabilities(compatibility_probabilities)
  compatibility_distribution_HR <- plot_hr_contrast_distribution(compatibility_replicates)
  compatibility_plot_ARR <- plot_arr_compatibility_probabilities(compatibility_probabilities)
  compatibility_distribution_ARR <- plot_arr_contrast_distribution(compatibility_replicates)

  ggsave(
    filename = file.path(output_dir, "compatibility_probability_HR_contrast.png"),
    plot = compatibility_plot_HR,
    width = 10.8,
    height = 4.8,
    dpi = 320,
    device = ragg::agg_png
  )

  ggsave(
    filename = file.path(output_dir, "compatibility_HR_ratio_distribution.png"),
    plot = compatibility_distribution_HR,
    width = 10.8,
    height = 4.8,
    dpi = 320,
    device = ragg::agg_png
  )

  if (!is.null(compatibility_plot_ARR)) {
    ggsave(
      filename = file.path(output_dir, "compatibility_probability_ARR_contrast.png"),
      plot = compatibility_plot_ARR,
      width = 10.8,
      height = 4.8,
      dpi = 320,
      device = ragg::agg_png
    )
  }

  if (!is.null(compatibility_distribution_ARR)) {
    ggsave(
      filename = file.path(output_dir, "compatibility_ARR_distribution.png"),
      plot = compatibility_distribution_ARR,
      width = 9.0,
      height = 4.8,
      dpi = 320,
      device = ragg::agg_png
    )
  }

  write.csv(anchor_summary, file.path(output_dir, "anchor_summary.csv"), row.names = FALSE)
  write.csv(anchor_event_tables, file.path(output_dir, "anchor_event_tables.csv"), row.names = FALSE)
  write.csv(anchor_diagnostics_by_grs, file.path(output_dir, "anchor_diagnostics_by_grs.csv"), row.names = FALSE)
  write.csv(anchor_prs_decile_tables, file.path(output_dir, "anchor_prs_decile_tables.csv"), row.names = FALSE)
  write.csv(anchor_analytic_tables, file.path(output_dir, "anchor_analytic_group_risks.csv"), row.names = FALSE)
  write.csv(anchor_analytic_contrasts, file.path(output_dir, "anchor_analytic_high_low_contrasts.csv"), row.names = FALSE)
  write.csv(anchor_cox_tables_adjusted, file.path(output_dir, "anchor_cox_tables_adjusted_raw.csv"), row.names = FALSE)
  write.csv(anchor_continuous_interactions_adjusted, file.path(output_dir, "anchor_continuous_interactions_adjusted_raw.csv"), row.names = FALSE)
  write.csv(anchor_categorical_interactions_adjusted, file.path(output_dir, "anchor_categorical_interactions_adjusted_raw.csv"), row.names = FALSE)
  write.csv(published_trial_contrasts, file.path(output_dir, "published_trial_contrasts_used_for_compatibility.csv"), row.names = FALSE)
  write.csv(compatibility_replicates, file.path(output_dir, "compatibility_replicates.csv"), row.names = FALSE)
  write.csv(compatibility_probabilities, file.path(output_dir, "compatibility_probabilities.csv"), row.names = FALSE)
  write.csv(compatibility_probability_table, file.path(output_dir, "compatibility_probability_table_pretty.csv"), row.names = FALSE)
}

if (run_beta_AU_sensitivity) {
  beta_AU_sensitivity <- run_beta_AU_grid(n_reps = n_reps)
  beta_AU_sensitivity_adjusted <- run_beta_AU_grid_adjusted(n_reps = n_reps)

  beta_AU_compatibility_table <- make_beta_AU_compatibility_table_pretty(beta_AU_sensitivity)
  beta_AU_adjusted_compatibility_table <- make_beta_AU_adjusted_compatibility_table_pretty(beta_AU_sensitivity_adjusted)

  print(beta_AU_sensitivity)
  print(beta_AU_sensitivity_adjusted)
  print(beta_AU_compatibility_table)
  print(beta_AU_adjusted_compatibility_table)

  if (write_plots) {
    if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

    sensitivity_plot <- plot_beta_AU_sensitivity(beta_AU_sensitivity)
    sensitivity_plot_adjusted <- plot_beta_AU_sensitivity_adjusted(beta_AU_sensitivity_adjusted)
    sensitivity_plot_hr_ratio <- plot_beta_AU_HR_ratio_sensitivity(beta_AU_sensitivity)

    ggsave(
      filename = file.path(output_dir, "beta_AU_sensitivity_compatibility_probability.png"),
      plot = sensitivity_plot,
      width = 12,
      height = 6,
      dpi = 300,
      device = ragg::agg_png,
      bg = fps$white
    )

    ggsave(
      filename = file.path(output_dir, "beta_AU_sensitivity_adjusted_compatibility_probability.png"),
      plot = sensitivity_plot_adjusted,
      width = 10.8,
      height = 5.1,
      dpi = 320,
      device = ragg::agg_png
    )

    ggsave(
      filename = file.path(output_dir, "beta_AU_sensitivity_HR_contrast.png"),
      plot = sensitivity_plot_hr_ratio,
      width = 10.8,
      height = 5.1,
      dpi = 320,
      device = ragg::agg_png
    )

    # Backwards-compatible filenames for downstream notebooks.
    ggsave(
      filename = file.path(output_dir, "beta_AU_sensitivity_overview.png"),
      plot = sensitivity_plot,
      width = 12,
      height = 6,
      dpi = 300,
      device = ragg::agg_png,
      bg = fps$white
    )

    ggsave(
      filename = file.path(output_dir, "beta_AU_sensitivity_adjusted_overview.png"),
      plot = sensitivity_plot_adjusted,
      width = 10.8,
      height = 5.1,
      dpi = 320,
      device = ragg::agg_png
    )

    write.csv(beta_AU_sensitivity, file.path(output_dir, "beta_AU_sensitivity.csv"), row.names = FALSE)
    write.csv(beta_AU_sensitivity_adjusted, file.path(output_dir, "beta_AU_sensitivity_adjusted_models.csv"), row.names = FALSE)
    write.csv(beta_AU_compatibility_table, file.path(output_dir, "beta_AU_compatibility_probability_table_pretty.csv"), row.names = FALSE)
    write.csv(beta_AU_adjusted_compatibility_table, file.path(output_dir, "beta_AU_adjusted_compatibility_probability_table_pretty.csv"), row.names = FALSE)
  }
}
