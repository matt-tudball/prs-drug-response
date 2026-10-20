# Polygenic prediction of treatment response for precision medicine is vulnerable to collider bias

Code to reproduce the simulations, UK Biobank analyses and figures for:

> Matthew J Tudball, Elizabeth T Cirulli, Michael V Holmes, George Davey Smith, J Brent Richards. **Polygenic prediction of treatment response for precision medicine is vulnerable to collider bias.** *medRxiv* (2026). [TODO: DOI]

**Preprint:** [PLACEHOLDER: medRxiv link to be added once the preprint is posted]

The paper argues that apparent treatment-by-polygenic-risk-score (PRS) interactions in drug-response studies can arise without any true genetic modification of drug effect. Two designs are studied:

1. **Randomised trials with risk-based eligibility.** Trials such as FOURIER, ODYSSEY OUTCOMES and WOSCOPS enrol people selected on elevated LDL-C and presence or absence of pre-existing CAD. Both criteria depend on the PRS (G) and on non-genetic risk (U). Selection on their common effect induces a G–U correlation among participants. If U modifies the treatment effect, a spurious treatment × G interaction appears.
2. **Observational biobank/EHR studies of LDL-C response.** Adjusting for, stratifying on, or computing change scores against baseline LDL-C conditions on a collider of G and U. This can produce a spurious statin × LDL-PRS interaction.

In every simulation the true data-generating model contains **no treatment × G term**. Any estimated interaction is therefore artefactual.

---

## Repository layout

```
.
├── simulation/
│   ├── clinical_trial_one_shot_sim_refactor.r   # Trial simulations (FOURIER/ODYSSEY/WOSCOPS-like)
│   ├── log_normal_prs_statin_interaction.r      # Observational statin–LDL-C simulations
│   ├── simulation-outputs/                      # Figures (PNG, tracked) and tables (CSV, untracked)
│   └── old-code/                                # Superseded earlier versions; not used for the paper
├── ukb/                                         # UK Biobank analysis notebooks (run on the DNAnexus RAP)
│   ├── format_ukb_covariates.ipynb              # Python: extract covariates, PRS, biomarkers
│   ├── format_inpatient_data.ipynb              # PySpark: acute coronary events, PCI/CABG
│   ├── format_gp_data.ipynb                     # R: primary-care landmark cohort (main design)
│   ├── format_gp_data_from_assessment.ipynb     # R: cohort anchored at UKB assessment (sensitivity)
│   ├── run_analysis.ipynb                       # R: full-cohort interaction models (LDL-C and CAD)
│   ├── run_analysis_from_assessment.ipynb       # R: same, assessment-anchored cohort
│   ├── run_stratified_analysis.ipynb            # R: interaction within baseline LDL-C tertiles
│   ├── run_initiator_only_analysis.ipynb        # R: change-score analyses in statin initiators
│   └── *.png                                    # UKB-only figures
└── figures/
    ├── fps_style.r                              # Shared ggplot2 theme and colour palette
    ├── collider_example.r                       # Schematic collider-bias figure
    ├── ukb_hrn_meta_analysis.r                  # UKB + HRN fixed-effect meta-analysis and figures
    ├── ukb_hrn_interaction_estimates.csv        # Cohort-level stratified interaction estimates (input)
    └── ukb_hrn_change_score_estimates.csv       # Cohort-level change-score estimates (input)
```

---

## Mapping from paper to code

[TODO: replace the "Paper item" column with the final figure/table numbers.]

| Paper item | Script | Output file |
|---|---|---|
| Collider-bias schematic | `figures/collider_example.r` | `figures/collider_example.png` |
| Trial simulations: treatment HR by PRS group vs published estimates | `simulation/clinical_trial_one_shot_sim_refactor.r` | `simulation-outputs/risk_reduction_by_grs.png` |
| Trial simulations: non-genetic risk by PRS decile, population vs trial | same | `simulation-outputs/nongenetic_risk_by_prs_decile.png` |
| Trial simulations: compatibility with published HR contrasts | same | `simulation-outputs/compatibility_probability_HR_contrast.png`, `compatibility_HR_ratio_distribution.png` |
| Trial simulations: compatibility with published ARR contrasts | same | `simulation-outputs/compatibility_probability_ARR_contrast.png`, `compatibility_ARR_distribution.png` |
| Trial simulations: sensitivity to β<sub>AU</sub> | same | `simulation-outputs/beta_AU_sensitivity_*.png` |
| Observational simulation: interaction by adjustment strategy | `simulation/log_normal_prs_statin_interaction.r` | `simulation-outputs/mc_boxplot_interaction.png`, `mc_hist_interaction.png` |
| Observational simulation: interaction within baseline LDL-C tertiles | same | `simulation-outputs/mc_boxplot_interaction_strata_tertiles.png` |
| Observational simulation: tertiles with randomised treatment | same | `simulation-outputs/mc_boxplot_interaction_rand_tertiles.png` |
| Observational simulation: change-score designs | same | `simulation-outputs/mc_boxplot_change_score_designs.png` |
| UKB + HRN: interaction within baseline LDL-C tertiles | `figures/ukb_hrn_meta_analysis.r` | `figures/statin_ldl_stratified_forest_ukb_hrn_meta.png` |
| UKB + HRN: change-score designs | same | `figures/statin_ldl_change_score_forest_ukb_hrn_meta.png` |
| UKB-only results and diagnostics | `ukb/*.ipynb` | `ukb/*.png` |

---

## Software requirements

The simulations and figures were last run with **R 4.4.1** and these package versions:

| Package | Version | Package | Version |
|---|---|---|---|
| dplyr | 1.1.4 | survival | 3.2.13 |
| tidyr | 1.3.1 | stringr | 1.5.1 |
| purrr | 1.0.2 | ragg | 1.3.2 |
| ggplot2 | 3.5.2 | systemfonts | 1.1.0 |
| broom | 1.0.6 | this.path | 2.8.0 |

Install them with:

```r
install.packages(c("dplyr", "tidyr", "purrr", "ggplot2", "broom", "survival",
                   "stringr", "ragg", "systemfonts", "this.path"))
```

The UK Biobank notebooks also need `data.table` (1.15.4), `glmnet` (4.1.8), `lubridate` (1.9.3), `readxl`, `stringi`, `scales`, `patchwork`, `bit64` and `R.utils`. Each R notebook installs any missing packages when it starts. The Python notebooks need `dxpy`, `dxdata`, `pandas`, `numpy` and `pyspark`, all of which come preinstalled on the DNAnexus Research Analysis Platform (RAP).

**Fonts.** `figures/fps_style.r` uses Raleway if it is installed, then Source Sans Pro, then DejaVu Sans. Results do not depend on the font, but figure appearance does.

---

## Reproducing the simulations

The simulations need no external data. They can be fully reproduced from this repository.

### Paths

Every R script sets its working directory to its own folder with `setwd(this.path::here())`. Paths inside the scripts are relative to that folder, so the repository can be cloned anywhere and the scripts run from any directory:

```bash
git clone https://github.com/matt-tudball/prs-drug-response
cd prs-drug-response
```

`this.path::here()` needs the script to be run with `Rscript` or `source()`. It fails if the code is pasted line by line into an interactive session.

### 1. Trial simulations

```bash
Rscript simulation/clinical_trial_one_shot_sim_refactor.r
```

**Design.** A source population is simulated with the following structure:

```
G → LDL ← U          (G: CAD PRS; U: non-genetic risk; both N(0,1), independent)
G → B ← U, LDL → B   (B: latent atherosclerotic burden)
eligibility = f(B, LDL), trial-specific
A randomised among eligible participants
Y ~ A + G + U + LDL + B + A×U      (Cox-type event model; no A×G term)
```

Trial-specific inputs (in `trial_specs`):

| Trial | Source N | Trial N (treat / control) | Eligibility | Placebo event rate | PRS grouping |
|---|---|---|---|---|---|
| FOURIER | 600,000 | 14,298 (7,163 / 7,135) | top 5% of B and LDL ≥ 130 mg/dL | 11.5% | bottom 20% / middle / top 20% |
| ODYSSEY | 1,500,000 | 11,953 (5,967 / 5,986) | top 2% of B and LDL ≥ 140 mg/dL | 12.0% | top 10% vs rest |
| WOSCOPS | 300,000 | 4,892 (2,452 / 2,440) | LDL ≥ 174 mg/dL, excluding top 10% of B | 14.2% | top 20% vs rest |

For FOURIER and ODYSSEY, the LDL thresholds are back-calculated from the on-statin eligibility criteria to untreated-equivalent values. Population parameters are shared across trials and defined in `pars_anchor`, for example β<sub>A</sub> = log(0.80) and β<sub>AU</sub> = 0.30. The published subgroup HRs used for the compatibility analysis are in `published_trial_contrasts`.

Each simulation calibrates the baseline cumulative hazard to the target placebo event rate. It then fits:

- Cox models within each PRS group.
- Cox models with a treatment × PRS interaction, both unadjusted and adjusted for baseline LDL-C and a noisy proxy for B.

A simulation is counted as *compatible* with a published trial if its high-vs-lower PRS HR ratio falls within the published 95% CI for that ratio. That CI is derived from the published subgroup CIs by the delta method (`ci_hr_ratio()`).

**Settings** (top of script):

| Variable | Default | Meaning |
|---|---|---|
| `n_reps` | `1000` | Simulations per trial and per β<sub>AU</sub> grid point. Use 5 for a quick check. |
| `write_plots` | `TRUE` | Write PNGs and CSVs to `output_dir`. |
| `run_beta_AU_sensitivity` | `TRUE` | Repeat over β<sub>AU</sub> ∈ {0, 0.1, 0.2, 0.3, 0.4}, unadjusted and adjusted. |

**Seeds.** `set.seed(20260425)` is set globally. Each simulation then uses its own fixed seed: `20260425 + 10000·trial_index + r` for the main runs, `41000 + …` for the β<sub>AU</sub> grid and `51000 + …` for the adjusted grid. Results should therefore be identical for the same `n_reps`, R version and package versions.

**Runtime.** With `n_reps = 1000` and the sensitivity grid on, the script runs about 33,000 simulations, each with a source population of up to 1.5 million. Expect this to take many hours on one core. Run it in the background, for example with `nohup Rscript … &`.

**Outputs** (`simulation-outputs/`): the PNGs listed above, and CSV tables including `anchor_summary.csv`, `anchor_event_tables.csv`, `anchor_prs_decile_tables.csv`, `compatibility_probabilities.csv`, `compatibility_probability_table_pretty.csv`, `beta_AU_sensitivity.csv` and `beta_AU_sensitivity_adjusted_models.csv`. CSVs are git-ignored and must be regenerated.

### 2. Observational statin–LDL-C simulations

```bash
Rscript simulation/log_normal_prs_statin_interaction.r
```

**Design** (`simulate_1period()`):

```
log(LDL0) = 0.12·G + 0.15·U + ε,  ε ~ N(0, 0.25²)
P(A = 1)  = logit⁻¹(−2.5 + 0.03·G + log(LDL0))         # statin initiation
log(LDL1) = log(0.70)·A + log(LDL0) + 0.20·U − 0.25·A·U + ε,  ε ~ N(0, 0.20²)
```

There is no A×G term, so the true interaction is zero. Each simulation draws n = 78,145 individuals, and each analysis uses 1,000 simulations. The script estimates the A×G coefficient under:

- **Adjustment strategies:** no adjustment for LDL0; adjustment for log(LDL0); and a marginal structural model fitted by IPTW with the true propensity scores (stabilised, untruncated).
- **Baseline LDL-C tertiles:** IPTW within each tertile of log(LDL0). Cut points come from a pilot simulation of 500,000.
- **Tertiles with randomised treatment:** as above, with treatment independent of G and LDL0.
- **Change scores in initiators only (A = 1):** LDL1 − LDL0 ~ G; log(LDL1/LDL0) ~ G; and log(LDL1/LDL0) ~ G + log(LDL0).

The seed is `set.seed(20260114)`. The script runs in minutes to tens of minutes.

### 3. Collider schematic

```bash
Rscript figures/collider_example.r
```

This uses `set.seed(20260428)`, n = 1,200 and a selection threshold at the 88th percentile of G + U.

---

## Reproducing the UK Biobank analyses

### Data access

These analyses use individual-level UK Biobank data, which this repository does not include. To reproduce them you need:

1. An approved UK Biobank application covering primary-care (GP) records, Hospital Episode Statistics, death registry, standard PRS (Category 301) and biochemistry. The paper's analyses were conducted under application 27449.
2. A DNAnexus RAP project with the dispensed dataset.
3. The UK Biobank primary-care files `gp_clinical.txt`, `gp_scripts.txt` and `gp_registrations.txt`, plus the code-mapping workbook `all_lkps_maps_v4.xlsx` (UKB Resource 592), in the project root.

The notebooks contain placeholders for the DNAnexus project ID (`project-XXXXXXXXXXXXXXXXXXXXXXXX`), the Spark database name in `format_inpatient_data.ipynb` (`database_XXXXXXXXXXXXXXXXXXXXXXXX__appXXXXX_XXXXXXXXXXXXXX`) and the scratch folder for intermediate files (`/path/to/scratch/`). Replace them with your own before running.

The notebooks read and write files in the working directory of the Jupyter session, so they do not use `this.path::here()`.

### Run order

Run on the RAP in this order. Kernels: Python for 1, PySpark for 2, R for the rest.

| Step | Notebook | Produces |
|---|---|---|
| 1 | `format_ukb_covariates.ipynb` | `ukb_covariates_formatted.csv.gz` |
| 2 | `format_inpatient_data.ipynb` | `ukb_acute_coronary_events_inpatient_plus_death.csv.gz`, `ukb_pci_cabg_events.csv.gz` |
| 3 | `format_gp_data.ipynb` | `ukb_gp_formatted.csv.gz` |
| 4 | `run_analysis.ipynb` | `ukb_gp_for_analysis.rds`, full-cohort results, IPTW diagnostics |
| 5 | `run_stratified_analysis.ipynb` | UKB estimates within baseline LDL-C tertiles |
| 6 | `run_initiator_only_analysis.ipynb` | UKB change-score estimates |
| (S) | `format_gp_data_from_assessment.ipynb` → `run_analysis_from_assessment.ipynb` | Sensitivity analysis anchored at UKB assessment |

Steps 5 and 6 read the `.rds` file written in step 4.

### Key definitions

**Covariates (step 1).** UKB showcase fields 31, 34, 52, 53, 904, 1558, 2335, 2443, 4079, 4080, 6138, 6177, 20116, 21001, 21022, 22006, 22009 (PCs 1–40), 22189, 23099 and 30500–30880 (biochemistry). The standard PRS are 26223 (CVD), 26227 (CAD) and 26250 (LDL-C).

**Outcome events (step 2).**
- Acute coronary events: hospital admissions with a primary ICD-10 diagnosis of I21–I24, plus deaths with I21–I24 recorded anywhere among the causes. Events are de-duplicated by person and day.
- Revascularisation: OPCS-4 codes K49, K50 and K75 (PCI) and K40–K44 (CABG).

**Landmark cohort (step 3).**
- Time zero (t0) is a GP LDL-C measurement (Read code `44P6.`, 1–8 mmol/L) on or after 1990-01-01. One index date is used per person.
- Eligibility: prescription records observed for at least 365 days before and after t0, and no statin prescription in the 365 days before t0.
- Statins are identified by BNF chapter 2.12.1 or by drug name.
- Exposure: A = 1 if a statin was prescribed within 365 days after t0. The landmark is t0 + 365 days.
- Outcome: the first LDL-C measurement (LDL1) within 365 days after the landmark.
- Covariates: prior CAD (Read v2/v3 codes mapped to ICD-10 I20–I25), prior statin use, and GP utilisation (prescription and clinical contact days) in the year before t0.

**Models (steps 4–6).**
- Analyses are restricted to participants of white British ancestry (field 22006) with complete data.
- Outcome: log(LDL1). Exposure: A. Modifier: LDL PRS (`A_x_ldl_prs`).
- Three estimators:
  - Unadjusted.
  - Adjusted for log(LDL0), age, sex, GP utilisation, both PRS, prior statin use and prior CAD.
  - MSM/IPTW: the propensity model is a lasso logistic regression (`cv.glmnet`, `lambda.min`) on the same covariates. Weights are stabilised and truncated at the 5th and 95th percentiles. The outcome model includes age and sex only.
- The CAD analysis in `run_analysis.ipynb` uses incident acute coronary events within 10 years of the landmark, with an A × CAD PRS interaction in logistic models.
- Stratified analysis: the MSM/IPTW model, refitted within tertiles of LDL0. Weights are estimated in the full cohort.
- Change scores, in initiators only: LDL1 − LDL0; log(LDL1/LDL0); and log(LDL1/LDL0) adjusted for log(LDL0). Each is regressed on the LDL PRS with the adjustment covariates.

---

## Reproducing the UKB + HRN meta-analysis

```bash
Rscript figures/ukb_hrn_meta_analysis.r
```

This step needs only the two summary-statistic files in `figures/`. Each file holds one estimate, its 95% CI and its SE per cohort and analysis:

- `ukb_hrn_interaction_estimates.csv`: A × LDL-PRS interaction in the full cohort and within each baseline LDL-C tertile.
- `ukb_hrn_change_score_estimates.csv`: LDL-PRS coefficient for each change-score outcome, with sample size.

The UKB rows come from steps 5 and 6 above. The HRN rows come from the same analyses, run by the Helix team in the Helix Research Network (HRN) cohorts under HRN application 20250641. HRN individual-level data are not included in this repository. [TODO: HRN data access route and contact.]

The script pools estimates by fixed-effect inverse-variance meta-analysis. It reports Cochran's Q and I², writes `ukb_hrn_meta_analysis_results.csv` and `ukb_hrn_change_score_meta_analysis_results.csv`, and draws the two forest plots.

---

## Ethics and data availability

- **Helix / HRN.** The Helix cohorts were reviewed and approved by Salus IRB (reliance on Salus for all sites; approval 21143), WCG IRB (Western Institutional Review Board, WIRB Copernicus Group; approval 20224919), the MUSC Institutional Review Board for Human Research (approval Pro00129083) and the University of Nevada, Reno Institutional Review Board (approval 7701703417). All participants gave informed, written consent before taking part. All data used for research were de-identified. The HRN analyses were conducted under approved HRN application 20250641.
- **UK Biobank.** The UK Biobank analyses were conducted under application 27449. Researchers can apply for access at <https://www.ukbiobank.ac.uk/enable-your-research>.
- **Simulations.** The simulation analyses involved no human participants and required no additional ethics approval.

Only summary statistics (`figures/ukb_hrn_*_estimates.csv`) are distributed with this repository.

---

## Notes and known limitations

- `simulation/old-code/` holds earlier versions of the trial and observational simulations. They are kept for provenance and are not used for any reported result.
- Simulation CSV outputs are git-ignored by the `*.csv` rule. Only the two meta-analysis input CSVs are tracked.
- The UKB notebooks were run interactively on the RAP. Cell outputs have been cleared so that no individual-level data are distributed, and some exploratory cells remain.
- Exact UKB results depend on the data release and on participant withdrawals. [TODO: state the dataset release / dispensal date used.]

## Citation

[TODO: BibTeX entry once published.]

## Contact

Matt Tudball, 5 Prime Sciences: matt.tudball@5primesciences.com

## License

The code is released under the MIT License. See [LICENSE](LICENSE).
