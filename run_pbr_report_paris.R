##
## Launcher for the Paris Expert Workshop illustrative scenario -- a
## single fictional wind farm ("Wind Farm A"), generic "Bat model"
## species, 100% simulated Section 7 data. No real client data is read
## anywhere in this path (contrast run_pbr_report_BSH_DGY.R, which
## requires data-raw/*.xlsx).
##
## Source this file (Ctrl+Shift+S in RStudio, or
## source("run_pbr_report_paris.R") in the console) to render the
## workshop document.
##

pbr_settings_file <- "pbrSettings_paris.R"
output_file <- "outputs/pbr_paris_workshop.docx"
template    <- "report/pbr_paris_workshop_template.Rmd"

source("R/pbr_analysis.R")
source("R/simulate_paris_mortality_data.R")
source("R/paris_operational_analysis.R")
source("R/paris_turbine_validation.R")
source("R/paris_turbine_classification.R")
source("R/pbr_report.R")
source(file.path("inputs", pbr_settings_file))

pbr_params <- run_pbr_analysis(fig_dir = "outputs/figures")

sim <- simulate_paris_mortality_data()
operational_params <- run_paris_operational_analysis(fig_dir = "outputs/figures", sim = sim)
validation_params  <- run_paris_turbine_validation(fig_dir = "outputs/figures", sim = sim)
classification_params <- run_paris_turbine_classification(fig_dir = "outputs/figures", sim = sim)

report_params <- c(
  list(
    project_ref = project_ref,
    species_common_name = species_common_name,
    threshold = pbr_thresholds$threshold,
    n_turbines_total = n_turbines_total,
    n_turbines_curtailed_target = n_turbines_curtailed,
    correction_factor = genest_correction_factor_paris,
    nmin_assumed = pbr_params$nmin_assumed,
    pbr_corrected_fr = pbr_params$pbr_corrected_fr
  ),
  operational_params,
  validation_params,
  classification_params
)

build_pbr_report(
  output_file    = output_file,
  report_params  = report_params,
  template       = template,
  reference_docx = reference_docx_path
)
