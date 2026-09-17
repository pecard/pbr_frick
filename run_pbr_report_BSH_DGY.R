##
## Launcher for the BSH & DGY Parti-coloured bat PBR technical note.
## Source A ESTE FICHEIRO (Ctrl+Shift+S no RStudio, ou
## source("run_pbr_report_BSH_DGY.R") na consola) para gerar o relatorio.
##
## Para uma nova especie/facility/threshold (um novo cenario), copiar
## inputs/pbrSettings_BSH_DGY.R para inputs/pbrSettings_<cenario>.R,
## editar os literais la' dentro, e apontar pbr_settings_file abaixo para
## esse ficheiro -- nao e' preciso tocar em R/pbr_analysis.R,
## R/pbr_report.R nem no template.
##
## Mais pequeno que o equivalente run_monthly_report_BSH.R/
## IDF_monthly_report.R (sem etapa de importacao/cache de dados brutos:
## esta analise e' pura simulacao a partir dos literais em
## inputs/pbrSettings_*.R), por isso o launcher trata tambem do que la'
## fica no IDF_monthly_report.R intermedio.
##

pbr_settings_file <- "pbrSettings_BSH_DGY.R"
output_file        <- "outputs/pbr_frick_technical_note_BSH_DGY.docx"
lender_output_file <- "outputs/pbr_lender_summary_BSH_DGY.docx"
lender_template     <- "report/pbr_lender_summary_template.Rmd"

source("R/pbr_analysis.R")
source("R/adaptive_management_analysis.R")
source("R/load_real_mortality_data.R")
source("R/adaptive_management_real.R")
source("R/turbine_selection_validation.R")
source("R/curtailment_year_effectiveness.R")
source("R/global_pbr_uncertainty.R")
source("R/candidate_pbr_reference.R")
source("R/turbine_2027_classification.R")
source("R/pbr_report.R")
source(file.path("inputs", pbr_settings_file))

bash_data_path <- "data-raw/BashWPP_Weekly_PCFM_PBR.xlsx"
djangeldy_data_path <- "data-raw/DjangeldyWPP_Weekly_PCFM_PBR.xlsx"
official_bash_path <- "data-raw/official_turbine_selection_bash.csv"
official_djangeldy_path <- "data-raw/official_turbine_selection_djangeldy.csv"

if (!file.exists(bash_data_path) || !file.exists(djangeldy_data_path)) {
  message(
    "Real PCFM data not found at '", bash_data_path, "' / '", djangeldy_data_path, "' (gitignored, local-only) -- ",
    "section 7 will need them to render. See R/adaptive_management_demo.R for the synthetic version if only that is needed."
  )
  stop("Missing the weekly PCFM workbooks: place Paulo's BashWPP_Weekly_PCFM_PBR.xlsx and DjangeldyWPP_Weekly_PCFM_PBR.xlsx in data-raw/ before rendering the report.")
}
if (!file.exists(official_bash_path) || !file.exists(official_djangeldy_path)) {
  stop(
    "Official turbine selection reference tables not found in data-raw/ (gitignored, local-only). ",
    "See R/turbine_selection_validation.R."
  )
}

pbr_params <- run_pbr_analysis(fig_dir = "outputs/figures")
adaptive_params <- run_adaptive_management_real(fig_dir = "outputs/figures", bash_path = bash_data_path, djangeldy_path = djangeldy_data_path)
turbine_validation_params <- run_turbine_validation(
  fig_dir = "outputs/figures", bash_path = bash_data_path, djangeldy_path = djangeldy_data_path,
  official_bash_path = official_bash_path, official_djangeldy_path = official_djangeldy_path
)
curtailment_year_params <- run_curtailment_year_effectiveness(fig_dir = "outputs/figures", bash_path = bash_data_path, djangeldy_path = djangeldy_data_path)
global_pbr_params <- run_global_pbr_uncertainty(
  fig_dir = "outputs/figures", nmin_assumed = pbr_params$nmin_assumed, old_pbr_quantiles = pbr_params$pbr_quantiles
)
candidate_pbr_params <- run_candidate_pbr_reference(
  fig_dir = "outputs/figures", nmin_assumed = pbr_params$nmin_assumed,
  global_pbr_quantiles_by_project = global_pbr_params$global_pbr_quantiles_by_project,
  real_turbine_hist = adaptive_params$real_turbine_hist,
  curtailment_year_full_counterfactual = curtailment_year_params$curtailment_year_full_counterfactual
)
turbine_2027_params <- run_turbine_2027_classification(
  fig_dir = "outputs/figures", bash_path = bash_data_path, djangeldy_path = djangeldy_data_path,
  official_bash_path = official_bash_path, official_djangeldy_path = official_djangeldy_path
)

report_params <- c(
  pbr_params,
  adaptive_params,
  turbine_validation_params,
  curtailment_year_params,
  global_pbr_params,
  candidate_pbr_params,
  turbine_2027_params
)

build_pbr_report(
  output_file    = output_file,
  report_params  = report_params,
  reference_docx = reference_docx_path
)

build_pbr_report(
  output_file    = lender_output_file,
  report_params  = report_params,
  template       = lender_template,
  reference_docx = reference_docx_path
)
