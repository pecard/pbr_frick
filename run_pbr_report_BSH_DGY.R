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

source("R/pbr_analysis.R")
source("R/adaptive_management_analysis.R")
source("R/load_real_mortality_data.R")
source("R/adaptive_management_real.R")
source("R/pbr_report.R")
source(file.path("inputs", pbr_settings_file))

real_data_path <- "data-raw/pcfm_bat_summary.xlsx"
adaptive_params <- if (file.exists(real_data_path)) {
  run_adaptive_management_real(fig_dir = "outputs/figures", xlsx_path = real_data_path)
} else {
  message(
    "Real PCFM data not found at '", real_data_path, "' (gitignored, local-only) -- ",
    "section 7 will need it to render. See R/adaptive_management_demo.R for the synthetic version if only that is needed."
  )
  stop("Missing '", real_data_path, "': place Paulo's pcfm_bat_summary.xlsx there before rendering the report.")
}

report_params <- c(
  run_pbr_analysis(fig_dir = "outputs/figures"),
  adaptive_params
)

build_pbr_report(
  output_file    = output_file,
  report_params  = report_params,
  reference_docx = reference_docx_path
)
