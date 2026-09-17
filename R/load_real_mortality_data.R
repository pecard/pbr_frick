##
## Loads real weekly PCFM carcass data for V. murinus from the two
## per-project weekly workbooks Paulo receives (2026-09) --
## data-raw/BashWPP_Weekly_PCFM_PBR.xlsx and
## data-raw/DjangeldyWPP_Weekly_PCFM_PBR.xlsx, both gitignored, never
## committed; this script only reads them. Both files are meant to be
## overwritten in place each week with the latest export (generic names,
## no week number), so this loader always reads "the current file", not a
## fixed historical snapshot.
##
## Each workbook's "Carcass_ID" sheet holds every carcass search record
## (all taxa, all search types), not just the usable bat records for this
## analysis -- three filters select the right subset:
##   Group == "Bat": excludes birds and other taxa recorded in the same log;
##   Scheduled_search == "Yes": excludes carcasses found incidentally
##     outside a scheduled PCFM search, which are not part of the
##     standardised search effort GenEst-style correction assumes;
##   Inside_search_plot == "Yes": excludes carcasses found outside the
##     turbine's defined search plot, for the same reason.
## (Paulo, 2026-09.) The species filter to V. murinus itself is applied
## after these three, on top of the same species-name-column difference
## between files already handled below.
##
## Project_Name in the raw data is the real facility name ("Bash",
## "Djangeldy"); mapped here to the same facility_labels/threshold order
## already established in inputs/pbrSettings_BSH_DGY.R (Bash = BSH =
## 112 bats/yr = "Project 1"; Djangeldy = DGY = 88 bats/yr = "Project 2")
## so nothing downstream of this file, including the report, ever sees the
## real names -- consistent with the anonymisation already done for the
## workshop presentation.
##
## Raw carcass counts are corrected for detection probability using
## Paulo's stated multiplicative rule (GenEst-style unbiased estimate =
## 4 x raw carcasses found) -- a project-level rule of thumb, not a
## per-search-effort GenEst model fit; flagged as such wherever it's used.
##

suppressPackageStartupMessages({ library(readxl); library(dplyr); library(lubridate) })

genest_correction_factor <- 4
curtailment_start_date <- as.Date("2026-05-05")

load_real_mortality_data <- function(bash_path = "data-raw/BashWPP_Weekly_PCFM_PBR.xlsx",
                                      djangeldy_path = "data-raw/DjangeldyWPP_Weekly_PCFM_PBR.xlsx") {
  if (!file.exists(bash_path)) {
    stop("Real PCFM data not found at '", bash_path, "'. This file is gitignored and local-only -- ",
         "place Paulo's BashWPP_Weekly_PCFM_PBR.xlsx there to run the real operational analysis.")
  }
  if (!file.exists(djangeldy_path)) {
    stop("Real PCFM data not found at '", djangeldy_path, "'. This file is gitignored and local-only -- ",
         "place Paulo's DjangeldyWPP_Weekly_PCFM_PBR.xlsx there to run the real operational analysis.")
  }

  read_project <- function(xlsx_path, species_col, project_label) {
    read_excel(xlsx_path, sheet = "Carcass_ID") %>%
      rename(species = !!species_col) %>%
      filter(Group == "Bat", Scheduled_search == "Yes", Inside_search_plot == "Yes",
             species == "Vespertilio murinus") %>%
      transmute(
        project = project_label,
        turbine = Turbine,
        date = as.Date(Date_found),
        year = lubridate::year(date),
        iso_week = lubridate::isoweek(date)
      )
  }

  bind_rows(
    read_project(bash_path, "Species_name", facility_labels[1]),
    read_project(djangeldy_path, "Species_name (lat.)", facility_labels[2])
  ) %>%
    mutate(
      period = case_when(
        year < 2026 ~ "2025 (pre-curtailment baseline)",
        date < curtailment_start_date ~ "2026, pre-curtailment",
        TRUE ~ "2026, post-curtailment"
      )
    )
}
