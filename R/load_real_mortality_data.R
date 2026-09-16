##
## Loads real weekly PCFM carcass data for V. murinus from
## data-raw/pcfm_bat_summary.xlsx (Paulo, 2026-09) -- gitignored, never
## committed; this script only reads it. Only the Carcass_BSH and
## Carcass_DGY sheets are used, per Paulo's instruction (Bash/Djangeldy/
## BCKP_* are weekly all-species/activity summaries and backups, not the
## per-carcass record this analysis needs).
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

load_real_mortality_data <- function(xlsx_path = "data-raw/pcfm_bat_summary.xlsx") {
  if (!file.exists(xlsx_path)) {
    stop("Real PCFM data not found at '", xlsx_path, "'. This file is gitignored and local-only -- ",
         "place Paulo's pcfm_bat_summary.xlsx there to run the real operational analysis.")
  }

  read_project <- function(sheet, species_col, real_name, project_label) {
    read_excel(xlsx_path, sheet = sheet) %>%
      rename(species = !!species_col) %>%
      filter(species == "Vespertilio murinus") %>%
      transmute(
        project = project_label,
        turbine = Turbine,
        date = as.Date(Date_found),
        year = lubridate::year(date),
        iso_week = lubridate::isoweek(date)
      )
  }

  bind_rows(
    read_project("Carcass_BSH", "Species_name", "Bash", facility_labels[1]),
    read_project("Carcass_DGY", "Species_name (lat.)", "Djangeldy", facility_labels[2])
  ) %>%
    mutate(
      period = case_when(
        year < 2026 ~ "2025 (pre-curtailment baseline)",
        date < curtailment_start_date ~ "2026, pre-curtailment",
        TRUE ~ "2026, post-curtailment"
      )
    )
}
