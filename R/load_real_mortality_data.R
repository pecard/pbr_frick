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
## Raw carcass counts are corrected for detection probability using a
## project-specific multiplicative factor -- NOT a per-search-effort
## GenEst model fit in its own right, but no longer a single flat rule
## either (Paulo, 2026-09): each value is the aggregate, mortality-
## weighted ratio between the real GenEst-estimated mortality per turbine
## (Blanket Curtailment Plan, Table 3/Table 4) and the matching raw
## carcass count over the same baseline season, computed and validated in
## R/turbine_selection_validation.R. A single shared factor -- the
## originally-assumed flat 4x, or even the combined per-turbine median
## (~4.3x) -- was checked against this and rejected: Project 2's real
## ratio (6.21x) is meaningfully higher than Project 1's (4.64x), so a
## shared factor would have understated Project 2's true corrected
## mortality by roughly a third. Per-turbine variability remains even
## within a project (Project 1: 2.9x-9.75x; Project 2: 3.9x-12.9x) --
## these project-level aggregates are the best available working
## estimate given the real GenEst output already published for the
## baseline season, not a claim of turbine-level precision; flagged as
## such wherever this is used.
##

suppressPackageStartupMessages({ library(readxl); library(dplyr); library(lubridate) })

genest_correction_factor <- c("Project 1" = 4.64, "Project 2" = 6.21)
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

  d <- bind_rows(
    read_project(bash_path, "Species_name", facility_labels[1]),
    read_project(djangeldy_path, "Species_name (lat.)", facility_labels[2])
  )

  ## Data-quality guard: a Date_found in the future is impossible -- a
  ## carcass cannot be found before it happens. Caught in practice
  ## (2026-09): two Bash rows dated 25/27 September 2026 whose own ISO
  ## week column read 35 (~24-30 August), i.e. a month-digit typo in the
  ## source workbook (08 entered as 09), silently pulling the pipeline's
  ## "current" checkpoint 2-3 weeks into the future. Excluded here with a
  ## loud warning rather than guessed-and-corrected, since only Paulo can
  ## fix the source row; re-run after he does to pick the record back up
  ## with its real date.
  future_records <- d %>% filter(date > Sys.Date())
  if (nrow(future_records) > 0) {
    warning(
      "Excluded ", nrow(future_records), " record(s) with Date_found in the future (impossible -- ",
      "likely a data-entry typo in the source workbook): ",
      paste(sprintf("%s/%s/%s", future_records$project, future_records$turbine, future_records$date), collapse = "; "),
      ". Fix the source row(s) and re-run to include the real record."
    )
    d <- d %>% filter(date <= Sys.Date())
  }

  d %>%
    mutate(
      period = case_when(
        year < 2026 ~ "2025 (pre-curtailment baseline)",
        date < curtailment_start_date ~ "2026, pre-curtailment",
        TRUE ~ "2026, post-curtailment"
      )
    )
}
