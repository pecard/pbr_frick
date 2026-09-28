##
## Support for writing Section "What this shows, and what still needs to
## happen" (Paulo, 2026-09): that closing synthesis is hand-written prose,
## not computed from params -- deliberately, since it is an editorial
## judgement call across several signals, not a number to report. But
## several of its sentences assert something about THIS run's specific
## results ("directionally sound", "facility-differentiated signals",
## "only one untreated year", "part of the year uncovered") that nothing
## in the code actually re-checks -- if a future run's data told a
## different story, the prose would not know.
##
## This does not generate that prose. It computes the underlying facts
## each claim rests on, as a plain checklist, so whoever writes/revises
## the section (Paulo, or an LLM assisting him) has the current numbers
## in front of them instead of re-deriving them by eye or trusting the
## last version's wording. Not rendered in the client-facing docx --
## printed to the console and written to a side file for drafting use.
##

suppressPackageStartupMessages({ library(dplyr) })

check_section_interpretation_claims <- function(turbine_validation_overlap,
                                                 curtailment_year_reduction_so_far,
                                                 curtailment_year_checkpoint_week,
                                                 bash_path = "data-raw/BashWPP_Weekly_PCFM_PBR.xlsx",
                                                 djangeldy_path = "data-raw/DjangeldyWPP_Weekly_PCFM_PBR.xlsx",
                                                 facility_diff_threshold_pp = 10,
                                                 validation_sound_threshold_pct = 60,
                                                 validation_weak_threshold_pct = 30,
                                                 full_year_weeks = 52) {

  d <- suppressWarnings(load_real_mortality_data(bash_path, djangeldy_path))

  ## ---- Claim 1: "facility-differentiated signals" ---------------------------
  reductions <- curtailment_year_reduction_so_far %>% arrange(project) %>% pull(pct_reduction)
  facility_diff_pp <- if (length(reductions) == 2) abs(diff(reductions)) else NA_real_
  claim_facility_diff <- tibble::tibble(
    claim = "Facilities show differentiated (not near-identical) reduction signals",
    metric = "abs difference in pct_reduction between projects (pp)",
    value = round(facility_diff_pp, 1),
    status = case_when(
      is.na(facility_diff_pp) ~ "CANNOT CHECK -- need exactly 2 projects",
      facility_diff_pp >= facility_diff_threshold_pp ~ "HOLDS -- keep 'facility-differentiated' framing",
      TRUE ~ paste0("REVIEW -- projects within ", facility_diff_threshold_pp, "pp of each other; soften or drop the claim")
    )
  )

  ## ---- Claim 2: "only one untreated year, none will exist while
  ## curtailment continues" -----------------------------------------------------
  n_untreated_years <- d %>%
    filter(period == "2025 (pre-curtailment baseline)") %>%
    summarise(n = dplyr::n_distinct(lubridate::year(date))) %>%
    pull(n)
  claim_one_untreated_year <- tibble::tibble(
    claim = "Seasonal risk profile rests on a single untreated year, unimprovable while curtailment continues",
    metric = "distinct pre-curtailment-baseline years in the data",
    value = n_untreated_years,
    status = if (n_untreated_years <= 1) {
      "HOLDS -- still true, keep as written"
    } else {
      paste0("STALE -- ", n_untreated_years, " untreated years now available; revisit both the seasonal-profile ",
             "claim and the DiD-not-pursued reasoning below it")
    }
  )

  ## ---- Claim 3: turbine validation "directionally sound" --------------------
  claim_validation_sound <- turbine_validation_overlap %>%
    transmute(
      claim = paste0("Independent turbine-selection replication is directionally sound (", project, ")"),
      metric = "% of officially-selected turbines also recovered by the independent replication",
      value = round(pct_of_official_recovered, 0),
      status = case_when(
        pct_of_official_recovered >= validation_sound_threshold_pct ~ "HOLDS -- keep 'directionally sound'",
        pct_of_official_recovered >= validation_weak_threshold_pct ~
          paste0("MIXED -- ", round(pct_of_official_recovered), "%; soften to 'partial' or 'mixed' agreement"),
        TRUE ~ paste0("WEAK -- only ", round(pct_of_official_recovered), "%; do not describe as sound without caveats")
      )
    )

  ## ---- Claim 4: "current curtailment window leaves part of the
  ## operational year uncovered" ------------------------------------------------
  claim_year_coverage <- tibble::tibble(
    claim = "Operational year is only partially covered by real data so far",
    metric = paste0("weeks of real data out of ", full_year_weeks),
    value = curtailment_year_checkpoint_week,
    status = if (curtailment_year_checkpoint_week >= full_year_weeks) {
      "STALE -- full operational year now covered; drop the 'part of the year uncovered' framing"
    } else {
      paste0("HOLDS -- ", full_year_weeks - curtailment_year_checkpoint_week, " week(s) still outstanding")
    }
  )

  checklist <- bind_rows(claim_facility_diff, claim_one_untreated_year, claim_validation_sound, claim_year_coverage)

  list(
    section_interpretation_checklist = checklist
  )
}
