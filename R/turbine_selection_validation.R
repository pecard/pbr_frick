##
## Turbine-level retrospective validation (Paulo, 2026-09): how well would
## an independently-derived turbine ranking -- built only from the raw
## per-carcass PCFM records and the project-specific working GenEst
## correction factor -- have reproduced the REAL turbine selection
## actually agreed with the lenders and documented in the Blanket
## Curtailment Plan technical note (May 2026)? (This script is itself the
## source of those project-specific factors, R/load_real_mortality_data.R
## -- see the note on that below.)
##
## Real methodology (from that document, Steps 1-5): rank turbines by
## GenEst-estimated mortality over a defined baseline season, cumulative-
## sum in ranked order, select turbines up to 50% of the season's total
## estimated mortality. Baseline seasons: Bash WPP, April-August 2025;
## Djangeldy WPP, March-August 2025. This produced 17 selected turbines at
## Bash (of 68) and 7 at Djangeldy (of 27), cumulative to exactly 50%.
##
## This script reproduces the SAME ranking procedure using only what an
## independent, project-level analysis could compute from the raw carcass
## log: per-turbine V. murinus carcass counts x the project-specific
## working factor (genest_correction_factor, R/load_real_mortality_data.R
## -- itself derived FROM the comparison this script performs, see the
## note below), over the same baseline windows. A per-project SCALAR
## multiplier does not change turbine rank order within a project (it
## scales every turbine's count equally), so which turbines this
## replication selects is unaffected by that circularity -- only the
## absolute mortality-value columns reported alongside the ranking are.
## The two real per-turbine reference tables (official rank, real GenEst estimate,
## real turbine code, real selection flag) are the actual Table 3/Table 4
## content from that document -- real turbine identifiers, kept local-only
## in data-raw/ (gitignored) alongside the raw carcass data, never
## committed. Anything derived here that reaches the report is anonymised
## via the OFFICIAL rank (turbine_anon = "T-<official rank, zero-padded>"),
## so the report can compare "my rank" against "their rank" for the same
## anonymised turbine without ever naming it.
##
## Known data-coverage caveat: the raw carcass file's earliest record may
## fall after the nominal window start (this has been the case for both
## projects at some point -- Bash's records starting partway into its
## nominal April 2025 window, Djangeldy's missing essentially all of its
## nominal March 2025 window). The actual data-start date used is reported
## explicitly alongside the comparison below, not assumed -- it is one of
## the two likely reasons for any turbine-level mismatch between the two
## rankings, the other being genuine per-turbine variability that even a
## validated PROJECT-level correction factor cannot capture, addressed
## further down.
##

suppressPackageStartupMessages({ library(dplyr); library(ggplot2) })

run_turbine_validation <- function(fig_dir = "outputs/figures",
                                    bash_path = "data-raw/BashWPP_Weekly_PCFM_PBR.xlsx",
                                    djangeldy_path = "data-raw/DjangeldyWPP_Weekly_PCFM_PBR.xlsx",
                                    official_bash_path = "data-raw/official_turbine_selection_bash.csv",
                                    official_djangeldy_path = "data-raw/official_turbine_selection_djangeldy.csv") {

  dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)
  fig_dir <- normalizePath(fig_dir, winslash = "/")

  if (!file.exists(official_bash_path) || !file.exists(official_djangeldy_path)) {
    stop("Official turbine selection reference tables not found in data-raw/. ",
         "These are extracted from the Blanket Curtailment Plan technical note and, like the raw ",
         "carcass data, are local-only and gitignored -- see R/turbine_selection_validation.R.")
  }

  ## Zero-pad the numeric suffix of a turbine code (raw data: "BSH1"; official
  ## table / this analysis: "BSH01") so the two sources join correctly.
  pad_bsh <- function(x) sprintf("BSH%02d", as.integer(sub("^BSH", "", x)))
  pad_dzh <- function(x) x  # Djangeldy raw codes are already zero-padded (DZH01, DZH03, ...)

  ## ---- My own replication: raw carcasses x4, same baseline windows ---------
  ## Reuses load_real_mortality_data() (R/load_real_mortality_data.R), which
  ## already applies the Group == "Bat" / Scheduled_search == "Yes" /
  ## Inside_search_plot == "Yes" / species == "V. murinus" filters -- rather
  ## than re-implementing that filtering logic here.
  d_all <- load_real_mortality_data(bash_path, djangeldy_path)
  bash_raw <- d_all %>% filter(project == facility_labels[1]) %>% select(turbine, date) %>% mutate(turbine = pad_bsh(turbine))
  dzh_raw  <- d_all %>% filter(project == facility_labels[2]) %>% select(turbine, date) %>% mutate(turbine = pad_dzh(turbine))

  bash_window <- c(as.Date("2025-04-01"), as.Date("2025-08-31"))
  dzh_window  <- c(as.Date("2025-03-01"), as.Date("2025-08-31"))
  bash_data_start <- min(bash_raw$date)
  dzh_data_start  <- min(dzh_raw$date)

  rank_and_select <- function(raw_dt, window, project_label) {
    raw_dt %>%
      filter(date >= window[1], date <= window[2]) %>%
      count(turbine, name = "raw") %>%
      mutate(my_estimated_mortality = raw * genest_correction_factor[project_label]) %>%
      arrange(desc(my_estimated_mortality)) %>%
      mutate(
        my_rank = row_number(),
        my_cum_pct = cumsum(my_estimated_mortality) / sum(my_estimated_mortality),
        my_selected = my_cum_pct <= 0.50 | dplyr::lag(my_cum_pct, default = 0) < 0.50
      )
  }

  bash_mine <- rank_and_select(bash_raw, bash_window, facility_labels[1])
  dzh_mine  <- rank_and_select(dzh_raw, dzh_window, facility_labels[2])

  ## ---- Official reference (real Table 3 / Table 4 content) ------------------
  official_bash <- read.csv(official_bash_path, stringsAsFactors = FALSE) %>%
    rename(turbine = turbine_real) %>%
    mutate(official_selected = as.logical(official_selected))
  official_dzh  <- read.csv(official_djangeldy_path, stringsAsFactors = FALSE) %>%
    rename(turbine = turbine_real) %>%
    mutate(official_selected = as.logical(official_selected))

  ## ---- Join, anonymise by OFFICIAL rank --------------------------------------
  build_comparison <- function(official, mine, label) {
    official %>%
      full_join(mine, by = "turbine") %>%
      mutate(
        official_rank = tidyr::replace_na(official_rank, NA_integer_),
        official_selected = tidyr::replace_na(official_selected, FALSE),
        my_selected = tidyr::replace_na(my_selected, FALSE),
        agreement = case_when(
          official_selected & my_selected   ~ "Selected by both",
          official_selected & !my_selected  ~ "Official only",
          !official_selected & my_selected  ~ "Replication only",
          TRUE ~ "Neither"
        ),
        turbine_anon = ifelse(!is.na(official_rank), sprintf("T-%02d", official_rank),
                               sprintf("X-%02d", my_rank)),
        project = label
      ) %>%
      arrange(official_rank)
  }

  comparison_bash <- build_comparison(official_bash, bash_mine, facility_labels[1])
  comparison_dzh  <- build_comparison(official_dzh, dzh_mine, facility_labels[2])
  comparison_all  <- bind_rows(comparison_bash, comparison_dzh)

  ## ---- Overlap statistics -----------------------------------------------------
  overlap_stats <- comparison_all %>%
    group_by(project) %>%
    summarise(
      n_official_selected = sum(official_selected, na.rm = TRUE),
      n_my_selected = sum(my_selected, na.rm = TRUE),
      n_both = sum(official_selected & my_selected, na.rm = TRUE),
      n_union = sum(official_selected | my_selected, na.rm = TRUE),
      jaccard = n_both / n_union,
      pct_of_official_recovered = 100 * n_both / n_official_selected,
      official_mortality_covered_by_my_selection = 100 *
        sum(official_estimated_mortality[official_selected & my_selected], na.rm = TRUE) /
        sum(official_estimated_mortality[official_selected], na.rm = TRUE),
      .groups = "drop"
    )

  ## ---- Where a flat rule would most over/under-shoot the real GenEst
  ## estimate, and the project-specific aggregate ratio this note now uses
  ## as genest_correction_factor (R/load_real_mortality_data.R) ------------
  ## (only meaningful where both a real and a replicated estimate exist)
  correction_check <- comparison_all %>%
    filter(!is.na(official_estimated_mortality), !is.na(my_estimated_mortality)) %>%
    mutate(implied_correction_ratio = official_estimated_mortality / pmax(raw, 0.5)) %>%
    select(project, turbine_anon, official_rank, my_rank, raw,
           official_estimated_mortality, my_estimated_mortality, implied_correction_ratio)
  correction_ratio_range <- range(correction_check$implied_correction_ratio)
  correction_ratio_median <- median(correction_check$implied_correction_ratio)
  correction_ratio_by_project <- correction_check %>%
    group_by(project) %>%
    summarise(
      median_ratio = median(implied_correction_ratio),
      min_ratio = min(implied_correction_ratio),
      max_ratio = max(implied_correction_ratio),
      aggregate_ratio = sum(official_estimated_mortality) / sum(raw),
      .groups = "drop"
    )

  ## ---- Figure: official rank vs. replicated rank, coloured by agreement ----
  rank_plot_dt <- comparison_all %>%
    filter(!is.na(official_rank) | !is.na(my_rank)) %>%
    mutate(
      official_rank_plot = tidyr::replace_na(official_rank, max(official_rank, na.rm = TRUE) + 3),
      my_rank_plot = tidyr::replace_na(my_rank, max(my_rank, na.rm = TRUE) + 3)
    )
  agreement_colours <- c(
    "Selected by both" = "forestgreen", "Official only" = "firebrick",
    "Replication only" = "steelblue", "Neither" = "grey70"
  )
  cutoff_dt <- overlap_stats %>% select(project, n_official_selected, n_my_selected)

  fig_turbine_validation <- file.path(fig_dir, "turbine_selection_validation.png")
  ggsave(fig_turbine_validation, width = 10, height = 5.5, dpi = 150, plot = {
    ggplot(rank_plot_dt, aes(x = official_rank_plot, y = my_rank_plot, colour = agreement)) +
      geom_vline(data = cutoff_dt, aes(xintercept = n_official_selected + 0.5), linetype = "dashed", colour = "grey40", linewidth = 0.4) +
      geom_hline(data = cutoff_dt, aes(yintercept = n_my_selected + 0.5), linetype = "dashed", colour = "grey40", linewidth = 0.4) +
      geom_abline(slope = 1, intercept = 0, linetype = "dotted", colour = "grey60") +
      geom_point(size = 1.8, alpha = 0.8) +
      facet_wrap(~project, scales = "free") +
      scale_colour_manual(name = "Agreement", values = agreement_colours) +
      labs(
        x = "Official rank (Blanket Curtailment Plan, Table 3/4)",
        y = "Independent replication rank (raw carcasses x project-specific working factor, same baseline window)",
        title = "Turbine-level retrospective validation: official vs. independently-replicated selection",
        subtitle = paste0(
          "Dashed lines mark each method's own 50%-cumulative-mortality cutoff; points below/left of both lines\n",
          "are selected by both. Turbines are anonymised (T-xx = official rank; X-xx = replication-only, no official match)."
        )
      ) +
      theme_minimal() +
      theme(plot.subtitle = element_text(size = 11), strip.text = element_text(face = "bold"),
            legend.position = "bottom")
  })

  list(
    turbine_validation_correction_ratio_range = correction_ratio_range,
    turbine_validation_correction_ratio_median = correction_ratio_median,
    turbine_validation_correction_ratio_by_project = correction_ratio_by_project,
    fig_turbine_validation = fig_turbine_validation,
    turbine_validation_bash_window = paste(format(bash_window, "%d %B %Y"), collapse = " - "),
    turbine_validation_dzh_window = paste(format(dzh_window, "%d %B %Y"), collapse = " - "),
    turbine_validation_bash_data_start = format(bash_data_start, "%d %B %Y"),
    turbine_validation_dzh_data_start = format(dzh_data_start, "%d %B %Y"),
    turbine_validation_comparison = comparison_all,
    turbine_validation_overlap = overlap_stats,
    turbine_validation_correction_check = correction_check
  )
}
