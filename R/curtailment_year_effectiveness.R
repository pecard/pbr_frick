##
## Curtailment-year effectiveness (Paulo, 2026-09): builds a full
## "operational year" running from the curtailment start date (5 May 2026)
## through the same date the following year, to evaluate the May-August
## 2026 emergency response on its own terms and set up the question of how
## it can be refined into a longer-term strategy.
##
## Two things this script does NOT claim:
## 1. It does not forecast what 2026/2027 mortality WILL be under the
##    current curtailment regime once it is reintroduced or continued --
##    that forward, current-regime forecast already exists in
##    R/adaptive_management_real.R (project_eoy()). This script instead
##    builds the complementary NO-CURTAILMENT counterfactual explicitly
##    requested: what mortality would have accumulated across the same
##    calendar weeks WITHOUT curtailment, using 2025 (the only full
##    pre-curtailment year available) as the reference, shifted forward
##    exactly one year week-by-week from the curtailment start date.
## 2. It does not claim the real, actually-recorded reduction extends
##    through August: this checkpoint covers however much of the response
##    window the current weekly workbooks (data-raw/*_Weekly_PCFM_PBR.xlsx)
##    happen to contain when this is run, honestly reported as such rather
##    than assumed to be the full May-August period Paulo described.
##

suppressPackageStartupMessages({ library(dplyr); library(lubridate); library(ggplot2) })

run_curtailment_year_effectiveness <- function(fig_dir = "outputs/figures",
                                                bash_path = "data-raw/BashWPP_Weekly_PCFM_PBR.xlsx",
                                                djangeldy_path = "data-raw/DjangeldyWPP_Weekly_PCFM_PBR.xlsx") {

  dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)
  fig_dir <- normalizePath(fig_dir, winslash = "/")

  d <- load_real_mortality_data(bash_path, djangeldy_path)
  project_thresholds <- tibble::tibble(project = facility_labels, threshold = pbr_thresholds$threshold)

  n_weeks <- 53
  op_week_start <- curtailment_start_date + lubridate::weeks(0:(n_weeks - 1))
  op_week_end   <- op_week_start + lubridate::days(6)
  ref_week_start <- op_week_start %m-% lubridate::years(1)
  ref_week_end   <- op_week_end %m-% lubridate::years(1)

  data_max_date <- max(d$date)

  week_grid <- tibble::tibble(
    week = 1:n_weeks,
    start_date = op_week_start, end_date = op_week_end,
    ref_start = ref_week_start, ref_end = ref_week_end
  )

  count_in_window <- function(proj, w_start, w_end) {
    sum(d$project == proj & d$date >= w_start & d$date <= w_end)
  }

  operational_year <- tidyr::expand_grid(project = project_thresholds$project, week = 1:n_weeks) %>%
    left_join(week_grid, by = "week") %>%
    rowwise() %>%
    mutate(
      actual_raw = count_in_window(project, start_date, pmin(end_date, data_max_date)),
      counterfactual_raw = count_in_window(project, ref_start, ref_end)
    ) %>%
    ungroup() %>%
    mutate(
      actual_corrected = ifelse(start_date <= data_max_date, actual_raw * genest_correction_factor, NA_real_),
      counterfactual_corrected = counterfactual_raw * genest_correction_factor
    ) %>%
    left_join(project_thresholds, by = "project")

  checkpoint_op_week <- max(operational_year$week[!is.na(operational_year$actual_corrected)])

  ## ---- Reduction achieved so far (real data window only) -------------------
  reduction_so_far <- operational_year %>%
    filter(week <= checkpoint_op_week) %>%
    group_by(project, threshold) %>%
    summarise(
      weeks_of_real_data = checkpoint_op_week,
      window_end = max(end_date[!is.na(actual_corrected)]),
      actual_corrected_total = sum(actual_corrected, na.rm = TRUE),
      counterfactual_corrected_total = sum(counterfactual_corrected),
      avoided_fatalities = counterfactual_corrected_total - actual_corrected_total,
      pct_reduction = 100 * avoided_fatalities / counterfactual_corrected_total,
      .groups = "drop"
    )

  ## ---- Full operational-year counterfactual (no curtailment at all) --------
  full_year_counterfactual <- operational_year %>%
    group_by(project, threshold) %>%
    summarise(
      counterfactual_annual_total = sum(counterfactual_corrected),
      pct_of_threshold = 100 * counterfactual_annual_total / threshold[1],
      .groups = "drop"
    )

  ## ---- Cumulative trajectories for the figure -------------------------------
  cumulative_dt <- operational_year %>%
    arrange(project, week) %>%
    group_by(project) %>%
    mutate(
      actual_cum = cumsum(coalesce(actual_corrected, 0)),
      actual_cum = ifelse(week <= checkpoint_op_week, actual_cum, NA_real_),
      counterfactual_cum = cumsum(counterfactual_corrected)
    ) %>%
    ungroup()

  fig_curtailment_year <- file.path(fig_dir, "curtailment_year_effectiveness.png")
  ggsave(fig_curtailment_year, width = 11, height = 5.5, dpi = 150, plot = {
    ggplot(cumulative_dt, aes(x = week)) +
      geom_line(aes(y = counterfactual_cum, linetype = "No-curtailment counterfactual (2025 pattern, +1 yr)"), colour = "firebrick", linewidth = 0.8) +
      geom_line(aes(y = actual_cum, linetype = "Actual (2026, with curtailment)"), colour = "forestgreen", linewidth = 1) +
      geom_hline(aes(yintercept = threshold), linetype = "dotted", colour = "grey30") +
      geom_vline(xintercept = checkpoint_op_week, linetype = "dotted", colour = "steelblue") +
      facet_wrap(~project, nrow = 1) +
      scale_linetype_manual(name = NULL, values = c("Actual (2026, with curtailment)" = "solid",
                                                      "No-curtailment counterfactual (2025 pattern, +1 yr)" = "dashed")) +
      labs(
        x = paste0("Week of the operational year (week 1 = ", format(curtailment_start_date, "%d %b %Y"), ")"),
        y = "Cumulative corrected fatalities (V. murinus)",
        title = "Operational year (May 2026-May 2027): actual vs. a no-curtailment counterfactual",
        subtitle = paste0(
          "Dotted horizontal: annual threshold. Dotted vertical: last week with real 2026 data (through ",
          format(data_max_date, "%d %b %Y"), ").\n",
          "Counterfactual: 2025's weekly pattern shifted one year forward -- what would have happened had curtailment never started."
        )
      ) +
      theme_minimal() +
      theme(plot.subtitle = element_text(size = 8), strip.text = element_text(face = "bold"),
            legend.position = "bottom")
  })

  list(
    curtailment_year_data_max_date = format(data_max_date, "%d %B %Y"),
    curtailment_year_checkpoint_week = checkpoint_op_week,
    curtailment_year_reduction_so_far = reduction_so_far,
    curtailment_year_full_counterfactual = full_year_counterfactual,
    fig_curtailment_year = fig_curtailment_year
  )
}
