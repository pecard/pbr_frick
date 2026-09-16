##
## Real operational example: applies the adaptive-management engine
## (R/adaptive_management_analysis.R) to the actual 2025-2026 PCFM carcass
## data (Paulo, 2026-09; R/load_real_mortality_data.R). Vespertilio murinus
## only, GenEst-corrected (x4 multiplicative rule on raw carcasses found,
## per Paulo -- a project-level rule of thumb, not a fitted detection model).
##
## Framing: curtailment started 5 May 2026. Months before that in 2026 (and
## all of 2025) are historical reference, not part of the "how is the
## current regime doing" question -- but they DO count toward the year's
## total against the annual PBR-derived reference, which the data show was
## already exceeded before curtailment even began. Both facts are reported:
## the already-exhausted annual budget (a real, factual finding), and the
## forward-looking question that is actually still actionable -- is the
## post-curtailment regime, on its own terms, trending toward a materially
## better trajectory, and where does mortality concentrate.
##
## The empirical seasonal risk profile (q_t) is estimated from 2025 -- the
## only full pre-curtailment year available -- pooled across both projects
## for a less noisy weekly estimate (consistent with treating this as
## potentially one regional/migratory population, discussed earlier in the
## report). One year of data cannot fully separate true seasonal risk from
## that year's particular weather/timing -- flagged wherever this profile
## is used, and the first thing that should update once a second
## pre-curtailment year (or more) becomes available.
##

suppressPackageStartupMessages({ library(dplyr); library(tidyr); library(ggplot2); library(lubridate) })

run_adaptive_management_real <- function(fig_dir, xlsx_path = "data-raw/pcfm_bat_summary.xlsx") {

  dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)
  fig_dir <- normalizePath(fig_dir, winslash = "/")

  d <- load_real_mortality_data(xlsx_path)
  project_thresholds <- tibble::tibble(project = facility_labels, threshold = pbr_thresholds$threshold)
  weeks_full <- tibble::tibble(iso_week = 1:53)
  curtailment_week <- lubridate::isoweek(curtailment_start_date)

  ## ---- Empirical seasonal risk profile from 2025 (pooled) -------------------
  weekly_2025 <- d %>% filter(period == "2025 (pre-curtailment baseline)") %>% count(iso_week, name = "n")
  seasonal_2025 <- weeks_full %>% left_join(weekly_2025, by = "iso_week") %>% mutate(n = tidyr::replace_na(n, 0))
  seasonal_weight_real <- seasonal_2025$n / sum(seasonal_2025$n)

  ## ---- 2026 weekly corrected counts per project ------------------------------
  weekly_2026 <- d %>% filter(year == 2026) %>% count(project, iso_week, name = "raw") %>%
    tidyr::complete(project = project_thresholds$project, iso_week = 1:53, fill = list(raw = 0)) %>%
    mutate(corrected = raw * genest_correction_factor) %>%
    left_join(project_thresholds, by = "project")

  checkpoint_week <- d %>% filter(year == 2026) %>% summarise(m = max(iso_week)) %>% pull(m)

  ## ---- Factual headline: full-year-to-date vs. annual reference -------------
  overage_summary <- weekly_2026 %>%
    filter(iso_week <= checkpoint_week) %>%
    group_by(project, threshold) %>%
    summarise(raw = sum(raw), corrected = sum(corrected), .groups = "drop") %>%
    mutate(pct_of_threshold = 100 * corrected / threshold)

  ## ---- Year-over-year comparison, same checkpoint week -----------------------
  weekly_2025_by_project <- d %>% filter(year == 2025) %>% count(project, iso_week, name = "raw") %>%
    tidyr::complete(project = project_thresholds$project, iso_week = 1:53, fill = list(raw = 0)) %>%
    mutate(corrected = raw * genest_correction_factor)
  yoy_summary <- bind_rows(
    weekly_2025_by_project %>% filter(iso_week <= checkpoint_week) %>% mutate(year = 2025),
    weekly_2026 %>% filter(iso_week <= checkpoint_week) %>% mutate(year = 2026)
  ) %>%
    group_by(project, year) %>% summarise(corrected = sum(corrected), .groups = "drop") %>%
    tidyr::pivot_wider(names_from = year, values_from = corrected, names_prefix = "corrected_") %>%
    mutate(pct_change = 100 * (corrected_2026 - corrected_2025) / corrected_2025)

  ## ---- MRI trajectory across 2026 to date, with real trigger status --------
  mri_data <- lapply(seq_len(nrow(project_thresholds)), function(i) {
    proj <- project_thresholds$project[i]; thr <- project_thresholds$threshold[i]
    build_control_trajectory(thr, seasonal_weight_real) %>%
      left_join(
        weekly_2026 %>% filter(project == proj) %>% select(iso_week, corrected) %>% rename(week = iso_week, observed_weekly = corrected),
        by = "week"
      ) %>%
      mutate(observed_weekly = tidyr::replace_na(observed_weekly, 0), observed_cum = cumsum(observed_weekly),
             project = proj, threshold = thr) %>%
      filter(week <= checkpoint_week)
  }) %>%
    bind_rows() %>%
    mutate(trigger = classify_trigger(observed_cum, expected_cum, ucl_cum)) %>%
    left_join(trigger_response, by = "trigger")

  trigger_colours <- c(GREEN = "forestgreen", AMBER = "orange", RED = "firebrick")
  fig_real_dashboard <- file.path(fig_dir, "adaptive_management_real.png")
  ggsave(fig_real_dashboard, width = 11, height = 6.5, dpi = 150, plot = {
    ggplot(mri_data, aes(x = week)) +
      geom_ribbon(aes(ymin = lcl_cum, ymax = ucl_cum), fill = "grey70", alpha = 0.3) +
      geom_line(aes(y = expected_cum), linetype = "dashed", colour = "grey30", linewidth = 0.6) +
      geom_vline(xintercept = curtailment_week, linetype = "dotted", colour = "steelblue", linewidth = 0.7) +
      geom_line(aes(y = observed_cum, colour = trigger, group = project), linewidth = 1) +
      geom_point(aes(y = observed_cum, colour = trigger), size = 1.5) +
      facet_wrap(~project, nrow = 1) +
      scale_colour_manual(name = "Trigger status", values = trigger_colours) +
      labs(
        x = "ISO week, 2026", y = "Cumulative corrected fatalities (V. murinus)",
        title = "Real 2026 cumulative mortality vs. the annual PBR reference",
        subtitle = paste0(
          "Dashed line: expected pace from the 2025 empirical seasonal profile; grey band: 90% Poisson control limits;\n",
          "blue dotted line: curtailment start (5 May, ISO week ", curtailment_week, "). GenEst-corrected (x4 raw carcasses)."
        )
      ) +
      theme_minimal() +
      theme(plot.subtitle = element_text(size = 9.5), strip.text = element_text(face = "bold"))
  })

  ## ---- Forward forecast under the CURRENT (post-curtailment) regime --------
  post_curt_lookback <- max(1, checkpoint_week - curtailment_week + 1)
  forecast_summary <- lapply(seq_len(nrow(project_thresholds)), function(i) {
    proj <- project_thresholds$project[i]; thr <- project_thresholds$threshold[i]
    obs_full <- weeks_full %>%
      left_join(weekly_2026 %>% filter(project == proj) %>% select(iso_week, corrected), by = "iso_week") %>%
      mutate(corrected = tidyr::replace_na(corrected, 0)) %>% pull(corrected)
    eoy <- project_eoy(obs_full, checkpoint_week, seasonal_weight_real, thr, lookback = post_curt_lookback)
    bind_cols(tibble::tibble(project = proj, threshold = thr, week = checkpoint_week, post_curt_weeks_observed = post_curt_lookback), eoy)
  }) %>% bind_rows()

  ## ---- Historical turbine concentration (pre-curtailment, both years) ------
  ## Turbine IDs are anonymised to a rank-based label (T-01, T-02, ...) --
  ## the raw codes (e.g. "BSH36") embed the real facility name, which must
  ## not appear in report-facing figures given the Project 1/Project 2
  ## anonymisation already applied throughout this note.
  turbine_hist <- d %>%
    filter(period != "2026, post-curtailment") %>%
    count(project, turbine, name = "raw") %>%
    mutate(corrected = raw * genest_correction_factor) %>%
    group_by(project) %>%
    arrange(desc(corrected), .by_group = TRUE) %>%
    mutate(rank = row_number(), cum_pct = cumsum(corrected) / sum(corrected),
           turbine_anon = sprintf("T-%02d", rank)) %>%
    ungroup()

  focus_project <- project_thresholds$project[which.max(overage_summary$corrected)]
  turbine_focus <- turbine_hist %>% filter(project == focus_project)
  n_turbines_focus <- nrow(turbine_focus)
  n_turbines_top80 <- min(which(turbine_focus$cum_pct >= 0.80))

  fig_real_turbine_pareto <- file.path(fig_dir, "turbine_pareto_real.png")
  ggsave(fig_real_turbine_pareto, width = 10, height = 5.5, dpi = 150, plot = {
    ggplot(turbine_focus, aes(x = reorder(turbine_anon, rank))) +
      geom_col(aes(y = corrected), fill = "firebrick", alpha = 0.8) +
      geom_line(aes(y = cum_pct * max(corrected), group = 1), colour = "grey20", linewidth = 0.6) +
      geom_point(aes(y = cum_pct * max(corrected)), colour = "grey20", size = 1) +
      scale_y_continuous(
        name = "GenEst-corrected fatalities (pre-curtailment period)",
        sec.axis = sec_axis(~ . / max(turbine_focus$corrected), name = "Cumulative % of mortality", labels = scales::percent)
      ) +
      labs(
        x = "Turbine (ranked by pre-curtailment corrected mortality)",
        title = paste0("Where pre-curtailment V. murinus mortality concentrated (", focus_project, ")"),
        subtitle = paste0(
          "2025 + 2026 pre-curtailment records only, GenEst-corrected (x4). ", n_turbines_top80, " of ", n_turbines_focus,
          " turbines (", round(100 * n_turbines_top80 / n_turbines_focus), "%) account for 80% of pre-curtailment mortality."
        )
      ) +
      theme_minimal() +
      theme(axis.text.x = element_text(angle = 90, vjust = 0.5, size = 6), plot.subtitle = element_text(size = 8.5))
  })

  list(
    real_correction_factor = genest_correction_factor,
    real_curtailment_start = format(curtailment_start_date, "%d %B %Y"),
    real_curtailment_week = curtailment_week,
    real_checkpoint_week = checkpoint_week,
    real_overage_summary = overage_summary,
    real_yoy_summary = yoy_summary,
    real_forecast_summary = forecast_summary,
    real_focus_project = focus_project,
    real_n_turbines_focus = n_turbines_focus,
    real_n_turbines_top80 = n_turbines_top80,
    fig_real_dashboard = fig_real_dashboard,
    fig_real_turbine_pareto = fig_real_turbine_pareto
  )
}
