##
## Adaptive mortality management framework -- PROPOSED methodology, wired
## into the report as its own analysis function (mirrors run_pbr_analysis()'s
## "R computes, Rmd only formats" convention). Ported from the exploratory
## R/adaptive_management_dashboard.R (kept standalone for reference).
##
## Framing (Paulo, 2026-09): keep the PBR-derived annual reference
## relatively stable, and make the RESPONSE to monitoring dynamic instead
## -- three layers: (1) stable biological reference (the PBR/Niel-Lebreton
## analysis above, revisited only on materially new evidence), (2) dynamic
## mortality assessment (this section: weekly, from PCFM), (3) dynamic,
## turbine-level curtailment (the operational output of (2)).
##
## This file holds the SHARED ENGINE (control trajectory, trigger tiers,
## EOY forecast, required reduction, turbine allocation) as top-level
## functions, reused by both this illustrative/synthetic demo and
## R/adaptive_management_real.R (the real 2025-2026 PCFM analysis). The
## demo below is entirely synthetic -- see R/adaptive_management_real.R
## for the real analysis the report actually presents.
##

suppressPackageStartupMessages({ library(dplyr); library(tidyr); library(ggplot2) })

## ---- Expected cumulative trajectory + Poisson control limits -------------
build_control_trajectory <- function(annual_threshold, seasonal_weight, ci_level = 0.90) {
  alpha <- 1 - ci_level
  expected_cum <- annual_threshold * cumsum(seasonal_weight)
  tibble::tibble(
    week = seq_along(seasonal_weight),
    expected_weekly = annual_threshold * seasonal_weight,
    expected_cum = expected_cum,
    ucl_cum = qpois(1 - alpha / 2, pmax(expected_cum, 1e-6)),
    lcl_cum = qpois(alpha / 2, pmax(expected_cum, 1e-6))
  )
}

## ---- Trigger tiers + curtailment escalation ladder ------------------------
classify_trigger <- function(observed_cum, expected_cum, ucl_cum) {
  dplyr::case_when(
    observed_cum <= expected_cum ~ "GREEN",
    observed_cum <= ucl_cum ~ "AMBER",
    TRUE ~ "RED"
  )
}
trigger_response <- tibble::tibble(
  trigger = c("GREEN", "AMBER", "RED"),
  curtailment_action = c(
    "Maintain current prescription on critical turbines; scheduled review at next checkpoint.",
    "Increase curtailment intensity on critical turbines (lower cut-in speed by one defined step, or extend curtailed hours); reassess weekly.",
    "Escalate: widen curtailment to the next turbine tier and/or move to the most conservative cut-in speed on critical turbines for the remainder of the high-risk window; trigger a structured review; document for lenders."
  )
)

## ---- Stochastic end-of-season forecast under the CURRENT regime ----------
project_eoy <- function(observed_weekly_vec, current_week, seasonal_weight, annual_threshold,
                         lookback = 6, n_sim = 4000) {
  expected_weekly <- annual_threshold * seasonal_weight
  recent_weeks <- max(1, current_week - lookback + 1):current_week
  recent_expected <- sum(expected_weekly[recent_weeks])
  recent_observed <- sum(observed_weekly_vec[recent_weeks])
  intensity <- if (recent_expected > 0.05) recent_observed / recent_expected else 1
  future_weeks <- seq(current_week + 1, length(seasonal_weight))
  future_expected <- pmax(expected_weekly[future_weeks] * intensity, 0)
  sims <- if (length(future_weeks) > 0) {
    replicate(n_sim, sum(rpois(length(future_weeks), future_expected)))
  } else rep(0, n_sim)
  m_t <- sum(observed_weekly_vec[1:current_week])
  eoy_dist <- m_t + sims
  tibble::tibble(
    m_t = m_t, regime_intensity = intensity,
    eoy_median = median(eoy_dist), eoy_p90 = quantile(eoy_dist, 0.90, names = FALSE),
    p_exceed = mean(eoy_dist > annual_threshold),
    remaining_forecast_median = median(sims), remaining_forecast_p90 = quantile(sims, 0.90, names = FALSE)
  )
}

## ---- Required reduction in the REMAINING season ---------------------------
required_reduction <- function(m_t, threshold, remaining_forecast) {
  remaining_budget <- max(0, threshold - m_t)
  pct_reduction_needed <- if (remaining_forecast > 0) max(0, 1 - remaining_budget / remaining_forecast) else 0
  tibble::tibble(remaining_budget = remaining_budget, remaining_forecast = remaining_forecast,
                 required_abs_reduction = max(0, remaining_forecast - remaining_budget),
                 pct_reduction_needed = pct_reduction_needed)
}

## ---- Turbine-level allocation ---------------------------------------------
allocate_curtailment <- function(turbines, required_abs_reduction, effectiveness = 0.70) {
  turbines %>%
    mutate(
      avoided_if_curtailed = expected_remaining * effectiveness,
      cum_avoided = cumsum(avoided_if_curtailed),
      curtail = cum_avoided < required_abs_reduction | lag(cum_avoided, default = 0) < required_abs_reduction
    )
}

## ---- Illustrative synthetic demo (kept for reference / fallback) ---------
run_adaptive_management_demo <- function(fig_dir) {

  dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)
  fig_dir <- normalizePath(fig_dir, winslash = "/")

  n_weeks <- 52
  background_weight <- 0.3
  spring_weeks <- 14:26
  autumn_weeks <- 35:44
  raw_weight <- rep(background_weight, n_weeks)
  raw_weight[spring_weeks] <- 10
  raw_weight[autumn_weeks] <- 8
  seasonal_weight <- raw_weight / sum(raw_weight)

  weekly_fatalities_synthetic <- function(annual_threshold, curtailment_start_week = 31) {
    tibble::tibble(week = 1:n_weeks) %>%
      mutate(
        in_high_risk = week %in% c(spring_weeks, autumn_weeks),
        pre_curtailment_rate = ifelse(in_high_risk, annual_threshold * 0.9 / length(c(spring_weeks, autumn_weeks)), 0.05),
        post_curtailment_rate = ifelse(in_high_risk, 0.4, 0.05),
        rate = ifelse(week < curtailment_start_week, pre_curtailment_rate, post_curtailment_rate),
        observed_weekly = rpois(n(), rate)
      ) %>%
      mutate(observed_cum = cumsum(observed_weekly)) %>%
      select(week, observed_weekly, observed_cum)
  }

  project_thresholds <- tibble::tibble(project = facility_labels, threshold = pbr_thresholds$threshold)

  set.seed(1)
  dashboard_data <- lapply(seq_len(nrow(project_thresholds)), function(i) {
    proj <- project_thresholds$project[i]
    thr <- project_thresholds$threshold[i]
    build_control_trajectory(thr, seasonal_weight) %>%
      left_join(weekly_fatalities_synthetic(thr), by = "week") %>%
      mutate(project = proj, threshold = thr)
  }) %>%
    bind_rows() %>%
    mutate(trigger = classify_trigger(observed_cum, expected_cum, ucl_cum)) %>%
    left_join(trigger_response, by = "trigger")

  trigger_colours <- c(GREEN = "forestgreen", AMBER = "orange", RED = "firebrick")

  fig_adaptive_dashboard <- file.path(fig_dir, "adaptive_management_dashboard_demo.png")
  ggsave(fig_adaptive_dashboard, width = 11, height = 6.5, dpi = 150, plot = {
    ggplot(dashboard_data, aes(x = week)) +
      geom_ribbon(aes(ymin = lcl_cum, ymax = ucl_cum), fill = "grey70", alpha = 0.3) +
      geom_line(aes(y = expected_cum), linetype = "dashed", colour = "grey30", linewidth = 0.6) +
      geom_line(aes(y = observed_cum, colour = trigger, group = project), linewidth = 1) +
      geom_point(aes(y = observed_cum, colour = trigger), size = 1.2) +
      facet_wrap(~project, nrow = 1) +
      scale_colour_manual(name = "Trigger status", values = trigger_colours) +
      labs(
        x = "Week of year", y = "Cumulative fatalities",
        title = "Adaptive management control chart (illustrative synthetic data)",
        subtitle = paste0(
          "Dashed line: expected cumulative pace toward the annual threshold, given the seasonal risk profile;\n",
          "grey band: 90% Poisson control limits. NOT real PCFM counts -- demonstrates the mechanism only."
        )
      ) +
      theme_minimal() +
      theme(plot.subtitle = element_text(size = 10), strip.text = element_text(face = "bold"))
  })

  set.seed(3)
  checkpoint_week <- 22

  forecast_summary <- lapply(seq_len(nrow(project_thresholds)), function(i) {
    proj <- project_thresholds$project[i]; thr <- project_thresholds$threshold[i]
    obs <- dashboard_data %>% filter(project == proj) %>% arrange(week) %>% pull(observed_weekly)
    eoy <- project_eoy(obs, checkpoint_week, seasonal_weight, thr)
    red <- required_reduction(eoy$m_t, thr, eoy$remaining_forecast_median)
    bind_cols(tibble::tibble(project = proj, threshold = thr, week = checkpoint_week), eoy, red)
  }) %>% bind_rows()

  n_turbines <- 25
  set.seed(3)
  turbine_weights <- sort(rgamma(n_turbines, shape = 0.6, rate = 1), decreasing = TRUE)
  focus_project <- project_thresholds$project[1]
  remaining_forecast_focus <- forecast_summary$remaining_forecast_median[forecast_summary$project == focus_project]
  turbines <- tibble::tibble(
    turbine_id = sprintf("WTG_%02d", seq_len(n_turbines)),
    expected_remaining = turbine_weights / sum(turbine_weights) * remaining_forecast_focus
  ) %>%
    arrange(desc(expected_remaining)) %>%
    mutate(rank = row_number(), cum_expected = cumsum(expected_remaining),
           cum_pct = cum_expected / sum(expected_remaining))

  required_abs_focus <- forecast_summary$required_abs_reduction[forecast_summary$project == focus_project]
  turbines_allocated <- allocate_curtailment(turbines, required_abs_focus)
  n_curtailed <- sum(turbines_allocated$curtail)

  fig_turbine_pareto <- file.path(fig_dir, "turbine_pareto_allocation_demo.png")
  ggsave(fig_turbine_pareto, width = 9, height = 5.5, dpi = 150, plot = {
    ggplot(turbines_allocated, aes(x = reorder(turbine_id, -expected_remaining))) +
      geom_col(aes(y = expected_remaining, fill = curtail)) +
      geom_line(aes(y = cum_pct * max(expected_remaining), group = 1), colour = "grey20", linewidth = 0.6) +
      geom_point(aes(y = cum_pct * max(expected_remaining)), colour = "grey20", size = 1.2) +
      scale_y_continuous(
        name = "Expected remaining fatalities (turbine)",
        sec.axis = sec_axis(~ . / max(turbines_allocated$expected_remaining), name = "Cumulative % of remaining mortality", labels = scales::percent)
      ) +
      scale_fill_manual(name = "Curtailment\nprescribed", values = c(`TRUE` = "firebrick", `FALSE` = "grey70"),
                         labels = c(`TRUE` = "Yes", `FALSE` = "No")) +
      labs(
        x = "Turbine (ranked by expected remaining mortality)",
        title = paste0("Turbine-level allocation to meet the required reduction (", focus_project, ", illustrative)"),
        subtitle = "Synthetic data, demonstrates the allocation logic only."
      ) +
      theme_minimal() +
      theme(axis.text.x = element_text(angle = 90, vjust = 0.5, size = 7), plot.subtitle = element_text(size = 8.5))
  })

  list(
    demo_checkpoint_week = checkpoint_week,
    demo_forecast_summary = forecast_summary,
    demo_focus_project = focus_project,
    demo_n_turbines = n_turbines,
    demo_n_curtailed = n_curtailed,
    demo_required_abs_reduction = required_abs_focus,
    fig_adaptive_dashboard_demo = fig_adaptive_dashboard,
    fig_turbine_pareto_demo = fig_turbine_pareto
  )
}
