##
## Illustrative Section 7 (operational adaptive management) for the Paris
## Expert Workshop scenario -- single fictional wind farm, "Bat model"
## species, entirely simulated weekly mortality (R/simulate_paris_mortality_data.R).
## Mirrors the real report's methodology (turbine concentration -> 50%-
## cumulative selection rule; actual-vs-counterfactual effectiveness) but
## on generic (period, week) time instead of real calendar dates.
##

suppressPackageStartupMessages({ library(dplyr); library(ggplot2) })

run_paris_operational_analysis <- function(fig_dir = "outputs/figures", sim = NULL) {
  if (is.null(sim)) sim <- simulate_paris_mortality_data()
  dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)
  fig_dir <- normalizePath(fig_dir, winslash = "/")

  factor <- genest_correction_factor_paris

  ## ---- Turbine concentration (Baseline period) -----------------------------
  turbine_hist <- sim$turbine_totals_baseline %>%
    transmute(turbine, raw = total, corrected = raw * factor) %>%
    arrange(desc(corrected)) %>%
    mutate(rank = row_number(), cum_pct = cumsum(corrected) / sum(corrected),
           turbine_anon = sprintf("T-%02d", rank),
           curtailed = turbine %in% sim$curtailed_turbines)
  n_selected <- sum(turbine_hist$curtailed)

  fig_turbine_pareto <- file.path(fig_dir, "paris_turbine_pareto.png")
  ggsave(fig_turbine_pareto, width = 9, height = 5, dpi = 150, plot = {
    ggplot(turbine_hist, aes(x = reorder(turbine_anon, rank))) +
      geom_col(aes(y = corrected, fill = curtailed), alpha = 0.85) +
      geom_line(aes(y = cum_pct * max(corrected), group = 1), colour = "grey20", linewidth = 0.6) +
      geom_point(aes(y = cum_pct * max(corrected)), colour = "grey20", size = 1) +
      geom_hline(yintercept = 0.5 * max(turbine_hist$corrected), linetype = "dotted", colour = "steelblue") +
      scale_fill_manual(values = c(`TRUE` = "firebrick", `FALSE` = "grey70"), name = "Curtailed", labels = c("No", "Yes")) +
      scale_y_continuous(
        name = "Corrected fatalities (Baseline period)",
        sec.axis = sec_axis(~ . / max(turbine_hist$corrected), name = "Cumulative % of mortality", labels = scales::percent)
      ) +
      labs(
        x = "Turbine (ranked by Baseline corrected mortality)",
        title = "Where Baseline mortality concentrated",
        subtitle = paste0(
          n_selected, " of ", nrow(turbine_hist), " turbines account for ",
          round(100 * max(turbine_hist$cum_pct[turbine_hist$curtailed])), "% of Baseline mortality ",
          "(minimum set covering >=50%)."
        )
      ) +
      theme_minimal() +
      theme(axis.text.x = element_text(angle = 90, vjust = 0.5, size = 7), plot.subtitle = element_text(size = 10))
  })

  ## ---- Effectiveness: actual (Response) vs Baseline-pattern counterfactual --
  weekly_cum <- sim$weekly %>%
    group_by(period, week) %>%
    summarise(raw = sum(raw_count), .groups = "drop") %>%
    arrange(period, week) %>%
    group_by(period) %>%
    mutate(cum_raw = cumsum(raw), cum_corrected = cum_raw * factor) %>%
    ungroup()

  baseline_cum <- weekly_cum %>% filter(period == "Baseline") %>% select(week, counterfactual_cum = cum_corrected)
  response_cum <- weekly_cum %>% filter(period == "Response") %>% select(week, actual_cum = cum_corrected)
  effectiveness <- full_join(baseline_cum, response_cum, by = "week") %>% arrange(week)

  fig_effectiveness <- file.path(fig_dir, "paris_effectiveness.png")
  ggsave(fig_effectiveness, width = 8, height = 4.8, dpi = 150, plot = {
    ggplot(effectiveness, aes(x = week)) +
      geom_line(aes(y = counterfactual_cum, linetype = "No-response counterfactual (Baseline pattern)"), colour = "firebrick", linewidth = 0.9) +
      geom_line(aes(y = actual_cum, linetype = "Actual (Response, with curtailment)"), colour = "forestgreen", linewidth = 1.1) +
      geom_hline(yintercept = pbr_thresholds$threshold, linetype = "dotted", colour = "grey30") +
      geom_vline(xintercept = sim$response_start_week, linetype = "dotted", colour = "steelblue") +
      scale_linetype_manual(name = NULL, values = c("Actual (Response, with curtailment)" = "solid",
                                                      "No-response counterfactual (Baseline pattern)" = "dashed")) +
      labs(
        x = "Week", y = "Cumulative corrected fatalities (Bat model)",
        title = "Response period: actual vs. a no-response counterfactual",
        subtitle = "Dotted horizontal: annual reference. Dotted vertical: curtailment start."
      ) +
      theme_minimal() +
      theme(plot.subtitle = element_text(size = 10), legend.position = "bottom")
  })

  raw_baseline_total <- sum(sim$weekly$raw_count[sim$weekly$period == "Baseline"])
  raw_response_total <- sum(sim$weekly$raw_count[sim$weekly$period == "Response"])
  actual_post <- sum(sim$weekly$raw_count[sim$weekly$period == "Response" & sim$weekly$week >= sim$response_start_week])
  cf_post <- sum(sim$weekly$raw_count[sim$weekly$period == "Baseline" & sim$weekly$week >= sim$response_start_week])

  summary_tbl <- tibble::tibble(
    metric = c("Baseline corrected total", "% of annual reference (Baseline)",
               "Response corrected total (full period)", "Actual corrected (from response start)",
               "Counterfactual corrected (from response start)", "Reduction vs. counterfactual",
               "Turbines curtailed"),
    value = c(
      round(raw_baseline_total * factor), paste0(round(100 * raw_baseline_total * factor / pbr_thresholds$threshold), "%"),
      round(raw_response_total * factor), round(actual_post * factor), round(cf_post * factor),
      paste0(round(100 * (1 - actual_post / cf_post)), "%"),
      paste0(n_selected, " of ", nrow(turbine_hist))
    )
  )

  list(
    turbine_hist = turbine_hist, effectiveness = effectiveness, summary_tbl = summary_tbl,
    fig_turbine_pareto = fig_turbine_pareto, fig_effectiveness = fig_effectiveness
  )
}
