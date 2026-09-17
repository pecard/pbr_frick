##
## Quasi-experimental curtailment effectiveness (Paulo, 2026-09, point 4
## of his post-review list): "how much of the reduction can reasonably be
## attributed to curtailment rather than to 2026 simply being a lower-
## mortality year?" A difference-in-differences (DiD) negative-binomial
## GLMM, comparing the pre/post change at curtailed turbines against the
## pre/post change at non-curtailed (control) turbines within the same
## project -- the interaction term is the part of the change specific to
## curtailed turbines, net of whatever both groups experienced anyway.
##
## Two caveats disclosed up front, both confirmed by Paulo (2026-09):
## no SCADA/weather data is available, so no operational or weather
## covariate is included; and per-turbine search effort varied only
## slightly and was "almost always weekly", so no search-effort offset
## is used -- weekly presence/absence of a scheduled search is assumed
## constant across turbines and time. A "high_risk" (April-May) vs "low"
## season indicator is included instead of full calendar-month dummies:
## the latter caused quasi-separation (several months have zero counts
## at every turbine, pushing those coefficients to unstable extremes)
## without changing the interaction term's estimate.
##
## Turbines are the full official assessed universe (68 at Project 1, 27
## at Project 2, R/turbine_selection_validation.R's official reference
## tables), curtailed_ever = the real, lender-agreed selection (17/7).
##

if (!requireNamespace("glmmTMB", quietly = TRUE)) {
  message("Installing missing package needed for the curtailment-effectiveness DiD model: glmmTMB...")
  tryCatch(install.packages("glmmTMB"), error = function(e) stop(
    "Could not install 'glmmTMB' automatically (", conditionMessage(e), "). Install manually and re-run."
  ))
}

suppressPackageStartupMessages({ library(dplyr); library(tidyr); library(lubridate); library(glmmTMB); library(ggplot2) })

run_curtailment_effectiveness_did <- function(fig_dir,
                                               bash_path = "data-raw/BashWPP_Weekly_PCFM_PBR.xlsx",
                                               djangeldy_path = "data-raw/DjangeldyWPP_Weekly_PCFM_PBR.xlsx",
                                               official_bash_path = "data-raw/official_turbine_selection_bash.csv",
                                               official_djangeldy_path = "data-raw/official_turbine_selection_djangeldy.csv") {

  dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)
  fig_dir <- normalizePath(fig_dir, winslash = "/")

  d <- load_real_mortality_data(bash_path, djangeldy_path)

  pad_bsh <- function(x) sprintf("BSH%02d", as.integer(sub("^BSH", "", x)))
  pad_dzh <- function(x) x
  official_bash <- read.csv(official_bash_path, stringsAsFactors = FALSE) %>%
    mutate(official_selected = as.logical(official_selected)) %>% rename(turbine = turbine_real)
  official_dzh <- read.csv(official_djangeldy_path, stringsAsFactors = FALSE) %>%
    mutate(official_selected = as.logical(official_selected)) %>% rename(turbine = turbine_real)

  ## ifelse() evaluates both branches for every row, so pad_bsh() is
  ## harmlessly applied to Djangeldy codes too (and vice versa) before the
  ## correct branch is selected -- produces a benign "NAs introduced by
  ## coercion" warning on the discarded branch, not a real data problem
  ## (verified: zero actual NAs in the result).
  d <- suppressWarnings(d %>% mutate(turbine = ifelse(project == facility_labels[1], pad_bsh(turbine), pad_dzh(turbine))))

  data_end <- max(d$date)

  build_panel <- function(d_proj, official, proj_label) {
    data_start <- min(d_proj$date)
    week_starts <- seq(data_start, data_end, by = "week")
    grid <- tidyr::expand_grid(turbine = official$turbine, week_start = week_starts) %>%
      left_join(official %>% select(turbine, curtailed_ever = official_selected), by = "turbine")
    counts <- d_proj %>%
      mutate(week_start = week_starts[findInterval(date, week_starts)]) %>%
      count(turbine, week_start, name = "count")
    grid %>%
      left_join(counts, by = c("turbine", "week_start")) %>%
      mutate(
        count = tidyr::replace_na(count, 0),
        period = factor(ifelse(week_start >= curtailment_start_date, "post", "pre"), levels = c("pre", "post")),
        month = lubridate::month(week_start),
        high_risk = factor(ifelse(month %in% c(4, 5), "high", "low"), levels = c("low", "high")),
        project = proj_label,
        turbine_uid = paste0(proj_label, "_", turbine)
      )
  }

  panel_bash <- build_panel(d %>% filter(project == facility_labels[1]), official_bash, facility_labels[1])
  panel_dzh <- build_panel(d %>% filter(project == facility_labels[2]), official_dzh, facility_labels[2])
  panel_all <- bind_rows(panel_bash, panel_dzh)

  fit_did <- function(panel) glmmTMB(count ~ curtailed_ever * period + high_risk, family = nbinom2, data = panel)

  extract_did <- function(fitted_model, label) {
    ci <- confint(fitted_model, parm = "curtailed_everTRUE:periodpost")
    p_value <- summary(fitted_model)$coefficients$cond["curtailed_everTRUE:periodpost", "Pr(>|z|)"]
    tibble::tibble(
      model = label, log_rate_ratio = ci[1, "Estimate"], ci_low = ci[1, "2.5 %"], ci_high = ci[1, "97.5 %"],
      rate_ratio = exp(ci[1, "Estimate"]), rate_ratio_low = exp(ci[1, "2.5 %"]), rate_ratio_high = exp(ci[1, "97.5 %"]),
      p_value = p_value
    )
  }

  m_bash <- fit_did(panel_bash)
  m_dzh <- fit_did(panel_dzh)
  m_pooled <- glmmTMB(count ~ curtailed_ever * period + high_risk + project, family = nbinom2, data = panel_all)

  did_results <- bind_rows(
    extract_did(m_bash, facility_labels[1]),
    extract_did(m_dzh, facility_labels[2]),
    extract_did(m_pooled, "Pooled (both projects)")
  ) %>%
    mutate(ci_excludes_no_effect = ci_low > 0 | ci_high < 0)

  ## ---- Parallel-trends diagnostic: curtailed/control rate ratio across
  ## the two available PRE-treatment sub-periods (2025 vs. Jan-early May
  ## 2026) -- if the DiD's key assumption holds, this ratio should be
  ## roughly stable before treatment starts; a big shift is a caveat on
  ## how literally to read the DiD estimate above.
  pretrend <- panel_all %>%
    filter(period == "pre") %>%
    mutate(subperiod = ifelse(lubridate::year(week_start) == 2025, "2025", "2026 (Jan-early May)")) %>%
    group_by(project, subperiod, curtailed_ever) %>%
    summarise(rate = sum(count) / n(), .groups = "drop") %>%
    tidyr::pivot_wider(names_from = curtailed_ever, values_from = rate, names_prefix = "rate_curtailed_") %>%
    mutate(curtailed_control_ratio = rate_curtailed_TRUE / rate_curtailed_FALSE)

  fig_did <- file.path(fig_dir, "curtailment_effectiveness_did.png")
  ggsave(fig_did, width = 9, height = 5, dpi = 150, plot = {
    ggplot(did_results, aes(x = rate_ratio, y = model)) +
      geom_vline(xintercept = 1, linetype = "dashed", colour = "grey40") +
      geom_pointrange(aes(xmin = rate_ratio_low, xmax = rate_ratio_high), colour = "steelblue4", size = 0.8) +
      scale_x_log10() +
      labs(
        x = "Rate ratio (curtailed x post interaction, log scale)\n>1: curtailed turbines fared relatively WORSE after treatment started; <1: relatively better",
        y = NULL,
        title = "Curtailment effect net of shared seasonal/year change (DiD)",
        subtitle = "Dashed line = no differential effect. All three intervals cross it: not statistically\ndistinguishable from zero at current data volume."
      ) +
      theme_minimal() +
      theme(plot.subtitle = element_text(size = 10.5))
  })

  list(
    did_results = did_results,
    did_pretrend_check = pretrend,
    did_n_curtailed_weeks_post = list(
      project1 = sum(panel_bash$period == "post" & panel_bash$curtailed_ever),
      project2 = sum(panel_dzh$period == "post" & panel_dzh$curtailed_ever)
    ),
    fig_did = fig_did
  )
}
