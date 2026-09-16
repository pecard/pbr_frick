##
## Adaptive management control chart: turns the annual PBR-derived
## threshold into a WEEKLY monitoring rule that can actually be acted on
## during the season, instead of a single pass/fail check at year-end.
##
## Core idea (standard statistical process control / CUSUM logic, applied
## to fatality monitoring): mortality risk for this species is not spread
## evenly across the year -- it concentrates in the migration windows
## (Apr-Jun, Sep-Oct, per Paulo's account). So the annual threshold has to
## be converted into an EXPECTED CUMULATIVE TRAJECTORY that follows the
## same seasonal shape, with statistical control limits around it. Compare
## the actual cumulative PCFM count against that trajectory every week,
## not against the flat annual number -- this is what lets curtailment be
## adjusted mid-season, while there is still a season left to act in.
##
## ILLUSTRATIVE ONLY: the seasonal risk profile and the weekly fatality
## series below are synthetic, built only to demonstrate the mechanism and
## match the qualitative pattern Paulo described (pre-curtailment: dozens
## /month in the high-risk windows; post-curtailment: 1-2/month). Replace
## weekly_fatalities_synthetic with the real weekly PCFM series once
## available -- everything downstream (control limits, trigger status,
## the figure) recomputes automatically from that one input.
##

suppressPackageStartupMessages({ library(dplyr); library(tidyr); library(ggplot2) })

fig_dir <- "outputs/figures"
dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)

# ---- 1. Seasonal risk profile (illustrative) -------------------------------
# Weight per ISO week (1-52): background risk everywhere, elevated in the
# two migration windows. Normalised to sum to 1 over the year.
n_weeks <- 52
background_weight <- 0.3
spring_weeks <- 14:26   # ~ April-June
autumn_weeks <- 35:44   # ~ September-October
raw_weight <- rep(background_weight, n_weeks)
raw_weight[spring_weeks] <- 10
raw_weight[autumn_weeks] <- 8
seasonal_weight <- raw_weight / sum(raw_weight)

# ---- 2. Expected cumulative trajectory + Poisson control limits -----------
# For an annual threshold T, the expected cumulative count by week t is
# C_t = T * cumsum(seasonal_weight)[t]. If mortality were exactly on the
# threshold's own pace, the observed cumulative count by week t is
# approximately Poisson(C_t) -- gives statistically grounded, not
# arbitrary, control limits.
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

# ---- 3. Trigger tiers -------------------------------------------------------
# GREEN: on or under pace. AMBER: ahead of pace but within statistical
# noise of the threshold's own trajectory (90% Poisson band). RED:
# cumulative count statistically exceeds what the annual threshold's own
# pace would produce -- not just "a bad week", a genuine early-warning
# signal with a full season's worth of curtailment response time left.
classify_trigger <- function(observed_cum, expected_cum, ucl_cum) {
  dplyr::case_when(
    observed_cum <= expected_cum ~ "GREEN",
    observed_cum <= ucl_cum ~ "AMBER",
    TRUE ~ "RED"
  )
}

# Curtailment response tied to each tier -- an escalation ladder, not a
# continuously optimised dial: with only one before/after curtailment data
# point so far (critical turbines on vs off), a smooth dose-response curve
# (fatality rate vs. cut-in speed or % capacity curtailed) cannot yet be
# fitted. That becomes possible once curtailment intensity itself has been
# varied (not just on/off) across enough weeks to estimate a response --
# see note at the end of this script.
trigger_response <- tibble::tibble(
  trigger = c("GREEN", "AMBER", "RED"),
  curtailment_action = c(
    "Maintain current prescription on critical turbines; scheduled review at next checkpoint.",
    "Increase curtailment intensity on critical turbines (lower cut-in speed by one defined step, or extend curtailed hours); reassess weekly.",
    "Escalate: widen curtailment to the next turbine tier and/or move to the most conservative cut-in speed on critical turbines for the remainder of the high-risk window; trigger the structured review (magnitude, persistence, PCFM uncertainty) from the report's closing section; document for lenders."
  )
)

# ---- 4. Synthetic illustrative weekly series -------------------------------
# Built only to demonstrate the mechanism end-to-end: elevated mortality
# in the migration windows before curtailment (weeks 1-30), then curtailed
# from week 31 (illustrative "curtailment start"), matching Paulo's
# qualitative pre/post description. NOT real PCFM data.
set.seed(1)
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

# ---- 5. Assemble the dashboard for each project ----------------------------
project_thresholds <- tibble::tibble(project = c("Project 1", "Project 2"), threshold = c(144, 120))

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

current_status <- dashboard_data %>%
  group_by(project) %>%
  filter(week == max(week)) %>%
  ungroup() %>%
  select(project, week, observed_cum, expected_cum, ucl_cum, trigger, curtailment_action)

cat("=== Current trigger status (illustrative synthetic data, last week in series) ===\n")
print(as.data.frame(current_status))

# ---- 6. Dashboard figure ----------------------------------------------------
trigger_colours <- c(GREEN = "forestgreen", AMBER = "orange", RED = "firebrick")

fig_adaptive_dashboard <- file.path(fig_dir, "adaptive_management_dashboard.png")
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
cat("\nWrote", fig_adaptive_dashboard, "\n")

# ---- 7. Note on moving from an escalation ladder to a fitted dose-response -
cat(
  "\nNote: the curtailment_action column above is a fixed escalation ladder\n",
  "(GREEN/AMBER/RED -> a bigger or smaller curtailment step), not a continuously\n",
  "optimised rule, because only one on/off curtailment data point exists so far\n",
  "(critical turbines curtailed vs not). Once curtailment INTENSITY itself has\n",
  "been varied across enough weeks (e.g. different cut-in speed thresholds on\n",
  "different nights, a standard experimental design for curtailment studies),\n",
  "a real dose-response curve (fatality rate ~ cut-in speed, or ~ % capacity\n",
  "curtailed) could be fitted directly from the weekly PCFM series, turning the\n",
  "ladder into a numeric rule: 'to bring the trajectory from RED back to GREEN,\n",
  "raise cut-in speed by X m/s' -- not yet possible with only on/off data.\n",
  sep = ""
)

##
## ---- PART 2: MRI / EOY forecast / exceedance probability / turbine ------
##      allocation (Paulo, 2026-09) -- upgrades the tier-based control
##      chart above with the richer layer he proposed: a Mortality
##      Reference Index (observed / expected-to-date, not just a raw
##      count), a stochastic end-of-year forecast under the CURRENT
##      regime (not the annual threshold's own pace), the resulting
##      exceedance probability, the reduction still required in the
##      remaining season specifically (not the whole year -- fatalities
##      already incurred cannot be recovered), and a turbine-level
##      Pareto/allocation step translating that reduction into which
##      turbines would need tighter curtailment.
##
## STILL ILLUSTRATIVE: uses the same synthetic weekly series as Part 1,
## and a synthetic turbine-level fatality distribution loosely following
## the "mortality is spatially concentrated" pattern. The forecast
## uncertainty (M_50/M_90 for the SEASON REMAINING) is genuinely computed
## by forward simulation, not fabricated -- but the cumulative-to-date
## estimate itself is still treated as a fixed point estimate here,
## because no real PCFM/GenEst detection-probability uncertainty exists
## in this session to draw from. Swap in real weekly (ideally per-turbine)
## PCFM output and this becomes a real tool, not a demonstration.
##

# ---- MRI: observed-to-date vs. expected-to-date under the annual reference
mri_t <- dashboard_data %>%
  mutate(MRI = observed_cum / expected_cum) %>%
  select(project, threshold, week, observed_cum, expected_cum, MRI)

# ---- Stochastic end-of-year forecast under the CURRENT regime -------------
# Estimates the recent realised rate relative to what the threshold's own
# pace would imply over a lookback window, then projects the remaining
# season forward stochastically at that same relative intensity -- i.e.
# "if the current mitigation regime continues unchanged, where do we
# finish, and how uncertain is that?" This is a forecast under status quo,
# not a claim about the true underlying mortality process.
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

# ---- Required reduction in the REMAINING season, not the whole year -------
required_reduction <- function(m_t, threshold, remaining_forecast) {
  remaining_budget <- max(0, threshold - m_t)
  pct_reduction_needed <- if (remaining_forecast > 0) max(0, 1 - remaining_budget / remaining_forecast) else 0
  tibble::tibble(remaining_budget = remaining_budget, remaining_forecast = remaining_forecast,
                 required_abs_reduction = max(0, remaining_forecast - remaining_budget),
                 pct_reduction_needed = pct_reduction_needed)
}

set.seed(3)
checkpoint_week <- 22  # illustrative mid-spring checkpoint, well before curtailment starts in
                        # the synthetic series -- chosen to demonstrate the case that matters
                        # (catching an off-track trajectory WHILE there is still season left to
                        # act in), not the already-controlled end-of-season state

forecast_summary <- lapply(seq_len(nrow(project_thresholds)), function(i) {
  proj <- project_thresholds$project[i]; thr <- project_thresholds$threshold[i]
  obs <- dashboard_data %>% filter(project == proj) %>% arrange(week) %>% pull(observed_weekly)
  eoy <- project_eoy(obs, checkpoint_week, seasonal_weight, thr)
  red <- required_reduction(eoy$m_t, thr, eoy$remaining_forecast_median)
  bind_cols(tibble::tibble(project = proj, threshold = thr, week = checkpoint_week), eoy, red)
}) %>% bind_rows()

cat("\n=== Six-number dashboard header, illustrative checkpoint at week", checkpoint_week, "===\n")
print(as.data.frame(forecast_summary %>%
  select(project, threshold, m_t, eoy_median, eoy_p90, p_exceed, required_abs_reduction, pct_reduction_needed)))

# ---- Turbine-level Pareto + allocation (Project 1 only, illustrative) -----
n_turbines <- 25
turbine_weights <- sort(rgamma(n_turbines, shape = 0.6, rate = 1), decreasing = TRUE)
remaining_forecast_p1 <- forecast_summary$remaining_forecast_median[forecast_summary$project == "Project 1"]
turbines <- tibble::tibble(
  turbine_id = sprintf("WTG_%02d", seq_len(n_turbines)),
  expected_remaining = turbine_weights / sum(turbine_weights) * remaining_forecast_p1
) %>%
  arrange(desc(expected_remaining)) %>%
  mutate(rank = row_number(), cum_expected = cumsum(expected_remaining),
         cum_pct = cum_expected / sum(expected_remaining))

allocate_curtailment <- function(turbines, required_abs_reduction, effectiveness = 0.70) {
  turbines %>%
    mutate(
      avoided_if_curtailed = expected_remaining * effectiveness,
      cum_avoided = cumsum(avoided_if_curtailed),
      curtail = cum_avoided < required_abs_reduction | lag(cum_avoided, default = 0) < required_abs_reduction
    )
}

required_abs_p1 <- forecast_summary$required_abs_reduction[forecast_summary$project == "Project 1"]
turbines_allocated <- allocate_curtailment(turbines, required_abs_p1)
n_curtailed <- sum(turbines_allocated$curtail)

cat(sprintf(
  "\nProject 1: required reduction in remaining-season mortality = %.1f fatalities.\n",
  required_abs_p1
))
cat(sprintf(
  "Illustrative allocation (70%% assumed curtailment effectiveness): curtailing the top %d of %d turbines\n",
  n_curtailed, n_turbines
))
cat(sprintf("(%.0f%% of turbines, accounting for %.0f%% of the project's remaining expected mortality) meets it.\n",
            100 * n_curtailed / n_turbines, 100 * turbines_allocated$cum_pct[n_curtailed]))

fig_turbine_pareto <- file.path(fig_dir, "turbine_pareto_allocation.png")
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
      title = "Turbine-level allocation to meet the required reduction (Project 1, illustrative)",
      subtitle = paste0(
        "Bars: expected remaining fatalities per turbine; line: cumulative %. Red bars: turbines selected to meet the ",
        round(required_abs_p1, 0), "-fatality reduction needed\n(70% assumed curtailment effectiveness) -- synthetic data, demonstrates the allocation logic only."
      )
    ) +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, size = 7), plot.subtitle = element_text(size = 8.5))
})
cat("\nWrote", fig_turbine_pareto, "\n")
