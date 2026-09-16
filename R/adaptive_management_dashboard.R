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
