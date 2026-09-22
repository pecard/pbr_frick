##
## Synthetic weekly mortality data for the Paris Expert Workshop
## illustrative scenario -- generates fresh, fictional records rather
## than scaling the real BSH/DGY carcass log. Scaling real per-turbine,
## per-week counts by a fixed factor would keep the real pattern's exact
## fingerprint (which weeks spike, which turbines concentrate risk,
## down to non-integer artefacts of the scaling itself); simulating from
## a generative model avoids that entirely, at the cost of the numbers
## no longer being literally real -- the trade-off Paulo asked for.
##
## Output uses (period, week) instead of real calendar dates -- "Baseline"
## / "Response" per Paulo (2026-09), not real years or months. No Date
## object is ever constructed, so there is no real calendar to leak and
## no locale-dependent formatting risk (unlike R/load_real_mortality_data.R,
## which this deliberately does not reuse).
##
## Turbine-level concentration: risk weights follow a power law
## (weight_rank = rank^-0.71 after random ID shuffling), calibrated so
## that the same 50%-cumulative-mortality selection rule used in the real
## report lands at ~10 of 50 turbines organically (Paulo's target count),
## rather than hard-coding which turbines are "curtailed" -- the
## selection emerges from the same method the real report uses.
##

suppressPackageStartupMessages({ library(dplyr); library(tidyr) })

## ---- Seasonal risk shape (generic bat-migration pattern, not
## project-specific -- same shape already used illustratively in
## R/adaptive_management_dashboard.R) --------------------------------------
paris_seasonal_weight <- function(n_weeks = 52) {
  background_weight <- 0.3
  spring_weeks <- 14:26
  autumn_weeks <- 35:44
  raw <- rep(background_weight, n_weeks)
  raw[spring_weeks] <- 10
  raw[autumn_weeks] <- 8
  raw / sum(raw)
}

## ---- Turbine risk weights (power law, calibrated for ~50%/top-10-of-50) --
paris_turbine_weights <- function(n_turbines, pareto_b = 0.71, seed = 101) {
  set.seed(seed)
  rank_weight <- (seq_len(n_turbines))^(-pareto_b)
  shuffled <- sample(rank_weight)  # so "highest risk" isn't turbine 1 by construction
  turbine_id <- sprintf("T-%02d", seq_len(n_turbines))
  tibble::tibble(turbine = turbine_id, risk_weight = shuffled / sum(shuffled))
}

## ---- Main simulation -------------------------------------------------------
## baseline_total_raw: expected total raw carcass count over the Baseline
## period, summed across all turbines/weeks (Poisson noise around this).
## reduction_at_curtailed: multiplicative rate cut applied at curtailed
## turbines from response_start_week onward in the Response period
## (e.g. 0.25 = a 75% reduction in their weekly rate).
simulate_paris_mortality_data <- function(n_turbines = 50,
                                           n_curtailed_target = 10,
                                           baseline_total_raw = 60,
                                           response_start_week = 18,
                                           reduction_at_curtailed = 0.02,
                                           n_weeks = 52,
                                           seed = 515) {
  seasonal <- paris_seasonal_weight(n_weeks)
  turbines <- paris_turbine_weights(n_turbines)

  ## Baseline period: full year, no curtailment anywhere.
  set.seed(seed)
  baseline_grid <- tidyr::expand_grid(turbine = turbines$turbine, week = seq_len(n_weeks)) %>%
    left_join(turbines, by = "turbine") %>%
    mutate(
      seasonal_w = seasonal[week],
      lambda = baseline_total_raw * risk_weight * seasonal_w,
      raw_count = rpois(n(), lambda),
      period = "Baseline"
    )

  ## Curtailed set: rank turbines by Baseline total mortality, select the
  ## minimum set covering >=50% of cumulative mortality -- same rule the
  ## real report applies, not a hard-coded turbine list. With only ~50 raw
  ## records spread over 50 turbines, several turbines often tie on the
  ## same integer count right at the cutoff; arrange()'s stable sort would
  ## silently break the tie by turbine ID, curtailing some and not others
  ## among turbines with IDENTICAL observed mortality -- indefensible, and
  ## visible on the pareto plot as equal-height bars split red/grey.
  ## Fixed by never splitting a tied group: if the 50% cutoff falls inside
  ## a run of equal counts, the whole run is included.
  turbine_totals <- baseline_grid %>%
    group_by(turbine) %>%
    summarise(total = sum(raw_count), .groups = "drop") %>%
    left_join(turbines, by = "turbine") %>%
    arrange(desc(total), desc(risk_weight)) %>%
    mutate(cum_frac = cumsum(total) / sum(total), rank = row_number())
  n_raw_cutoff <- min(which(turbine_totals$cum_frac >= 0.5))
  boundary_value <- turbine_totals$total[n_raw_cutoff]
  n_selected <- max(which(turbine_totals$total == boundary_value))  # extend through the tied group
  curtailed_turbines <- turbine_totals$turbine[seq_len(n_selected)]

  ## Response period: same generative process, except curtailed turbines'
  ## rate is cut from response_start_week onward.
  response_grid <- tidyr::expand_grid(turbine = turbines$turbine, week = seq_len(n_weeks)) %>%
    left_join(turbines, by = "turbine") %>%
    mutate(
      seasonal_w = seasonal[week],
      is_curtailed = turbine %in% curtailed_turbines,
      is_active_curtailment = is_curtailed & week >= response_start_week,
      rate_multiplier = ifelse(is_active_curtailment, reduction_at_curtailed, 1),
      lambda = baseline_total_raw * risk_weight * seasonal_w * rate_multiplier,
      raw_count = rpois(n(), lambda),
      period = "Response"
    )

  weekly <- bind_rows(
    baseline_grid %>% select(turbine, period, week, raw_count),
    response_grid %>% select(turbine, period, week, raw_count)
  ) %>%
    mutate(project = facility_labels[1]) %>%
    relocate(project)

  list(
    weekly = weekly,
    turbine_weights = turbines,
    turbine_totals_baseline = turbine_totals,
    curtailed_turbines = curtailed_turbines,
    n_curtailed = length(curtailed_turbines),
    n_curtailed_target = n_curtailed_target,
    response_start_week = response_start_week,
    reduction_at_curtailed = reduction_at_curtailed
  )
}
