##
## PVA-lite robustness (Paulo, 2026-09, point 6 of his post-review list):
## the PVA-lite everywhere else in this note (Section 4, Section 7.5)
## fixes four assumptions that plausibly matter most for the conclusion,
## and Paulo asked for a targeted sensitivity on exactly those four,
## rather than a general re-derivation of the model:
##   1. N0 fixed at nmin_assumed rather than drawn from the same
##      project-footprint x density distribution used in Section 7.5;
##   2. mortality split evenly by sex/stage rather than concentrated in
##      adult females (the reproductive engine of the population -- a
##      standard conservation-biology asymmetry);
##   3. environmental stochasticity at a 10% vital-rate CV, a working
##      assumption with no inter-annual variance estimate located for
##      this species;
##   4. no density dependence -- the "no additional mortality" baseline
##      grows the population unchecked for 25 years, which is not
##      biologically credible over that horizon.
##
## Sensitivities 1 and 3 are directly computable from data already used
## elsewhere in this note. Sensitivities 2 and 4 are NOT: this project's
## carcass records have no sex field, and no independent carrying-
## capacity estimate exists for this species/region -- both are run as
## explicit, labelled BOUNDING scenarios (a worst-case assumption), not
## as data-calibrated corrections, and are reported as such throughout.
## Density dependence's ceiling K is set to the same distribution's own
## upper bound (project_footprint_km2 x density_proxy_high) used
## generatively in Section 7.5, rather than an unrelated round number --
## the most defensible ceiling available without a real abundance
## estimate, not a precise one.
##
## Each sensitivity changes ONE assumption at a time from the Section 4 /
## 7.5 baseline, at each project's own imposed threshold as the annual
## removal (the operationally relevant scenario, not a candidate level).
##

if (!requireNamespace("ggrepel", quietly = TRUE)) {
  message("Installing missing package needed for the PVA robustness figure: ggrepel...")
  tryCatch(install.packages("ggrepel"), error = function(e) stop(
    "Could not install 'ggrepel' automatically (", conditionMessage(e), "). Install manually and re-run."
  ))
}

suppressPackageStartupMessages({ library(dplyr); library(tidyr); library(ggplot2); library(ggrepel) })

run_pva_robustness <- function(fig_dir, nmin_assumed) {

  dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)
  fig_dir <- normalizePath(fig_dir, winslash = "/")

  build_dekker_stage_matrix <- function(n_stages, s_juv, s_adult, fecundity) {
    A <- matrix(0, n_stages, n_stages)
    A[1, n_stages] <- fecundity
    for (i in 1:(n_stages - 1)) A[i + 1, i] <- s_juv
    A[n_stages, n_stages] <- s_adult
    A
  }
  leslie_F <- leslie_p_breed * leslie_litter * leslie_sex_ratio * leslie_s_juv
  A_leslie_validated <- build_dekker_stage_matrix(2, leslie_s_juv, leslie_s_adult_f, leslie_F)
  stable_dist <- Re(eigen(A_leslie_validated)$vectors[, 1]); stable_dist <- stable_dist / sum(stable_dist)

  beta_params <- function(mean, cv) {
    var <- (mean * cv)^2
    common <- mean * (1 - mean) / var - 1
    list(shape1 = mean * common, shape2 = (1 - mean) * common)
  }

  nmin_ceiling_K <- density_proxy_high * project_footprint_km2

  ## female_bias: 0.5 = neutral (current baseline elsewhere), 1.0 = a
  ## bounding scenario where all removed individuals are adult females.
  ## density_dependence: NULL = none (baseline); a numeric K applies a
  ## logistic-style compensation to recruitment as N approaches K.
  simulate_trajectory <- function(n0_juv, n0_adult, n_years, annual_removal, cv, female_bias = 0.5, K = NULL) {
    bp_s_adult <- beta_params(leslie_s_adult_f, cv)
    bp_s_juv <- beta_params(leslie_s_juv, cv)
    bp_p_breed <- beta_params(leslie_p_breed, cv)

    n_juv <- as.integer(round(n0_juv)); n_adult <- as.integer(round(n0_adult))
    traj <- numeric(n_years + 1); traj[1] <- 2 * (n_juv + n_adult)
    for (yr in seq_len(n_years)) {
      s_adult_t <- rbeta(1, bp_s_adult$shape1, bp_s_adult$shape2)
      s_juv_t <- rbeta(1, bp_s_juv$shape1, bp_s_juv$shape2)
      p_breed_t <- rbeta(1, bp_p_breed$shape1, bp_p_breed$shape2)

      density_factor <- if (!is.null(K)) max(0, 1 - (2 * (n_juv + n_adult)) / K) else 1

      n_breeders <- rbinom(1, n_adult, p_breed_t)
      n_pups_total <- rpois(1, n_breeders * leslie_litter * density_factor)
      n_female_pups <- rbinom(1, n_pups_total, leslie_sex_ratio)
      n_juv_recruits <- rbinom(1, n_female_pups, s_juv_t)
      n_juv_survive <- rbinom(1, n_juv, s_juv_t)
      n_adult_survive <- rbinom(1, n_adult, s_adult_t)
      n_juv_next <- n_juv_recruits
      n_adult_next <- n_juv_survive + n_adult_survive

      ## female_bias = 0.5 (neutral): same removal logic as elsewhere in
      ## this note (half the removal is female, split across stages by
      ## relative abundance). female_bias = 1.0: ALL removal falls on
      ## adult females specifically (the bounding scenario).
      if (female_bias >= 0.999) {
        n_adult_next <- max(0, n_adult_next - annual_removal)
      } else {
        female_removal <- annual_removal * female_bias
        total_female <- n_juv_next + n_adult_next
        if (total_female > 0 && female_removal > 0) {
          juv_share <- round(female_removal * n_juv_next / total_female)
          adult_share <- round(female_removal * n_adult_next / total_female)
          n_juv_next <- max(0, n_juv_next - juv_share)
          n_adult_next <- max(0, n_adult_next - adult_share)
        }
      }
      n_juv <- n_juv_next; n_adult <- n_adult_next
      traj[yr + 1] <- 2 * (n_juv + n_adult)
    }
    traj
  }

  n0_female <- nmin_assumed / 2
  n0_juv_fixed <- round(n0_female * stable_dist[1]); n0_adult_fixed <- n0_female - n0_juv_fixed
  quasi_ext_threshold <- pva_quasi_extinction_fraction * nmin_assumed

  ## variable_N0 = TRUE: redraw N0 (and its juv/adult split) each
  ## replicate from the same project-footprint x density distribution
  ## used generatively in Section 7.5, instead of a single fixed value.
  run_risk <- function(annual_removal, cv = pva_vital_rate_cv, female_bias = 0.5, K = NULL, variable_N0 = FALSE) {
    set.seed(pva_seed)
    reps <- replicate(pva_n_reps, {
      if (variable_N0) {
        nmin_draw <- runif(1, density_proxy_low, density_proxy_high) * project_footprint_km2
        nf <- nmin_draw / 2
        nj <- round(nf * stable_dist[1]); na <- nf - nj
      } else {
        nj <- n0_juv_fixed; na <- n0_adult_fixed
      }
      simulate_trajectory(nj, na, pva_n_years, annual_removal, cv, female_bias, K)
    })
    final_N <- reps[pva_n_years + 1, ]
    ever_below <- apply(reps, 2, function(traj) any(traj < quasi_ext_threshold))
    tibble::tibble(
      p_decline = mean(final_N < nmin_assumed) * 100,
      p_quasi_extinction = mean(ever_below) * 100
    )
  }

  project_thresholds <- tibble::tibble(project = facility_labels, threshold = pbr_thresholds$threshold)

  sensitivities <- tibble::tibble(
    sensitivity = c(
      "Baseline (Section 4/7.5 assumptions)",
      "Variable N0 (area x density, not fixed)",
      "Female-biased removal (bounding scenario: 100% adult female)",
      "Higher environmental stochasticity (25% vital-rate CV)",
      "Density dependence (bounding scenario: logistic ceiling)"
    ),
    is_bounding_scenario = c(FALSE, FALSE, TRUE, FALSE, TRUE)
  )

  pva_robustness_table <- tidyr::expand_grid(project_thresholds, sensitivities) %>%
    rowwise() %>%
    mutate(
      risk = list(switch(sensitivity,
        "Baseline (Section 4/7.5 assumptions)" = run_risk(threshold),
        "Variable N0 (area x density, not fixed)" = run_risk(threshold, variable_N0 = TRUE),
        "Female-biased removal (bounding scenario: 100% adult female)" = run_risk(threshold, female_bias = 1.0),
        "Higher environmental stochasticity (25% vital-rate CV)" = run_risk(threshold, cv = 0.25),
        "Density dependence (bounding scenario: logistic ceiling)" = run_risk(threshold, K = nmin_ceiling_K)
      ))
    ) %>%
    ungroup() %>%
    tidyr::unnest_wider(risk)

  fig_pva_robustness <- file.path(fig_dir, "pva_robustness.png")
  ggsave(fig_pva_robustness, width = 10.5, height = 5.5, dpi = 150, plot = {
    ggplot(pva_robustness_table, aes(x = reorder(sensitivity, p_decline), y = p_decline, fill = is_bounding_scenario)) +
      geom_col() +
      geom_text(aes(label = sprintf("%.0f%%", p_decline)), hjust = -0.15, size = 3) +
      coord_flip(clip = "off") +
      facet_wrap(~project) +
      scale_y_continuous(limits = c(0, 105), expand = expansion(mult = c(0, 0.15))) +
      scale_fill_manual(name = NULL, values = c(`FALSE` = "steelblue", `TRUE` = "firebrick"),
                         labels = c(`FALSE` = "Directly computable", `TRUE` = "Bounding scenario (no data to calibrate)")) +
      labs(
        x = NULL, y = "25-year PVA-lite: P(decline below N0), %",
        title = "PVA-lite robustness: one assumption changed at a time",
        subtitle = paste0(
          "Annual removal = each project's own imposed threshold. Blue: computable from data used elsewhere in this note.\n",
          "Red: an explicit worst-case bound, not a calibrated estimate (no sex or abundance data available)."
        )
      ) +
      theme_minimal() +
      theme(plot.subtitle = element_text(size = 10.5), strip.text = element_text(face = "bold"), legend.position = "bottom")
  })

  list(
    pva_robustness_table = pva_robustness_table,
    pva_robustness_K = nmin_ceiling_K,
    fig_pva_robustness = fig_pva_robustness
  )
}
