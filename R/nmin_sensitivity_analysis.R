##
## Nmin sensitivity: how much does doubling the assumed local population
## (Nmin = 6044 instead of 3022, e.g. if the region hosts a migratory
## passage roughly doubling the resident count) change each of the 3
## approaches already explored (PBR/Niel-Lebreton, the validated Leslie
## matrix, PVA-lite)? Nmin = 3022 is not a literature figure -- it is the
## median-fallback value the four threshold x lambda_max combinations
## back-solve to under the analysis's own original-IUCN-guess assumption
## (see R/pbr_analysis.R, "Reverse-engineering the imposed thresholds");
## this script treats it as one scenario among others, not a fixed truth.
##

suppressPackageStartupMessages({ library(dplyr); library(tidyr); library(ggplot2) })

source("R/pbr_functions.R")
source("R/leslie_dekker_limpens.R")     # S_j, S_af, p_breed, litter, sex_ratio, F_recruitment,
                                         # build_dekker_stage_matrix(), leslie_lambda(), lambda_female
source("inputs/pbrSettings_BSH_DGY.R")  # s_range, alpha_range, fr_corrected inputs, pbr_thresholds, mc_*, leslie_*, pva_*

fr_corrected <- fr_from_iucn(iucn_status_corrected)
nmin_scenarios <- c(3022, 6044)
fig_dir <- "outputs/figures"
dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)

# ---- 1. PBR / Niel-Lebreton -----------------------------------------------
# Fixed-benchmark PBR at each Nmin (Fr corrected = 0.5; lambda_max = 1.20/1.24).
pbr_fixed_by_nmin <- tidyr::expand_grid(nmin = nmin_scenarios, lambda_max_fixed = lambda_max_benchmarks) %>%
  mutate(PBR = pbr_from_components(nmin, fr_corrected, lambda_max_fixed))

cat("=== PBR at fixed lambda_max benchmarks, by Nmin (Fr = ", fr_corrected, ") ===\n", sep = "")
print(pbr_fixed_by_nmin)

# Monte Carlo: same s/alpha ranges as the main analysis, only Nmin changes --
# isolates the effect of Nmin alone on where the imposed thresholds sit.
mc_for_nmin <- function(nmin) {
  set.seed(mc_seed)
  sim <- tibble::tibble(
    s = runif(mc_n_sim, s_range[1], s_range[2]),
    alpha = runif(mc_n_sim, alpha_range[1], alpha_range[2])
  ) %>%
    mutate(lambda_max = lambda_max_niel(s, alpha), PBR = pbr_from_components(nmin, fr_corrected, lambda_max))
  list(
    sim = sim,
    quantiles = quantile(sim$PBR, c(0.05, 0.25, 0.50, 0.75, 0.95)),
    threshold_percentiles = pbr_thresholds %>%
      mutate(percentile_rank = sapply(threshold, function(x) mean(sim$PBR < x) * 100))
  )
}
mc_results <- setNames(lapply(nmin_scenarios, mc_for_nmin), paste0("nmin_", nmin_scenarios))

cat("\n=== Monte Carlo PBR quantiles, by Nmin ===\n")
for (nm in names(mc_results)) { cat(nm, ":\n"); print(round(mc_results[[nm]]$quantiles, 1)) }

cat("\n=== Where the imposed thresholds sit (percentile of the simulated PBR distribution), by Nmin ===\n")
for (nm in names(mc_results)) { cat(nm, ":\n"); print(mc_results[[nm]]$threshold_percentiles) }

# ---- 2. Leslie matrix: implied maximum sustainable yield -------------------
# The Leslie matrix has no Nmin input (lambda is a rate, density-independent
# in this structure) -- but its OWN lambda (empirically parametrised, not a
# theoretical max) implies a deterministic surplus-production-style ceiling
# on removals: yield ~ (lambda - 1) * N, the population growth in absolute
# numbers each year at that lambda. This is a stricter, more conservative
# comparator than PBR's Rmax*Nmin term (lambda_observed < lambda_max here),
# not a substitute for PBR's own logic -- included for direct comparability
# with the two other methods at the same Nmin values.
leslie_yield_by_nmin <- tibble::tibble(nmin = nmin_scenarios) %>%
  mutate(lambda_leslie = lambda_female, sustainable_yield_leslie = (lambda_female - 1) * nmin)

cat("\n=== Leslie-matrix-implied sustainable yield (lambda_observed - 1) * N, by Nmin ===\n")
print(leslie_yield_by_nmin)

# Imposed thresholds as a fraction of Nmin, for context.
threshold_fraction_by_nmin <- tidyr::expand_grid(pbr_thresholds, nmin = nmin_scenarios) %>%
  mutate(pct_of_nmin = 100 * threshold / nmin)
cat("\n=== Imposed thresholds as a percentage of Nmin ===\n")
print(threshold_fraction_by_nmin)

# ---- 3. PVA-lite at each Nmin ----------------------------------------------
leslie_F <- leslie_p_breed * leslie_litter * leslie_sex_ratio * leslie_s_juv
A_leslie_validated <- build_dekker_stage_matrix(2, leslie_s_juv, leslie_s_adult_f, leslie_F)
stable_dist <- Re(eigen(A_leslie_validated)$vectors[, 1]); stable_dist <- stable_dist / sum(stable_dist)

beta_params <- function(mean, cv) {
  var <- (mean * cv)^2
  common <- mean * (1 - mean) / var - 1
  list(shape1 = mean * common, shape2 = (1 - mean) * common)
}
bp_s_adult <- beta_params(leslie_s_adult_f, pva_vital_rate_cv)
bp_s_juv <- beta_params(leslie_s_juv, pva_vital_rate_cv)
bp_p_breed <- beta_params(leslie_p_breed, pva_vital_rate_cv)

simulate_pva_trajectory <- function(n0_juv, n0_adult, n_years, annual_removal) {
  n_juv <- as.integer(round(n0_juv)); n_adult <- as.integer(round(n0_adult))
  traj <- numeric(n_years + 1); traj[1] <- 2 * (n_juv + n_adult)
  for (yr in seq_len(n_years)) {
    s_adult_t <- rbeta(1, bp_s_adult$shape1, bp_s_adult$shape2)
    s_juv_t <- rbeta(1, bp_s_juv$shape1, bp_s_juv$shape2)
    p_breed_t <- rbeta(1, bp_p_breed$shape1, bp_p_breed$shape2)

    n_breeders <- rbinom(1, n_adult, p_breed_t)
    n_pups_total <- rpois(1, n_breeders * leslie_litter)
    n_female_pups <- rbinom(1, n_pups_total, leslie_sex_ratio)
    n_juv_recruits <- rbinom(1, n_female_pups, s_juv_t)

    n_juv_survive <- rbinom(1, n_juv, s_juv_t)
    n_adult_survive <- rbinom(1, n_adult, s_adult_t)

    n_juv_next <- n_juv_recruits
    n_adult_next <- n_juv_survive + n_adult_survive

    female_removal <- annual_removal / 2
    total_female <- n_juv_next + n_adult_next
    if (total_female > 0 && female_removal > 0) {
      juv_share <- round(female_removal * n_juv_next / total_female)
      adult_share <- round(female_removal * n_adult_next / total_female)
      n_juv_next <- max(0, n_juv_next - juv_share)
      n_adult_next <- max(0, n_adult_next - adult_share)
    }
    n_juv <- n_juv_next; n_adult <- n_adult_next
    traj[yr + 1] <- 2 * (n_juv + n_adult)
  }
  traj
}

run_pva_for_n0 <- function(n0_total) {
  n0_female <- n0_total / 2
  n0_juv <- round(n0_female * stable_dist[1]); n0_adult <- n0_female - n0_juv

  scenario_labels <- sprintf("%s threshold (%s/yr)", pbr_thresholds$facility, pbr_thresholds$threshold)
  scenarios <- bind_rows(
    tibble::tibble(scenario = "No additional mortality", threshold = 0),
    tibble::tibble(scenario = scenario_labels, threshold = pbr_thresholds$threshold)
  ) %>% mutate(scenario = factor(scenario, levels = c("No additional mortality", scenario_labels)))

  set.seed(pva_seed)
  trajectories <- lapply(seq_len(nrow(scenarios)), function(i) {
    reps <- replicate(pva_n_reps, simulate_pva_trajectory(n0_juv, n0_adult, pva_n_years, scenarios$threshold[i]))
    tibble::tibble(
      scenario = scenarios$scenario[i], year = rep(0:pva_n_years, pva_n_reps),
      rep = rep(seq_len(pva_n_reps), each = pva_n_years + 1), N = as.vector(reps)
    )
  }) %>% bind_rows()

  quasi_ext_threshold <- pva_quasi_extinction_fraction * n0_total
  risk_summary <- trajectories %>%
    group_by(scenario, rep) %>%
    summarise(final_N = N[year == pva_n_years], ever_below_threshold = any(N < quasi_ext_threshold), .groups = "drop") %>%
    group_by(scenario) %>%
    summarise(
      median_final_N = median(final_N), q05_final_N = quantile(final_N, 0.05), q95_final_N = quantile(final_N, 0.95),
      p_decline = mean(final_N < n0_total) * 100, p_quasi_extinction = mean(ever_below_threshold) * 100, .groups = "drop"
    )
  list(n0_total = n0_total, trajectories = trajectories, risk_summary = risk_summary, quasi_ext_threshold = quasi_ext_threshold)
}

pva_by_nmin <- setNames(lapply(nmin_scenarios, run_pva_for_n0), paste0("nmin_", nmin_scenarios))

cat("\n=== PVA-lite risk summary, by Nmin (N0) ===\n")
for (nm in names(pva_by_nmin)) { cat(nm, ":\n"); print(pva_by_nmin[[nm]]$risk_summary) }

# ---- 4. Combined comparison table ------------------------------------------
comparison_by_nmin <- tibble::tibble(nmin = nmin_scenarios) %>%
  rowwise() %>%
  mutate(
    pbr_range = sprintf(
      "%.0f - %.0f",
      pbr_fixed_by_nmin$PBR[pbr_fixed_by_nmin$nmin == nmin & pbr_fixed_by_nmin$lambda_max_fixed == 1.20],
      pbr_fixed_by_nmin$PBR[pbr_fixed_by_nmin$nmin == nmin & pbr_fixed_by_nmin$lambda_max_fixed == 1.24]
    ),
    leslie_sustainable_yield = round((lambda_female - 1) * nmin),
    pva_p_decline_range = {
      rs <- pva_by_nmin[[paste0("nmin_", nmin)]]$risk_summary
      rs <- rs[rs$scenario != "No additional mortality", ]
      sprintf("%.0f - %.0f%%", min(rs$p_decline), max(rs$p_decline))
    },
    pva_p_quasi_ext_range = {
      rs <- pva_by_nmin[[paste0("nmin_", nmin)]]$risk_summary
      rs <- rs[rs$scenario != "No additional mortality", ]
      sprintf("%.0f - %.0f%%", min(rs$p_quasi_extinction), max(rs$p_quasi_extinction))
    }
  ) %>%
  ungroup()

cat("\n=== Summary across the 3 approaches, by Nmin ===\n")
print(as.data.frame(comparison_by_nmin))

# What Nmin (or N0) would each imposed threshold represent as a % of, for
# context against the hypothesis of ~1000 combined fatalities/year: purely
# illustrative cross-check, not a 4th model.
cat("\n=== Reality-check: a combined 1000 fatalities/year against each Nmin scenario ===\n")
print(tibble::tibble(nmin = nmin_scenarios, pct_of_nmin_if_1000 = 100 * 1000 / nmin_scenarios))

# ---- 5. Figure: PVA trajectories, both Nmin scenarios side by side --------
pva_summary_by_year <- bind_rows(lapply(names(pva_by_nmin), function(nm) {
  pva_by_nmin[[nm]]$trajectories %>%
    group_by(scenario, year) %>%
    summarise(median_N = median(N), q05 = quantile(N, 0.05), q25 = quantile(N, 0.25),
              q75 = quantile(N, 0.75), q95 = quantile(N, 0.95), .groups = "drop") %>%
    mutate(nmin_label = paste0("N0 = ", pva_by_nmin[[nm]]$n0_total))
}))

fig_pva_nmin_sensitivity <- file.path(fig_dir, "pva_nmin_sensitivity.png")
ggsave(fig_pva_nmin_sensitivity, width = 11, height = 6.5, dpi = 150, plot = {
  ggplot(pva_summary_by_year, aes(x = year)) +
    geom_ribbon(aes(ymin = q05, ymax = q95, fill = scenario), alpha = 0.15) +
    geom_ribbon(aes(ymin = q25, ymax = q75, fill = scenario), alpha = 0.3) +
    geom_line(aes(y = median_N, colour = scenario), linewidth = 0.9) +
    facet_grid(nmin_label ~ scenario, scales = "free_y") +
    scale_colour_viridis_d(option = "C", end = 0.8, guide = "none") +
    scale_fill_viridis_d(option = "C", end = 0.8, guide = "none") +
    labs(
      x = "Year", y = "Population (both sexes)",
      title = "PVA-lite: sensitivity to the assumed starting population (N0 = Nmin)",
      subtitle = "Same absolute annual removal (112/88 bats/yr) applied to a 2x larger starting population"
    ) +
    theme_minimal() +
    theme(strip.text = element_text(face = "bold", size = 8), plot.subtitle = element_text(size = 8))
})
cat("\nWrote", fig_pva_nmin_sensitivity, "\n")
