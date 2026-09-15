##
## Stochastic population projection ("PVA-lite") built on the validated
## Dekker & Limpens (2024) V. murinus stage matrix (R/leslie_dekker_limpens.R,
## lambda = 1.047 matching the published value exactly). Adds environmental
## stochasticity (vital rates redrawn each year) and demographic
## stochasticity (binomial/Poisson realisations of those rates applied to
## the actual, finite population), then projects 25 years under:
##   (1) no additional mortality -- the model's own baseline trend;
##   (2) the imposed BSH threshold (144 bats/year) as extra removal;
##   (3) the imposed DGY threshold (120 bats/year) as extra removal.
## This is the direct, practical question behind the whole PBR review:
## given a validated, species-specific demographic model, do the
## currently-imposed fatality thresholds look sustainable over a
## generation-scale horizon?
##
## Starting population: N0 = 4000 (both sexes), the Nmin assumed per
## facility throughout the PBR analysis (see R/pbr_analysis.R /
## inputs/pbrSettings_BSH_DGY.R) -- halved here to female-only, since the
## validated matrix (and hence this projection) tracks females only, the
## male block never feeding back into it (see leslie_dekker_limpens.R).
## Turbine removal is split evenly across sexes (no strong evidence of an
## extreme sex bias was located for this analysis) and, within females,
## proportionally across the juvenile/adult stages by their relative
## abundance that year -- both simplifications, flagged rather than hidden.
##
## Environmental stochasticity: no inter-annual variance estimate for
## V. murinus vital rates was available this session; a CV of 10% of the
## mean is used for each rate as a working assumption (moment-matched to
## Beta distributions, which respect the [0, 1] bound), not a
## literature-derived figure -- a sensitivity target for follow-up.
##

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(ggplot2)
})

source("R/leslie_dekker_limpens.R")

# ---- Beta distribution shape parameters from mean + CV --------------------
beta_params <- function(mean, cv) {
  var <- (mean * cv)^2
  stopifnot(var < mean * (1 - mean))  # CV small enough to stay a valid Beta
  common <- mean * (1 - mean) / var - 1
  list(shape1 = mean * common, shape2 = (1 - mean) * common)
}

cv_vital_rates <- 0.10
bp_s_adult <- beta_params(S_af, cv_vital_rates)
bp_s_juv <- beta_params(S_j, cv_vital_rates)
bp_p_breed <- beta_params(p_breed, cv_vital_rates)

# ---- One stochastic 25-year trajectory -------------------------------------
simulate_trajectory <- function(n0_juv, n0_adult, n_years = 25, annual_removal = 0, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  n_juv <- n0_juv
  n_adult <- n0_adult
  traj <- numeric(n_years + 1)
  traj[1] <- 2 * (n_juv + n_adult)  # both sexes

  for (yr in seq_len(n_years)) {
    s_adult_t <- rbeta(1, bp_s_adult$shape1, bp_s_adult$shape2)
    s_juv_t <- rbeta(1, bp_s_juv$shape1, bp_s_juv$shape2)
    p_breed_t <- rbeta(1, bp_p_breed$shape1, bp_p_breed$shape2)

    # Reproduction (demographic stochasticity: binomial/Poisson realisations)
    n_breeders <- rbinom(1, n_adult, p_breed_t)
    n_pups_total <- rpois(1, n_breeders * litter)
    n_female_pups <- rbinom(1, n_pups_total, sex_ratio)
    n_juv_recruits <- rbinom(1, n_female_pups, s_juv_t)

    # Survival transitions
    n_juv_survive <- rbinom(1, n_juv, s_juv_t)     # this year's juveniles -> next year's adults
    n_adult_survive <- rbinom(1, n_adult, s_adult_t)

    n_juv_next <- n_juv_recruits
    n_adult_next <- n_juv_survive + n_adult_survive

    # Turbine removal: split across sexes evenly, and within females across
    # stages proportional to relative abundance
    female_removal <- annual_removal / 2
    total_female <- n_juv_next + n_adult_next
    if (total_female > 0 && female_removal > 0) {
      juv_share <- round(female_removal * n_juv_next / total_female)
      adult_share <- round(female_removal * n_adult_next / total_female)
      n_juv_next <- max(0, n_juv_next - juv_share)
      n_adult_next <- max(0, n_adult_next - adult_share)
    }

    n_juv <- n_juv_next
    n_adult <- n_adult_next
    traj[yr + 1] <- 2 * (n_juv + n_adult)
  }
  traj
}

# ---- Monte Carlo: n_reps trajectories per scenario -------------------------
n0_total <- 4000            # Nmin per facility, established in the PBR analysis
n0_female <- n0_total / 2
n0_juv <- round(n0_female * 0.35)   # initial stage split ~ stable stage distribution (see below)
n0_adult <- n0_female - n0_juv

# Use the validated matrix's own stable stage distribution for a sensible starting split
stable_dist <- Re(eigen(A_female)$vectors[, 1]); stable_dist <- stable_dist / sum(stable_dist)
n0_juv <- round(n0_female * stable_dist[1])
n0_adult <- n0_female - n0_juv

n_reps <- 1000
n_years <- 25
scenarios <- tibble::tibble(
  scenario = c("No additional mortality", "BSH threshold (144/yr)", "DGY threshold (120/yr)"),
  annual_removal = c(0, 144, 120)
)

set.seed(42)
all_trajectories <- lapply(seq_len(nrow(scenarios)), function(i) {
  reps <- replicate(n_reps, simulate_trajectory(n0_juv, n0_adult, n_years, scenarios$annual_removal[i]))
  tibble::tibble(
    scenario = scenarios$scenario[i],
    year = rep(0:n_years, n_reps),
    rep = rep(seq_len(n_reps), each = n_years + 1),
    N = as.vector(reps)
  )
}) %>% bind_rows() %>%
  mutate(scenario = factor(scenario, levels = scenarios$scenario))

# ---- Summaries --------------------------------------------------------------
summary_by_year <- all_trajectories %>%
  group_by(scenario, year) %>%
  summarise(
    median_N = median(N), q05 = quantile(N, 0.05), q25 = quantile(N, 0.25),
    q75 = quantile(N, 0.75), q95 = quantile(N, 0.95), .groups = "drop"
  )

final_year <- all_trajectories %>% filter(year == n_years)
quasi_ext_threshold <- 0.10 * (2 * n0_female)  # 10% of starting N, a common PVA convention

risk_summary <- all_trajectories %>%
  group_by(scenario, rep) %>%
  summarise(
    final_N = N[year == n_years],
    ever_below_threshold = any(N < quasi_ext_threshold),
    .groups = "drop"
  ) %>%
  group_by(scenario) %>%
  summarise(
    median_final_N = median(final_N),
    q05_final_N = quantile(final_N, 0.05),
    q95_final_N = quantile(final_N, 0.95),
    p_decline = mean(final_N < 2 * n0_female),
    p_quasi_extinction = mean(ever_below_threshold),
    .groups = "drop"
  )

cat(sprintf("Starting population (both sexes): %d (female: %d = %d juv + %d adult)\n",
            2 * n0_female, n0_female, n0_juv, n0_adult))
cat(sprintf("Quasi-extinction threshold: %d individuals (10%% of N0)\n\n", round(quasi_ext_threshold)))
cat("Risk summary after 25 years, ", n_reps, " replicate trajectories per scenario:\n", sep = "")
print(risk_summary)

# ---- Figure -----------------------------------------------------------------
fig_dir <- "outputs/figures"
dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)

p <- ggplot(summary_by_year, aes(x = year)) +
  geom_ribbon(aes(ymin = q05, ymax = q95, fill = scenario), alpha = 0.15) +
  geom_ribbon(aes(ymin = q25, ymax = q75, fill = scenario), alpha = 0.3) +
  geom_line(aes(y = median_N, colour = scenario), linewidth = 0.9) +
  geom_hline(yintercept = 2 * n0_female, linetype = "dotted", colour = "grey40") +
  facet_wrap(~scenario, nrow = 1) +
  scale_colour_viridis_d(option = "C", end = 0.8, guide = "none") +
  scale_fill_viridis_d(option = "C", end = 0.8, guide = "none") +
  labs(
    x = "Year", y = "Population (both sexes)",
    title = "25-year stochastic projection, validated V. murinus demographic model",
    subtitle = paste0(
      "N0 = ", 2 * n0_female, " (Nmin per facility); shaded: 50%/90% of ", n_reps, " trajectories; dotted: N0"
    )
  ) +
  theme_minimal() +
  theme(strip.text = element_text(face = "bold"), plot.subtitle = element_text(size = 8))
ggsave(file.path(fig_dir, "leslie_pva_25yr_projection.png"), p, width = 11, height = 4.5, dpi = 150)

cat("\nWrote", file.path(fig_dir, "leslie_pva_25yr_projection.png"), "\n")
