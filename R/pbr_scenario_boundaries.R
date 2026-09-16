##
## Conservative vs. relaxed PBR scenario boundaries, for an informed
## lender/promoter discussion -- NOT an argument to raise the imposed
## thresholds (112 Project 1 / 88 Project 2), which stay fixed throughout. Instead:
## (a) makes explicit which knobs PBR's own uncertainty lives in (Fr's
## dependence on an IUCN-status judgement call, the choice of lambda_max
## benchmark, Nmin's dependence on an unresolved local-vs-regional/
## migratory-passage population size -- see R/nmin_sensitivity_analysis.R);
## (b) cross-checks the resulting boundary against BOTH the imposed
## thresholds and the actual mortality observed since curtailment was
## applied to the most critical turbines, which is the strongest empirical
## anchor available (measured outcome, not a modelled one).
##
## Curtailment context (Paulo, 2026-09): mortality during the high-risk
## migration months (Apr-Jun, Sep-Oct) fell from "dozens of bats/month"
## pre-curtailment to ~1-2 individuals/month post-curtailment, on the
## turbines curtailed. Exact monthly counts were not given in this
## session, so the post-curtailment annual estimate below is built from
## explicit, labelled brackets (illustrative, to be replaced with the
## actual monthly series once available) -- not treated as a precise
## figure anywhere it is used.
##

suppressPackageStartupMessages({ library(dplyr); library(tidyr); library(ggplot2) })

source("R/pbr_functions.R")
source("R/leslie_dekker_limpens.R")
source("inputs/pbrSettings_BSH_DGY.R")

fig_dir <- "outputs/figures"
dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)

# ---- 1. The three knobs PBR's own uncertainty lives in --------------------
# Fr: the regulatory-judgement-call knob (which IUCN reading applies).
fr_scenarios_boundary <- tibble::tibble(
  fr_label = c("Conservative (Fr=0.3, Near Threatened)",
               "Central (Fr=0.5, Least Concern, declining/unknown)",
               "Relaxed (Fr=1.0, Least Concern, stable/increasing)"),
  fr = c(0.3, 0.5, 1.0)
)
# lambda_max: the paper's own two fixed benchmarks (not widened further here).
lambda_max_boundary <- lambda_max_benchmarks  # 1.24, 1.20
# Nmin: the population-size knob (see R/nmin_sensitivity_analysis.R).
nmin_boundary <- c(3022, 6044)

boundary_grid <- tidyr::expand_grid(fr_scenarios_boundary, lambda_max = lambda_max_boundary, nmin = nmin_boundary) %>%
  mutate(PBR = pbr_from_components(nmin, fr, lambda_max))

cat("=== Full PBR boundary grid (Fr x lambda_max x Nmin) ===\n")
print(as.data.frame(boundary_grid))

boundary_summary <- boundary_grid %>%
  group_by(nmin) %>%
  summarise(PBR_min = min(PBR), PBR_max = max(PBR), .groups = "drop")
cat("\n=== Overall PBR boundary (most conservative to most relaxed combination), by Nmin ===\n")
print(boundary_summary)

# ---- 2. Post-curtailment actual mortality (illustrative bracket) ----------
# Pre-curtailment: "dozens/month" (bracketed 15-40/month) across the 5
# highest-risk months (Apr, May, Jun, Sep, Oct); negligible outside them
# (no evidence otherwise this session). Post-curtailment: 1-2 individuals/
# month over the same 5 months, ~0 assumed outside them.
high_risk_months <- 5
pre_curtailment_annual  <- c(low = 15, high = 40) * high_risk_months
post_curtailment_annual <- c(low = 1, high = 2) * high_risk_months

cat(sprintf(
  "\nIllustrative annual mortality bracket -- pre-curtailment: %d-%d; post-curtailment: %d-%d (both sexes, combined months, SINGLE facility's critical turbines as described; NOT yet split by facility or confirmed against exact monthly counts)\n",
  pre_curtailment_annual["low"], pre_curtailment_annual["high"],
  post_curtailment_annual["low"], post_curtailment_annual["high"]
))

# ---- 3. Where do the imposed thresholds and curtailment outcome sit? ------
reference_points <- tibble::tibble(
  label = c("Project 1 imposed threshold", "Project 2 imposed threshold",
            "Pre-curtailment (illustrative)", "Post-curtailment (illustrative)"),
  value_low  = c(112, 88, pre_curtailment_annual["low"], post_curtailment_annual["low"]),
  value_high = c(112, 88, pre_curtailment_annual["high"], post_curtailment_annual["high"]),
  type = c("Imposed threshold", "Imposed threshold", "Observed mortality", "Observed mortality")
)

cat("\n=== Reference points against the boundary ===\n")
print(as.data.frame(reference_points))

# ---- 4. Figure: full boundary vs thresholds vs curtailment outcome --------
plot_grid <- boundary_grid %>%
  mutate(
    nmin_label = paste0("Nmin = ", format(nmin, big.mark = ",")),
    scenario_label = paste0(fr_label, "\nlambda_max=", lambda_max)
  )

fig_pbr_boundaries <- file.path(fig_dir, "pbr_scenario_boundaries.png")
ggsave(fig_pbr_boundaries, width = 11, height = 6.5, dpi = 150, plot = {
  ggplot(plot_grid, aes(x = reorder(fr_label, fr), y = PBR, colour = factor(lambda_max))) +
    geom_rect(
      data = tibble::tibble(nmin_label = unique(plot_grid$nmin_label)),
      aes(xmin = -Inf, xmax = Inf, ymin = post_curtailment_annual["low"], ymax = post_curtailment_annual["high"]),
      fill = "forestgreen", alpha = 0.25, inherit.aes = FALSE
    ) +
    geom_rect(
      data = tibble::tibble(nmin_label = unique(plot_grid$nmin_label)),
      aes(xmin = -Inf, xmax = Inf, ymin = pre_curtailment_annual["low"], ymax = pre_curtailment_annual["high"]),
      fill = "firebrick", alpha = 0.12, inherit.aes = FALSE
    ) +
    geom_point(position = position_dodge(width = 0.4), size = 3) +
    geom_line(aes(group = factor(lambda_max)), position = position_dodge(width = 0.4), linewidth = 0.4) +
    geom_hline(data = pbr_thresholds, aes(yintercept = threshold, linetype = facility), colour = "grey30") +
    facet_wrap(~nmin_label) +
    scale_y_log10(limits = c(4, 1200)) +
    scale_colour_manual(name = "lambda_max benchmark", values = c("1.2" = "cyan4", "1.24" = "chartreuse4")) +
    scale_linetype_manual(name = "Imposed threshold",
                           values = setNames(c("dashed", "dotted"), facility_labels)) +
    labs(
      x = "Recovery factor (Fr) scenario", y = "PBR (bats/year), log scale",
      title = "PBR scenario boundaries: conservative to relaxed, vs. imposed thresholds and curtailment outcome",
      subtitle = paste0(
        "Green band: illustrative post-curtailment mortality (", post_curtailment_annual["low"], "-",
        post_curtailment_annual["high"], " bats/yr); red band: illustrative pre-curtailment (",
        pre_curtailment_annual["low"], "-", pre_curtailment_annual["high"], " bats/yr) -- both high-risk months only, pending exact counts"
      )
    ) +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 20, hjust = 1), plot.subtitle = element_text(size = 8),
          strip.text = element_text(face = "bold"))
})
cat("\nWrote", fig_pbr_boundaries, "\n")

# ---- 5. PVA-lite: imposed threshold vs. actual post-curtailment removal ---
# Same validated structure/settings as R/nmin_sensitivity_analysis.R, run at
# Nmin = 3022 (median-fallback, status quo), comparing the imposed-threshold
# removal against the illustrative actual post-curtailment removal rate.
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

n0_total <- 3022
n0_female <- n0_total / 2
n0_juv <- round(n0_female * stable_dist[1]); n0_adult <- n0_female - n0_juv

curtailment_scenarios <- tibble::tibble(
  scenario = c("No additional mortality", "Project 1 imposed threshold (112/yr)",
               "Pre-curtailment (illustrative, ~140/yr)", "Post-curtailment (illustrative, ~7.5/yr)"),
  annual_removal = c(0, 112, mean(pre_curtailment_annual), mean(post_curtailment_annual))
) %>%
  mutate(scenario = factor(scenario, levels = scenario))

set.seed(pva_seed)
curtailment_trajectories <- lapply(seq_len(nrow(curtailment_scenarios)), function(i) {
  reps <- replicate(pva_n_reps, simulate_pva_trajectory(n0_juv, n0_adult, pva_n_years, curtailment_scenarios$annual_removal[i]))
  tibble::tibble(
    scenario = curtailment_scenarios$scenario[i], year = rep(0:pva_n_years, pva_n_reps),
    rep = rep(seq_len(pva_n_reps), each = pva_n_years + 1), N = as.vector(reps)
  )
}) %>% bind_rows()

quasi_ext_threshold_curt <- pva_quasi_extinction_fraction * n0_total
curtailment_risk_summary <- curtailment_trajectories %>%
  group_by(scenario, rep) %>%
  summarise(final_N = N[year == pva_n_years], ever_below_threshold = any(N < quasi_ext_threshold_curt), .groups = "drop") %>%
  group_by(scenario) %>%
  summarise(
    median_final_N = median(final_N), q05_final_N = quantile(final_N, 0.05), q95_final_N = quantile(final_N, 0.95),
    p_decline = mean(final_N < n0_total) * 100, p_quasi_extinction = mean(ever_below_threshold) * 100, .groups = "drop"
  )

cat("\n=== PVA-lite: imposed threshold vs. illustrative pre/post-curtailment removal (Nmin=3022) ===\n")
print(as.data.frame(curtailment_risk_summary))

curtailment_summary_by_year <- curtailment_trajectories %>%
  group_by(scenario, year) %>%
  summarise(median_N = median(N), q05 = quantile(N, 0.05), q25 = quantile(N, 0.25),
            q75 = quantile(N, 0.75), q95 = quantile(N, 0.95), .groups = "drop")

fig_pva_curtailment <- file.path(fig_dir, "pva_curtailment_effect.png")
ggsave(fig_pva_curtailment, width = 11, height = 4.5, dpi = 150, plot = {
  ggplot(curtailment_summary_by_year, aes(x = year)) +
    geom_ribbon(aes(ymin = q05, ymax = q95, fill = scenario), alpha = 0.15) +
    geom_ribbon(aes(ymin = q25, ymax = q75, fill = scenario), alpha = 0.3) +
    geom_line(aes(y = median_N, colour = scenario), linewidth = 0.9) +
    geom_hline(yintercept = n0_total, linetype = "dotted", colour = "grey40") +
    facet_wrap(~scenario, nrow = 1) +
    scale_colour_viridis_d(option = "C", end = 0.8, guide = "none") +
    scale_fill_viridis_d(option = "C", end = 0.8, guide = "none") +
    labs(
      x = "Year", y = "Population (both sexes)",
      title = "PVA-lite: the effect of curtailment on projected trend (Nmin = 3022)",
      subtitle = "Pre/post-curtailment removal rates are illustrative brackets, not confirmed monthly counts"
    ) +
    theme_minimal() +
    theme(strip.text = element_text(face = "bold", size = 8), plot.subtitle = element_text(size = 8))
})
cat("\nWrote", fig_pva_curtailment, "\n")
