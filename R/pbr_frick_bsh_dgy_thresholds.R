# ============================================================
# Reverse-engineering the BSH/DGY Parti-coloured bat (Vespertilio murinus)
# fatality thresholds against the Frick et al. (2026) PBR framework, and
# testing their robustness to demographic parameter uncertainty.
#
# Context: BSH and DGY were each given a fixed annual fatality threshold
# for the Parti-coloured bat (144 and 120 bats/year respectively), with no
# documented derivation available to us. This script (1) finds the most
# parsimonious Nmin/lambda_max combination that reproduces both thresholds
# exactly under Frick et al.'s PBR equation, then (2) replaces the paper's
# single fixed lambda_max benchmarks with a plausible range for adult
# survival (s) and age at first breeding (alpha), to show how much the
# resulting PBR varies once demographic uncertainty is propagated instead
# of pinned to one point value.
# ============================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
})

# Mirrors fr_from_iucn() / lambda_max_niel() in pbr_frick_main.R; redefined
# here (not sourced) so this script doesn't also re-run that file's
# exploratory plotting code as a side effect.
fr_from_iucn <- function(iucn_status) {
  status <- tolower(iucn_status)
  dplyr::case_when(
    status %in% c("least concern (stable or increasing)") ~ 1.0,
    status %in% c("least concern", "lc", "least concern (declining or unknown)") ~ 0.5,
    status %in% c("near threatened", "nt") ~ 0.3,
    status %in% c("vulnerable", "vu") ~ 0.1,
    status %in% c("endangered", "en") ~ 0.0,
    status %in% c("critically endangered", "cr") ~ 0.0,
    status %in% c("data deficient", "dd") ~ 0.2,
    TRUE ~ NA_real_
  )
}

lambda_max_niel <- function(s, alpha) {
  disc <- (s - s * alpha - alpha - 1)^2 - 4 * s * alpha^2
  ((s * alpha - s + alpha + 1) + sqrt(pmax(disc, 0))) / (2 * alpha)
}

pbr_from_components <- function(Nmin, Fr, lambda_max) {
  Rmax <- lambda_max - 1
  0.5 * Rmax * Fr * Nmin
}

# ---- 1. Known/imposed thresholds ---------------------------------------
thresholds <- tibble::tibble(
  facility  = c("BSH", "DGY"),
  threshold = c(144, 120)  # bats/year
)

# Assumed status used throughout the prior analysis for this species here;
# to be confirmed against the actual national/local assessment.
iucn_status <- "Near Threatened"
Fr <- fr_from_iucn(iucn_status)  # 0.3, per paper's Table 1 and the Appendix
                                  # S1 RecoveryFactors sheet

# ---- 2. Backward solve: Nmin required to hit each threshold, for the ------
#         paper's two fixed lambda_max benchmarks (1.24 "general"; 1.20 for
#         alpha=2 taxa such as Pteropodidae) ------------------------------
required_nmin <- function(pbr_target, Fr, lambda_max) {
  Rmax <- lambda_max - 1
  pbr_target / (0.5 * Rmax * Fr)
}

backsolve_grid <- tidyr::expand_grid(
  thresholds,
  lambda_max_fixed = c(1.24, 1.20)
) %>%
  mutate(Nmin_required = required_nmin(threshold, Fr, lambda_max_fixed))

print(backsolve_grid)
# BSH/144 paired with lambda_max = 1.24, and DGY/120 paired with
# lambda_max = 1.20, both back-solve to the SAME Nmin_required = 4000.
# No other pairing of threshold x lambda_max gives a shared, round Nmin.
# This is the most parsimonious explanation available for how 144 and 120
# were obtained: a single assumed local population of ~4000 individuals,
# combined with the paper's two lambda_max benchmarks applied per facility.

Nmin_assumed <- 4000

# ---- 3. Deconstruction: replace fixed lambda_max with plausible (s, alpha) -
#         ranges and propagate the resulting PBR distribution -------------
# Ranges as used in the earlier sensitivity/elasticity analysis: adult
# survival 0.70-0.95; age at first breeding 2-4 years. Note alpha = 1.5
# (needed to obtain the paper's own lambda_max = 1.24 benchmark) sits
# outside this range -- the report's own text treats alpha = 2 as the
# "defensible and conservative" lower bound, with 3-4 years "very
# precautionary". The imposed thresholds therefore rest on a more
# optimistic maturation assumption than the range considered plausible.

set.seed(1)
n_sim <- 20000

sim <- tibble::tibble(
  s     = runif(n_sim, 0.70, 0.95),
  alpha = runif(n_sim, 2, 4)
) %>%
  mutate(
    lambda_max = lambda_max_niel(s, alpha),
    PBR = pbr_from_components(Nmin = Nmin_assumed, Fr = Fr, lambda_max = lambda_max)
  )

pbr_quantiles <- quantile(sim$PBR, c(0.05, 0.25, 0.50, 0.75, 0.95))
print(round(pbr_quantiles, 1))

# Where do the imposed thresholds fall within this simulated distribution?
threshold_percentiles <- thresholds %>%
  mutate(percentile_rank = sapply(threshold, function(x) mean(sim$PBR < x) * 100))
print(threshold_percentiles)
# BSH (144) and DGY (120) sit at roughly the 85th and 63rd percentile of the
# simulated distribution: i.e. both imposed thresholds are on the
# permissive (higher-PBR) side of what a plausible demographic range
# produces, not at its median or conservative tail.

# ---- 4. Plot: simulated PBR distribution vs. imposed thresholds ----------
p_pbr_dist <- ggplot(sim, aes(x = PBR)) +
  geom_density(fill = "grey80", colour = "grey40") +
  geom_vline(
    data = thresholds, aes(xintercept = threshold, colour = facility),
    linewidth = 0.8, linetype = "dashed"
  ) +
  labs(
    title = "Simulated PBR distribution for the Parti-coloured bat (Fr = 0.3, Nmin = 4000)",
    subtitle = "Adult survival ~ U(0.70, 0.95); age at first breeding ~ U(2, 4) years",
    x = "PBR (bats/year)", y = "Density", colour = "Imposed threshold"
  ) +
  theme_minimal()

p_pbr_dist
