# ============================================================
# Reverse-engineering the BSH/DGY Parti-coloured bat (Vespertilio murinus)
# fatality thresholds against the Frick et al. (2026) PBR framework, and
# testing their robustness to demographic parameter uncertainty.
#
# Context: BSH and DGY were each given a fixed annual fatality threshold
# for the Parti-coloured bat (144 and 120 bats/year respectively), with no
# documented derivation available to us, and no way to obtain an
# independent local population estimate for the region (population size
# data are simply not available -- there is no "Simon parametrisation" to
# fall back on). This script:
#  (1) finds the most parsimonious Nmin/lambda_max/Fr combination that
#      reproduces both thresholds under Frick et al.'s PBR equation;
#  (2) sanity-checks whether the implied Nmin is a plausible local
#      population size for the region, via the interaction area it would
#      imply under the paper's own default density proxy;
#  (3) corrects Fr from 0.3 (Near Threatened) to 0.5 (Least Concern) --
#      the species is not a threatened category, so Fr=0.3 is not
#      defensible -- and checks how far the resulting PBR sits from 144
#      and 120 with everything else unchanged;
#  (4) replaces the paper's single fixed lambda_max benchmarks with a
#      plausible range for adult survival (s) and age at first breeding
#      (alpha), building PBR response surfaces and a Monte Carlo
#      distribution, to test whether a flexible, range-based
#      parametrisation would justify materially higher thresholds than
#      the fixed 144/120 currently imposed.
# ============================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(ggplot2)
})

# Mirrors fr_from_iucn() / lambda_max_niel() in pbr_frick_main.R; redefined
# here (not sourced) so this script doesn't also re-run that file's
# exploratory plotting code as a side effect.
fr_from_iucn <- function(iucn_status) {
  status <- tolower(iucn_status)
  dplyr::case_when(
    status == "least concern (stable or increasing)" ~ 1.0,
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

nmin_from_bounds <- function(NL, NU, ci_level = 0.60) {
  z_lower <- qnorm((1 - ci_level) / 2)
  z_upper <- qnorm(1 - (1 - ci_level) / 2)
  Nhat <- sqrt(NL * NU)
  CV <- sqrt(exp((log(NU / NL) / (2 * z_upper))^2) - 1)
  Nhat * exp(z_lower * CV)
}

nmin_from_area_density <- function(area_km2, density_low, density_high,
                                    scaling_factor = 2, ci_level = 0.60) {
  NL <- area_km2 * scaling_factor * density_low
  NU <- area_km2 * scaling_factor * density_high
  nmin_from_bounds(NL, NU, ci_level)
}

# ---- 1. Known/imposed thresholds ---------------------------------------
thresholds <- tibble::tibble(
  facility  = c("BSH", "DGY"),
  threshold = c(144, 120)  # bats/year
)

# ---- 2. Backward solve under the apparent original assumption -----------
#         (Fr = 0.3, "Near Threatened") -- reconstructs how 144/120 were
#         most likely obtained, even though this Fr is later corrected.
required_nmin <- function(pbr_target, Fr, lambda_max) {
  Rmax <- lambda_max - 1
  pbr_target / (0.5 * Rmax * Fr)
}

Fr_original_guess <- fr_from_iucn("Near Threatened")  # 0.3

backsolve_grid <- tidyr::expand_grid(
  thresholds,
  lambda_max_fixed = c(1.24, 1.20)
) %>%
  mutate(Nmin_required = required_nmin(threshold, Fr_original_guess, lambda_max_fixed))

print(backsolve_grid)
# BSH/144 paired with lambda_max = 1.24, and DGY/120 paired with
# lambda_max = 1.20, both back-solve to the SAME Nmin_required = 4000 under
# Fr = 0.3. No other pairing gives a shared, round Nmin. This is the most
# parsimonious explanation for how 144 and 120 were obtained: a single
# assumed local population of ~4000 individuals, the paper's two
# lambda_max benchmarks applied per facility, and Fr = 0.3.

Nmin_assumed <- 4000

# ---- 3. Is Nmin = 4000 a plausible local population for the region? -----
# There is no independent population estimate for Vespertilio murinus in
# this region (and no realistic prospect of obtaining one), so plausibility
# is checked the other way around: what interaction area would Nmin = 4000
# imply, under the paper's own default density proxy (5.5-18.2 bats/km^2,
# based on Pipistrellus pipistrellus, used in absence of species-specific
# density data)?
implied_area <- function(nmin_target, density_low = 5.5, density_high = 18.2,
                          scaling_factor = 2, ci_level = 0.60) {
  f <- function(area) {
    nmin_from_area_density(area, density_low, density_high, scaling_factor, ci_level) - nmin_target
  }
  uniroot(f, c(1, 1e7))$root
}

area_for_4000  <- implied_area(4000)
area_for_40000 <- implied_area(40000)

cat(sprintf(
  "Implied interaction area for Nmin=%s (default density 5.5-18.2/km2): %.0f km2\n",
  format(4000, big.mark = ","), area_for_4000
))
cat(sprintf(
  "Implied interaction area for Nmin=%s (default density 5.5-18.2/km2): %.0f km2\n",
  format(40000, big.mark = ","), area_for_40000
))
# Nmin = 4000 implies an interaction area of ~400 km^2 -- large for a
# single wind farm's footprint, but within the range of a broader regional
# "local population" buffer as the paper's project-scale proxy intends.
# Nmin = 40,000 (10x) would imply ~4,000 km^2 -- an order of magnitude
# larger, well into sub-regional/national scale, which is hard to square
# with a *local*, project-relevant population as the framework defines it.
# This does not prove 4000 is correct, but it is the more defensible order
# of magnitude of the two.

# ---- 4. Correct Fr: the species is not a threatened category -------------
# Fr = 0.3 corresponds to "Near Threatened" (Table 1 / RecoveryFactors
# sheet). Vespertilio murinus is not threatened, so the defensible value is
# Fr = 0.5 ("Least Concern", declining-or-unknown trend; see Appendix S1
# RecoveryFactors sheet, which further distinguishes LC-stable/increasing
# at Fr = 1 -- not used here for lack of trend data).
Fr <- fr_from_iucn("Least Concern")  # 0.5

pbr_corrected_fr <- tibble::tibble(lambda_max_fixed = c(1.24, 1.20)) %>%
  mutate(PBR = pbr_from_components(Nmin_assumed, Fr, lambda_max_fixed))

print(pbr_corrected_fr)
# With Nmin held at 4000 and only Fr corrected (0.3 -> 0.5), PBR moves from
# 144 to 240 (lambda_max = 1.24) and from 120 to 200 (lambda_max = 1.20):
# ~1.67x higher, exactly the Fr ratio, before any change to lambda_max
# itself is even considered.

# ---- 5. Deconstruction: replace fixed lambda_max with plausible (s, alpha)-
#         ranges and propagate the resulting PBR distribution, now under
#         the corrected Fr = 0.5 -------------------------------------------
# Ranges as used in the earlier sensitivity/elasticity analysis: adult
# survival 0.70-0.95; age at first breeding 2-4 years. alpha = 1.5 (needed
# to obtain the paper's lambda_max = 1.24 benchmark) sits outside this
# range and is treated by the report itself as an optimistic edge case.

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

threshold_percentiles <- thresholds %>%
  mutate(
    percentile_rank_corrected = sapply(threshold, function(x) mean(sim$PBR < x) * 100),
    median_ratio = median(sim$PBR) / threshold,
    p95_ratio = quantile(sim$PBR, 0.95) / threshold
  )
print(threshold_percentiles)
# Under the corrected Fr = 0.5, the currently imposed thresholds (144, 120)
# now sit at the ~16th and ~5th percentile of the simulated PBR
# distribution -- i.e. LOW relative to a defensible, non-optimistic
# parametrisation, not high. Median simulated PBR is ~1.3-1.5x the current
# thresholds; the 95th percentile reaches ~1.7-1.9x. This partially
# supports the expectation that a flexible parametrisation would justify a
# materially higher threshold, approaching (though not uniformly reaching)
# double the current fixed values.

# ---- 6. PBR response surfaces over (s, alpha) -----------------------------
# Reproduces the lambda_max(s, alpha) and PBR(s, alpha) surfaces sketched
# in the earlier Word report (Figures 3-4), now consistently using
# Nmin = 4000 and the corrected Fr = 0.5.
surface_grid <- tidyr::expand_grid(
  s     = seq(0.70, 0.95, by = 0.01),
  alpha = seq(2, 4, by = 0.05)
) %>%
  mutate(
    lambda_max = lambda_max_niel(s, alpha),
    PBR = pbr_from_components(Nmin = Nmin_assumed, Fr = Fr, lambda_max = lambda_max)
  )

p_lambda_surface <- ggplot(surface_grid, aes(x = s, y = alpha, fill = lambda_max)) +
  geom_raster(interpolate = TRUE) +
  geom_contour(aes(z = lambda_max), breaks = c(1.20, 1.24), colour = "white", linewidth = 0.4) +
  scale_fill_viridis_c(name = expression(lambda[max]), option = "C") +
  labs(
    x = "Adult survival (s)", y = "Age at first breeding (α)",
    title = expression("Maximum population growth rate " * lambda[max]),
    subtitle = "White contours: λ[max] = 1.20 and 1.24 (paper's fixed benchmarks)"
  ) +
  theme_minimal()

p_pbr_surface <- ggplot(surface_grid, aes(x = s, y = alpha, fill = PBR)) +
  geom_raster(interpolate = TRUE) +
  geom_contour(aes(z = PBR), breaks = c(120, 144), colour = "white", linewidth = 0.4) +
  scale_fill_viridis_c(name = "PBR", option = "C") +
  labs(
    x = "Adult survival (s)", y = "Age at first breeding (α)",
    title = "PBR response surface (Fr = 0.5, Nmin = 4000)",
    subtitle = "White contours: PBR = 120 (DGY) and 144 (BSH), current imposed thresholds"
  ) +
  theme_minimal()

p_lambda_surface
p_pbr_surface

# ---- 7. Plot: simulated PBR distribution vs. imposed thresholds ----------
p_pbr_dist <- ggplot(sim, aes(x = PBR)) +
  geom_density(fill = "grey80", colour = "grey40") +
  geom_vline(
    data = thresholds, aes(xintercept = threshold, colour = facility),
    linewidth = 0.8, linetype = "dashed"
  ) +
  labs(
    title = "Simulated PBR distribution for the Parti-coloured bat (Fr = 0.5, Nmin = 4000)",
    subtitle = "Adult survival ~ U(0.70, 0.95); age at first breeding ~ U(2, 4) years",
    x = "PBR (bats/year)", y = "Density", colour = "Imposed threshold"
  ) +
  theme_minimal()

p_pbr_dist
