# ============================================================
# Shared pure functions for the PBR / Frick et al. (2026) analysis.
# Sourced by the analysis scripts (pbr_frick_main.R,
# pbr_frick_bsh_dgy_thresholds.R, pbr_frick_response_surfaces.R) and by
# report/pbr_frick_technical_note.Rmd, so the formulas exist in exactly one
# place.
# ============================================================

# Recovery factor (Fr) by IUCN status -- Table 1 of Frick et al. (2026),
# with the Least Concern split into declining/unknown vs stable/increasing
# from the Appendix S1 RecoveryFactors sheet.
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

# lambda_max via the Niel & Lebreton (2005) demographic-invariant
# approximation from adult survival (s) and age at first breeding (alpha).
lambda_max_niel <- function(s, alpha) {
  disc <- (s - s * alpha - alpha - 1)^2 - 4 * s * alpha^2
  ((s * alpha - s + alpha + 1) + sqrt(pmax(disc, 0))) / (2 * alpha)
}

# PBR = 0.5 * Rmax * Fr * Nmin, with Rmax = lambda_max - 1.
pbr_from_components <- function(Nmin, Fr, lambda_max) {
  Rmax <- lambda_max - 1
  0.5 * Rmax * Fr * Nmin
}

# Conservative lower-bound local abundance (Nmin) from upper/lower bounds
# (NL, NU) via the lognormal formulation in Frick et al. (2026).
nmin_from_bounds <- function(NL, NU, ci_level = 0.60) {
  z_lower <- qnorm((1 - ci_level) / 2)
  z_upper <- qnorm(1 - (1 - ci_level) / 2)
  Nhat <- sqrt(NL * NU)
  CV <- sqrt(exp((log(NU / NL) / (2 * z_upper))^2) - 1)
  Nhat * exp(z_lower * CV)
}

# NL/NU from project area and a density proxy (bats/km^2), then Nmin.
nmin_from_area_density <- function(area_km2, density_low, density_high,
                                    scaling_factor = 2, ci_level = 0.60) {
  NL <- area_km2 * scaling_factor * density_low
  NU <- area_km2 * scaling_factor * density_high
  nmin_from_bounds(NL, NU, ci_level)
}

# Inverse of pbr_from_components(): Nmin required to reach a target PBR at
# a given Fr and lambda_max.
required_nmin <- function(pbr_target, Fr, lambda_max) {
  Rmax <- lambda_max - 1
  pbr_target / (0.5 * Rmax * Fr)
}

# Interaction area (km^2) that would produce a given Nmin, under a density
# proxy, found by inverting nmin_from_area_density() numerically.
implied_area <- function(nmin_target, density_low = 5.5, density_high = 18.2,
                          scaling_factor = 2, ci_level = 0.60) {
  f <- function(area) {
    nmin_from_area_density(area, density_low, density_high, scaling_factor, ci_level) - nmin_target
  }
  uniroot(f, c(1, 1e7))$root
}

# Density (bats/km^2) that would produce a given Nmin over a KNOWN area,
# holding the density_high:density_low ratio fixed at density_ratio (that
# ratio alone sets CV_Nhat in the lognormal formulation, so Nmin is exactly
# linear in density_low for a fixed area/ratio -- solved in closed form,
# not by root-finding). Use this when the interaction area is known (e.g.
# a project's real footprint) and the question is what density it implies,
# the mirror image of implied_area().
implied_density <- function(nmin_target, area_km2, density_ratio,
                             scaling_factor = 2, ci_level = 0.60) {
  z_lower <- qnorm((1 - ci_level) / 2)
  z_upper <- qnorm(1 - (1 - ci_level) / 2)
  CV <- sqrt(exp((log(density_ratio) / (2 * z_upper))^2) - 1)
  density_low <- nmin_target / (area_km2 * scaling_factor * sqrt(density_ratio) * exp(z_lower * CV))
  c(low = density_low, high = density_low * density_ratio)
}
