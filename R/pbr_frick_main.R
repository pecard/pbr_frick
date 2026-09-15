# ============================================================
# PBR for bats at wind facilities (Frick et al. 2026)
# link: https://besjournals.onlinelibrary.wiley.com/doi/10.1002/2688-8319.70189
# Implements:
#  (1) "Classical" PBR structure: 0.5 * (lambda_max - 1) * Fr * Nmin
#  (2) Paper's project-scale Nmin using NL/NU uncertainty formulation
#  (3) "Adjusted PBR" variant: same as (2) but with lambda_max forced
#      to the authors' recommended benchmark value (e.g., 1.24 or 1.20)
# ============================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
})

# ---- Recovery factor (Fr) mapping from IUCN status (Table 1) ----
# Least concern 0.5; Near threatened 0.3; Vulnerable 0.1;
# Endangered 0; Critically endangered 0; Data deficient 0.2
fr_from_iucn <- function(iucn_status) {
  status <- tolower(iucn_status)
  dplyr::case_when(
    status %in% c("least concern", "lc") ~ 0.5,
    status %in% c("near threatened", "nt") ~ 0.3,
    status %in% c("vulnerable", "vu") ~ 0.1,
    status %in% c("endangered", "en") ~ 0.0,
    status %in% c("critically endangered", "cr") ~ 0.0,
    status %in% c("data deficient", "dd") ~ 0.2,
    TRUE ~ NA_real_
  )
}

# ---- lambda_max approximation (Niel & Lebreton demographic invariant) ----
# lambda_max ≈ ((s*alpha - s + alpha + 1) + sqrt((s - s*alpha - alpha - 1)^2 - 4*s*alpha^2)) / (2*alpha)
lambda_max_niel <- function(s, alpha) {
  stopifnot(all(alpha > 0), all(s > 0), all(s < 1))
  disc <- (s - s * alpha - alpha - 1)^2 - 4 * s * alpha^2
  if (any(disc < 0)) {
    warning("Discriminant < 0 for some rows; returning NA for those.")
  }
  out <- ((s * alpha - s + alpha + 1) + sqrt(pmax(disc, 0))) / (2 * alpha)
  out
}

# ---- Convert lambda_max to Rmax ----
rmax_from_lambda <- function(lambda_max) {
  lambda_max - 1
}

# ---- Local abundance bounds from area and density proxy ----
# Nlocal-U = area_km2 * 2 * upper_density
# Nlocal-L = area_km2 * 2 * lower_density
n_bounds_from_area_density <- function(area_km2, density_low, density_high, scaling_factor = 2) {
  NL <- area_km2 * scaling_factor * density_low
  NU <- area_km2 * scaling_factor * density_high
  tibble::tibble(NL = NL, NU = NU)
}

# ---- CV_Nhat from NL/NU and CI level (paper uses 60% CI) ----
# CV_Nhat = sqrt( exp( ( ln(NU/NL) / (2 * Z_{1-alpha/2}) )^2 ) - 1 )
# where alpha is the *CI alpha*, not age-at-first-breeding.
cv_nhat_lognormal <- function(NL, NU, ci_level = 0.60) {
  stopifnot(ci_level > 0, ci_level < 1)
  z_upper <- qnorm(1 - (1 - ci_level) / 2)  # e.g., 0.8 for 60% CI
  sqrt(exp((log(NU / NL) / (2 * z_upper))^2) - 1)
}

# ---- Nhat and Nmin (lower bound) ----
# Nhat = sqrt(NL * NU)
# Nmin = Nhat * exp(Z_lower * CV_Nhat)
# For a 60% CI, lower tail = 0.2 => Z_lower = qnorm(0.2) (~ -0.8416)
nmin_from_bounds <- function(NL, NU, ci_level = 0.60) {
  z_lower <- qnorm((1 - ci_level) / 2)      # e.g., 0.2 for 60% CI
  Nhat <- sqrt(NL * NU)
  CV <- cv_nhat_lognormal(NL, NU, ci_level = ci_level)
  Nmin <- Nhat * exp(z_lower * CV)
  tibble::tibble(Nhat = Nhat, CV_Nhat = CV, z_lower = z_lower, Nmin = Nmin)
}

# ---- PBR calculator (core) ----
# PBR = 0.5 * Rmax * Fr * Nmin
pbr_from_components <- function(Nmin, Fr, lambda_max) {
  Rmax <- rmax_from_lambda(lambda_max)
  0.5 * Rmax * Fr * Nmin
}

# ---- Wrapper 1: "Classical" PBR as described using demographic-invariant lambda_max ----
# Inputs: s, alpha (age at first breeding), Fr, and either Nmin or (NL, NU)
pbr_classical <- function(Fr, s, alpha, NL, NU, ci_level = 0.60) {
  lambda_max <- lambda_max_niel(s = s, alpha = alpha)
  nmin_tbl <- nmin_from_bounds(NL = NL, NU = NU, ci_level = ci_level)
  PBR <- pbr_from_components(Nmin = nmin_tbl$Nmin, Fr = Fr, lambda_max = lambda_max)
  dplyr::bind_cols(
    tibble::tibble(lambda_max = lambda_max, Rmax = rmax_from_lambda(lambda_max), Fr = Fr),
    nmin_tbl,
    tibble::tibble(PBR = PBR)
  )
}

# ---- Wrapper 2: "Adjusted PBR" variant forcing lambda_max to authors' benchmark ----
# Authors suggest: lambda_max = 1.24 (alpha=1.5,s=0.9) and 1.20 if alpha=2 for some families.
pbr_adjusted_lambda <- function(Fr, NL, NU, lambda_max_fixed = 1.24, ci_level = 0.60) {
  nmin_tbl <- nmin_from_bounds(NL = NL, NU = NU, ci_level = ci_level)
  PBR <- pbr_from_components(Nmin = nmin_tbl$Nmin, Fr = Fr, lambda_max = lambda_max_fixed)
  dplyr::bind_cols(
    tibble::tibble(lambda_max = lambda_max_fixed, Rmax = rmax_from_lambda(lambda_max_fixed), Fr = Fr),
    nmin_tbl,
    tibble::tibble(PBR = PBR)
  )
}

# ============================================================
# Simulated dataset for a hypothetical bat species & wind facility
# ============================================================

set.seed(42)

sim_species <- tibble::tibble(
  species_id      = "Hypothetic bat",
  iucn_status     = "Near threatened",
  Fr              = fr_from_iucn(iucn_status),
  # facility / local scale inputs
  facility_id     = "WF_001",
  area_km2        = 120,
  density_low     = 3.0,
  density_high    = 12.0,
  scaling_factor  = 2,
  # demographic inputs for "classical" lambda_max approximation
  s               = 0.88,
  alpha_breed     = 1.5
) %>%
  dplyr::bind_cols(
    n_bounds_from_area_density(
      area_km2 = .$area_km2,
      density_low = .$density_low,
      density_high = .$density_high,
      scaling_factor = .$scaling_factor
    )
  )

sim_species
#> A tibble: 1 × ... (includes NL/NU)

# ---- Run both metrics on the simulated dataset ----
results <- sim_species %>%
  rowwise() %>%
  mutate(
    classical = list(pbr_classical(
      Fr = Fr, s = s, alpha = alpha_breed,
      NL = NL, NU = NU, ci_level = 0.60
    )),
    adjusted_124 = list(pbr_adjusted_lambda(
      Fr = Fr, NL = NL, NU = NU,
      lambda_max_fixed = 1.24, ci_level = 0.60
    )),
    adjusted_120 = list(pbr_adjusted_lambda(
      Fr = Fr, NL = NL, NU = NU,
      lambda_max_fixed = 1.20, ci_level = 0.60
    ))
  ) %>%
  ungroup() %>%
  select(species_id, facility_id, iucn_status, Fr, s, alpha_breed, area_km2, density_low, density_high, NL, NU,
         classical, adjusted_124, adjusted_120) %>%
  tidyr::unnest_wider(classical, names_sep = "_") %>%
  tidyr::unnest_wider(adjusted_124, names_sep = "_") %>%
  tidyr::unnest_wider(adjusted_120, names_sep = "_")

results %>%
  select(
    species_id, facility_id, iucn_status, Fr,
    NL, NU,
    classical_lambda_max, classical_PBR,
    adjusted_124_lambda_max, adjusted_124_PBR,
    adjusted_120_lambda_max, adjusted_120_PBR
  )


#" Testing elastiticy

# ---- lambda_max approximation from the paper (Niel & Lebreton style) ----
lambda_max_niel <- function(s, alpha) {
  disc <- (s - s * alpha - alpha - 1)^2 - 4 * s * alpha^2
  # guard numerical issues
  out <- ((s * alpha - s + alpha + 1) + sqrt(pmax(disc, 0))) / (2 * alpha)
  out[disc < 0] <- NA_real_
  out
}

# ---- numerical partial derivatives (central difference) ----
partials_lambda <- function(s, alpha, eps_s = 1e-5, eps_a = 1e-5) {
  lam <- lambda_max_niel(s, alpha)
  dlds <- (lambda_max_niel(s + eps_s, alpha) - lambda_max_niel(s - eps_s, alpha)) / (2 * eps_s)
  dlda <- (lambda_max_niel(s, alpha + eps_a) - lambda_max_niel(s, alpha - eps_a)) / (2 * eps_a)
  
  tibble(
    lambda_max = lam,
    d_lambda_ds = dlds,
    d_lambda_dalpha = dlda,
    E_s = dlds * s / lam,
    E_alpha = dlda * alpha / lam
  )
}

# ---- build a grid over your proposed ranges ----
s_seq <- seq(0.80, 0.95, by = 0.01)
a_seq <- seq(1.0,  3.0,  by = 0.5)

sens_grid <- tidyr::expand_grid(
  s = s_seq,
  alpha = a_seq
) %>%
  rowwise() %>%
  mutate(out = list(partials_lambda(s, alpha))) %>%
  ungroup() %>%
  tidyr::unnest(out)

# ---- summary across the whole range ----
sens_summary <- sens_grid %>%
  summarise(
    lambda_min = min(lambda_max, na.rm = TRUE),
    lambda_max = max(lambda_max, na.rm = TRUE),
    
    E_s_min = min(E_s, na.rm = TRUE),
    E_s_max = max(E_s, na.rm = TRUE),
    
    E_alpha_min = min(E_alpha, na.rm = TRUE),
    E_alpha_max = max(E_alpha, na.rm = TRUE),
    
    # optional: typical magnitudes (medians)
    E_s_median = median(E_s, na.rm = TRUE),
    E_alpha_median = median(E_alpha, na.rm = TRUE)
  )

sens_summary

library(dplyr)
library(tidyr)

# Summarise elasticity ranges across the (s, alpha) grid
# Elasticity = proportional change in 𝜆max for a 1% change in the parameter, 
# holding the other parameter constant.


elasticity_ranges <- sens_grid %>%
  summarise(
    E_s_min        = min(E_s, na.rm = TRUE),
    E_s_max        = max(E_s, na.rm = TRUE),
    E_s_median     = median(E_s, na.rm = TRUE),
    E_alpha_min    = min(E_alpha, na.rm = TRUE),
    E_alpha_max    = max(E_alpha, na.rm = TRUE),
    E_alpha_median = median(E_alpha, na.rm = TRUE)
  ) %>%
  tidyr::pivot_longer(
    cols = everything(),
    names_to = c("Parameter", "Statistic"),
    names_pattern = "E_(s|alpha)_(min|max|median)",
    values_to = "Value"
  ) %>%
  mutate(
    Parameter = recode(Parameter,
                       s = "Adult survival (s)",
                       alpha = "Age at first breeding (α)")
  ) %>%
  tidyr::pivot_wider(
    names_from = Statistic,
    values_from = Value
  ) %>%
  select(Parameter, min, max, median) %>%
  arrange(match(Parameter, c("Adult survival (s)", "Age at first breeding (α)")))

elasticity_ranges


# PLOTS

# Heatmap of λmax(s, α)
library(ggplot2)

ggplot(sens_grid, aes(x = s, y = alpha)) +
  geom_raster(
    aes(fill = lambda_max),
    interpolate = TRUE
  ) +
  geom_contour(
    aes(z = lambda_max),
    breaks = c(1.2, 1.24),
    colour = "white",
    linewidth = 0.4
  ) +
  scale_fill_viridis_c(
    name = expression(lambda[max]),
    option = "C"
  ) +
  labs(
    x = "Adult survival (s)",
    y = "Age at first breeding (α)",
    title = expression("Maximum population growth rate " * lambda[max]),
    subtitle = "Contours show λmax = 1.20 and 1.24"
  ) +
  theme_minimal()


# Heatmap of elasticity to survival Es
ggplot(sens_grid, aes(x = s, y = alpha, fill = E_s)) +
  geom_raster(interpolate = TRUE) +
  scale_fill_viridis_c(
    name = expression(E[s]),
    option = "B",
    direction = -1
  ) +
  labs(
    x = "Adult survival (s)",
    y = "Age at first breeding (α)",
    title = expression("Elasticity of " * lambda[max] * " to adult survival")
  ) +
  theme_minimal()

# Heatmap of elasticity to age at first breeding
ggplot(sens_grid, aes(x = s, y = alpha, fill = E_alpha)) +
  geom_raster(interpolate = TRUE) +
  scale_fill_viridis_c(
    name = expression(E[alpha]),
    option = "A",
    direction = -1
  ) +
  labs(
    x = "Adult survival (s)",
    y = "Age at first breeding (α)",
    title = expression("Elasticity of " * lambda[max] * " to age at first breeding")
  ) +
  theme_minimal()
