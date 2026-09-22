##
## PBR analysis settings -- Paris Expert Workshop illustrative scenario.
##
## A single fictional wind farm, generic model species ("Bat model"),
## and independently-chosen scale parameters -- NOT derived from the real
## BSH/DGY numbers by any shared multiplier (no "real x 0.5" anywhere
## below). The point of this file is that nothing in it can be inverted
## back to the real client numbers even if both settings files are seen
## side by side. Threshold/turbine-count magnitudes per Paulo (2026-09):
## threshold 50 bats/yr, 50 WTG, 10 curtailed -- kept close to the real
## orders of magnitude, everything else (correction factor, footprint,
## mortality levels) chosen fresh.
##
## Literature-sourced values (lambda_max benchmarks, density proxy,
## Leslie vital rates, demographic ranges) are NOT client data and are
## kept as in pbrSettings_BSH_DGY.R -- they describe a real, published
## vesper-bat dataset, reused here under a generic "Bat model" label
## rather than replaced, per Paulo's instruction. The report text must
## not name the source species or cite the original papers (Safi 2006,
## Dekker & Limpens 2024, Zhigalin & Moskvitina 2017) -- those citations
## stay in code comments only, never in rendered output, so the
## parametrisation cannot be traced back to V. murinus specifically.
##

## Single fictional facility -- no real name, no real geography.
facility_labels       <- "Wind Farm A"
project_ref           <- facility_labels
species_name          <- "Bat model"
species_common_name   <- "Bat model"

## Imposed fatality threshold under review (bats/year). Chosen directly
## for this illustrative scenario -- not backsolved from an undocumented
## real value, so the "reverse-engineering the imposed threshold" framing
## from the real report does not apply here and is dropped from this
## scenario's document.
pbr_thresholds <- tibble::tibble(
  facility  = facility_labels,
  threshold = 50
)

## Turbine universe for the illustrative operational scenario (Section 7
## equivalent). Defined here (not backed by any real official selection
## file) since the whole Section 7 dataset for this scenario is
## simulated -- see R/simulate_paris_mortality_data.R.
n_turbines_total     <- 50
n_turbines_curtailed <- 10   # target -- the simulation calibrates turbine
                              # -level risk concentration so the same
                              # 50%-cumulative-mortality selection rule
                              # used in the real report lands close to
                              # this count, rather than hard-coding it

## IUCN status used for the corrected Fr throughout. run_pbr_analysis()
## still computes a backsolve step against an "original guess" status --
## kept only so that code path runs unchanged; the trimmed workshop
## document does not narrate a reverse-engineering story (there is no
## real, undocumented threshold to reverse-engineer here since 50 was
## chosen directly), so this value never surfaces in rendered text.
iucn_status_original_guess <- "Near Threatened"
iucn_status_corrected      <- "Least Concern"

## Fixed lambda_max benchmarks (Frick et al. 2026) -- literature values,
## not client data, unchanged from the real scenario.
lambda_max_benchmarks <- c(1.24, 1.20)

## Illustrative geography for the Nmin plausibility check -- chosen
## independently (not a fraction of the real 370 km2 / 7000 km2 figures).
project_footprint_km2 <- 240
combined_region_km2   <- 3000
interproject_distance_km <- NA_real_  # single facility -- not meaningful;
                                       # kept only so run_pbr_analysis()'s
                                       # return list is unchanged, not
                                       # referenced by the trimmed template

## Demographic parameter ranges -- literature-based, unchanged.
s_range     <- c(0.70, 0.95)
alpha_range <- c(1.5, 3.5)
alpha_range_prior <- c(2, 4)  # kept for the before/after comparison table's
                               # code path; not narrated in the trimmed doc

## Density proxy (bats/km2) for the Nmin plausibility check -- literature
## value (Pipistrellus pipistrellus, per Frick et al. 2026), unchanged.
density_proxy_low  <- 5.5
density_proxy_high <- 18.2

## Fr scenarios for the response-surface comparison figure.
fr_scenarios <- c(0.1, 0.3, 0.5, 1.0)

## Monte Carlo simulation controls.
mc_n_sim <- 20000
mc_seed  <- 1

## Reference document for Word styling only (no content) -- reused as-is.
reference_docx_path <- "report/reference_template.docx"

##
## Leslie-matrix / PVA-lite cross-check -- literature-sourced vital rates
## for a real vesper-bat dataset, reused under the generic "Bat model"
## label (see file header). Values unchanged from pbrSettings_BSH_DGY.R.
##
leslie_s_juv     <- 0.62
leslie_s_adult_f <- 0.76
leslie_s_adult_m <- 0.42
leslie_p_breed   <- 0.87
leslie_litter    <- 1.8
leslie_sex_ratio <- 0.5

leslie_stage_range <- 2:5

pva_n_years <- 25
pva_n_reps  <- 1000
pva_vital_rate_cv <- 0.10
pva_seed <- 42
pva_quasi_extinction_fraction <- 0.10

leslie_boundary_range_s_adult <- s_range
leslie_boundary_range_s_juv   <- c(0.40, 0.80)
leslie_boundary_range_p_breed <- c(0.70, 0.98)
leslie_boundary_range_litter  <- c(1.3, 2.9)
leslie_boundary_grid_resolution <- 12

##
## Illustrative operational-scenario parameters (Section 7 equivalent) --
## entirely new, not derived from the real correction factors (4.64x /
## 6.21x) or from real mortality counts by any shared multiplier.
##
genest_correction_factor_paris <- 1.3
baseline_period_label <- "Baseline"
response_period_label <- "Response"
