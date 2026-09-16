##
## PBR analysis settings -- BSH & DGY, Parti-coloured bat (Vespertilio murinus)
##
## Sourced by run_pbr_report_BSH_DGY.R via R/pbr_analysis_pipeline.R.
## Holds every literal the analysis needs (R/pbr_analysis.R) so a new
## scenario -- a different species, facility set, threshold values, or
## parameter ranges -- can be run by copying this file to a new
## inputs/pbrSettings_<scenario>.R and editing the values below, without
## touching the analysis or report code.
##

project_ref          <- "BSH & DGY"
species_name          <- "Vespertilio murinus"
species_common_name   <- "Parti-coloured bat"

## Imposed fatality thresholds under review -- one row per facility
## (bats/year). No documented derivation is available for these values;
## the analysis works backwards from them (see R/pbr_analysis.R).
pbr_thresholds <- tibble::tibble(
  facility  = c("BSH", "DGY"),
  threshold = c(144, 120)
)

## IUCN status assumed to back-solve the ORIGINAL (likely, undocumented)
## parametrisation behind the imposed thresholds.
iucn_status_original_guess <- "Near Threatened"

## Corrected/defensible IUCN status for this species -- drives Fr for
## every downstream analysis step after the reverse-engineering backsolve.
## The species is not a threatened category, so Least Concern (declining
## or unknown trend) applies, not Near Threatened.
iucn_status_corrected <- "Least Concern"

## Fixed lambda_max benchmarks from Frick et al. (2026) -- used only for
## the reverse-engineering backsolve step (section "Reverse-engineering
## the imposed thresholds"), not for the Monte Carlo / response surfaces,
## which treat lambda_max as derived from the ranges below.
lambda_max_benchmarks <- c(1.24, 1.20)

## Real, known project geography (Paulo, 2026-09) -- used for the Nmin
## plausibility check in place of a back-solved area: BSH and DGY are
## ~100 km apart; each project's own footprint is ~370 km2; the combined
## region spanning both (a landscape-scale reading of "local population")
## is ~7000 km2.
project_footprint_km2    <- 370
combined_region_km2      <- 7000
interproject_distance_km <- 100

## Demographic parameter ranges for the Monte Carlo simulation, response
## surfaces and elasticity analysis.
s_range     <- c(0.70, 0.95)  # adult survival
alpha_range <- c(1.5, 3.5)    # age at first breeding (years) -- literature
                               # -supported range; see
                               # references/alpha_first_breeding_evidence.md
alpha_range_prior <- c(2, 4)  # earlier, wider range -- kept only for the
                               # before/after comparison table

## Default density proxy (bats/km2) for the Nmin plausibility check --
## Pipistrellus pipistrellus, per Frick et al. (2026), used in the absence
## of species-specific density data. The check compares this against the
## density that Nmin_assumed would imply over the real project areas
## above (project_footprint_km2, combined_region_km2), not the other way
## around -- see R/pbr_analysis.R.
density_proxy_low  <- 5.5
density_proxy_high <- 18.2

## Fr scenarios shown in the response-surface comparison figure (one panel
## per value): Vulnerable, Near Threatened, Least Concern
## (declining/unknown), Least Concern (stable/increasing).
fr_scenarios <- c(0.1, 0.3, 0.5, 1.0)

## Monte Carlo simulation controls.
mc_n_sim <- 20000
mc_seed  <- 1

## Reference document whose styles (headings, tables, header/footer, page
## numbering) the rendered report reuses -- see R/pbr_report.R.
reference_docx_path <- "report/reference_template.docx"

##
## Leslie-matrix / PVA-lite cross-check (R/pbr_analysis.R section 5b-5c;
## see references/leslie_matrix_parametrisation.md for the literature and
## R/leslie_dekker_limpens.R for the validation this reuses).
##

## Species-specific vital rates (Safi 2006, as used by Dekker & Limpens
## 2024's published V. murinus model, which this analysis independently
## reproduces -- lambda = 1.047, matching exactly).
leslie_s_juv    <- 0.62   # juvenile survival, both sexes
leslie_s_adult_f <- 0.76  # adult female survival
leslie_s_adult_m <- 0.42  # adult male survival (not used in the female
                          # -only projection, kept for reference)
leslie_p_breed  <- 0.87   # fraction of adult females reproducing
leslie_litter   <- 1.8    # litter size, suburban colonies (Zhigalin & Moskvitina 2017)
leslie_sex_ratio <- 0.5

## Maturation-delay sensitivity: number of pre-reproductive stages tested
## (2 = the validated published structure), same adult-survival range as
## the PBR s x alpha analysis for direct comparability.
leslie_stage_range <- 2:5

## PVA-lite: 25-year stochastic projection controls.
pva_n_years <- 25
pva_n_reps  <- 1000
pva_vital_rate_cv <- 0.10   # working assumption -- no inter-annual
                            # variance estimate for V. murinus was found
pva_seed <- 42
pva_quasi_extinction_fraction <- 0.10  # of N0, a common PVA convention

## Inverse demographic boundary analysis ("stress test"): plausible ranges
## for the validated Leslie matrix's four vital rates, used to find which
## combinations reach the PBR benchmarks (1.20/1.24) and how they compare
## with Safi's (2006) empirical baseline. S_adult reuses s_range (same
## range used throughout the PBR/Niel-Lebreton analysis, for direct
## comparability). S_juv, p_breed and litter have no second
## species-specific estimate located this session; their ranges are
## plausibility bounds -- not literature point estimates -- documented in
## references/leslie_matrix_parametrisation.md. p_breed is capped below 1
## (a proportion); litter's upper bound is Zhigalin & Moskvitina's (2017)
## urban-colony estimate (2.7-2.9), the higher of the two figures reported.
leslie_boundary_range_s_adult <- s_range
leslie_boundary_range_s_juv   <- c(0.40, 0.80)
leslie_boundary_range_p_breed <- c(0.70, 0.98)
leslie_boundary_range_litter  <- c(1.3, 2.9)
leslie_boundary_grid_resolution <- 12  # per dimension, full 4-way factorial
