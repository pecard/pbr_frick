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
