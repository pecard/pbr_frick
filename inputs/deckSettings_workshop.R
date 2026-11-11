##
## Settings for the European-specialist workshop deck (Paulo, 2026-11):
## "Application of a new PBR framework and how it can be translated into
## operational decisions for bat curtailment". Mirrors the
## inputs/pbrSettings_*.R convention -- every literal the deck's figures
## need lives here, so re-styling or re-numbering never touches
## R/workshop_deck_figures.R or report/pbr_workshop_slides.Rmd.
##
## All facility identifiers and the host country are anonymised, and every
## derived/uncertain quantity is rounded DOWN (floor); only the imposed
## thresholds themselves are left exact (Paulo, 2026-11: "um ceiling
## apenas dos thresholds, nada mais" -- they are already clean integers,
## so ceiling has no visible effect here). Turbine counts and % reduction
## are reported as exact operational facts, not rounded, since flooring a
## small integer count would distort the ratio it is part of.
##

## ---- Anonymised facility labels & real operational facts (source:
## pbr_frick_technical_note_BSH_DGY.docx, Revision 17). Fixed mapping
## (Paulo, 2026-11): Project 1 = Bash WPP = 112 bats/yr; Project 2 =
## Djangeldy WPP = 88 bats/yr -- used consistently across every output
## (deck, Shiny simulator), independent of inputs/pbrSettings_BSH_DGY.R's
## own facility_labels, which Paulo toggles to real names for the actual
## client/lender technical note. --------------------------------------
deck_facility_labels <- c("Project 1", "Project 2")

deck_thresholds <- tibble::tibble(
  site      = deck_facility_labels,
  threshold = c(112, 88)                 # bats/year, imposed -- exact, no rounding
)

deck_turbines <- tibble::tibble(
  site       = deck_facility_labels,
  assessed   = c(68, 27),
  curtailed  = c(17, 7)
)

deck_cut_in_speed_ms <- 6.5

deck_reduction_pct <- tibble::tibble(
  site = deck_facility_labels,
  pct_reduction_vs_counterfactual = c(88, 53)   # % vs 2025 no-curtailment reference
)

## ---- Chapter 7 rationale (source: pbr_frick_technical_note_BSH_DGY.docx,
## Revision 17, Sections 7.1-7.4) -- replaces the old slides 5-6 with the
## fuller operational story: why an annual limit alone can't drive an
## in-season response, the 2026 mortality position read two ways (year-
## to-date vs. year-over-year), whether a simple turbine ranking recovers
## the real lender-agreed selection, and why that check is what drove this
## note to abandon a single shared GenEst correction factor for two
## project-specific ones. Large mortality/carcass counts floored to the
## nearest 10 (Paulo's rounding rule); percentages, ratios (floored to 1
## decimal) and small turbine counts left exact, as elsewhere in this file.

## 2026 year-to-date V. murinus mortality (through ISO week 39) against
## the annual PBR reference (Section 7.2's first table).
deck_ytd_status <- tibble::tibble(
  site             = deck_facility_labels,
  annual_reference = c(112, 88),
  raw_carcasses    = c(80, 30),
  genest_corrected = c(380, 240),
  pct_of_reference = c(344, 275)
)

## Year-over-year comparison, GenEst-corrected mortality through the same
## ISO week (39) in each year (Section 7.2's second table).
deck_yoy <- tibble::tibble(
  site        = deck_facility_labels,
  y2025       = c(940, 250),
  y2026       = c(380, 240),
  pct_change  = c(-59, -5)
)

## Turbine-level retrospective validation: independent replication (raw
## carcasses + this note's own working correction factor) against the
## official, lender-agreed selection (Section 7.3's first table). Counts
## left exact -- small integers, a 50%-cumulative-mortality selection rule.
deck_turbine_validation <- tibble::tibble(
  site                  = deck_facility_labels,
  officially_selected   = c(17, 7),
  independently_selected = c(18, 8),
  selected_by_both      = c(14, 5),
  jaccard               = c(0.67, 0.50),
  pct_official_recovered = c(82, 71),
  pct_mortality_captured = c(87, 76)
)

## Implied GenEst correction ratio (official GenEst estimate / matching raw
## carcass count), per turbine vs. the aggregate (mortality-weighted)
## ratio this note actually uses -- Section 7.3's second table, the
## finding that drove abandoning a single flat/combined factor. Ratios
## floored to 1 decimal.
deck_genest_ratio <- tibble::tibble(
  site            = deck_facility_labels,
  median_per_turbine = c(4.0, 5.7),
  min_per_turbine     = c(2.8, 3.9),
  max_per_turbine     = c(9.7, 12.9),
  aggregate_ratio     = c(4.6, 6.2)
)
deck_genest_rejected_flat   <- 4.0   # the flat factor rejected after this check
deck_genest_combined_median <- 4.3   # combined per-turbine median across both projects

## ---- Frick/PBR reverse-engineering facts (same source) -----------------
deck_lambda_max_benchmarks <- c(1.24, 1.20)

## Elasticity of lambda_max to each Niel & Lebreton parameter, at the
## explored range's lower/upper end (s in [0.70,0.95], alpha in [1.5,3.5]
## years) -- % change in lambda_max per 1% change in the parameter.
deck_elasticity <- tibble::tibble(
  parameter = c("Adult survival (s)", "Age at first breeding (α)"),
  low       = c(-0.13, -0.06),
  high      = c(1.35, 0.22)
)

## Nmin required to reproduce each imposed threshold, under Fr=0.3 and
## each fixed lambda_max benchmark -- floored to the nearest 100.
deck_nmin_range <- tibble::tibble(
  site                 = rep(deck_facility_labels, each = 2),
  lambda_max           = rep(deck_lambda_max_benchmarks, times = 2),
  nmin_required_floor  = c(3700, 3100, 2900, 2400)
)

## ---- Exploratory carrying-capacity comparison (this session's
## prototype, R/carrying_capacity_reference.R -- NOT part of the reviewed
## technical note; flag as exploratory wherever shown) --------------------
deck_cc_comparison <- tibble::tibble(
  scenario = c("Realistic demographic rate", "Literature-favourable survival ceiling"),
  H        = c(20, 100)   # floored
)

## ---- Plot controls -- the knobs Paulo asked for, so re-styling never
## means editing a plotting function. Every function in
## R/workshop_deck_figures.R reads ONLY from this list. --------------------
deck_plot_controls <- list(
  base_size     = 13,     # ggplot theme base font size
  title_size    = 15,     # in-figure title (not the slide title)
  subtitle_show = TRUE,   # flip to FALSE to drop every figure subtitle at once
  caption_show  = TRUE,
  point_size    = 2.6,
  line_width    = 1.0,
  bar_width     = 0.62,
  label_size    = 3.6,    # geom_text/data-label size (ggplot units, ~ pt/2.845)
  axis_text_size_mult = 1.0,  # multiplier on base_size for axis text

  colour_primary   = "#1C7293",  # teal -- main series
  colour_secondary = "#6CA6C1",  # sky -- secondary series
  colour_dark      = "#0B2545",  # navy -- "current policy" / emphasis
  colour_accent    = "#E8A33D",  # amber -- single highlight accent
  colour_muted     = "#4C4C4C",  # captions / reference lines

  dpi         = 300,
  fig_width   = 8.0,
  fig_height  = 5.0,
  bg          = "white"
)
