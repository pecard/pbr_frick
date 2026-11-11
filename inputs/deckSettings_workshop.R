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
## pbr_frick_technical_note_BSH_DGY.docx, Revision 17) --------------------
deck_facility_labels <- c("Site A", "Site B")

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
