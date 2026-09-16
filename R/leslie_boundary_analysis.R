##
## Inverse demographic boundary analysis ("stress test" of the Leslie
## matrix): where do the PBR benchmarks (lambda_max = 1.20 / 1.24) sit
## relative to the space of biologically plausible V. murinus vital-rate
## combinations, using the SAME validated Dekker & Limpens (2024) 2-stage
## structure (R/leslie_dekker_limpens.R, lambda = 1.047 matching the
## published value exactly)?
##
## Framing (important, per Paulo's caution): lambda from this matrix is
## driven by OBSERVED/empirical vital rates (Safi 2006's mark-recapture
## estimates) -- it is lambda_observed, not lambda_max. PBR's lambda_max is
## a theoretical ceiling under favourable, non-resource-limited conditions.
## The two are not directly comparable, so "the matrix gives 1.047" does
## NOT by itself show that 1.20/1.24 is biologically impossible as a
## lambda_max. What this script does instead: characterise WHAT demographic
## combinations (adult survival, juvenile survival, breeding fraction,
## litter size) the validated structure would need to reach 1.20/1.24, and
## how those combinations sit relative to Safi's empirical baseline and to
## plausible literature ranges for this species -- i.e. how favourable /
## extreme "favourable" would have to be. No calibration to the published
## elasticities (0.530/0.235/0.235) is attempted here; the lambda match
## alone (1.0474 vs 1.047) is what validates the reconstruction.
##

source("R/leslie_dekker_limpens.R")

suppressPackageStartupMessages({ library(dplyr); library(tidyr); library(ggplot2) })

fig_dir <- "outputs/figures"
dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)

# ---- 1. Plausible ranges for each vital rate ------------------------------
# Baseline = Safi (2006) via Dekker & Limpens (2024): S_j=0.62, S_af=0.76,
# p_breed=0.87, litter=1.8 (suburban colonies).
# S_adult: same range used throughout the PBR analysis (s_range).
# p_breed, litter: bounded above by biology (p_breed <= 1; litter's upper
# literature figure for this species is the higher, urban-colony estimate
# from Zhigalin & Moskvitina 2017, ~2.7-2.9 pups) -- see
# references/leslie_matrix_parametrisation.md.
# S_juv: no second species-specific estimate was located; the range below
# is a plausibility bound (juvenile survival in vespertilionids is
# consistently below adult survival), not a literature point estimate --
# flagged as such throughout.
range_s_adult <- c(0.70, 0.95)
range_s_juv   <- c(0.40, 0.80)
range_p_breed <- c(0.70, 0.98)
range_litter  <- c(1.3, 2.9)

baseline <- list(s_juv = S_j, s_adult = S_af, p_breed = p_breed, litter = litter)

lambda_of <- function(s_juv, s_adult, p_breed, litter) {
  F_rec <- p_breed * litter * sex_ratio * s_juv
  leslie_lambda(build_dekker_stage_matrix(2, s_juv, s_adult, F_rec))
}

stopifnot(abs(lambda_of(baseline$s_juv, baseline$s_adult, baseline$p_breed, baseline$litter) - lambda_female) < 1e-9)

benchmarks <- c(1.20, 1.24)
benchmark_colours <- c("cyan", "chartreuse")

# ---- 2. Univariate breakeven: vary ONE parameter, others at baseline ------
solve_breakeven <- function(param, target, search_upper) {
  f <- function(x) {
    args <- baseline
    args[[param]] <- x
    do.call(lambda_of, args) - target
  }
  lo <- baseline[[param]]
  if (f(lo) >= 0) return(lo)          # baseline already meets target (not the case here)
  if (f(search_upper) < 0) return(NA) # unreachable even at the search ceiling
  uniroot(f, c(lo, search_upper))$root
}

search_ceiling <- list(s_adult = 0.999, s_juv = 0.99, p_breed = 1.0, litter = 6)
param_labels <- c(s_adult = "Adult female survival", s_juv = "Juvenile survival",
                   p_breed = "Breeding fraction", litter = "Litter size")

plausible_ranges <- list(s_adult = range_s_adult, s_juv = range_s_juv,
                          p_breed = range_p_breed, litter = range_litter)

breakeven_tbl <- tidyr::expand_grid(param = names(baseline), target = benchmarks) %>%
  rowwise() %>%
  mutate(
    baseline_value = baseline[[param]],
    required_value = solve_breakeven(param, target, search_ceiling[[param]]),
    plausible_lower = plausible_ranges[[param]][1],
    plausible_upper = plausible_ranges[[param]][2],
    pct_of_plausible_range = 100 * (required_value - plausible_lower) / (plausible_upper - plausible_lower),
    within_plausible_range = !is.na(required_value) & required_value <= plausible_upper
  ) %>%
  ungroup() %>%
  mutate(parameter = param_labels[param]) %>%
  select(parameter, target, baseline_value, required_value, plausible_upper, pct_of_plausible_range, within_plausible_range)

cat("=== Univariate breakeven: value needed to reach lambda target, holding all other vital rates at Safi's baseline ===\n")
print(as.data.frame(breakeven_tbl), digits = 3)

# ---- 3. Joint plausible-space grid: where does the baseline lambda sit? ---
# Full factorial over all 4 vital rates across their plausible ranges
# (n_stages fixed at 2, the validated -- and fastest-maturing, hence most
# favourable to reaching high lambda -- structure).
grid_res <- 12
full_grid <- tidyr::expand_grid(
  s_adult = seq(range_s_adult[1], range_s_adult[2], length.out = grid_res),
  s_juv   = seq(range_s_juv[1], range_s_juv[2], length.out = grid_res),
  p_breed = seq(range_p_breed[1], range_p_breed[2], length.out = grid_res),
  litter  = seq(range_litter[1], range_litter[2], length.out = grid_res)
) %>%
  mutate(
    lambda = mapply(lambda_of, s_juv, s_adult, p_breed, litter),
    band = cut(lambda, breaks = c(-Inf, 1.20, 1.24, Inf),
               labels = c("< 1.20", "1.20-1.24", ">= 1.24"), right = FALSE)
  )

band_shares <- full_grid %>% count(band) %>% mutate(pct = 100 * n / sum(n))
pct_below_baseline <- 100 * mean(full_grid$lambda < lambda_female)

cat(sprintf(
  "\n=== Joint plausible-space grid (%d combinations, S_adult x S_juv x p_breed x litter, all within literature-plausible ranges) ===\n",
  nrow(full_grid)
))
print(band_shares)
cat(sprintf(
  "\nSafi's empirical baseline (lambda = %.3f) exceeds %.1f%% of this plausible-space grid.\n",
  lambda_female, pct_below_baseline
))
cat(sprintf(
  "Share of the plausible space reaching lambda >= 1.20: %.2f%%; reaching lambda >= 1.24: %.2f%%.\n",
  100 * mean(full_grid$lambda >= 1.20), 100 * mean(full_grid$lambda >= 1.24)
))

# ---- 4. Two-panel surface: baseline vitals vs favourable vitals -----------
# Panel A: S_juv, p_breed held at Safi's baseline (only S_adult, litter free)
# Panel B: S_juv, p_breed held at their plausible UPPER bound (a best-case
# assumption on the two parameters not shown), to see how much the boundary
# moves when the "other" vital rates are also pushed favourably.
surface_grid <- tidyr::expand_grid(
  s_adult = seq(range_s_adult[1], range_s_adult[2], length.out = 80),
  litter  = seq(range_litter[1], range_litter[2], length.out = 80),
  scenario = c("Baseline S_juv & breeding fraction", "Favourable S_juv & breeding fraction")
) %>%
  mutate(
    s_juv_use   = ifelse(scenario == "Baseline S_juv & breeding fraction", baseline$s_juv, range_s_juv[2]),
    p_breed_use = ifelse(scenario == "Baseline S_juv & breeding fraction", baseline$p_breed, range_p_breed[2]),
    lambda = mapply(lambda_of, s_juv_use, s_adult, p_breed_use, litter)
  )

baseline_point <- tibble::tibble(s_adult = baseline$s_adult, litter = baseline$litter,
                                  scenario = "Baseline S_juv & breeding fraction")

fig_boundary_surface <- file.path(fig_dir, "leslie_boundary_surface.png")
ggsave(fig_boundary_surface, width = 10, height = 4.8, dpi = 150, plot = {
  ggplot(surface_grid, aes(x = s_adult, y = litter, z = lambda)) +
    geom_raster(aes(fill = lambda)) +
    geom_contour(breaks = benchmarks, aes(colour = after_stat(factor(level))), linewidth = 0.8) +
    scale_colour_manual(name = "lambda contour", values = setNames(benchmark_colours, benchmarks)) +
    scale_fill_viridis_c(option = "D", name = "lambda") +
    geom_point(data = baseline_point, aes(x = s_adult, y = litter), inherit.aes = FALSE,
               colour = "red", size = 2.5, shape = 17) +
    facet_wrap(~scenario) +
    labs(
      x = "Adult female survival", y = "Litter size",
      title = "Where would adult survival and litter size need to sit to reach the PBR benchmarks?",
      subtitle = paste0(
        "Red triangle: Safi's empirical baseline (S_adult=", baseline$s_adult, ", litter=", baseline$litter,
        "); contours: PBR benchmarks 1.20 (cyan) / 1.24 (green)"
      )
    ) +
    theme_minimal() +
    theme(plot.subtitle = element_text(size = 8), strip.text = element_text(face = "bold"))
})

cat("\nWrote", fig_boundary_surface, "\n")
