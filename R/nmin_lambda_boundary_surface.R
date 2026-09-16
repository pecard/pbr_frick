##
## Separates the two uncertainties that were partly conflated so far:
## demography (lambda_max, now well characterised via the validated Leslie
## matrix stress test, R/leslie_boundary_analysis.R) vs. exposed population
## size (Nmin, essentially unconstrained -- no population data exist for
## this region/species this session, only weekly mortality counts at both
## facilities). PBR = 0.5*(lambda-1)*Fr*Nmin is LINEAR in Nmin, so this is
## now plausibly the dominant source of uncertainty in the whole framework.
##
## Two outputs:
## 1. The PBR(N, lambda) surface itself, with PBR contours and reference
##    lambda values (Leslie/Safi's own 1.047; the paper's 1.20/1.24
##    benchmarks) against a range of Nmin.
## 2. A diagnostic (NOT an inverse population estimate -- explicitly not
##    that, see caveat below): observed mortality as a fraction of Nmin,
##    across a range of assumed Nmin, against reference "sustainable
##    fraction" lines from the Leslie matrix's own deterministic surplus
##    and from the PBR conservative/relaxed bounds. This does not, and
##    cannot, back out Nmin from mortality alone (that would require
##    P(encounter) and P(fatality | encounter), neither known) -- it only
##    shows, for whichever Nmin turns out to be defensible, how that
##    reading compares with what either model calls sustainable.
##

suppressPackageStartupMessages({ library(dplyr); library(tidyr); library(ggplot2) })

source("R/pbr_functions.R")
source("R/leslie_dekker_limpens.R")  # lambda_female (Leslie/Safi baseline, ~1.0474)
source("inputs/pbrSettings_BSH_DGY.R")

fig_dir <- "outputs/figures"
dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)

observed_mortality_annual <- 1000  # combined, both facilities, illustrative per Paulo -- see caveat above

# ---- 1. PBR(N, lambda) surface ---------------------------------------------
n_seq <- 10^seq(log10(2000), log10(100000), length.out = 200)
lambda_seq <- seq(1.00, 1.30, length.out = 200)
fr_panels <- c(0.3, 0.5)

surface_grid <- tidyr::expand_grid(N = n_seq, lambda = lambda_seq, Fr = fr_panels) %>%
  mutate(PBR = pbr_from_components(N, Fr, lambda), Fr_label = paste0("Fr = ", Fr))

pbr_breaks <- sort(unique(c(pbr_thresholds$threshold, 240, 500, 1000)))
lambda_refs <- tibble::tibble(
  lambda = c(lambda_female, 1.20, 1.24),
  label = c(sprintf("Leslie/Safi (%.3f)", lambda_female), "PBR 1.20", "PBR 1.24")
)

fig_nmin_lambda_surface <- file.path(fig_dir, "nmin_lambda_boundary_surface.png")
ggsave(fig_nmin_lambda_surface, width = 11, height = 6, dpi = 150, plot = {
  ggplot(surface_grid, aes(x = N, y = lambda)) +
    geom_raster(aes(fill = PBR), interpolate = TRUE) +
    geom_contour(aes(z = PBR, colour = after_stat(factor(level))), breaks = pbr_breaks, linewidth = 0.4) +
    geom_hline(data = lambda_refs, aes(yintercept = lambda), linetype = "dashed", colour = "white", linewidth = 0.3) +
    geom_text(data = lambda_refs, aes(x = 2000, y = lambda, label = label), colour = "white",
              hjust = 0, vjust = -0.4, size = 2.6) +
    geom_vline(xintercept = 3022, linetype = "dotted", colour = "grey90") +
    annotate("text", x = 3022, y = 1.29, label = "Nmin=3,022\n(median fallback)", colour = "grey90", size = 2.6, hjust = -0.05) +
    scale_x_log10(labels = scales::comma) +
    facet_wrap(~Fr_label) +
    scale_fill_viridis_c(name = "PBR\n(bats/yr)", option = "C") +
    scale_colour_manual(name = "PBR contour\n(bats/yr)",
                         values = setNames(c("cyan", "deepskyblue", "chartreuse", "yellow", "red")[seq_along(pbr_breaks)],
                                            as.character(pbr_breaks))) +
    labs(
      x = "Nmin (log scale)", y = "lambda",
      title = "PBR as a function of Nmin and lambda: which uncertainty dominates?",
      subtitle = "PBR = 0.5*(lambda-1)*Fr*Nmin is linear in Nmin -- contours show how fast PBR moves along each axis"
    ) +
    theme_minimal() +
    theme(plot.subtitle = element_text(size = 8), strip.text = element_text(face = "bold"))
})
cat("Wrote", fig_nmin_lambda_surface, "\n")

# ---- 2. Diagnostic: observed mortality as a fraction of assumed Nmin ------
n_diag <- c(3022, 6044, 10000, 15000, 20000, 30000, 40000, 60000, 80000, 100000)
diag_tbl <- tibble::tibble(N = n_diag) %>%
  mutate(pct_of_N = 100 * observed_mortality_annual / N)

# Reference "sustainable fraction" lines -- NOT inverse population estimates,
# just the fraction of N each model's own logic treats as sustainable.
leslie_surplus_fraction <- lambda_female - 1  # deterministic, no Fr haircut, no Rmax theoretical ceiling
pbr_fraction_bounds <- tidyr::expand_grid(Fr = c(0.3, 1.0), lambda = c(1.20, 1.24)) %>%
  mutate(fraction = 0.5 * Fr * (lambda - 1))

cat("=== Observed mortality (", observed_mortality_annual, "/yr, illustrative, combined facilities) as % of assumed Nmin ===\n", sep = "")
print(as.data.frame(diag_tbl))
cat(sprintf("\nLeslie/Safi deterministic surplus fraction: %.1f%%\n", 100 * leslie_surplus_fraction))
cat("PBR sustainable-fraction bounds (Fr x lambda_max):\n")
print(as.data.frame(pbr_fraction_bounds %>% mutate(fraction_pct = 100 * fraction)))

# Nmin each reference would require for the observed mortality to sit exactly
# at that sustainable fraction (closed form: N = mortality / fraction).
required_n_tbl <- bind_rows(
  tibble::tibble(reference = "Leslie/Safi deterministic surplus", fraction = leslie_surplus_fraction),
  pbr_fraction_bounds %>% mutate(reference = sprintf("PBR (Fr=%.1f, lambda=%.2f)", Fr, lambda)) %>% select(reference, fraction)
) %>%
  mutate(required_N = observed_mortality_annual / fraction) %>%
  arrange(required_N)

cat("\n=== Nmin required for the observed mortality to sit exactly at each reference's own sustainable fraction ===\n")
print(as.data.frame(required_n_tbl))

fig_fatality_fraction <- file.path(fig_dir, "fatality_fraction_of_nmin.png")
ggsave(fig_fatality_fraction, width = 9, height = 5.5, dpi = 150, plot = {
  ggplot(diag_tbl, aes(x = N, y = pct_of_N)) +
    geom_rect(aes(xmin = 2000, xmax = 120000, ymin = min(pbr_fraction_bounds$fraction) * 100,
                  ymax = max(pbr_fraction_bounds$fraction) * 100),
              fill = "steelblue", alpha = 0.12, inherit.aes = FALSE) +
    geom_hline(yintercept = leslie_surplus_fraction * 100, linetype = "dashed", colour = "firebrick") +
    annotate("text", x = 100000, y = leslie_surplus_fraction * 100, label = "Leslie/Safi surplus (4.7%)",
              colour = "firebrick", size = 3, hjust = 1, vjust = -0.5) +
    annotate("text", x = 100000, y = mean(range(pbr_fraction_bounds$fraction)) * 100,
              label = "PBR sustainable-fraction band\n(3-12%, across Fr/lambda_max)", colour = "steelblue4",
              size = 3, hjust = 1, vjust = -0.3) +
    geom_line(colour = "grey20", linewidth = 0.9) +
    geom_point(size = 2, colour = "grey20") +
    geom_point(data = diag_tbl %>% filter(N %in% c(3022, 6044)), colour = "black", fill = "orange",
               shape = 21, size = 3.5) +
    scale_x_log10(labels = scales::comma) +
    labs(
      x = "Assumed Nmin (log scale)", y = "Observed mortality as % of Nmin",
      title = "Diagnostic consistency check: observed mortality vs. assumed Nmin",
      subtitle = paste0(
        observed_mortality_annual, " bats/yr (illustrative, combined facilities) against a range of Nmin -- ",
        "NOT an inverse estimate of population size (P(encounter), P(fatality|encounter) unknown)"
      )
    ) +
    theme_minimal() +
    theme(plot.subtitle = element_text(size = 7.5))
})
cat("\nWrote", fig_fatality_fraction, "\n")
