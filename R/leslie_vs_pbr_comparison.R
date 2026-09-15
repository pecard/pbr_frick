# ============================================================
# Compares Frick et al.'s PBR/lambda_max (Niel & Lebreton demographic
# -invariant approximation) against a standard Leslie matrix, using the
# same adult survival (s) and age at first breeding (alpha) ranges
# already established for the Parti-coloured bat analysis (s in
# [0.70, 0.95], alpha in [1.5, 3.5] -- see
# inputs/pbrSettings_BSH_DGY.R and references/alpha_first_breeding_evidence.md).
# The Leslie matrix additionally needs juvenile survival and fecundity;
# see references/leslie_matrix_parametrisation.md for their literature
# basis (Myotis lucifugus proxy studies) and the caveat that this
# session's web access could not verify them beyond search-result
# snippets.
# ============================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(ggplot2)
})

source("R/pbr_functions.R")
source("R/leslie_matrix.R")

s_range <- c(0.70, 0.95)
alpha_range <- c(1.5, 3.5)

# Juvenile survival: no species-specific figure for Vespertilio murinus;
# range and central value from Myotis lucifugus (Frick et al. 2010;
# Schorr et al. 2021) -- see references/leslie_matrix_parametrisation.md.
s_juv_range <- c(0.23, 0.71)
s_juv_central <- 0.45

female_fecundity <- 1  # V. murinus: litter size 2, 1:1 sex ratio (AnAge)

# ---- 1. Point comparison at the PBR analysis's own reference values -----
s0 <- 0.88; alpha0 <- 1.5  # matches the earlier dummy-species point in pbr_frick_main.R
lam_nl_0 <- lambda_max_niel(s0, alpha0)
lam_leslie_0 <- leslie_lambda_from_params(s0, s_juv_central, alpha0, female_fecundity)

cat(sprintf(
  "At s=%.2f, alpha=%.2f, juvenile survival=%.2f: lambda_NielLebreton=%.4f, lambda_Leslie=%.4f (diff %.4f, ratio %.3f)\n",
  s0, alpha0, s_juv_central, lam_nl_0, lam_leslie_0, lam_nl_0 - lam_leslie_0, lam_leslie_0 / lam_nl_0
))

# ---- 2. Comparison across the full (s, alpha) grid, central juvenile ----
#         survival ---------------------------------------------------------
grid <- tidyr::expand_grid(
  s = seq(s_range[1], s_range[2], length.out = 40),
  alpha = seq(alpha_range[1], alpha_range[2], length.out = 40)
) %>%
  mutate(
    lambda_NielLebreton = lambda_max_niel(s, alpha),
    lambda_Leslie = mapply(leslie_lambda_from_params, s, s_juv_central, alpha,
                            MoreArgs = list(female_fecundity = female_fecundity)),
    diff = lambda_NielLebreton - lambda_Leslie,
    ratio = lambda_Leslie / lambda_NielLebreton
  )

cat("\nDiscrepancy (Niel-Lebreton minus Leslie) across the (s, alpha) grid:\n")
print(summary(grid$diff))
cat("\nRatio (Leslie / Niel-Lebreton) across the (s, alpha) grid:\n")
print(summary(grid$ratio))

fig_dir <- "outputs/figures"
dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)

p_compare <- ggplot(grid, aes(x = s, y = alpha, fill = diff)) +
  geom_raster(interpolate = TRUE) +
  scale_fill_distiller(name = "lambda_NL -\nlambda_Leslie", palette = "RdBu", direction = -1) +
  labs(
    x = "Adult survival (s)", y = "Age at first breeding (alpha)",
    title = "PBR (Niel-Lebreton) vs Leslie matrix: lambda discrepancy",
    subtitle = sprintf("Leslie: juvenile survival = %.2f (central estimate), female fecundity = %s", s_juv_central, female_fecundity)
  ) +
  theme_minimal()
ggsave(file.path(fig_dir, "leslie_vs_pbr_diff.png"), p_compare, width = 7, height = 4.5, dpi = 150)

# ---- 3. How much does the comparison depend on juvenile survival? -------
#         (the parameter the PBR method does not need, and the one this
#         session could least pin down in the literature) -----------------
juv_sens <- tidyr::expand_grid(
  s_juv = seq(s_juv_range[1], s_juv_range[2], length.out = 30),
  s = c(0.70, 0.80, 0.88, 0.95)
) %>%
  mutate(
    alpha = alpha0,
    lambda_Leslie = mapply(leslie_lambda_from_params, s, s_juv, alpha,
                            MoreArgs = list(female_fecundity = female_fecundity)),
    lambda_NielLebreton = lambda_max_niel(s, alpha)
  )

p_juv_sens <- ggplot(juv_sens, aes(x = s_juv, y = lambda_Leslie, colour = factor(s))) +
  geom_line(linewidth = 0.8) +
  geom_hline(
    data = distinct(juv_sens, s, lambda_NielLebreton),
    aes(yintercept = lambda_NielLebreton, colour = factor(s)),
    linetype = "dashed", linewidth = 0.5
  ) +
  scale_colour_viridis_d(name = "Adult survival (s)", option = "C") +
  labs(
    x = "Juvenile (first-year) survival", y = "lambda",
    title = "Leslie lambda vs juvenile survival, at alpha = 1.5",
    subtitle = "Dashed: Niel-Lebreton lambda_max at same s (doesn't use juvenile survival)"
  ) +
  theme_minimal() +
  theme(plot.subtitle = element_text(size = 9))
ggsave(file.path(fig_dir, "leslie_juvenile_sensitivity.png"), p_juv_sens, width = 7, height = 4.5, dpi = 150)

cat("\nWrote", file.path(fig_dir, "leslie_vs_pbr_diff.png"), "and",
    file.path(fig_dir, "leslie_juvenile_sensitivity.png"), "\n")
