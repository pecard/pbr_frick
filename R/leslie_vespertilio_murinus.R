# ============================================================
# Leslie matrix for Vespertilio murinus using species-specific parameters
# (Safi 2006 survival; Zhigalin & Moskvitina 2017 fecundity) instead of
# the Myotis lucifugus proxy used in leslie_vs_pbr_comparison.R. See
# references/leslie_matrix_parametrisation.md for sources, the important
# caveat on Safi's survival estimates (simple return rates, no modern
# mark-recapture model, no confidence intervals), and why three
# age-at-first-breeding structures are compared explicitly here.
# ============================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(ggplot2)
})

source("R/pbr_functions.R")
source("R/leslie_matrix.R")

# ---- Point parametrisation (Safi 2006; Zhigalin & Moskvitina 2017) ------
s_juv_point <- 0.62
s_adult_point <- 0.76
p_breed <- 0.87
litter_suburban <- 1.8
litter_urban_upper <- 2.9

fecundity_suburban <- p_breed * litter_suburban * 0.5       # 0.783
fecundity_urban <- p_breed * litter_urban_upper * 0.5       # 1.262

# Sensitivity ranges (Paulo's proposal -- scenarios, not published values)
s_juv_range <- c(0.45, 0.70)
s_adult_range <- c(0.65, 0.85)

alpha_structures <- c(early = 1, intermediate = 2, delayed = 3)
lambda_max_benchmarks <- c(1.24, 1.20)

# ---- 1. Point lambda for each alpha structure, both fecundity levels ----
point_grid <- tidyr::expand_grid(
  alpha_label = names(alpha_structures),
  fecundity_label = c("suburban", "urban")
) %>%
  mutate(
    alpha = alpha_structures[alpha_label],
    fecundity = ifelse(fecundity_label == "suburban", fecundity_suburban, fecundity_urban),
    lambda_Leslie = mapply(leslie_lambda_from_params, s_adult_point, s_juv_point, alpha, fecundity),
    lambda_vs_1.24 = lambda_Leslie - 1.24,
    lambda_vs_1.20 = lambda_Leslie - 1.20
  )

cat("Leslie lambda by age-at-first-breeding structure (S_juv=0.62, S_adult=0.76):\n")
print(point_grid %>% select(alpha_label, alpha, fecundity_label, lambda_Leslie, lambda_vs_1.24, lambda_vs_1.20))

# Also compare against the continuous alpha range already used elsewhere
# in this analysis (1.5-3.5), same point survival/fecundity
continuous_alpha <- tibble::tibble(alpha = seq(1.5, 3.5, by = 0.5)) %>%
  mutate(
    lambda_Leslie_suburban = mapply(leslie_lambda_from_params, s_adult_point, s_juv_point, alpha,
                                     MoreArgs = list(female_fecundity = fecundity_suburban)),
    lambda_Leslie_urban = mapply(leslie_lambda_from_params, s_adult_point, s_juv_point, alpha,
                                  MoreArgs = list(female_fecundity = fecundity_urban))
  )
cat("\nSame survival/fecundity, across the continuous alpha range used elsewhere (1.5-3.5):\n")
print(continuous_alpha)

# ---- 2. Sensitivity across S_juv, S_adult ranges, at each alpha structure -
sens_grid <- tidyr::expand_grid(
  alpha_label = names(alpha_structures),
  s_juv = seq(s_juv_range[1], s_juv_range[2], length.out = 25),
  s_adult = seq(s_adult_range[1], s_adult_range[2], length.out = 25)
) %>%
  mutate(
    alpha = alpha_structures[alpha_label],
    lambda_Leslie = mapply(leslie_lambda_from_params, s_adult, s_juv, alpha,
                            MoreArgs = list(female_fecundity = fecundity_suburban)),
    alpha_label = factor(alpha_label, levels = c("early", "intermediate", "delayed"))
  )

fig_dir <- "outputs/figures"
dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)

p_sens <- ggplot(sens_grid, aes(x = s_adult, y = s_juv, fill = lambda_Leslie)) +
  geom_raster(interpolate = TRUE) +
  geom_contour(aes(z = lambda_Leslie, colour = after_stat(factor(level))),
               breaks = lambda_max_benchmarks, linewidth = 0.5) +
  facet_wrap(~alpha_label, nrow = 1) +
  scale_fill_viridis_c(name = "lambda", option = "C") +
  scale_colour_manual(name = "PBR benchmark", values = c("1.2" = "cyan", "1.24" = "chartreuse")) +
  labs(
    x = "Adult survival", y = "Juvenile survival",
    title = "Vespertilio murinus Leslie matrix: lambda by structure",
    subtitle = paste0("Fecundity = ", round(fecundity_suburban, 3), " (suburban); point estimate S_juv=0.62, S_adult=0.76 marked by cross"),
  ) +
  geom_point(
    data = tibble::tibble(
      s_adult = s_adult_point, s_juv = s_juv_point,
      alpha_label = factor(c("early", "intermediate", "delayed"), levels = c("early", "intermediate", "delayed"))
    ),
    aes(x = s_adult, y = s_juv),
    inherit.aes = FALSE, shape = 4, size = 2.5, colour = "white", stroke = 1.2
  ) +
  theme_minimal() +
  theme(strip.text = element_text(face = "bold"), plot.subtitle = element_text(size = 8))

ggsave(file.path(fig_dir, "leslie_vespertilio_sensitivity.png"), p_sens, width = 11, height = 4, dpi = 150)

# ---- 3. Elasticity decomposition: does adult survival still dominate? ---
elasticity_leslie <- function(s_adult, s_juv, alpha, fecundity, eps = 1e-4) {
  lam0 <- leslie_lambda_from_params(s_adult, s_juv, alpha, fecundity)
  d_sa <- (leslie_lambda_from_params(s_adult + eps, s_juv, alpha, fecundity) -
             leslie_lambda_from_params(s_adult - eps, s_juv, alpha, fecundity)) / (2 * eps)
  d_sj <- (leslie_lambda_from_params(s_adult, s_juv + eps, alpha, fecundity) -
             leslie_lambda_from_params(s_adult, s_juv - eps, alpha, fecundity)) / (2 * eps)
  d_f <- (leslie_lambda_from_params(s_adult, s_juv, alpha, fecundity + eps) -
            leslie_lambda_from_params(s_adult, s_juv, alpha, fecundity - eps)) / (2 * eps)
  tibble::tibble(
    lambda = lam0,
    E_s_adult = d_sa * s_adult / lam0,
    E_s_juv = d_sj * s_juv / lam0,
    E_fecundity = d_f * fecundity / lam0
  )
}

elasticity_by_structure <- tidyr::expand_grid(alpha_label = names(alpha_structures)) %>%
  mutate(alpha = alpha_structures[alpha_label]) %>%
  rowwise() %>%
  mutate(out = list(elasticity_leslie(s_adult_point, s_juv_point, alpha, fecundity_suburban))) %>%
  ungroup() %>%
  tidyr::unnest(out)

cat("\nElasticity of Leslie lambda to each vital rate, by structure (point estimates):\n")
print(elasticity_by_structure %>% select(alpha_label, alpha, lambda, E_s_adult, E_s_juv, E_fecundity))
cat(
  "\n(Compare to the Niel-Lebreton/PBR elasticity analysis: E_s ~ -0.13 to -1.35, ",
  "E_alpha ~ -0.06 to -0.22 -- adult survival dominated there too, but with the opposite ",
  "sign, since higher NL 's' trades off against fecundity in a way this explicit-fecundity ",
  "Leslie matrix does not.)\n", sep = ""
)

cat("\nWrote", file.path(fig_dir, "leslie_vespertilio_sensitivity.png"), "\n")

# ---- 4. Bias-aware correction: Safi's rates as a LOWER bound -------------
# A naive return rate (marked in year t, recaptured in t+1) is mathematically
# the product of true survival and recapture probability (and site fidelity),
# not survival alone -- so it is a lower bound on true survival, never an
# overestimate. No direct naive-vs-CJS comparison or a more recent
# V. murinus-specific CJS estimate could be located this session (see
# references/leslie_matrix_parametrisation.md), so the upper end of the
# corrected range is anchored instead on CJS-based estimates for
# ecologically similar, related vespertilionines already used elsewhere in
# this analysis: Pipistrellus pipistrellus adult survival ~0.80 (+/-0.05)
# from Sendor & Simon (2003), and the 0.9 "reasonable higher-end" bat
# adult survival benchmark Frick et al. themselves adopt. This is a
# plausibility anchor, not a formal bias correction (no correction factor
# specific to Safi's study design could be derived from what was
# accessible this session).
s_adult_corrected_range <- c(0.76, 0.90)  # Safi's raw value as the floor
s_juv_corrected_range <- c(0.62, 0.75)    # Safi's raw value as the floor

corrected_grid <- tidyr::expand_grid(
  alpha_label = names(alpha_structures),
  scenario = c("Safi raw (uncorrected)", "Bias-aware upper bound"),
  fecundity_label = c("suburban", "urban")
) %>%
  mutate(
    alpha = alpha_structures[alpha_label],
    fecundity = ifelse(fecundity_label == "suburban", fecundity_suburban, fecundity_urban),
    s_adult = ifelse(scenario == "Safi raw (uncorrected)", s_adult_corrected_range[1], s_adult_corrected_range[2]),
    s_juv = ifelse(scenario == "Safi raw (uncorrected)", s_juv_corrected_range[1], s_juv_corrected_range[2]),
    lambda_Leslie = mapply(leslie_lambda_from_params, s_adult, s_juv, alpha, fecundity),
    lambda_vs_1.24 = lambda_Leslie - 1.24,
    lambda_vs_1.20 = lambda_Leslie - 1.20
  )

cat("\nEffect of correcting for the known naive-return-rate downward bias\n",
    "(raw Safi values vs a bias-aware upper bound anchored on Pipistrellus\n",
    "pipistrellus CJS survival and Frick et al.'s own bat benchmark):\n", sep = "")
print(corrected_grid %>%
        select(alpha_label, scenario, fecundity_label, s_adult, s_juv, lambda_Leslie, lambda_vs_1.24, lambda_vs_1.20) %>%
        arrange(alpha_label, fecundity_label, scenario))

# Isolate how much of the correction's effect comes from s_adult alone,
# given it has by far the largest elasticity (step 3 above)
adult_only_effect <- tibble::tibble(
  s_adult = s_adult_corrected_range,
  lambda_at_alpha1_suburban = mapply(leslie_lambda_from_params, s_adult, s_juv_point, 1, fecundity_suburban)
)
cat("\nHolding juvenile survival and fecundity at Safi's raw values (alpha=1, suburban),\n",
    "correcting only adult survival (0.76 -> 0.90):\n", sep = "")
print(adult_only_effect)
