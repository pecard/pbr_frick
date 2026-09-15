##
## Reconstructs and validates the exact Vespertilio murinus stage matrix
## from Dekker & Limpens (2024), "Population dynamics of Dutch bats
## relevant to wind energy" (Figure 3.1), then generalises it to test how
## the length of the pre-reproductive phase (an explicit stand-in for age
## at first breeding, comparable to the PBR analysis's alpha) changes
## lambda -- using the SAME juvenile survival (S_j) and recruitment
## formula validated against the published model, instead of re-guessing
## a matrix structure as R/leslie_vespertilio_murinus.R (kept for
## reference) did before this was available.
##
## Source parameters (Safi 2006 survival; report's own reproduction/litter
## figures): juvenile survival 0.62 (both sexes), adult female survival
## 0.76, adult male survival 0.42, 87% of adult females breeding, litter
## size 1.8, 1:1 sex ratio at birth.
##
## Recruitment term (confirmed to reproduce the published lambda exactly):
##   F = p_breed * litter * sex_ratio * S_j
## i.e. fecundity already discounts for the offspring's own survival to
## the next census, through the SAME S_j used for the juvenile->adult
## transition -- not a separate first-year survival rate.
##

source("R/leslie_matrix.R")

S_j  <- 0.62
S_af <- 0.76
S_am <- 0.42
p_breed <- 0.87
litter <- 1.8
sex_ratio <- 0.5

F_recruitment <- p_breed * litter * sex_ratio * S_j  # 0.4855

# ---- 1. Validate against the published model (Figure 3.1) ---------------
# Female-only submatrix: the male block is lower in the block-triangular
# structure (M_juv recruits from F_ad, M_ad from M_juv) and never feeds
# back into the female block, so it cannot affect the dominant eigenvalue
# -- confirmed numerically against the full 4-stage matrix below.
build_female_matrix <- function(s_juv, s_adult, fecundity) {
  matrix(c(0, fecundity, s_juv, s_adult), nrow = 2, byrow = TRUE)
}

build_four_stage_matrix <- function(s_juv, s_adult_f, s_adult_m, fecundity) {
  matrix(
    c(
      0,     fecundity, 0,     0,
      s_juv, s_adult_f, 0,     0,
      0,     fecundity, 0,     0,
      0,     0,         s_juv, s_adult_m
    ),
    nrow = 4, byrow = TRUE
  )
}

A_female <- build_female_matrix(S_j, S_af, F_recruitment)
A_full <- build_four_stage_matrix(S_j, S_af, S_am, F_recruitment)

lambda_female <- leslie_lambda(A_female)
lambda_full <- leslie_lambda(A_full)

cat(sprintf("lambda (female 2-stage) = %.4f\n", lambda_female))
cat(sprintf("lambda (full 4-stage)   = %.4f\n", lambda_full))
cat(sprintf("Published value: 1.047 (%s)\n",
            if (abs(lambda_female - 1.047) < 0.001) "MATCHES" else "DOES NOT MATCH"))

# Elasticities via the standard Caswell (2001) formula: e_ij = (a_ij / lambda) *
# (v_i * w_j) / (v . w), where w is the right eigenvector (stable stage
# distribution) and v the left eigenvector (reproductive value).
leslie_elasticities <- function(A) {
  eig <- eigen(A)
  lam <- Re(eig$values[1])
  w <- Re(eig$vectors[, 1]); w <- w / sum(w)
  v <- Re(eigen(t(A))$vectors[, 1])
  sens <- outer(v, w) / sum(v * w)
  (A / lam) * sens
}

elas <- leslie_elasticities(A_female)
cat("\nElasticity matrix (rows/cols: juvenile, adult):\n")
print(round(elas, 3))
cat(sprintf(
  "\nThis session's elasticities: adult survival = %.3f, juvenile survival = %.3f, reproduction = %.3f (sum = %.3f)\n",
  elas[2, 2], elas[2, 1], elas[1, 2], sum(elas)
))
cat("Published (Dekker & Limpens 2024): adult survival = 0.530, juvenile survival = 0.235, reproduction = 0.235.\n")
cat(
  "Close but not identical -- consistent with (and possibly related to) the lambda=1.047 vs 1.096\n",
  "inconsistency Paulo already flagged between two places in the same report; not resolved here.\n",
  sep = ""
)

# ============================================================
# 2. Generalise: how many years of juvenile survival before breeding?
# ============================================================
# The validated baseline already bakes in TWO S_j-survival intervals
# before an individual breeds: one inside F_recruitment (birth to first
# census) and one in the explicit juvenile -> adult transition -- i.e.
# it structurally represents age at first breeding ~2, not ~1 (as Paulo
# noted). This builds the (n_stages)-class generalisation: n_stages - 1
# juvenile/subadult stages, each surviving at S_j, then an adult stage
# breeding at S_a -- n_stages = 2 reproduces the validated baseline
# exactly; n_stages = 3, 4, ... represent progressively delayed maturation,
# using the SAME recruitment formula and S_j/S_a throughout (no
# re-guessing of structure, unlike the exploratory alpha = 1/2/3 comparison
# in R/leslie_vespertilio_murinus.R).
build_dekker_stage_matrix <- function(n_stages, s_juv, s_adult, fecundity) {
  stopifnot(n_stages >= 2)
  A <- matrix(0, n_stages, n_stages)
  A[1, n_stages] <- fecundity
  for (i in 1:(n_stages - 1)) A[i + 1, i] <- s_juv
  A[n_stages, n_stages] <- s_adult
  A
}

stopifnot(abs(leslie_lambda(build_dekker_stage_matrix(2, S_j, S_af, F_recruitment)) - lambda_female) < 1e-9)

# ---- Grid: adult survival (same range as the PBR s analysis) x maturation
#      delay, juvenile survival fixed at Safi's value throughout ---------
suppressPackageStartupMessages({ library(dplyr); library(tidyr); library(ggplot2) })

s_adult_range <- c(0.70, 0.95)  # matches inputs/pbrSettings_BSH_DGY.R
stage_range <- 2:5              # n_stages = 2 (validated baseline) .. 5

grid <- tidyr::expand_grid(
  s_adult = seq(s_adult_range[1], s_adult_range[2], length.out = 50),
  n_stages = stage_range
) %>%
  mutate(
    lambda_Leslie = mapply(function(s_a, n) leslie_lambda(build_dekker_stage_matrix(n, S_j, s_a, F_recruitment)),
                            s_adult, n_stages),
    n_stages_label = factor(paste0("n_stages = ", n_stages), levels = paste0("n_stages = ", stage_range))
  )

fig_dir <- "outputs/figures"
dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)

p <- ggplot(grid, aes(x = s_adult, y = lambda_Leslie, colour = n_stages_label)) +
  geom_line(linewidth = 0.9) +
  geom_hline(yintercept = c(1.20, 1.24), linetype = "dashed", colour = c("cyan", "chartreuse")) +
  geom_vline(xintercept = S_af, linetype = "dotted", colour = "grey40") +
  annotate("text", x = S_af, y = min(grid$lambda_Leslie), label = "Safi S_adult=0.76", angle = 90,
           vjust = -0.5, hjust = 0, size = 3, colour = "grey40") +
  scale_colour_viridis_d(name = "Pre-reproductive\nstages", option = "C") +
  labs(
    x = "Adult female survival", y = "lambda (Dekker & Limpens structure)",
    title = "V. murinus lambda vs adult survival, by maturation delay",
    subtitle = "S_juv = 0.62 (Safi) throughout; dashed: PBR benchmarks 1.20 (cyan)/1.24 (green)"
  ) +
  theme_minimal() +
  theme(plot.subtitle = element_text(size = 8))
ggsave(file.path(fig_dir, "leslie_dekker_maturation_delay.png"), p, width = 7.5, height = 4.5, dpi = 150)

cat("\nlambda at Safi's S_adult=0.76, by number of pre-reproductive stages:\n")
print(grid %>% filter(abs(s_adult - S_af) < 0.01) %>% select(n_stages, lambda_Leslie) %>% distinct())

cat("\nWrote", file.path(fig_dir, "leslie_dekker_maturation_delay.png"), "\n")
