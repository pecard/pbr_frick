# Choose a baseline point (e.g., your hypothetical bat)
s0 <- 0.88
a0 <- 1.5

# Pull the local elasticities at the closest grid cell
pt <- sens_grid %>%
  mutate(ds = abs(s - s0), da = abs(alpha - a0)) %>%
  arrange(ds, da) %>%
  slice(1)

E_s0 <- pt$E_s
E_a0 <- pt$E_alpha
lam0 <- pt$lambda_max

# Specify plausible *relative* uncertainties (example)
# e.g., survival ±3% (0.03), age at first breeding ±10% (0.10)
rel_s  <- 0.03
rel_a  <- 0.10

# Approx relative change in lambda_max (worst-case magnitude)
rel_lam_wc <- abs(E_s0)*rel_s + abs(E_a0)*rel_a

lam_lo <- lam0 * (1 - rel_lam_wc)
lam_hi <- lam0 * (1 + rel_lam_wc)

# Convert to PBR band (holding Fr and Nmin fixed)
Fr <- 0.3
Nmin <- 628.4641

pbr <- function(lam, Fr, Nmin) 0.5 * (lam - 1) * Fr * Nmin

pbr0  <- pbr(lam0,  Fr, Nmin)
pbr_lo <- pbr(lam_lo, Fr, Nmin)
pbr_hi <- pbr(lam_hi, Fr, Nmin)

tibble::tibble(
  s0 = s0, alpha0 = a0,
  lambda0 = lam0,
  lambda_lo = lam_lo, lambda_hi = lam_hi,
  PBR0 = pbr0, PBR_lo = pbr_lo, PBR_hi = pbr_hi
)

# Monte Carlo Propagation

set.seed(1)

lambda_max_niel <- function(s, alpha) {
  disc <- (s - s * alpha - alpha - 1)^2 - 4 * s * alpha^2
  out <- ((s * alpha - s + alpha + 1) + sqrt(pmax(disc, 0))) / (2 * alpha)
  out[disc < 0] <- NA_real_
  out
}

# Baseline + uncertainty model (edit as needed)
s0 <- 0.88
a0 <- 1.5

# Example: truncated normal for s in (0,1); lognormal for alpha > 0
n <- 20000
s_draw <- pmin(pmax(rnorm(n, mean = s0, sd = 0.03), 0.5), 0.99)
a_draw <- rlnorm(n, meanlog = log(a0), sdlog = 0.10)

lam_draw <- lambda_max_niel(s_draw, a_draw)

Fr <- 0.3
Nmin <- 628.4641

pbr_draw <- 0.5 * (lam_draw - 1) * Fr * Nmin

# Summarise uncertainty band
quant <- c(0.05, 0.5, 0.95)
out <- tibble::tibble(
  lambda_q05 = quantile(lam_draw, quant[1], na.rm = TRUE),
  lambda_q50 = quantile(lam_draw, quant[2], na.rm = TRUE),
  lambda_q95 = quantile(lam_draw, quant[3], na.rm = TRUE),
  PBR_q05    = quantile(pbr_draw, quant[1], na.rm = TRUE),
  PBR_q50    = quantile(pbr_draw, quant[2], na.rm = TRUE),
  PBR_q95    = quantile(pbr_draw, quant[3], na.rm = TRUE)
)

out

# How to write this up (one sentence):
# Uncertainty in the demographic inputs s and alpha was propagated to λmax and to 
# the classical PBR using simulation; the resulting PBR distribution
# (e.g., 5th–95th percentile) provides an uncertainty band reflecting plausible 
# demographic variability.
