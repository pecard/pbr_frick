##
## Leslie (age-structured, pre-breeding-census, female-only) matrix model,
## built to cross-check Frick et al.'s PBR/lambda_max approach against a
## standard demographic method using the same adult survival (s) and age
## at first breeding (alpha) inputs. See
## references/leslie_matrix_parametrisation.md for the literature behind
## the extra parameters a Leslie matrix needs that the PBR/Niel-Lebreton
## approximation does not (juvenile survival, an explicit fecundity).
##

# Fraction of age class `a` (individuals aged a to a+1 at census) that has
# already reached age at first breeding `alpha`. 0 below alpha, 1 once
# fully mature, pro-rated for the single age class alpha falls within
# (only relevant for non-integer alpha, e.g. 1.5).
frac_mature <- function(a, alpha) {
  ifelse(
    a >= ceiling(alpha), 1,
    ifelse(a == floor(alpha) & alpha != floor(alpha), ceiling(alpha) - alpha, 0)
  )
}

# Builds the (max_age + 1) x (max_age + 1) pre-breeding-census Leslie
# matrix for ages 0..max_age. female_fecundity = female pups per breeding
# female per year (before that pup's own survival to the next census,
# which is applied via s_juv per the pre-breeding-census convention).
build_leslie_matrix <- function(s_adult, s_juv, alpha, female_fecundity = 1, max_age = 15) {
  ages <- 0:max_age
  n <- length(ages)
  L <- matrix(0, n, n, dimnames = list(ages, ages))

  L[2, 1] <- s_juv
  if (n > 2) {
    for (i in 2:(n - 1)) L[i + 1, i] <- s_adult
  }

  L[1, ] <- female_fecundity * frac_mature(ages, alpha) * s_juv
  L
}

# Dominant eigenvalue (asymptotic population growth rate, comparable to
# lambda_max from the PBR/Niel-Lebreton approximation).
leslie_lambda <- function(L) {
  max(Mod(eigen(L, only.values = TRUE)$values))
}

# Convenience wrapper: lambda directly from the biological parameters.
leslie_lambda_from_params <- function(s_adult, s_juv, alpha, female_fecundity = 1, max_age = 15) {
  leslie_lambda(build_leslie_matrix(s_adult, s_juv, alpha, female_fecundity, max_age))
}
