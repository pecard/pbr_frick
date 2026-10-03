##
## PROTOTYPE -- carrying-capacity-referenced mortality threshold (Paulo,
## 2026-10, following Thierry Chambert's comment and the Epstein et al.
## 2015 FCS paper). Not wired into the launcher/report; a standalone
## script to look at before deciding whether/how to formalise it.
##
## Classic PBR anchors removal to Nmin (a pessimistic estimate of CURRENT
## abundance): "how much can we take without falling below a
## conservative floor." Thierry's alternative anchors removal to K
## (carrying capacity): "what constant removal keeps the population near
## K, not just away from extinction" -- closer to the EU Habitats
## Directive's ecological-viability/FCS standard than to MVP-style
## demographic viability (Epstein et al. 2015).
##
## The bridge is classic logistic/Schaefer surplus-production harvesting
## theory (Clark, Mathematical Bioeconomics): under dN/dt = rN(1-N/K) - H,
## a constant harvest H has a stable equilibrium at N* = alpha*K when
##   H(alpha) = r * K * alpha * (1 - alpha)
## alpha near 1 (stay close to K, the FCS-style standard) implies a much
## SMALLER sustainable H than alpha = 0.5 (classic MSY) -- confirms the
## intuition that an ecological-viability standard is more conservative
## than PBR's demographic-viability standard, once K is pinned down.
##
## The catch, made concrete here rather than just asserted: alpha(1-alpha)
## is symmetric, so every H below the MSY ceiling (rK/4) has TWO
## equilibria -- a stable one near K and an UNSTABLE one near extinction.
## A population pushed below the unstable one by ordinary stochastic
## variation collapses even though H was "sustainable" at the
## deterministic equilibrium. This is tested below by stochastic
## simulation (reusing the same Leslie-based demographic engine as
## R/pva_robustness.R), not just reported as a theoretical caveat.
##
## K itself is the same project-footprint x density-proxy proxy already
## used for the density-dependence PVA-lite sensitivity (no independent
## abundance/habitat-capacity estimate exists for this species/region --
## this is a bounding scenario, not a calibrated K).
##

if (!requireNamespace("ggrepel", quietly = TRUE)) {
  message("Installing missing package needed for the litter-sweep figure: ggrepel...")
  tryCatch(install.packages("ggrepel"), error = function(e) stop(
    "Could not install 'ggrepel' automatically (", conditionMessage(e), "). Install manually and re-run."
  ))
}

suppressPackageStartupMessages({ library(dplyr); library(tidyr); library(ggplot2); library(ggrepel) })

run_carrying_capacity_reference <- function(fig_dir, nmin_assumed,
                                             alpha_grid = c(0.5, 0.6, 0.7, 0.8, 0.9),
                                             juv_mortality_ratio = 1.3) {

  dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)
  fig_dir <- normalizePath(fig_dir, winslash = "/")

  project_thresholds <- tibble::tibble(project = facility_labels, threshold = pbr_thresholds$threshold)

  ## ---- 0. Realistic r from the validated Leslie matrix (not the PBR
  ## Niel-Lebreton Rmax benchmark) ---------------------------------------------
  ## lambda_max_benchmarks (1.20/1.24) is the demographic-invariant
  ## THEORETICAL MAXIMUM growth rate PBR uses to set Rmax -- an idealised
  ## ceiling, deliberately optimistic for threshold-setting purposes, not
  ## a claim about what this species actually achieves. The validated
  ## Leslie matrix (same real vital rates used throughout this note)
  ## gives a realistic lambda = 1.047 (r = ln(lambda) = 0.046) -- roughly
  ## 5x lower. Using the PBR Rmax inside a logistic growth formula and
  ## then testing the resulting H against a stochastic simulator built
  ## from the REAL vital rates is internally inconsistent: the formula
  ## assumes a growth capacity the simulated demography cannot deliver.
  ## Included here as its own r scenario specifically to make this
  ## mismatch visible instead of silently baking it in.
  build_dekker_stage_matrix_r <- function(n_stages, s_juv, s_adult, fecundity) {
    A <- matrix(0, n_stages, n_stages); A[1, n_stages] <- fecundity
    for (i in 1:(n_stages - 1)) A[i + 1, i] <- s_juv
    A[n_stages, n_stages] <- s_adult; A
  }
  leslie_F_cc <- leslie_p_breed * leslie_litter * leslie_sex_ratio * leslie_s_juv
  A_leslie_cc <- build_dekker_stage_matrix_r(2, leslie_s_juv, leslie_s_adult_f, leslie_F_cc)
  lambda_leslie_validated_cc <- Re(eigen(A_leslie_cc)$values[1])
  r_validated <- log(lambda_leslie_validated_cc)

  ## ---- 1. Deterministic H(alpha, K) across the same K and r uncertainty
  ## already used elsewhere in this note -------------------------------------
  K_scenarios <- tibble::tibble(
    K_label = c("K low (density_proxy_low x footprint)", "K high (density_proxy_high x footprint)"),
    K = c(density_proxy_low, density_proxy_high) * project_footprint_km2
  )
  r_scenarios <- tibble::tibble(
    r_label = c(paste0("PBR Rmax (theoretical max) at lambda_max = ", lambda_max_benchmarks),
                "Realistic r (validated Leslie matrix, lambda=1.047)"),
    Rmax = c(lambda_max_benchmarks - 1, r_validated)
  )

  cc_table <- tidyr::expand_grid(project_thresholds, K_scenarios, r_scenarios, alpha = alpha_grid) %>%
    mutate(
      N_target = alpha * K,
      H = Rmax * K * alpha * (1 - alpha),
      pct_of_current_threshold = 100 * H / threshold,
      is_MSY = alpha == 0.5
    )

  fig_H_curve <- file.path(fig_dir, "carrying_capacity_H_curve.png")
  ggsave(fig_H_curve, width = 9.5, height = 5, dpi = 150, bg = "white", plot = {
    curve_dt <- tidyr::expand_grid(project_thresholds, K_scenarios, r_scenarios,
                                    alpha = seq(0.5, 0.95, length.out = 100)) %>%
      mutate(H = Rmax * K * alpha * (1 - alpha))
    ggplot(curve_dt, aes(x = alpha, y = H, colour = r_label, linetype = K_label)) +
      geom_line(linewidth = 0.8) +
      geom_hline(data = project_thresholds, aes(yintercept = threshold), linetype = "dotted", colour = "grey30") +
      geom_vline(xintercept = 0.5, linetype = "dotted", colour = "grey60") +
      facet_wrap(~project, scales = "free_y") +
      labs(
        x = "Target alpha = N*/K (distance from carrying capacity)",
        y = "Sustainable constant removal H(alpha)",
        title = "Carrying-capacity-referenced removal vs. classic PBR threshold",
        subtitle = paste0(
          "Dotted horizontal: current imposed PBR threshold. Dotted vertical: alpha=0.5 (classic MSY).\n",
          "H falls as alpha moves toward 1 (closer to K) -- an ecological-viability/FCS-style standard is more conservative than PBR."
        )
      ) +
      theme_minimal() +
      theme(plot.subtitle = element_text(size = 9.5), strip.text = element_text(face = "bold"), legend.position = "bottom")
  })

  ## ---- 2. Stochastic check: does H(alpha) actually hold, and does it
  ## matter where the population starts? (the bistability risk made
  ## concrete, not just asserted) ---------------------------------------------
  build_dekker_stage_matrix <- function(n_stages, s_juv, s_adult, fecundity) {
    A <- matrix(0, n_stages, n_stages)
    A[1, n_stages] <- fecundity
    for (i in 1:(n_stages - 1)) A[i + 1, i] <- s_juv
    A[n_stages, n_stages] <- s_adult
    A
  }
  leslie_F <- leslie_p_breed * leslie_litter * leslie_sex_ratio * leslie_s_juv
  A_leslie <- build_dekker_stage_matrix(2, leslie_s_juv, leslie_s_adult_f, leslie_F)
  stable_dist <- Re(eigen(A_leslie)$vectors[, 1]); stable_dist <- stable_dist / sum(stable_dist)

  ## A Beta(mean, var) distribution only exists for var < mean*(1-mean); the
  ## fixed cv=10% assumption breaches that ceiling once mean gets close to 0
  ## or 1 (e.g. the s_adult~0.997 breakeven value Section 5 finds is needed
  ## to reach lambda=1.24 alone -- common would go negative and rbeta()
  ## would silently return NaN downstream). Clamping var just below the
  ## ceiling is a no-op for every mean this note actually uses elsewhere
  ## (0.4-0.95 at cv=0.10, nowhere near the ceiling) and only engages for
  ## these extreme, already-flagged-as-implausible elasticity points.
  beta_params <- function(mean, cv) {
    var <- min((mean * cv)^2, mean * (1 - mean) * 0.98)
    common <- mean * (1 - mean) / var - 1
    list(shape1 = mean * common, shape2 = (1 - mean) * common)
  }

  ## Same discrete-time, density-dependent, stochastic demographic engine
  ## as R/pva_robustness.R's simulate_trajectory() (duplicated here, not
  ## imported, since this is a standalone prototype) -- constant annual
  ## removal H, logistic compensation in recruitment as N approaches K.
  ## juv_mortality_ratio: per-capita collision-mortality risk for the
  ## juvenile-to-adult cohort, relative to adults (1 = equal risk; >1 =
  ## juveniles somewhat more likely to be hit -- first-year bats are
  ## less experienced flyers, consistent with the hoary bat literature's
  ## sex/age-class differences in collision risk). H is split between
  ## the two cohorts that form next year's population (maturing
  ## juveniles and surviving adults) in proportion to their
  ## risk-weighted abundance, not dumped entirely on adults as before.
  ## Turbine mortality is assumed to fall only on already-volant
  ## individuals (this year's standing juv/adult classes), never on
  ## newborn pups not yet flying, so the differential applies here and
  ## not to n_juv_recruits.
  simulate_constant_harvest <- function(n0_juv, n0_adult, n_years, H, K, cv = pva_vital_rate_cv,
                                         juv_mortality_ratio = 1) {
    bp_s_adult <- beta_params(leslie_s_adult_f, cv)
    bp_s_juv <- beta_params(leslie_s_juv, cv)
    bp_p_breed <- beta_params(leslie_p_breed, cv)
    n_juv <- as.integer(round(n0_juv)); n_adult <- as.integer(round(n0_adult))
    traj <- numeric(n_years + 1); traj[1] <- 2 * (n_juv + n_adult)
    for (yr in seq_len(n_years)) {
      s_adult_t <- rbeta(1, bp_s_adult$shape1, bp_s_adult$shape2)
      s_juv_t <- rbeta(1, bp_s_juv$shape1, bp_s_juv$shape2)
      p_breed_t <- rbeta(1, bp_p_breed$shape1, bp_p_breed$shape2)
      density_factor <- max(0, 1 - (2 * (n_juv + n_adult)) / K)
      n_breeders <- rbinom(1, n_adult, p_breed_t)
      n_pups_total <- rpois(1, n_breeders * leslie_litter * density_factor)
      n_female_pups <- rbinom(1, n_pups_total, leslie_sex_ratio)
      n_juv_recruits <- rbinom(1, n_female_pups, s_juv_t)
      n_juv_survive <- rbinom(1, n_juv, s_juv_t)
      n_adult_survive <- rbinom(1, n_adult, s_adult_t)

      juv_weight <- n_juv_survive * juv_mortality_ratio
      adult_weight <- n_adult_survive
      total_weight <- juv_weight + adult_weight
      if (total_weight > 0 && H > 0) {
        H_juv <- H * juv_weight / total_weight
        H_adult <- H * adult_weight / total_weight
      } else {
        H_juv <- 0; H_adult <- 0
      }
      ## H is a continuous removal rate (not necessarily an integer
      ## number of animals); round after subtracting it, otherwise
      ## n_adult drifts non-integer and the next iteration's rbinom()
      ## silently returns NaN (rbinom requires an integer `size`).
      juv_after_H <- as.integer(round(max(0, n_juv_survive - H_juv)))
      adult_after_H <- as.integer(round(max(0, n_adult_survive - H_adult)))
      n_adult_next <- juv_after_H + adult_after_H
      n_juv <- n_juv_recruits; n_adult <- n_adult_next
      traj[yr + 1] <- 2 * (n_juv + n_adult)
    }
    traj
  }

  run_from_start <- function(N0, H, K, n_years = pva_n_years, n_reps = pva_n_reps, seed = pva_seed,
                              juv_mortality_ratio = 1) {
    set.seed(seed)
    nf <- N0 / 2; nj <- round(nf * stable_dist[1]); na <- nf - nj
    reps <- replicate(n_reps, simulate_constant_harvest(nj, na, n_years, H, K, juv_mortality_ratio = juv_mortality_ratio))
    list(
      trajectories = reps,
      median_traj = apply(reps, 1, median),
      p10_traj = apply(reps, 1, quantile, probs = 0.10),
      p90_traj = apply(reps, 1, quantile, probs = 0.90),
      p_collapse = mean(reps[n_years + 1, ] < 0.1 * K) * 100
    )
  }

  ## ---- 2b. Empirical search for the actually-sustainable H at a given
  ## target, by simulation -- not the closed-form Schaefer shortcut,
  ## which Section 2 already showed is too optimistic once realistic age
  ## structure and stochasticity are in play. Binary search over H:
  ## collapse risk is monotonically increasing in H, so bisection finds
  ## the largest H keeping P(falls below quasi-extinction in n_years)
  ## at or under an explicit, named acceptable_risk -- not left implicit.
  ## Uses fewer reps during the search (directional signal only) and
  ## reruns the converged H at full pva_n_reps for a reportable number.
  find_sustainable_H <- function(N0, K, acceptable_risk, juv_mortality_ratio = 1.3,
                                  H_upper = NULL, tol = 0.5, max_iter = 25, search_reps = 300) {
    if (is.null(H_upper)) H_upper <- K * 0.1
    H_lo <- 0; H_hi <- H_upper
    for (iter in seq_len(max_iter)) {
      H_mid <- (H_lo + H_hi) / 2
      risk <- run_from_start(N0, H_mid, K, n_reps = search_reps, juv_mortality_ratio = juv_mortality_ratio)$p_collapse / 100
      if (risk > acceptable_risk) H_hi <- H_mid else H_lo <- H_mid
      if ((H_hi - H_lo) < tol) break
    }
    H_lo
  }

  acceptable_risk_grid <- c(0.05, 0.10, 0.20)
  realistic_r <- r_validated
  sustainable_H_table <- project_thresholds %>%
    rowwise() %>%
    reframe(
      project = project, threshold = threshold,
      K = density_proxy_high * project_footprint_km2,
      N_target = 0.8 * K,
      acceptable_risk = acceptable_risk_grid
    ) %>%
    rowwise() %>%
    mutate(
      H_sustainable = find_sustainable_H(N_target, K, acceptable_risk, juv_mortality_ratio = juv_mortality_ratio),
      ## Confirm at full rep count -- the search uses fewer reps for speed.
      p_collapse_confirmed = run_from_start(N_target, H_sustainable, K, juv_mortality_ratio = juv_mortality_ratio)$p_collapse,
      pct_of_current_threshold = 100 * H_sustainable / threshold,
      pct_of_closedform_H = 100 * H_sustainable / (realistic_r * K * 0.8 * 0.2)
    ) %>%
    ungroup()

  ## Risk-vs-H sweep (one project -- K/threshold formulae are identical
  ## across projects here, only the imposed threshold used for the %
  ## comparison differs) for a visual of where the three acceptable-risk
  ## points above actually sit relative to the closed-form H and the
  ## current imposed threshold.
  sweep_K <- sustainable_H_table$K[1]
  sweep_N_target <- 0.8 * sweep_K
  sweep_H_max <- realistic_r * sweep_K * 0.8 * 0.2
  sweep_dt <- tibble::tibble(H = seq(0, sweep_H_max * 1.2, length.out = 18)) %>%
    rowwise() %>%
    mutate(p_collapse = run_from_start(sweep_N_target, H, sweep_K, n_reps = 400,
                                        juv_mortality_ratio = juv_mortality_ratio)$p_collapse) %>%
    ungroup()

  fig_risk_sweep <- file.path(fig_dir, "carrying_capacity_risk_sweep.png")
  ggsave(fig_risk_sweep, width = 8.5, height = 5, dpi = 150, bg = "white", plot = {
    ggplot(sweep_dt, aes(x = H, y = p_collapse)) +
      geom_line(colour = "steelblue", linewidth = 0.9) +
      geom_point(colour = "steelblue", size = 1.8) +
      geom_vline(xintercept = sweep_H_max, linetype = "dashed", colour = "firebrick") +
      geom_hline(data = tibble::tibble(acceptable_risk = acceptable_risk_grid * 100),
                 aes(yintercept = acceptable_risk), linetype = "dotted", colour = "grey40") +
      annotate("text", x = sweep_H_max, y = 102, label = "Closed-form\nSchaefer H", colour = "firebrick",
                size = 3, hjust = 1.05) +
      labs(
        x = "Constant annual removal H", y = "P(falls below quasi-extinction in 25 yr), %",
        title = "Empirically-found risk curve at alpha=0.8 (realistic r, simulation, not closed-form)",
        subtitle = paste0(
          "Dotted horizontals: acceptable-risk conventions (5/10/20%). Dashed red: what the Schaefer\n",
          "formula alone would have called 'sustainable' -- already well into high-risk territory here."
        )
      ) +
      theme_minimal() +
      theme(plot.subtitle = element_text(size = 9.5))
  })

  ## Representative scenario for the stochastic check: alpha = 0.8
  ## (a plausible "near K" ecological-viability target), K high, crossed
  ## with BOTH r choices -- the PBR Rmax (theoretical max, what the
  ## closed-form formula used originally) and the realistic validated-
  ## Leslie r -- to show directly whether the earlier 100% collapse was
  ## about the removal-distribution assumption or about r itself.
  check_alpha <- 0.8
  r_check_scenarios <- tibble::tibble(
    r_label = c("PBR Rmax (theoretical max)", "Realistic r (validated Leslie)"),
    Rmax = c(lambda_max_benchmarks[1] - 1, r_validated)
  )
  stoch_check <- tidyr::expand_grid(project_thresholds, r_check_scenarios) %>%
    mutate(
      K = density_proxy_high * project_footprint_km2,
      H = Rmax * K * check_alpha * (1 - check_alpha),
      N_target = check_alpha * K
    )

  ## Compare three removal-distribution assumptions: the original
  ## prototype's "all on adults" (ratio=0, i.e. zero weight on
  ## juveniles), an equal-risk split (ratio=1), and Paulo's assumption
  ## of somewhat higher juvenile risk (ratio=juv_mortality_ratio).
  ratio_scenarios <- tibble::tibble(
    ratio_label = c("All removal on adults (original)", "Equal risk, juv & adult",
                    paste0("Juveniles ", juv_mortality_ratio, "x adult risk")),
    ratio_value = c(0, 1, juv_mortality_ratio)
  )

  stoch_results <- list()
  traj_plot_rows <- list()
  for (i in seq_len(nrow(stoch_check))) {
    proj <- stoch_check$project[i]; rl_r <- stoch_check$r_label[i]
    K <- stoch_check$K[i]; H <- stoch_check$H[i]; N_target <- stoch_check$N_target[i]
    for (j in seq_len(nrow(ratio_scenarios))) {
      rl <- ratio_scenarios$ratio_label[j]; rv <- ratio_scenarios$ratio_value[j]
      at_target <- run_from_start(N_target, H, K, juv_mortality_ratio = rv)
      at_nmin <- run_from_start(nmin_assumed, H, K, juv_mortality_ratio = rv)
      key <- paste(proj, rl_r, rl)
      stoch_results[[key]] <- list(
        project = proj, r_label = rl_r, ratio_label = rl,
        p_collapse_starting_at_target = at_target$p_collapse,
        p_collapse_starting_at_nmin = at_nmin$p_collapse
      )
      traj_plot_rows[[paste(key, "target")]] <- tibble::tibble(
        project = proj, r_label = rl_r, ratio_label = rl, year = 0:pva_n_years, start = "Starting at N* = alpha*K",
        median = at_target$median_traj, p10 = at_target$p10_traj, p90 = at_target$p90_traj
      )
      traj_plot_rows[[paste(key, "nmin")]] <- tibble::tibble(
        project = proj, r_label = rl_r, ratio_label = rl, year = 0:pva_n_years, start = "Starting at Nmin (today's plausible status)",
        median = at_nmin$median_traj, p10 = at_nmin$p10_traj, p90 = at_nmin$p90_traj
      )
    }
  }
  traj_plot_dt <- bind_rows(traj_plot_rows) %>%
    mutate(ratio_label = factor(ratio_label, levels = ratio_scenarios$ratio_label))
  stoch_summary <- bind_rows(lapply(stoch_results, as.data.frame)) %>%
    tibble::as_tibble() %>%
    mutate(ratio_label = factor(ratio_label, levels = ratio_scenarios$ratio_label)) %>%
    left_join(stoch_check %>% select(project, r_label, K, H, N_target), by = c("project", "r_label")) %>%
    arrange(project, r_label, ratio_label)

  ## Trajectory figure uses the realistic-r scenario only (the
  ## internally-consistent test); the PBR-Rmax scenario's numbers are in
  ## stoch_summary for comparison but are not re-plotted in full --
  ## already shown to collapse regardless of distribution assumption.
  fig_stoch_check <- file.path(fig_dir, "carrying_capacity_stochastic_check.png")
  ggsave(fig_stoch_check, width = 11, height = 7.5, dpi = 150, bg = "white", plot = {
    realistic_r_label <- "Realistic r (validated Leslie)"
    hline_dt <- stoch_check %>% filter(r_label == realistic_r_label) %>% select(project, N_target, K) %>%
      tidyr::crossing(ratio_label = ratio_scenarios$ratio_label) %>%
      mutate(ratio_label = factor(ratio_label, levels = ratio_scenarios$ratio_label))
    ggplot(traj_plot_dt %>% filter(r_label == realistic_r_label), aes(x = year, y = median, colour = start, fill = start)) +
      geom_ribbon(aes(ymin = p10, ymax = p90), alpha = 0.15, colour = NA) +
      geom_line(linewidth = 1) +
      geom_hline(data = hline_dt, aes(yintercept = N_target), linetype = "dotted", colour = "forestgreen") +
      geom_hline(data = hline_dt, aes(yintercept = 0.1 * K), linetype = "dotted", colour = "firebrick") +
      facet_grid(ratio_label ~ project, scales = "free_y") +
      labs(
        x = "Year", y = "Population size (median, 10-90th percentile band)",
        title = paste0("H from the REALISTIC r (alpha=", check_alpha, " target), three removal-distribution assumptions"),
        subtitle = "Green dotted: N* = alpha*K target. Red dotted: quasi-extinction (10% of K)."
      ) +
      theme_minimal() +
      theme(plot.subtitle = element_text(size = 9.5), strip.text = element_text(size = 8.5, face = "bold"), legend.position = "bottom")
  })

  ## ---- 4. Putting the realistic-r choice in perspective: sweep r across
  ## the plausibility range THIS NOTE has already established (not a new
  ## invented range) and find H_sustainable at each point -- shows
  ## directly how much of a biologically plausible fecundity/structure
  ## change it would take for the carrying-capacity standard to converge
  ## with the current PBR threshold (Paulo, 2026-10: "onde conseguimos
  ## chegar a um cenario relativamente real").
  ##
  ## Lever used: litter size. Safi (2006) + Zhigalin & Moskvitina (2017)
  ## give suburban colonies 1.8 pups/female (the value used throughout
  ## this note) vs urban colonies 2.7-2.9 (upper bound), independently
  ## corroborated (references/leslie_matrix_parametrisation.md). Holding
  ## the early-maturity (alpha=1) 2-stage structure fixed (the one this
  ## simulator already uses) and only moving litter size across that
  ## established empirical range moves lambda from 1.05 to ~1.17 (see
  ## R/leslie_vespertilio_murinus.R) -- entirely within evidence already
  ## in this analysis, not a fresh assumption.
  ##
  ## Separately, that same reference note flags that Safi's survival
  ## rates are naive return rates (marked-in-year-t, recaptured-in-t+1),
  ## not corrected for imperfect detection via an open-population model
  ## (e.g. Cormack-Jolly-Seber) -- a well-documented source of downward
  ## bias in survival (and hence lambda) estimates. This means even the
  ## "urban" end of the swept range below is plausibly still an
  ## UNDERESTIMATE of the true rate, not an optimistic ceiling -- noted
  ## here, not quantified (no correction factor for this bias was
  ## located for this species).
  simulate_constant_harvest_litter <- function(n0_juv, n0_adult, n_years, H, K, litter, cv = pva_vital_rate_cv,
                                                juv_mortality_ratio = 1) {
    bp_s_adult <- beta_params(leslie_s_adult_f, cv)
    bp_s_juv <- beta_params(leslie_s_juv, cv)
    bp_p_breed <- beta_params(leslie_p_breed, cv)
    n_juv <- as.integer(round(n0_juv)); n_adult <- as.integer(round(n0_adult))
    traj <- numeric(n_years + 1); traj[1] <- 2 * (n_juv + n_adult)
    for (yr in seq_len(n_years)) {
      s_adult_t <- rbeta(1, bp_s_adult$shape1, bp_s_adult$shape2)
      s_juv_t <- rbeta(1, bp_s_juv$shape1, bp_s_juv$shape2)
      p_breed_t <- rbeta(1, bp_p_breed$shape1, bp_p_breed$shape2)
      density_factor <- max(0, 1 - (2 * (n_juv + n_adult)) / K)
      n_breeders <- rbinom(1, n_adult, p_breed_t)
      n_pups_total <- rpois(1, n_breeders * litter * density_factor)
      n_female_pups <- rbinom(1, n_pups_total, leslie_sex_ratio)
      n_juv_recruits <- rbinom(1, n_female_pups, s_juv_t)
      n_juv_survive <- rbinom(1, n_juv, s_juv_t)
      n_adult_survive <- rbinom(1, n_adult, s_adult_t)
      juv_weight <- n_juv_survive * juv_mortality_ratio
      adult_weight <- n_adult_survive
      total_weight <- juv_weight + adult_weight
      if (total_weight > 0 && H > 0) {
        H_juv <- H * juv_weight / total_weight; H_adult <- H * adult_weight / total_weight
      } else { H_juv <- 0; H_adult <- 0 }
      juv_after_H <- as.integer(round(max(0, n_juv_survive - H_juv)))
      adult_after_H <- as.integer(round(max(0, n_adult_survive - H_adult)))
      n_adult_next <- juv_after_H + adult_after_H
      n_juv <- n_juv_recruits; n_adult <- n_adult_next
      traj[yr + 1] <- 2 * (n_juv + n_adult)
    }
    traj
  }

  find_sustainable_H_litter <- function(litter, N0, K, acceptable_risk = 0.10,
                                         juv_mortality_ratio = 1.3, search_reps = 300,
                                         tol = 0.5, max_iter = 25) {
    F_lit <- leslie_p_breed * litter * leslie_sex_ratio * leslie_s_juv
    A_lit <- build_dekker_stage_matrix_r(2, leslie_s_juv, leslie_s_adult_f, F_lit)
    lambda_lit <- Re(eigen(A_lit)$values[1])
    sd_lit <- Re(eigen(A_lit)$vectors[, 1]); sd_lit <- sd_lit / sum(sd_lit)
    nf <- N0 / 2; nj <- round(nf * sd_lit[1]); na <- nf - nj
    run_risk_lit <- function(H, n_reps) {
      set.seed(pva_seed)
      reps <- replicate(n_reps, simulate_constant_harvest_litter(nj, na, pva_n_years, H, K, litter,
                                                                   juv_mortality_ratio = juv_mortality_ratio))
      mean(reps[pva_n_years + 1, ] < 0.1 * K) * 100 / 100
    }
    H_lo <- 0; H_hi <- K * 0.15
    for (iter in seq_len(max_iter)) {
      H_mid <- (H_lo + H_hi) / 2
      risk <- run_risk_lit(H_mid, search_reps)
      if (risk > acceptable_risk) H_hi <- H_mid else H_lo <- H_mid
      if ((H_hi - H_lo) < tol) break
    }
    tibble::tibble(litter = litter, lambda = lambda_lit, r = log(lambda_lit), H_sustainable = H_lo,
                   p_collapse_confirmed = run_risk_lit(H_lo, pva_n_reps) * 100)
  }

  litter_scenarios <- tibble::tibble(
    litter = c(1.8, 2.1, 2.4, 2.7, 2.9),
    litter_label = c("Suburban (1.8, this note's value)", "", "", "", "Urban upper bound (2.9)")
  )
  K_sweep <- density_proxy_high * project_footprint_km2
  N_target_sweep <- 0.8 * K_sweep
  litter_sweep <- lapply(litter_scenarios$litter, function(lit) {
    find_sustainable_H_litter(lit, N_target_sweep, K_sweep, acceptable_risk = 0.10,
                               juv_mortality_ratio = juv_mortality_ratio)
  }) %>% bind_rows() %>%
    left_join(litter_scenarios, by = "litter") %>%
    tidyr::crossing(project_thresholds) %>%
    mutate(pct_of_current_threshold = 100 * H_sustainable / threshold)

  fig_litter_sweep <- file.path(fig_dir, "carrying_capacity_litter_sweep.png")
  ggsave(fig_litter_sweep, width = 9, height = 5.5, dpi = 150, bg = "white", plot = {
    ggplot(litter_sweep, aes(x = litter, y = H_sustainable)) +
      geom_line(colour = "darkorange", linewidth = 0.9) +
      geom_point(colour = "darkorange", size = 2.2) +
      geom_hline(data = project_thresholds, aes(yintercept = threshold), linetype = "dotted", colour = "grey30") +
      ggrepel::geom_text_repel(data = litter_sweep %>% filter(litter_label != ""),
                                aes(label = litter_label), size = 2.8, nudge_y = 8, seed = 1) +
      facet_wrap(~project) +
      labs(
        x = "Litter size (pups/breeding female/year)", y = "Empirically-sustainable H (alpha=0.8, 10% risk)",
        title = "H sustainable vs. fecundity, across this note's own established plausibility range",
        subtitle = paste0(
          "Dotted: current imposed PBR threshold. Suburban-to-urban litter range (Zhigalin & Moskvitina 2017)\n",
          "is the only lever swept here; Safi survival rates likely underestimate true survival on top of this (not quantified)."
        )
      ) +
      theme_minimal() +
      theme(plot.subtitle = element_text(size = 9))
  })

  ## ---- 5. Elasticity: how far would each vital rate need to move to reach
  ## Frick's lambda_max benchmarks (1.20/1.24), and what carrying-capacity H
  ## does that combination actually sustain once it is run through the real
  ## age-structured stochastic engine (not the closed-form Schaefer formula
  ## used in Section 1)? Paulo, 2026-10: "se subirmos o lambda para 1.2 ou
  ## 1.24 como Frick propoe, como fica o nosso H ... falta-nos um teste de
  ## elasticidade dos parametros do carrying capacity".
  ##
  ## Reuses the SAME baseline vital rates, plausible ranges and univariate-
  ## breakeven method as R/leslie_boundary_analysis.R (that script's own
  ## range_s_adult/range_s_juv/range_p_breed/range_litter, now promoted to
  ## inputs/pbrSettings_BSH_DGY.R as leslie_boundary_range_* / s_range) --
  ## not a fresh set of assumptions. What is new here: each breakeven vital-
  ## rate combination is fed into this script's own stochastic H-finder
  ## (same engine as Sections 2-4), so "lambda=1.2/1.24" is translated into
  ## an actual simulated population trajectory and a reportable H, instead
  ## of being plugged into r*K*alpha*(1-alpha) as if the age structure and
  ## maturation lag did not matter (Section 0 already showed that shortcut
  ## is what drove the earlier 100% collapse finding).
  build_leslie_general <- function(s_juv, s_adult, p_breed, litter) {
    F_rec <- p_breed * litter * leslie_sex_ratio * s_juv
    build_dekker_stage_matrix_r(2, s_juv, s_adult, F_rec)
  }
  lambda_of_vitals <- function(s_juv, s_adult, p_breed, litter) {
    Re(eigen(build_leslie_general(s_juv, s_adult, p_breed, litter))$values[1])
  }
  stable_dist_of_vitals <- function(s_juv, s_adult, p_breed, litter) {
    v <- Re(eigen(build_leslie_general(s_juv, s_adult, p_breed, litter))$vectors[, 1])
    v / sum(v)
  }
  baseline_vitals <- list(s_juv = leslie_s_juv, s_adult = leslie_s_adult_f,
                           p_breed = leslie_p_breed, litter = leslie_litter)
  stopifnot(abs(do.call(lambda_of_vitals, baseline_vitals) - lambda_leslie_validated_cc) < 1e-9)

  ## Same discrete-time engine as simulate_constant_harvest()/
  ## simulate_constant_harvest_litter() above, generalised to take all four
  ## vital rates as arguments rather than reading them (or just litter)
  ## from the global baseline -- so any plausible combination, not only a
  ## single-lever sweep, can be run through it.
  simulate_constant_harvest_general <- function(n0_juv, n0_adult, n_years, H, K,
                                                  s_juv, s_adult, p_breed, litter,
                                                  cv = pva_vital_rate_cv, juv_mortality_ratio = 1) {
    bp_s_adult <- beta_params(s_adult, cv)
    bp_s_juv <- beta_params(s_juv, cv)
    bp_p_breed <- beta_params(p_breed, cv)
    n_juv <- as.integer(round(n0_juv)); n_adult <- as.integer(round(n0_adult))
    traj <- numeric(n_years + 1); traj[1] <- 2 * (n_juv + n_adult)
    for (yr in seq_len(n_years)) {
      s_adult_t <- rbeta(1, bp_s_adult$shape1, bp_s_adult$shape2)
      s_juv_t <- rbeta(1, bp_s_juv$shape1, bp_s_juv$shape2)
      p_breed_t <- rbeta(1, bp_p_breed$shape1, bp_p_breed$shape2)
      density_factor <- max(0, 1 - (2 * (n_juv + n_adult)) / K)
      n_breeders <- rbinom(1, n_adult, p_breed_t)
      n_pups_total <- rpois(1, n_breeders * litter * density_factor)
      n_female_pups <- rbinom(1, n_pups_total, leslie_sex_ratio)
      n_juv_recruits <- rbinom(1, n_female_pups, s_juv_t)
      n_juv_survive <- rbinom(1, n_juv, s_juv_t)
      n_adult_survive <- rbinom(1, n_adult, s_adult_t)
      juv_weight <- n_juv_survive * juv_mortality_ratio
      adult_weight <- n_adult_survive
      total_weight <- juv_weight + adult_weight
      if (total_weight > 0 && H > 0) {
        H_juv <- H * juv_weight / total_weight; H_adult <- H * adult_weight / total_weight
      } else { H_juv <- 0; H_adult <- 0 }
      juv_after_H <- as.integer(round(max(0, n_juv_survive - H_juv)))
      adult_after_H <- as.integer(round(max(0, n_adult_survive - H_adult)))
      n_juv <- n_juv_recruits; n_adult <- juv_after_H + adult_after_H
      traj[yr + 1] <- 2 * (n_juv + n_adult)
    }
    traj
  }

  find_sustainable_H_general <- function(N0, K, s_juv, s_adult, p_breed, litter,
                                          acceptable_risk = 0.10, juv_mortality_ratio = 1.3,
                                          search_reps = 300, confirm_reps = pva_n_reps,
                                          H_upper = NULL, tol = 0.5, max_iter = 25) {
    sd_v <- stable_dist_of_vitals(s_juv, s_adult, p_breed, litter)
    nf <- N0 / 2; nj <- round(nf * sd_v[1]); na <- nf - nj
    run_risk <- function(H, n_reps) {
      set.seed(pva_seed)
      reps <- replicate(n_reps, simulate_constant_harvest_general(
        nj, na, pva_n_years, H, K, s_juv, s_adult, p_breed, litter,
        juv_mortality_ratio = juv_mortality_ratio))
      mean(reps[pva_n_years + 1, ] < 0.1 * K)
    }
    if (is.null(H_upper)) H_upper <- K * 0.2
    H_lo <- 0; H_hi <- H_upper
    for (iter in seq_len(max_iter)) {
      H_mid <- (H_lo + H_hi) / 2
      risk <- run_risk(H_mid, search_reps)
      if (risk > acceptable_risk) H_hi <- H_mid else H_lo <- H_mid
      if ((H_hi - H_lo) < tol) break
    }
    tibble::tibble(H_sustainable = H_lo, p_collapse_confirmed = run_risk(H_lo, confirm_reps) * 100)
  }

  ## 5a. Univariate breakeven: value each vital rate alone would need to
  ## reach lambda = 1.20 / 1.24, holding the other three at Safi's (2006)
  ## baseline -- same method as R/leslie_boundary_analysis.R's
  ## solve_breakeven(), not re-derived differently here.
  vital_param_labels <- c(s_adult = "Adult female survival", s_juv = "Juvenile survival",
                           p_breed = "Breeding fraction", litter = "Litter size")
  vital_plausible_ranges <- list(s_adult = s_range, s_juv = leslie_boundary_range_s_juv,
                                  p_breed = leslie_boundary_range_p_breed, litter = leslie_boundary_range_litter)
  vital_search_ceiling <- list(s_adult = 0.999, s_juv = 0.99, p_breed = 1.0, litter = 6)

  solve_breakeven_vital <- function(param, target, search_upper) {
    f <- function(x) {
      args <- baseline_vitals; args[[param]] <- x
      do.call(lambda_of_vitals, args) - target
    }
    lo <- baseline_vitals[[param]]
    if (f(lo) >= 0) return(lo)
    if (f(search_upper) < 0) return(NA_real_)
    uniroot(f, c(lo, search_upper))$root
  }

  K_elastic <- density_proxy_high * project_footprint_km2
  N_target_elastic <- check_alpha * K_elastic  # same alpha=0.8 scenario as Sections 2-4

  breakeven_H_table <- tidyr::expand_grid(param = names(baseline_vitals), lambda_target = lambda_max_benchmarks) %>%
    rowwise() %>%
    mutate(
      parameter = vital_param_labels[param],
      baseline_value = baseline_vitals[[param]],
      required_value = solve_breakeven_vital(param, lambda_target, vital_search_ceiling[[param]]),
      plausible_lower = vital_plausible_ranges[[param]][1],
      plausible_upper = vital_plausible_ranges[[param]][2],
      within_plausible_range = !is.na(required_value) & required_value <= plausible_upper
    ) %>%
    ungroup()

  run_breakeven_row <- function(param, required_value) {
    if (is.na(required_value)) return(tibble::tibble(H_sustainable = NA_real_, p_collapse_confirmed = NA_real_))
    v <- baseline_vitals; v[[param]] <- required_value
    find_sustainable_H_general(N_target_elastic, K_elastic, s_juv = v$s_juv, s_adult = v$s_adult,
                                p_breed = v$p_breed, litter = v$litter,
                                acceptable_risk = 0.10, juv_mortality_ratio = juv_mortality_ratio)
  }
  breakeven_H_table <- breakeven_H_table %>%
    rowwise() %>%
    mutate(H_result = list(run_breakeven_row(param, required_value))) %>%
    ungroup() %>%
    tidyr::unnest(H_result) %>%
    tidyr::crossing(project_thresholds) %>%
    mutate(pct_of_current_threshold = round(100 * H_sustainable / threshold))

  ## 5b. Full one-at-a-time sweep across each vital rate's plausible range
  ## (6 points, others held at baseline) -- the elasticity of simulated H,
  ## not just the single breakeven point.
  run_vital_sweep <- function(param) {
    rng <- vital_plausible_ranges[[param]]
    values <- sort(unique(c(seq(rng[1], rng[2], length.out = 6), baseline_vitals[[param]])))
    sweep_rows <- lapply(values, function(val) {
      v <- baseline_vitals; v[[param]] <- val
      lambda_v <- do.call(lambda_of_vitals, v)
      H_res <- find_sustainable_H_general(N_target_elastic, K_elastic, s_juv = v$s_juv, s_adult = v$s_adult,
                                           p_breed = v$p_breed, litter = v$litter,
                                           acceptable_risk = 0.10, juv_mortality_ratio = juv_mortality_ratio,
                                           confirm_reps = 500)
      tibble::tibble(param = param, value = val, lambda = lambda_v, H_sustainable = H_res$H_sustainable)
    })
    bind_rows(sweep_rows)
  }
  elasticity_sweep <- bind_rows(lapply(names(baseline_vitals), run_vital_sweep)) %>%
    mutate(parameter = vital_param_labels[param], baseline_value = unlist(baseline_vitals[param]))

  fig_elasticity_vitals <- file.path(fig_dir, "carrying_capacity_elasticity_vitals.png")
  ggsave(fig_elasticity_vitals, width = 11, height = 5, dpi = 150, bg = "white", plot = {
    hline_dt <- tidyr::crossing(project_thresholds, parameter = unique(elasticity_sweep$parameter))
    ggplot(elasticity_sweep, aes(x = value, y = H_sustainable)) +
      geom_hline(data = hline_dt, aes(yintercept = threshold, colour = project), linetype = "dotted", linewidth = 0.6) +
      geom_vline(aes(xintercept = baseline_value), linetype = "dashed", colour = "grey50") +
      geom_line(colour = "darkorange", linewidth = 0.9) +
      geom_point(aes(fill = lambda), shape = 21, size = 2.6, colour = "black") +
      scale_fill_viridis_c(option = "D", name = "lambda") +
      scale_colour_manual(name = "Current imposed\nPBR threshold", values = c("grey20", "grey50")) +
      facet_wrap(~parameter, scales = "free_x", nrow = 1) +
      labs(
        x = "Vital-rate value swept (others held at Safi 2006 baseline)",
        y = "Simulated sustainable H (alpha=0.8, 10% risk)",
        title = "Elasticity of simulated carrying-capacity H to each vital rate",
        subtitle = paste0(
          "Dashed grey vertical: Safi (2006) baseline value. Dotted horizontals: current imposed PBR thresholds.\n",
          "Point colour: resulting lambda. Each panel varies one rate across its literature-plausible range; others fixed."
        )
      ) +
      theme_minimal() +
      theme(plot.subtitle = element_text(size = 9), strip.text = element_text(face = "bold"), legend.position = "bottom")
  })

  ## 5c. Tornado: rank ALL the parameters this note has treated as
  ## uncertain -- the four vital rates AND the carrying-capacity-specific
  ## choices (alpha, acceptable_risk, K, juv_mortality_ratio) -- by how much
  ## moving each one alone, low vs. high across the SAME ranges already
  ## used elsewhere in this script (not new bounds invented for this
  ## figure), shifts simulated H away from the alpha=0.8/K-high/10%-risk
  ## reference point. Vital rates have a literature-grounded plausible
  ## range (flagged "literature range"); the other four are this analysis'
  ## own structural/policy choices, not biological estimates (flagged
  ## "modelling choice") -- the two kinds of uncertainty are not
  ## interchangeable and the figure keeps them visually distinct.
  tornado_reference <- find_sustainable_H_general(
    N_target_elastic, K_elastic, s_juv = leslie_s_juv, s_adult = leslie_s_adult_f,
    p_breed = leslie_p_breed, litter = leslie_litter,
    acceptable_risk = 0.10, juv_mortality_ratio = juv_mortality_ratio
  )$H_sustainable

  tornado_vital <- function(param) {
    rng <- vital_plausible_ranges[[param]]
    H_each <- sapply(rng, function(val) {
      v <- baseline_vitals; v[[param]] <- val
      find_sustainable_H_general(N_target_elastic, K_elastic, s_juv = v$s_juv, s_adult = v$s_adult,
                                  p_breed = v$p_breed, litter = v$litter,
                                  acceptable_risk = 0.10, juv_mortality_ratio = juv_mortality_ratio,
                                  confirm_reps = 500)$H_sustainable
    })
    tibble::tibble(parameter = vital_param_labels[param], kind = "Literature range",
                    H_low = H_each[1], H_high = H_each[2])
  }
  tornado_cc <- tibble::tribble(
    ~parameter, ~H_low_args, ~H_high_args,
    "Target distance from K (alpha)", list(N0_mult = 0.5), list(N0_mult = 0.9),
    "Acceptable quasi-extinction risk", list(risk = min(acceptable_risk_grid)), list(risk = max(acceptable_risk_grid)),
    "Carrying capacity (K)", list(K_use = "low"), list(K_use = "high"),
    "Juvenile:adult mortality risk ratio", list(ratio = 0), list(ratio = max(ratio_scenarios$ratio_value))
  )
  tornado_cc_row <- function(parameter, args_low, args_high) {
    run_one <- function(a) {
      N0_mult <- if (!is.null(a$N0_mult)) a$N0_mult else check_alpha
      risk <- if (!is.null(a$risk)) a$risk else 0.10
      K_use <- if (!is.null(a$K_use) && a$K_use == "low") density_proxy_low * project_footprint_km2 else K_elastic
      ratio <- if (!is.null(a$ratio)) a$ratio else juv_mortality_ratio
      find_sustainable_H_general(N0_mult * K_use, K_use, s_juv = leslie_s_juv, s_adult = leslie_s_adult_f,
                                  p_breed = leslie_p_breed, litter = leslie_litter,
                                  acceptable_risk = risk, juv_mortality_ratio = ratio,
                                  confirm_reps = 500)$H_sustainable
    }
    tibble::tibble(parameter = parameter, kind = "Modelling choice", H_low = run_one(args_low), H_high = run_one(args_high))
  }
  tornado_table <- bind_rows(
    bind_rows(lapply(names(baseline_vitals), tornado_vital)),
    bind_rows(lapply(seq_len(nrow(tornado_cc)), function(i) {
      tornado_cc_row(tornado_cc$parameter[i], tornado_cc$H_low_args[[i]], tornado_cc$H_high_args[[i]])
    }))
  ) %>%
    mutate(
      H_reference = tornado_reference,
      range_width = abs(H_high - H_low)
    ) %>%
    arrange(desc(range_width)) %>%
    mutate(parameter = factor(parameter, levels = rev(parameter)))

  fig_tornado <- file.path(fig_dir, "carrying_capacity_tornado.png")
  ggsave(fig_tornado, width = 9, height = 5.5, dpi = 150, bg = "white", plot = {
    tornado_long <- tornado_table %>%
      tidyr::pivot_longer(c(H_low, H_high), names_to = "end", values_to = "H")
    ggplot(tornado_table) +
      geom_segment(aes(x = pmin(H_low, H_high), xend = pmax(H_low, H_high), y = parameter, yend = parameter, colour = kind),
                   linewidth = 5, alpha = 0.6) +
      geom_point(data = tornado_long, aes(x = H, y = parameter, colour = kind), size = 2.5) +
      geom_vline(aes(xintercept = H_reference), linetype = "dashed", colour = "grey30") +
      scale_colour_manual(name = NULL, values = c("Literature range" = "steelblue", "Modelling choice" = "firebrick")) +
      annotate("text", x = tornado_reference, y = Inf, label = "Reference\n(baseline vitals,\nalpha=0.8, 10% risk)",
               vjust = 1.3, hjust = -0.05, size = 2.9, colour = "grey30") +
      labs(
        x = "Simulated sustainable H, varying ONE parameter low to high (others at reference)",
        y = NULL,
        title = "Which parameter moves the carrying-capacity H the most?",
        subtitle = "Bar = low-to-high H across each parameter's own plausible range (same ranges used elsewhere in this script).\nBlue = literature-grounded vital rate; red = this analysis' own structural/policy choice, not a biological estimate."
      ) +
      theme_minimal() +
      theme(plot.subtitle = element_text(size = 9), legend.position = "bottom")
  })

  ## ---- 6. Back to the Section 3 risk-vs-H curve, now at the ONE
  ## single-parameter breakeven combination Section 5 found to be within a
  ## literature-plausible range (S_adult = breakeven value for lambda=1.20;
  ## S_juv, p_breed, litter left at Safi 2006 baseline) -- Paulo, 2026-10:
  ## "se voltarmos a simulacao do Prob quasi-extinc...com o cenario que nos
  ## deu plausivel, o que acontece a curva?". Same engine, same K, same
  ## alpha=0.8 target as the original realistic-r curve -- only the vital
  ## rates differ -- so the two curves are a direct, apples-to-apples
  ## before/after of what raising survival toward the top of its plausible
  ## range (not swapping in a bare Rmax) buys in removal headroom.
  plausible_lambda_target <- min(lambda_max_benchmarks)  # 1.20 -- the only breakeven that stayed in-range
  plausible_breakeven_row <- breakeven_H_table %>%
    filter(param == "s_adult", lambda_target == plausible_lambda_target) %>%
    slice(1)
  plausible_s_adult <- plausible_breakeven_row$required_value
  plausible_H_sustainable <- plausible_breakeven_row$H_sustainable
  plausible_vitals <- baseline_vitals
  plausible_vitals$s_adult <- plausible_s_adult

  run_risk_general <- function(N0, H, K, s_juv, s_adult, p_breed, litter, juv_mortality_ratio,
                                n_reps = 400, seed = pva_seed) {
    sd_v <- stable_dist_of_vitals(s_juv, s_adult, p_breed, litter)
    nf <- N0 / 2; nj <- round(nf * sd_v[1]); na <- nf - nj
    set.seed(seed)
    reps <- replicate(n_reps, simulate_constant_harvest_general(
      nj, na, pva_n_years, H, K, s_juv, s_adult, p_breed, litter, juv_mortality_ratio = juv_mortality_ratio))
    mean(reps[pva_n_years + 1, ] < 0.1 * K) * 100
  }

  sweep_H_max_plausible <- plausible_H_sustainable * 1.3
  sweep_dt_plausible <- tibble::tibble(H = seq(0, sweep_H_max_plausible, length.out = 18)) %>%
    rowwise() %>%
    mutate(p_collapse = run_risk_general(sweep_N_target, H, sweep_K, s_juv = plausible_vitals$s_juv,
                                          s_adult = plausible_vitals$s_adult, p_breed = plausible_vitals$p_breed,
                                          litter = plausible_vitals$litter, juv_mortality_ratio = juv_mortality_ratio,
                                          n_reps = 400)) %>%
    ungroup()

  scenario_realistic <- "Realistic r (validated Leslie, Safi 2006 baseline vitals)"
  scenario_plausible <- paste0("S_adult=", round(plausible_s_adult, 3), " (breakeven for lambda=1.20, within literature range)")
  sweep_dt_combined <- bind_rows(
    sweep_dt %>% mutate(scenario = scenario_realistic),
    sweep_dt_plausible %>% mutate(scenario = scenario_plausible)
  )

  fig_risk_sweep_plausible <- file.path(fig_dir, "carrying_capacity_risk_sweep_plausible.png")
  ggsave(fig_risk_sweep_plausible, width = 9.5, height = 5.5, dpi = 150, bg = "white", plot = {
    ggplot(sweep_dt_combined, aes(x = H, y = p_collapse, colour = scenario)) +
      geom_line(linewidth = 0.9) +
      geom_point(size = 1.8) +
      geom_vline(xintercept = plausible_H_sustainable, linetype = "dashed", colour = "darkorange") +
      geom_vline(data = project_thresholds, aes(xintercept = threshold), linetype = "dotted", colour = "grey40") +
      geom_hline(data = tibble::tibble(acceptable_risk = acceptable_risk_grid * 100),
                 aes(yintercept = acceptable_risk), linetype = "dotted", colour = "grey40") +
      scale_colour_manual(name = NULL, values = setNames(c("steelblue", "darkorange"), c(scenario_realistic, scenario_plausible))) +
      annotate("text", x = plausible_H_sustainable, y = 102, label = "H sustainable\n(10% risk) at\nS_adult=0.95",
               colour = "darkorange", size = 2.8, hjust = -0.05) +
      labs(
        x = "Constant annual removal H", y = "P(falls below quasi-extinction in 25 yr), %",
        title = "Risk-vs-H curve: realistic-r baseline vs. the one literature-plausible lambda=1.20 route",
        subtitle = paste0(
          "Dotted vertical grey: current imposed PBR thresholds (", paste(project_thresholds$project, project_thresholds$threshold, sep = "=", collapse = ", "), ").\n",
          "Dotted horizontal: acceptable-risk conventions (5/10/20%). Raising S_adult to its breakeven value shifts the\n",
          "whole curve right, but the 10%-risk H it buys (", round(plausible_H_sustainable), ") still sits ",
          round(100 * plausible_H_sustainable / min(project_thresholds$threshold)), "-",
          round(100 * plausible_H_sustainable / max(project_thresholds$threshold)), "% of the current thresholds."
        )
      ) +
      theme_minimal() +
      theme(plot.subtitle = element_text(size = 8.7), legend.position = "bottom")
  })

  ## ---- 7. Joint (not one-at-a-time) Monte Carlo over the family/order-
  ## informed biological envelope (references/leslie_matrix_parametrisation.md,
  ## "Family/order-level envelope" section) -- Paulo, 2026-10: "exigir que
  ## o PBR ou H seja sempre apresentado com intervalo plausivel ... nós
  ## poderiamos buscar parametros biologicos de especies proximas ... para
  ## sugerir a plausibilidade do envelope biologico". Section 5 moved ONE
  ## vital rate at a time, holding the other three at Safi's point values --
  ## useful for attributing sensitivity, but not a genuine uncertainty
  ## interval (the four rates are not independently pinned down in reality
  ## either, and a real population draws all four jointly, not one-at-a-
  ## time). This section draws all four TOGETHER, every draw run through
  ## the same stochastic engine as the rest of this script, and reports H
  ## as a distribution, not a point.
  ##
  ## Distributions: s_adult centred on Safi's 0.76 (independently confirmed
  ## close to the Lentini et al. 2015 order-wide meta-analytic mean of
  ## 0.774, 44 species/7 families) with spread matching that paper's own
  ## 95% CI [0.617, 0.890]; s_juv centred on Safi's 0.62 with spread
  ## matching the confamilial range found across *M. lucifugus*,
  ## *M. daubentonii* and *E. fuscus* (0.23-0.71); p_breed centred on
  ## Safi's 0.87 with spread matching *E. fuscus*'s first-time-breeder to
  ## experienced-breeder range (0.64-0.98, this note's own established
  ## ceiling); litter centred on 1.8 (suburban -- the Uzbekistan site is a
  ## natural/suburban setting, not urban, per Paulo) via a Beta rescaled
  ## onto V. murinus's own documented suburban-urban span (1.3-2.9), so
  ## some mass still reaches toward the urban figure without the mean
  ## being pulled up to it -- deliberately not widened using confamilial
  ## litter data either, which points toward smaller litters for most
  ## Vespertilionidae, not larger (see the references note).
  envelope_dist <- function(mean, half_width_95) {
    sd_eq <- half_width_95 / 1.96
    bp <- beta_params(mean, sd_eq / mean)
    function(n) rbeta(n, bp$shape1, bp$shape2)
  }
  draw_s_adult <- envelope_dist(leslie_s_adult_f, (0.890 - 0.617) / 2)
  draw_s_juv   <- envelope_dist(leslie_s_juv, (0.71 - 0.23) / 2)
  draw_p_breed <- envelope_dist(leslie_p_breed, (0.98 - 0.64) / 2)

  ## Litter: centred on 1.8 (suburban), not spread uniformly up to the
  ## urban figure (2.9) -- Paulo, 2026-10: the Uzbekistan site is clearly
  ## a natural/suburban setting, not urban, so weighting the distribution
  ## toward urban-level fecundity would be the LESS cautious assumption
  ## (it inflates lambda and hence H). Kept as a Beta rescaled onto the
  ## full literature-documented span [1.3, 2.9] -- not a narrower
  ## arbitrary band -- so some probability mass still reaches toward the
  ## urban figure (colony classification is not certain), but the mean
  ## and most of the mass sit at the suburban value this note uses
  ## throughout.
  litter_lo <- leslie_boundary_range_litter[1]; litter_hi <- leslie_boundary_range_litter[2]
  litter_p_mean <- (leslie_litter - litter_lo) / (litter_hi - litter_lo)
  litter_bp <- beta_params(litter_p_mean, 0.12 / litter_p_mean)
  draw_litter <- function(n) litter_lo + (litter_hi - litter_lo) * rbeta(n, litter_bp$shape1, litter_bp$shape2)

  n_envelope_draws <- 200
  set.seed(pva_seed)
  envelope_draws <- tibble::tibble(
    draw_id = seq_len(n_envelope_draws),
    s_adult = draw_s_adult(n_envelope_draws),
    s_juv   = draw_s_juv(n_envelope_draws),
    p_breed = draw_p_breed(n_envelope_draws),
    litter  = draw_litter(n_envelope_draws)
  ) %>%
    rowwise() %>%
    mutate(lambda = lambda_of_vitals(s_juv, s_adult, p_breed, litter)) %>%
    ungroup()

  ## 7a. Distribution of the empirically-sustainable H (10% risk, alpha=0.8)
  ## across the joint envelope -- cheaper search settings than Section 5's
  ## single breakeven points (this runs 200x), still reportable.
  envelope_draws <- envelope_draws %>%
    rowwise() %>%
    mutate(H_sustainable = find_sustainable_H_general(
      N_target_elastic, K_elastic, s_juv = s_juv, s_adult = s_adult, p_breed = p_breed, litter = litter,
      acceptable_risk = 0.10, juv_mortality_ratio = juv_mortality_ratio,
      search_reps = 120, confirm_reps = 1, max_iter = 12
    )$H_sustainable) %>%
    ungroup()

  H_envelope_percentiles <- tibble::tibble(
    percentile = c("5th", "25th", "median", "75th", "95th"),
    H_sustainable = quantile(envelope_draws$H_sustainable, probs = c(0.05, 0.25, 0.5, 0.75, 0.95))
  ) %>%
    tidyr::crossing(project_thresholds) %>%
    mutate(pct_of_current_threshold = round(100 * H_sustainable / threshold))

  fig_envelope_H <- file.path(fig_dir, "carrying_capacity_envelope_H.png")
  ggsave(fig_envelope_H, width = 9, height = 5, dpi = 150, bg = "white", plot = {
    ggplot(envelope_draws, aes(x = H_sustainable)) +
      geom_histogram(bins = 25, fill = "steelblue", alpha = 0.75, colour = "white") +
      geom_vline(data = project_thresholds, aes(xintercept = threshold, linetype = project), colour = "grey20") +
      geom_vline(xintercept = median(envelope_draws$H_sustainable), colour = "darkorange", linewidth = 0.9) +
      annotate("text", x = median(envelope_draws$H_sustainable), y = Inf,
               label = paste0("Median = ", round(median(envelope_draws$H_sustainable))),
               colour = "darkorange", vjust = 1.4, hjust = -0.05, size = 3) +
      scale_linetype_manual(name = "Current imposed\nPBR threshold", values = c("dashed", "dotted")) +
      labs(
        x = "Empirically-sustainable H (alpha=0.8, 10% risk), one joint draw per bar",
        y = paste0("Count (of ", n_envelope_draws, " joint draws)"),
        title = "Carrying-capacity H as a distribution, not a point: joint draws over the family/order-informed envelope",
        subtitle = paste0(
          "Each draw moves S_adult, S_juv, p_breed and litter TOGETHER (not one-at-a-time as in the breakeven table),\n",
          "within ranges grounded in Lentini et al. (2015), Sendor & Simon (2003), Frick et al. (2010), O'Shea et al. (2010)\n",
          "and this note's own V. murinus-specific sources (references/leslie_matrix_parametrisation.md)."
        )
      ) +
      theme_minimal() +
      theme(plot.subtitle = element_text(size = 8.3))
  })

  ## 7b. The flip side of the same question, and cheaper to compute: fixing
  ## H at the CURRENT imposed PBR threshold, what does quasi-extinction risk
  ## look like across the same joint envelope? Directly answers whether
  ## today's policy is safe once the biological uncertainty this note has
  ## now documented is accounted for, not just at the Safi point estimate.
  run_risk_at_H <- function(s_juv, s_adult, p_breed, litter, H, K, N0, n_reps) {
    sd_v <- stable_dist_of_vitals(s_juv, s_adult, p_breed, litter)
    nf <- N0 / 2; nj <- round(nf * sd_v[1]); na <- nf - nj
    reps <- replicate(n_reps, simulate_constant_harvest_general(
      nj, na, pva_n_years, H, K, s_juv, s_adult, p_breed, litter, juv_mortality_ratio = juv_mortality_ratio))
    mean(reps[pva_n_years + 1, ] < 0.1 * K) * 100
  }
  risk_at_pbr_dt <- tidyr::crossing(
    envelope_draws %>% select(draw_id, s_juv, s_adult, p_breed, litter), project_thresholds
  ) %>%
    rowwise() %>%
    mutate(p_collapse_at_pbr = run_risk_at_H(s_juv, s_adult, p_breed, litter, H = threshold, K = K_elastic,
                                              N0 = N_target_elastic, n_reps = 150)) %>%
    ungroup()

  risk_at_pbr_summary <- risk_at_pbr_dt %>%
    group_by(project, threshold) %>%
    summarise(
      pct_draws_above_10pct_risk = round(100 * mean(p_collapse_at_pbr > 10)),
      pct_draws_above_20pct_risk = round(100 * mean(p_collapse_at_pbr > 20)),
      median_p_collapse = round(median(p_collapse_at_pbr), 1),
      .groups = "drop"
    )

  fig_envelope_risk_at_pbr <- file.path(fig_dir, "carrying_capacity_envelope_risk_at_pbr.png")
  ggsave(fig_envelope_risk_at_pbr, width = 9, height = 5, dpi = 150, bg = "white", plot = {
    ggplot(risk_at_pbr_dt, aes(x = p_collapse_at_pbr)) +
      geom_histogram(bins = 25, fill = "firebrick", alpha = 0.75, colour = "white") +
      geom_vline(data = tibble::tibble(acceptable_risk = acceptable_risk_grid * 100),
                 aes(xintercept = acceptable_risk), linetype = "dotted", colour = "grey30") +
      facet_wrap(~project) +
      labs(
        x = "P(falls below quasi-extinction in 25 yr), %, removal fixed at the current PBR threshold",
        y = paste0("Count (of ", n_envelope_draws, " joint draws)"),
        title = "Is the CURRENT PBR threshold safe across the biological envelope, not just at Safi's point estimate?",
        subtitle = "Dotted verticals: acceptable-risk conventions (5/10/20%). Same joint draws as the H-distribution figure above."
      ) +
      theme_minimal() +
      theme(plot.subtitle = element_text(size = 9))
  })

  ## ---- 8. Back to PVA-lite itself (R/pva_robustness.R's Section 4/7.5
  ## baseline: N0 = nmin_assumed, NOT alpha*K; annual removal = each
  ## project's own imposed PBR threshold; NO density dependence in the
  ## baseline scenario -- that file's own comment flags unchecked 25-year
  ## growth as "not biologically credible", but keeps it as the baseline
  ## because Section 4/7.5 elsewhere in the real report use it that way).
  ## Paulo, 2026-10: "se voltassemos a PVA lite agora com essa
  ## parametrizacao, o que aconteceria a tendencia a 25 anos?"
  ##
  ## PVA-lite as it stands draws only environmental/demographic noise each
  ## year around FIXED point vital rates (Safi's values, 10% CV). This
  ## reruns its exact N0/removal/no-K setup, but every trajectory's vital
  ## rates are first drawn from the SAME joint family/order-informed
  ## envelope as Section 7 (litter centred on 1.8) and held fixed for that
  ## trajectory's 25 years, with year-to-year environmental noise layered
  ## on top as before -- so the resulting band pools parameter uncertainty
  ## AND environmental stochasticity, not environmental noise alone.
  ## K = Inf disables Section 7's density-dependent compensation (the CC
  ## engine's density_factor is 1 - 2N/K, which -> 1 as K -> Inf), matching
  ## PVA-lite's own "no density dependence" baseline exactly.
  stable_dist_baseline_pvalite <- stable_dist_of_vitals(leslie_s_juv, leslie_s_adult_f, leslie_p_breed, leslie_litter)
  n0_juv_pvalite <- round((nmin_assumed / 2) * stable_dist_baseline_pvalite[1])
  n0_adult_pvalite <- (nmin_assumed / 2) - n0_juv_pvalite
  quasi_ext_threshold_pvalite <- pva_quasi_extinction_fraction * nmin_assumed

  pvalite_scenario_baseline <- "PVA-lite baseline (fixed Safi point vitals)"
  pvalite_scenario_envelope <- "+ joint biological envelope (this note's broader uncertainty)"

  pvalite_baseline_traj <- tidyr::crossing(project_thresholds, rep_id = seq_len(pva_n_reps)) %>%
    rowwise() %>%
    mutate(traj = list(simulate_constant_harvest_general(
      n0_juv_pvalite, n0_adult_pvalite, pva_n_years, H = threshold, K = Inf,
      s_juv = leslie_s_juv, s_adult = leslie_s_adult_f, p_breed = leslie_p_breed, litter = leslie_litter,
      juv_mortality_ratio = juv_mortality_ratio
    ))) %>%
    ungroup() %>%
    mutate(scenario = pvalite_scenario_baseline)

  reps_per_envelope_draw <- 10
  pvalite_envelope_traj <- tidyr::crossing(
    envelope_draws %>% select(draw_id, s_juv, s_adult, p_breed, litter),
    project_thresholds,
    rep_id = seq_len(reps_per_envelope_draw)
  ) %>%
    rowwise() %>%
    mutate(traj = list(simulate_constant_harvest_general(
      n0_juv_pvalite, n0_adult_pvalite, pva_n_years, H = threshold, K = Inf,
      s_juv = s_juv, s_adult = s_adult, p_breed = p_breed, litter = litter,
      juv_mortality_ratio = juv_mortality_ratio
    ))) %>%
    ungroup() %>%
    mutate(scenario = pvalite_scenario_envelope)

  pvalite_combined_traj <- bind_rows(
    pvalite_baseline_traj %>% select(project, threshold, rep_id, traj, scenario),
    pvalite_envelope_traj %>% select(project, threshold, rep_id, traj, scenario)
  )

  pvalite_trend_summary <- pvalite_combined_traj %>%
    rowwise() %>%
    mutate(year = list(0:pva_n_years)) %>%
    ungroup() %>%
    tidyr::unnest(c(traj, year)) %>%
    rename(N = traj) %>%
    group_by(scenario, project, year) %>%
    summarise(median = median(N), p10 = quantile(N, 0.10), p90 = quantile(N, 0.90), .groups = "drop") %>%
    mutate(scenario = factor(scenario, levels = c(pvalite_scenario_baseline, pvalite_scenario_envelope)))

  pvalite_trend_risk <- pvalite_combined_traj %>%
    rowwise() %>%
    mutate(final_N = traj[length(traj)], ever_below_quasi_ext = any(traj < quasi_ext_threshold_pvalite)) %>%
    ungroup() %>%
    group_by(scenario, project) %>%
    summarise(
      n_trajectories = dplyr::n(),
      p_decline = round(100 * mean(final_N < nmin_assumed)),
      p_quasi_extinction = round(100 * mean(ever_below_quasi_ext)),
      median_final_N = round(median(final_N)),
      .groups = "drop"
    ) %>%
    mutate(scenario = factor(scenario, levels = c(pvalite_scenario_baseline, pvalite_scenario_envelope)))

  fig_pvalite_envelope <- file.path(fig_dir, "carrying_capacity_pvalite_envelope.png")
  ggsave(fig_pvalite_envelope, width = 10.5, height = 6.5, dpi = 150, bg = "white", plot = {
    ggplot(pvalite_trend_summary, aes(x = year, y = median, colour = scenario, fill = scenario)) +
      geom_ribbon(aes(ymin = p10, ymax = p90), alpha = 0.18, colour = NA) +
      geom_line(linewidth = 1) +
      geom_hline(yintercept = nmin_assumed, linetype = "dotted", colour = "grey30") +
      geom_hline(yintercept = quasi_ext_threshold_pvalite, linetype = "dotted", colour = "firebrick") +
      facet_wrap(~project) +
      scale_colour_manual(name = NULL, values = setNames(c("steelblue", "darkorange"),
                                                           c(pvalite_scenario_baseline, pvalite_scenario_envelope))) +
      scale_fill_manual(name = NULL, values = setNames(c("steelblue", "darkorange"),
                                                         c(pvalite_scenario_baseline, pvalite_scenario_envelope))) +
      labs(
        x = "Year", y = "Population size (median, 10-90th percentile band)",
        title = "PVA-lite's 25-year trend, fixed Safi point vitals vs. the joint biological envelope",
        subtitle = paste0(
          "N0 = Nmin = ", format(nmin_assumed, big.mark = ","), " (dotted grey); removal = each project's own imposed PBR\n",
          "threshold; no density dependence, matching PVA-lite's own baseline. Dotted red: quasi-extinction (",
          round(100 * pva_quasi_extinction_fraction), "% of N0)."
        )
      ) +
      theme_minimal() +
      theme(plot.subtitle = element_text(size = 9), legend.position = "bottom")
  })

  ## ---- 9. Updating Section 6's risk-vs-H comparison with the literature-
  ## grounded survival CEILING, not a breakeven-forced value. Paulo,
  ## 2026-10: "como e que ficaria a Curva P(quasi-extincao) vs H: baseline
  ## (r realista) vs. cenario plausivel (S_adult e lambda)" -- now that the
  ## references update (Section 5 of leslie_matrix_parametrisation.md)
  ## shows the Lentini et al. (2015) order-wide meta-analysis' own 95% CI
  ## upper bound for adult female survival (0.890, 44 species/7 families)
  ## sits BELOW the S_adult=0.949 Section 6 used -- a value forced to hit
  ## lambda=1.20 exactly, which no literature source actually reports, not
  ## a real estimate. The honest question is not "what survival would it
  ## take to reach Frick's benchmark" but "what does the best cross-
  ## species-supported survival estimate alone actually buy", holding the
  ## other three vital rates at Safi's baseline as Section 6 did.
  literature_ceiling_s_adult <- 0.890  # Lentini et al. 2015 order-wide 95% CI upper bound
  literature_ceiling_vitals <- baseline_vitals
  literature_ceiling_vitals$s_adult <- literature_ceiling_s_adult
  literature_ceiling_lambda <- do.call(lambda_of_vitals, literature_ceiling_vitals)
  literature_ceiling_H <- find_sustainable_H_general(
    N_target_elastic, K_elastic, s_juv = literature_ceiling_vitals$s_juv, s_adult = literature_ceiling_vitals$s_adult,
    p_breed = literature_ceiling_vitals$p_breed, litter = literature_ceiling_vitals$litter,
    acceptable_risk = 0.10, juv_mortality_ratio = juv_mortality_ratio
  )$H_sustainable

  sweep_H_max_litceiling <- literature_ceiling_H * 1.3
  sweep_dt_litceiling <- tibble::tibble(H = seq(0, sweep_H_max_litceiling, length.out = 18)) %>%
    rowwise() %>%
    mutate(p_collapse = run_risk_general(sweep_N_target, H, sweep_K, s_juv = literature_ceiling_vitals$s_juv,
                                          s_adult = literature_ceiling_vitals$s_adult, p_breed = literature_ceiling_vitals$p_breed,
                                          litter = literature_ceiling_vitals$litter, juv_mortality_ratio = juv_mortality_ratio,
                                          n_reps = 400)) %>%
    ungroup()

  scenario_litceiling <- paste0("S_adult=", literature_ceiling_s_adult,
                                 " (Lentini et al. 2015 order-wide 95% CI upper bound, lambda=",
                                 round(literature_ceiling_lambda, 3), ")")
  sweep_dt_litceiling_combined <- bind_rows(
    sweep_dt %>% mutate(scenario = scenario_realistic),
    sweep_dt_litceiling %>% mutate(scenario = scenario_litceiling)
  )

  fig_risk_sweep_litceiling <- file.path(fig_dir, "carrying_capacity_risk_sweep_litceiling.png")
  ggsave(fig_risk_sweep_litceiling, width = 9.5, height = 5.5, dpi = 150, bg = "white", plot = {
    ggplot(sweep_dt_litceiling_combined, aes(x = H, y = p_collapse, colour = scenario)) +
      geom_line(linewidth = 0.9) +
      geom_point(size = 1.8) +
      geom_vline(xintercept = literature_ceiling_H, linetype = "dashed", colour = "darkorange") +
      geom_vline(data = project_thresholds, aes(xintercept = threshold), linetype = "dotted", colour = "grey40") +
      geom_hline(data = tibble::tibble(acceptable_risk = acceptable_risk_grid * 100),
                 aes(yintercept = acceptable_risk), linetype = "dotted", colour = "grey40") +
      scale_colour_manual(name = NULL, values = setNames(c("steelblue", "darkorange"), c(scenario_realistic, scenario_litceiling))) +
      annotate("text", x = literature_ceiling_H, y = 102, label = paste0("H sustainable\n(10% risk) = ", round(literature_ceiling_H)),
               colour = "darkorange", size = 2.8, hjust = -0.05) +
      labs(
        x = "Constant annual removal H", y = "P(falls below quasi-extinction in 25 yr), %",
        title = "Risk-vs-H curve: realistic-r baseline vs. the literature-grounded survival ceiling (not a breakeven-forced value)",
        subtitle = paste0(
          "S_adult capped at the Lentini et al. (2015) order-wide meta-analysis' own 95% CI upper bound (0.890), not the\n",
          "S_adult=0.949 forced to hit lambda=1.20 exactly -- that value sat ABOVE every literature estimate found, this one doesn't.\n",
          "Resulting lambda=", round(literature_ceiling_lambda, 3), ", still below both PBR Rmax benchmarks (1.20/1.24)."
        )
      ) +
      theme_minimal() +
      theme(plot.subtitle = element_text(size = 8.3), legend.position = "bottom")
  })

  list(
    cc_table = cc_table,
    stoch_summary = stoch_summary,
    sustainable_H_table = sustainable_H_table,
    litter_sweep = litter_sweep,
    breakeven_H_table = breakeven_H_table,
    elasticity_sweep = elasticity_sweep,
    tornado_table = tornado_table,
    sweep_dt_combined = sweep_dt_combined,
    envelope_draws = envelope_draws,
    H_envelope_percentiles = H_envelope_percentiles,
    risk_at_pbr_summary = risk_at_pbr_summary,
    pvalite_trend_summary = pvalite_trend_summary,
    pvalite_trend_risk = pvalite_trend_risk,
    fig_H_curve = fig_H_curve,
    fig_stoch_check = fig_stoch_check,
    fig_risk_sweep = fig_risk_sweep,
    fig_risk_sweep_plausible = fig_risk_sweep_plausible,
    fig_litter_sweep = fig_litter_sweep,
    fig_elasticity_vitals = fig_elasticity_vitals,
    fig_tornado = fig_tornado,
    fig_envelope_H = fig_envelope_H,
    fig_envelope_risk_at_pbr = fig_envelope_risk_at_pbr,
    fig_pvalite_envelope = fig_pvalite_envelope,
    literature_ceiling_lambda = literature_ceiling_lambda,
    literature_ceiling_H = literature_ceiling_H,
    fig_risk_sweep_litceiling = fig_risk_sweep_litceiling
  )
}
