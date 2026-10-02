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

  beta_params <- function(mean, cv) {
    var <- (mean * cv)^2
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

  list(
    cc_table = cc_table,
    stoch_summary = stoch_summary,
    sustainable_H_table = sustainable_H_table,
    litter_sweep = litter_sweep,
    fig_H_curve = fig_H_curve,
    fig_stoch_check = fig_stoch_check,
    fig_risk_sweep = fig_risk_sweep,
    fig_litter_sweep = fig_litter_sweep
  )
}
