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

suppressPackageStartupMessages({ library(dplyr); library(tidyr); library(ggplot2) })

run_carrying_capacity_reference <- function(fig_dir, nmin_assumed,
                                             alpha_grid = c(0.5, 0.6, 0.7, 0.8, 0.9)) {

  dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)
  fig_dir <- normalizePath(fig_dir, winslash = "/")

  project_thresholds <- tibble::tibble(project = facility_labels, threshold = pbr_thresholds$threshold)

  ## ---- 1. Deterministic H(alpha, K) across the same K and r uncertainty
  ## already used elsewhere in this note -------------------------------------
  K_scenarios <- tibble::tibble(
    K_label = c("K low (density_proxy_low x footprint)", "K high (density_proxy_high x footprint)"),
    K = c(density_proxy_low, density_proxy_high) * project_footprint_km2
  )
  r_scenarios <- tibble::tibble(
    r_label = paste0("Rmax at lambda_max = ", lambda_max_benchmarks),
    Rmax = lambda_max_benchmarks - 1
  )

  cc_table <- tidyr::expand_grid(project_thresholds, K_scenarios, r_scenarios, alpha = alpha_grid) %>%
    mutate(
      N_target = alpha * K,
      H = Rmax * K * alpha * (1 - alpha),
      pct_of_current_threshold = 100 * H / threshold,
      is_MSY = alpha == 0.5
    )

  fig_H_curve <- file.path(fig_dir, "carrying_capacity_H_curve.png")
  ggsave(fig_H_curve, width = 9.5, height = 5, dpi = 150, plot = {
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
  simulate_constant_harvest <- function(n0_juv, n0_adult, n_years, H, K, cv = pva_vital_rate_cv) {
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
      ## H is a continuous removal rate (not necessarily an integer
      ## number of animals); round after subtracting it, otherwise
      ## n_adult drifts non-integer and the next iteration's rbinom()
      ## silently returns NaN (rbinom requires an integer `size`).
      n_adult_next <- as.integer(round(max(0, n_juv_survive + n_adult_survive - H)))
      n_juv <- n_juv_recruits; n_adult <- n_adult_next
      traj[yr + 1] <- 2 * (n_juv + n_adult)
    }
    traj
  }

  run_from_start <- function(N0, H, K, n_years = pva_n_years, n_reps = pva_n_reps, seed = pva_seed) {
    set.seed(seed)
    nf <- N0 / 2; nj <- round(nf * stable_dist[1]); na <- nf - nj
    reps <- replicate(n_reps, simulate_constant_harvest(nj, na, n_years, H, K))
    list(
      trajectories = reps,
      median_traj = apply(reps, 1, median),
      p10_traj = apply(reps, 1, quantile, probs = 0.10),
      p90_traj = apply(reps, 1, quantile, probs = 0.90),
      p_collapse = mean(reps[n_years + 1, ] < 0.1 * K) * 100
    )
  }

  ## Representative scenario for the stochastic check: alpha = 0.8
  ## (a plausible "near K" ecological-viability target), K high.
  check_alpha <- 0.8
  stoch_check <- project_thresholds %>%
    rowwise() %>%
    mutate(
      K = density_proxy_high * project_footprint_km2,
      Rmax = lambda_max_benchmarks[1] - 1,
      H = Rmax * K * check_alpha * (1 - check_alpha),
      N_target = check_alpha * K
    ) %>%
    ungroup()

  stoch_results <- list()
  traj_plot_rows <- list()
  for (i in seq_len(nrow(stoch_check))) {
    proj <- stoch_check$project[i]; K <- stoch_check$K[i]; H <- stoch_check$H[i]; N_target <- stoch_check$N_target[i]
    at_target <- run_from_start(N_target, H, K)
    at_nmin <- run_from_start(nmin_assumed, H, K)
    stoch_results[[proj]] <- list(
      p_collapse_starting_at_target = at_target$p_collapse,
      p_collapse_starting_at_nmin = at_nmin$p_collapse
    )
    traj_plot_rows[[paste0(proj, "_target")]] <- tibble::tibble(
      project = proj, year = 0:pva_n_years, start = "Starting at N* = alpha*K",
      median = at_target$median_traj, p10 = at_target$p10_traj, p90 = at_target$p90_traj
    )
    traj_plot_rows[[paste0(proj, "_nmin")]] <- tibble::tibble(
      project = proj, year = 0:pva_n_years, start = "Starting at Nmin (today's plausible status)",
      median = at_nmin$median_traj, p10 = at_nmin$p10_traj, p90 = at_nmin$p90_traj
    )
  }
  traj_plot_dt <- bind_rows(traj_plot_rows)
  stoch_summary <- stoch_check %>%
    mutate(
      p_collapse_starting_at_target = sapply(project, function(p) stoch_results[[p]]$p_collapse_starting_at_target),
      p_collapse_starting_at_nmin = sapply(project, function(p) stoch_results[[p]]$p_collapse_starting_at_nmin)
    )

  fig_stoch_check <- file.path(fig_dir, "carrying_capacity_stochastic_check.png")
  ggsave(fig_stoch_check, width = 10, height = 5.5, dpi = 150, plot = {
    ggplot(traj_plot_dt, aes(x = year, y = median, colour = start, fill = start)) +
      geom_ribbon(aes(ymin = p10, ymax = p90), alpha = 0.15, colour = NA) +
      geom_line(linewidth = 1) +
      geom_hline(data = stoch_check, aes(yintercept = N_target), linetype = "dotted", colour = "forestgreen") +
      geom_hline(data = stoch_check, aes(yintercept = 0.1 * K), linetype = "dotted", colour = "firebrick") +
      facet_wrap(~project, scales = "free_y") +
      labs(
        x = "Year", y = "Population size (median, 10-90th percentile band)",
        title = paste0("Same constant removal H (alpha=", check_alpha, " target), two starting points"),
        subtitle = "Green dotted: N* = alpha*K target. Red dotted: quasi-extinction (10% of K). Same H can be safe from one start and not the other."
      ) +
      theme_minimal() +
      theme(plot.subtitle = element_text(size = 9.5), strip.text = element_text(face = "bold"), legend.position = "bottom")
  })

  list(
    cc_table = cc_table,
    stoch_summary = stoch_summary,
    fig_H_curve = fig_H_curve,
    fig_stoch_check = fig_stoch_check
  )
}
