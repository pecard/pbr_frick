##
## Candidate PBR / management-reference analysis (Paulo, 2026-09), step 1
## of the sequence he set out after reviewing R/global_pbr_uncertainty.R
## (step 2): compares the current imposed thresholds against several
## quantiles of the global (demographic + project-specific Nmin)
## uncertainty distribution -- NOT a proposal to change the compliance
## reference, but a decision-support table making explicit what each
## candidate level would imply biologically and operationally, exactly as
## Paulo specified: "what changes biologically and operationally if the
## management reference moves from 112/88 towards the centre of the
## plausible PBR range?"
##
## For each project and each candidate level (current threshold, Q05,
## Q25, median, Q75 of that project's global PBR distribution), reports:
##   - where the level sits in the distribution (percentile);
##   - the Nmin that level implies at that project's own median simulated
##     lambda_max (a biological "what would have to be true" read, not a
##     new estimate);
##   - PVA-lite decline/quasi-extinction risk under that level applied as
##     a constant annual removal (same validated structure, same N0 =
##     nmin_assumed, as Section 4 -- so directly comparable to the numbers
##     already reported there);
##   - the mortality-budgeting consequence: given the no-curtailment
##     counterfactual annual total (R/curtailment_year_effectiveness.R),
##     the reduction required to sit at that level, and how many turbines
##     (ranked by pre-curtailment historical contribution, same ranking as
##     Section 7.2) a uniform 70% curtailment-effectiveness assumption
##     would need to cover to deliver it.
##
## All five rows per project use the SAME methodology (including the
## "current threshold" row), so they are directly comparable to each
## other -- this is deliberately NOT the same turbine count as the real,
## already-executed selection (17/7 turbines, Section 7.2-7.3), which
## used a different method (50% cumulative mortality in a single baseline
## season) and a different, empirically-observed effectiveness (Section
## 7.4). That real-world figure is the actual decision already taken;
## this table is a like-for-like comparison ACROSS candidate levels, not
## a replacement for it -- cross-referenced explicitly in the report text.
##

## Installs ggrepel automatically if missing (used below for non-
## overlapping labels on the candidate-level figure) -- same pattern as
## the optional 3D-figure packages in R/pbr_analysis.R.
if (!requireNamespace("ggrepel", quietly = TRUE)) {
  message("Installing missing package needed for the candidate-PBR figure labels: ggrepel...")
  tryCatch(
    install.packages("ggrepel"),
    error = function(e) stop(
      "Could not install 'ggrepel' automatically (", conditionMessage(e),
      "). Install it manually (install.packages(\"ggrepel\")) and re-run."
    )
  )
}

suppressPackageStartupMessages({ library(dplyr); library(tidyr); library(ggplot2); library(ggrepel) })

## Per-turbine curtailment effectiveness default: 0.63, the published
## average reduction in bat mortality from operational curtailment at
## North American wind farms --
##   Adams E.M., Gulka J., Williams K.A. (2021). A review of the
##   effectiveness of operational curtailment for reducing bat fatalities
##   at terrestrial wind farms in North America. PLoS ONE 16(9): e0256382.
##   https://doi.org/10.1371/journal.pone.0256382
## -- the same benchmark already cited in the Blanket Curtailment Plan
## technical note to justify the 6.5 m/s cut-in speed threshold used
## operationally at both facilities. NOT the same thing as the real,
## empirically-observed 88%/53% facility-level reduction reported in
## Section 7.4, which is an aggregate outcome of curtailing a specific
## 17/7-turbine subset for part of a season, not a per-turbine rate;
## using it here would conflate that realised, selection-specific outcome
## with a general per-turbine assumption and overstate how effective
## curtailing ADDITIONAL, untested turbines would be.
run_candidate_pbr_reference <- function(fig_dir, nmin_assumed, global_pbr_quantiles_by_project,
                                         real_turbine_hist, curtailment_year_full_counterfactual,
                                         curtailment_effectiveness = 0.63) {

  dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)
  fig_dir <- normalizePath(fig_dir, winslash = "/")

  fr_corrected <- fr_from_iucn(iucn_status_corrected)
  project_thresholds <- tibble::tibble(project = facility_labels, threshold = pbr_thresholds$threshold)

  ## ---- Self-contained PVA-lite runner (same validated structure/settings
  ## as Section 4; a constant annual removal, not a weekly control chart) --
  ## build_dekker_stage_matrix() is duplicated locally here rather than
  ## sourced from R/leslie_dekker_limpens.R, matching the same pattern
  ## R/pbr_analysis.R itself already uses (see its own PVA-lite section).
  build_dekker_stage_matrix <- function(n_stages, s_juv, s_adult, fecundity) {
    A <- matrix(0, n_stages, n_stages)
    A[1, n_stages] <- fecundity
    for (i in 1:(n_stages - 1)) A[i + 1, i] <- s_juv
    A[n_stages, n_stages] <- s_adult
    A
  }
  leslie_F <- leslie_p_breed * leslie_litter * leslie_sex_ratio * leslie_s_juv
  A_leslie_validated <- build_dekker_stage_matrix(2, leslie_s_juv, leslie_s_adult_f, leslie_F)
  stable_dist <- Re(eigen(A_leslie_validated)$vectors[, 1]); stable_dist <- stable_dist / sum(stable_dist)

  beta_params <- function(mean, cv) {
    var <- (mean * cv)^2
    common <- mean * (1 - mean) / var - 1
    list(shape1 = mean * common, shape2 = (1 - mean) * common)
  }
  bp_s_adult <- beta_params(leslie_s_adult_f, pva_vital_rate_cv)
  bp_s_juv <- beta_params(leslie_s_juv, pva_vital_rate_cv)
  bp_p_breed <- beta_params(leslie_p_breed, pva_vital_rate_cv)

  simulate_pva_trajectory <- function(n0_juv, n0_adult, n_years, annual_removal) {
    n_juv <- as.integer(round(n0_juv)); n_adult <- as.integer(round(n0_adult))
    traj <- numeric(n_years + 1); traj[1] <- 2 * (n_juv + n_adult)
    for (yr in seq_len(n_years)) {
      s_adult_t <- rbeta(1, bp_s_adult$shape1, bp_s_adult$shape2)
      s_juv_t <- rbeta(1, bp_s_juv$shape1, bp_s_juv$shape2)
      p_breed_t <- rbeta(1, bp_p_breed$shape1, bp_p_breed$shape2)
      n_breeders <- rbinom(1, n_adult, p_breed_t)
      n_pups_total <- rpois(1, n_breeders * leslie_litter)
      n_female_pups <- rbinom(1, n_pups_total, leslie_sex_ratio)
      n_juv_recruits <- rbinom(1, n_female_pups, s_juv_t)
      n_juv_survive <- rbinom(1, n_juv, s_juv_t)
      n_adult_survive <- rbinom(1, n_adult, s_adult_t)
      n_juv_next <- n_juv_recruits
      n_adult_next <- n_juv_survive + n_adult_survive
      female_removal <- annual_removal / 2
      total_female <- n_juv_next + n_adult_next
      if (total_female > 0 && female_removal > 0) {
        juv_share <- round(female_removal * n_juv_next / total_female)
        adult_share <- round(female_removal * n_adult_next / total_female)
        n_juv_next <- max(0, n_juv_next - juv_share)
        n_adult_next <- max(0, n_adult_next - adult_share)
      }
      n_juv <- n_juv_next; n_adult <- n_adult_next
      traj[yr + 1] <- 2 * (n_juv + n_adult)
    }
    traj
  }

  n0_female <- nmin_assumed / 2
  n0_juv <- round(n0_female * stable_dist[1]); n0_adult <- n0_female - n0_juv
  quasi_ext_threshold <- pva_quasi_extinction_fraction * nmin_assumed

  pva_risk_for_level <- function(annual_removal) {
    set.seed(pva_seed)
    reps <- replicate(pva_n_reps, simulate_pva_trajectory(n0_juv, n0_adult, pva_n_years, annual_removal))
    final_N <- reps[pva_n_years + 1, ]
    ever_below <- apply(reps, 2, function(traj) any(traj < quasi_ext_threshold))
    tibble::tibble(
      pva_p_decline = mean(final_N < nmin_assumed) * 100,
      pva_p_quasi_extinction = mean(ever_below) * 100
    )
  }

  ## ---- Candidate levels per project: current threshold + Q05/Q25/median/Q75
  candidate_rows <- lapply(seq_len(nrow(project_thresholds)), function(i) {
    proj <- project_thresholds$project[i]; thr <- project_thresholds$threshold[i]
    q <- global_pbr_quantiles_by_project %>% filter(project == proj)
    median_lambda_max <- {
      # recover a representative lambda_max for this project's simulated
      # distribution: same relationship as pbr_from_components(), inverted
      # at the distribution's own median Nmin and median PBR.
      pbr_val <- q$median; nmin_val <- q$nmin_median
      1 + (2 * pbr_val) / (fr_corrected * nmin_val)
    }
    turb <- real_turbine_hist %>% filter(project == proj) %>% arrange(desc(corrected)) %>%
      rename(expected_remaining = corrected)
    counterfactual_total <- curtailment_year_full_counterfactual %>% filter(project == proj) %>% pull(counterfactual_annual_total)

    levels_tbl <- tibble::tibble(
      candidate = c("Current imposed threshold", "Q05 (conservative tail)", "Q25", "Median", "Q75 (relaxed tail)"),
      level = c(thr, q$q05, q$q25, q$median, q$q75),
      percentile_in_distribution = c(q$threshold_percentile, 5, 25, 50, 75)
    )

    levels_tbl %>%
      rowwise() %>%
      mutate(
        project = proj,
        implied_nmin_at_median_lambda = level / (0.5 * fr_corrected * (median_lambda_max - 1)),
        required_abs_reduction = max(0, counterfactual_total - level),
        pct_reduction_needed = if (counterfactual_total > 0) 100 * required_abs_reduction / counterfactual_total else 0,
        n_turbines_required = {
          alloc <- allocate_curtailment(turb, required_abs_reduction, effectiveness = curtailment_effectiveness)
          sum(alloc$curtail)
        },
        n_turbines_available = nrow(turb)
      ) %>%
      ungroup() %>%
      bind_cols(bind_rows(lapply(levels_tbl$level, pva_risk_for_level)))
  }) %>% bind_rows()

  candidate_pbr_table <- candidate_rows %>%
    select(project, candidate, level, percentile_in_distribution, implied_nmin_at_median_lambda,
           pva_p_decline, pva_p_quasi_extinction, required_abs_reduction, pct_reduction_needed,
           n_turbines_required, n_turbines_available)

  ## ---- Figure: level vs turbines-required vs PVA decline risk, per project -
  fig_candidate_pbr <- file.path(fig_dir, "candidate_pbr_reference.png")
  ggsave(fig_candidate_pbr, width = 11, height = 5.5, dpi = 150, plot = {
    ggplot(candidate_pbr_table, aes(x = n_turbines_required, y = level)) +
      geom_line(colour = "grey50", linewidth = 0.5) +
      geom_point(aes(colour = pva_p_decline, size = pct_reduction_needed)) +
      ggrepel::geom_text_repel(aes(label = candidate), size = 2.6, seed = 1, max.overlaps = Inf) +
      facet_wrap(~project, scales = "free") +
      scale_colour_viridis_c(option = "C", name = "PVA-lite:\nP(decline), %") +
      scale_size_continuous(name = "% reduction\nneeded", range = c(2, 7)) +
      labs(
        x = "Turbines required (ranked by pre-curtailment contribution, 70% assumed effectiveness)",
        y = "Candidate management reference (bats/year)",
        title = "2027 mitigation frontier: candidate reference vs. turbines required vs. biological risk",
        subtitle = "Each point is one candidate level (current threshold, and Q05/Q25/median/Q75 of the global Nmin+demographic PBR distribution)"
      ) +
      theme_minimal() +
      theme(plot.subtitle = element_text(size = 11), strip.text = element_text(face = "bold"))
  })

  list(
    candidate_pbr_table = candidate_pbr_table,
    candidate_pbr_effectiveness = curtailment_effectiveness,
    fig_candidate_pbr = fig_candidate_pbr
  )
}
