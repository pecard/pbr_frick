##
## Global PBR uncertainty (Paulo, 2026-09): the Monte Carlo in
## R/pbr_analysis.R (section 4) propagates uncertainty in s and alpha
## only, at a single, project-shared Nmin (nmin_assumed, the
## reverse-engineered value from "Reverse-engineering the imposed
## thresholds"). Section 5.4 already shows numerically that PBR is far
## more sensitive to Nmin than to lambda_max (doubling Nmin doubles PBR;
## the paper's two lambda_max benchmarks differ by only 20%) -- so a
## Monte Carlo that holds Nmin fixed is not actually propagating the
## uncertainty that matters most.
##
## This script extends it: Nmin is now itself drawn from a project-
## specific distribution -- that project's own footprint x a plausible
## density range -- rather than held fixed, and each project gets its own
## PBR distribution rather than sharing one. Projects are treated in
## ISOLATION (Paulo, 2026-09, explicit instruction): Nmin is drawn from
## project_footprint_km2 (each project's own ~370 km2 footprint) x
## Uniform(density_proxy_low, density_proxy_high) bats/km2 -- the same
## density proxy (Pipistrellus pipistrellus, Frick et al.) already used
## for the Nmin plausibility check in Section 2.4, now used generatively
## rather than only as a comparison point. The combined_region_km2
## (~7,000 km2) regional/migratory-population scenario discussed
## elsewhere in this note is deliberately NOT folded in here -- pooling
## the two projects into one shared population is a separate, still-
## unresolved question, not a free parameter of this particular analysis.
##
## Fr stays fixed at the corrected value (0.5, Least Concern) throughout --
## this is about the Nmin x demographic-capacity axes, not Fr, which is
## already treated as a discrete scenario elsewhere (Section 2.5-2.6, 5.2).
##

suppressPackageStartupMessages({ library(dplyr); library(tidyr); library(ggplot2) })

run_global_pbr_uncertainty <- function(fig_dir, nmin_assumed, old_pbr_quantiles) {

  dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)
  fig_dir <- normalizePath(fig_dir, winslash = "/")

  fr_corrected <- fr_from_iucn(iucn_status_corrected)
  project_thresholds <- tibble::tibble(project = facility_labels, threshold = pbr_thresholds$threshold)

  ## ---- Reproduce the OLD (demographic-only, Nmin-fixed) simulation ---------
  ## Same seed/ranges as R/pbr_analysis.R section 4 -- an exact
  ## reproduction of that simulation, not a new estimate, so the "before"
  ## comparison below is on identical footing to what the rest of this
  ## note already reports.
  set.seed(mc_seed)
  sim_old <- tibble::tibble(
    s     = runif(mc_n_sim, s_range[1], s_range[2]),
    alpha = runif(mc_n_sim, alpha_range[1], alpha_range[2])
  ) %>%
    mutate(lambda_max = lambda_max_niel(s, alpha), PBR = pbr_from_components(nmin_assumed, fr_corrected, lambda_max))

  ## ---- NEW: demographic + project-specific Nmin uncertainty ----------------
  set.seed(mc_seed)
  sim_global <- lapply(project_thresholds$project, function(proj) {
    tibble::tibble(
      project = proj,
      s       = runif(mc_n_sim, s_range[1], s_range[2]),
      alpha   = runif(mc_n_sim, alpha_range[1], alpha_range[2]),
      density = runif(mc_n_sim, density_proxy_low, density_proxy_high)
    ) %>%
      mutate(
        Nmin = density * project_footprint_km2,
        lambda_max = lambda_max_niel(s, alpha),
        PBR = pbr_from_components(Nmin, fr_corrected, lambda_max)
      )
  }) %>% bind_rows()

  quantiles_by_project <- sim_global %>%
    group_by(project) %>%
    summarise(
      nmin_q05 = quantile(Nmin, 0.05), nmin_median = quantile(Nmin, 0.50), nmin_q95 = quantile(Nmin, 0.95),
      q05 = quantile(PBR, 0.05), q25 = quantile(PBR, 0.25), median = quantile(PBR, 0.50),
      q75 = quantile(PBR, 0.75), q95 = quantile(PBR, 0.95),
      .groups = "drop"
    ) %>%
    left_join(project_thresholds, by = "project") %>%
    rowwise() %>%
    mutate(
      threshold_percentile = mean(sim_global$PBR[sim_global$project == project] < threshold) * 100,
      median_ratio = median / threshold
    ) %>%
    ungroup()

  old_new_comparison <- bind_rows(
    tibble::tibble(
      scenario = "Demographic uncertainty only (s, alpha; Nmin fixed)",
      q05 = old_pbr_quantiles[["5%"]], median = old_pbr_quantiles[["50%"]], q95 = old_pbr_quantiles[["95%"]]
    ),
    quantiles_by_project %>% transmute(scenario = paste0("+ project-specific Nmin (", project, ")"), q05, median, q95)
  )

  ## ---- Figure: old (shared) vs. new (per-project) PBR distributions --------
  plot_dt <- bind_rows(
    sim_old %>% transmute(project = "Demographic uncertainty only\n(shared, Nmin fixed)", PBR),
    sim_global %>% transmute(project, PBR)
  ) %>%
    mutate(project = factor(project, levels = c("Demographic uncertainty only\n(shared, Nmin fixed)", project_thresholds$project)))

  threshold_lines <- bind_rows(
    tibble::tibble(project = "Demographic uncertainty only\n(shared, Nmin fixed)", threshold = pbr_thresholds$threshold, facility = pbr_thresholds$facility),
    project_thresholds %>% rename(facility = project) %>% mutate(project = facility)
  ) %>% mutate(project = factor(project, levels = levels(plot_dt$project)))

  fig_global_pbr_uncertainty <- file.path(fig_dir, "global_pbr_uncertainty.png")
  ggsave(fig_global_pbr_uncertainty, width = 12, height = 4.5, dpi = 150, plot = {
    ggplot(plot_dt, aes(x = PBR)) +
      geom_density(fill = "steelblue", alpha = 0.4, colour = "steelblue4") +
      geom_vline(data = threshold_lines, aes(xintercept = threshold, colour = facility), linewidth = 0.8, linetype = "dashed") +
      facet_wrap(~project, nrow = 1, scales = "free") +
      scale_colour_manual(name = "Imposed threshold", values = setNames(c("grey20", "firebrick"), pbr_thresholds$facility)) +
      labs(
        title = "Effect of adding project-specific Nmin uncertainty to the PBR distribution",
        subtitle = paste0(
          "Left: original Monte Carlo (Nmin fixed at ", format(nmin_assumed, big.mark = ","), ", shared across both projects).\n",
          "Right: Nmin drawn from each project's own footprint (", project_footprint_km2, " km2) x Uniform(",
          density_proxy_low, ", ", density_proxy_high, ") bats/km2 -- projects treated in isolation."
        ),
        x = "PBR (bats/year)", y = "Density"
      ) +
      theme_minimal() +
      theme(plot.subtitle = element_text(size = 11), strip.text = element_text(face = "bold"))
  })

  list(
    global_nmin_area_km2 = project_footprint_km2,
    global_density_low = density_proxy_low,
    global_density_high = density_proxy_high,
    global_pbr_quantiles_by_project = quantiles_by_project,
    global_pbr_old_new_comparison = old_new_comparison,
    fig_global_pbr_uncertainty = fig_global_pbr_uncertainty
  )
}
