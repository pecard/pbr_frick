##
## Runs the full PBR reverse-engineering + demographic-uncertainty analysis
## for one scenario and returns every object the report template
## (report/pbr_technical_note_template.Rmd) needs as params. Mirrors the
## IDF report convention: this function computes, the Rmd only formats.
##
## Expects the scenario's settings (inputs/pbrSettings_*.R) to already be
## sourced into the calling environment -- see run_pbr_report_BSH_DGY.R.
## Figures are ggsave()'d to fig_dir and passed onward as file paths (same
## convention as R/report.R / IDF_monthly_report.R), not as live ggplot
## objects, so the report template only ever does
## knitr::include_graphics().
##

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(ggplot2)
})

source("R/pbr_functions.R")

run_pbr_analysis <- function(fig_dir) {

  dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)
  # Absolute path: rmarkdown::render() knits with the working directory set
  # to the .Rmd's own folder, so a path relative to the project root (as
  # passed in by the launcher) would not resolve inside the template.
  # winslash = "/" matters on Windows: normalizePath()'s default backslash
  # paths break the markdown image reference pandoc generates from
  # knitr::include_graphics(), silently falling back to the image's alt
  # -text description instead of embedding it.
  fig_dir <- normalizePath(fig_dir, winslash = "/")

  ## ---- 1. Reverse-engineer the imposed thresholds ------------------------
  fr_original <- fr_from_iucn(iucn_status_original_guess)
  fr_corrected <- fr_from_iucn(iucn_status_corrected)

  backsolve_grid <- tidyr::expand_grid(
    pbr_thresholds,
    lambda_max_fixed = lambda_max_benchmarks
  ) %>%
    mutate(Nmin_required = required_nmin(threshold, fr_original, lambda_max_fixed))

  # Most parsimonious Nmin: the value shared by at least two facility x
  # lambda_max pairings (i.e. some assignment of benchmarks to facilities
  # back-solves to the *same* local population). Falls back to the median
  # of all candidates, with a warning, if nothing repeats.
  find_shared_value <- function(candidates, tol = 1) {
    for (i in seq_along(candidates)) {
      dup <- which(abs(candidates - candidates[i]) < tol & seq_along(candidates) != i)
      if (length(dup) > 0) return(candidates[i])
    }
    warning("No shared Nmin found across threshold x lambda_max pairings; falling back to the median.")
    median(candidates)
  }
  nmin_assumed <- find_shared_value(backsolve_grid$Nmin_required)

  # Plausibility check against the REAL, known project geography (not a
  # back-solved area): what density would Nmin_assumed imply if it applies
  # to a single project's own footprint, versus if it applies to the whole
  # landscape spanning both projects, holding the density_high:density_low
  # ratio fixed at the paper's own default proxy?
  density_ratio <- density_proxy_high / density_proxy_low

  density_at_footprint <- implied_density(nmin_assumed, project_footprint_km2, density_ratio)
  density_at_region    <- implied_density(nmin_assumed, combined_region_km2, density_ratio)

  nmin_plausibility <- tibble::tibble(
    area_label = c(
      sprintf("Single project footprint (%s km2)", project_footprint_km2),
      sprintf("Combined region, both projects (%s km2)", combined_region_km2)
    ),
    area_km2 = c(project_footprint_km2, combined_region_km2),
    density_low = c(density_at_footprint["low"], density_at_region["low"]),
    density_high = c(density_at_footprint["high"], density_at_region["high"])
  )

  pbr_corrected_fr <- tibble::tibble(lambda_max_fixed = lambda_max_benchmarks) %>%
    mutate(PBR = pbr_from_components(nmin_assumed, fr_corrected, lambda_max_fixed))

  ## ---- 2. Sensitivity / elasticity ---------------------------------------
  s_seq <- seq(s_range[1], s_range[2], by = 0.01)
  alpha_seq <- seq(alpha_range[1], alpha_range[2], by = 0.1)

  partials_lambda <- function(s, alpha, eps_s = 1e-5, eps_a = 1e-5) {
    lam <- lambda_max_niel(s, alpha)
    dlds <- (lambda_max_niel(s + eps_s, alpha) - lambda_max_niel(s - eps_s, alpha)) / (2 * eps_s)
    dlda <- (lambda_max_niel(s, alpha + eps_a) - lambda_max_niel(s, alpha - eps_a)) / (2 * eps_a)
    tibble::tibble(
      lambda_max = lam, d_lambda_ds = dlds, d_lambda_dalpha = dlda,
      E_s = dlds * s / lam, E_alpha = dlda * alpha / lam
    )
  }

  sens_grid <- tidyr::expand_grid(s = s_seq, alpha = alpha_seq) %>%
    rowwise() %>%
    mutate(out = list(partials_lambda(s, alpha))) %>%
    ungroup() %>%
    tidyr::unnest(out)

  elasticity_ranges <- sens_grid %>%
    summarise(
      E_s_min = min(E_s), E_s_max = max(E_s), E_s_median = median(E_s),
      E_alpha_min = min(E_alpha), E_alpha_max = max(E_alpha), E_alpha_median = median(E_alpha)
    ) %>%
    tidyr::pivot_longer(everything(), names_to = c("Parameter", "Statistic"),
                         names_pattern = "E_(s|alpha)_(min|max|median)", values_to = "Value") %>%
    mutate(Parameter = recode(Parameter, s = "Adult survival (s)", alpha = "Age at first breeding (alpha)")) %>%
    tidyr::pivot_wider(names_from = Statistic, values_from = Value) %>%
    select(Parameter, min, max, median)

  # lambda_max is > 1 (population growing) everywhere in sens_grid -- an
  # inherent property of the Niel & Lebreton formula for s in (0,1), not a
  # feature of the chosen ranges: for any s, lambda_max only APPROACHES 1
  # asymptotically as alpha -> infinity, and never actually reaches it for
  # any finite, realistic alpha (e.g. still ~1.0001 at alpha = 10,000
  # years). A "lambda_max = 1" reference line therefore cannot be drawn
  # within any realistic plot range; instead, both elasticity plots below
  # state the actual lambda_max range spanned, so a reader does not read
  # "elasticity is negative" as "the population is declining" -- it means
  # lambda_max *decreases* (while staying well above 1) as survival/age
  # increase. See fig_lambda_max for the lambda_max values themselves.
  lambda_range_txt <- sprintf("lambda_max = %.2f-%.2f here (always > 1, growing)",
                               min(sens_grid$lambda_max), max(sens_grid$lambda_max))

  fig_elasticity_survival <- file.path(fig_dir, "elasticity_survival.png")
  ggsave(fig_elasticity_survival, width = 7, height = 4.5, dpi = 150, plot = {
    ggplot(sens_grid, aes(x = s, y = alpha, fill = E_s)) +
      geom_raster(interpolate = TRUE) +
      scale_fill_viridis_c(name = "E[s]", option = "B", direction = -1) +
      labs(x = "Adult survival (s)", y = "Age at first breeding (alpha)",
           title = "Elasticity of lambda_max to adult survival",
           subtitle = lambda_range_txt) +
      theme_minimal() +
      theme(plot.subtitle = element_text(size = 9))
  })

  fig_elasticity_alpha <- file.path(fig_dir, "elasticity_alpha.png")
  ggsave(fig_elasticity_alpha, width = 7, height = 4.5, dpi = 150, plot = {
    ggplot(sens_grid, aes(x = s, y = alpha, fill = E_alpha)) +
      geom_raster(interpolate = TRUE) +
      scale_fill_viridis_c(name = "E[alpha]", option = "A", direction = -1) +
      labs(x = "Adult survival (s)", y = "Age at first breeding (alpha)",
           title = "Elasticity of lambda_max to age at first breeding",
           subtitle = lambda_range_txt) +
      theme_minimal() +
      theme(plot.subtitle = element_text(size = 9))
  })

  # One colour per benchmark contour (named by its own value) instead of a
  # single colour for both lines -- two same-coloured, unlabelled isolines
  # are not distinguishable from each other.
  benchmark_labels <- as.character(lambda_max_benchmarks)
  # Avoid white: invisible as a legend key swatch against the (white)
  # legend background, even though it reads fine on the dark raster itself.
  benchmark_colours <- setNames(c("cyan", "chartreuse")[seq_along(lambda_max_benchmarks)], benchmark_labels)

  fig_lambda_max <- file.path(fig_dir, "lambda_max.png")
  ggsave(fig_lambda_max, width = 8, height = 4.5, dpi = 150, plot = {
    ggplot(sens_grid, aes(x = s, y = alpha)) +
      geom_raster(aes(fill = lambda_max), interpolate = TRUE) +
      geom_contour(aes(z = lambda_max, colour = after_stat(factor(level))),
                   breaks = lambda_max_benchmarks, linewidth = 0.5) +
      scale_fill_viridis_c(name = "lambda_max", option = "C") +
      scale_colour_manual(name = "Benchmark", values = benchmark_colours) +
      labs(
        x = "Adult survival (s)", y = "Age at first breeding (alpha)",
        title = "lambda_max derived from the demographic-invariant approximation",
        subtitle = paste0(
          "Contours: paper's fixed benchmarks; lambda_max > 1 (growing) throughout"
        )
      ) +
      theme_minimal() +
      theme(plot.subtitle = element_text(size = 9))
  })

  ## ---- 3. Response surfaces across s, alpha and Fr -----------------------
  surface_grid <- tidyr::expand_grid(
    s = seq(s_range[1], s_range[2], length.out = 60),
    alpha = seq(alpha_range[1], alpha_range[2], length.out = 60),
    Fr = fr_scenarios
  ) %>%
    mutate(
      lambda_max = lambda_max_niel(s, alpha),
      PBR = pbr_from_components(nmin_assumed, Fr, lambda_max),
      Fr_label = factor(paste0("Fr = ", Fr), levels = paste0("Fr = ", fr_scenarios))
    )

  # One colour per threshold contour (named by facility) -- same-coloured,
  # unlabelled isolines for different thresholds are not distinguishable.
  threshold_labels_chr <- as.character(pbr_thresholds$threshold)
  threshold_colours <- setNames(c("cyan", "chartreuse")[seq_len(nrow(pbr_thresholds))], threshold_labels_chr)
  threshold_display_labels <- setNames(
    sprintf("%s (%s)", pbr_thresholds$facility, pbr_thresholds$threshold), threshold_labels_chr
  )

  fig_response_surfaces <- file.path(fig_dir, "response_surfaces.png")
  ggsave(fig_response_surfaces, width = 10, height = 3.6, dpi = 150, plot = {
    ggplot(surface_grid, aes(x = s, y = alpha)) +
      geom_raster(aes(fill = PBR), interpolate = TRUE) +
      geom_contour(aes(z = PBR, colour = after_stat(factor(level))),
                   breaks = pbr_thresholds$threshold, linewidth = 0.4) +
      facet_wrap(~Fr_label, nrow = 1) +
      scale_fill_viridis_c(name = "PBR\n(bats/yr)", option = "C") +
      scale_colour_manual(name = "Imposed threshold", values = threshold_colours, labels = threshold_display_labels) +
      labs(
        x = "Adult survival (s)", y = "Age at first breeding (alpha, yrs)",
        title = "PBR response surfaces across s, alpha and Fr",
        subtitle = paste0("Nmin = ", format(nmin_assumed, big.mark = ","), " bats/year")
      ) +
      theme_minimal() +
      theme(strip.text = element_text(face = "bold"))
  })

  ## ---- 3b. Response surfaces, 3D (static print, matching the original ----
  ##          report's plot.ly figure) -- optional: needs the 'plotly' and
  ##          'kaleido' packages (install.packages(c("plotly", "kaleido")))
  ##          to export a static image; skipped with a message otherwise,
  ##          the 2D faceted version above always covers the same content.
  fig_response_surfaces_3d <- NULL
  has_3d_deps <- requireNamespace("plotly", quietly = TRUE) &&
    requireNamespace("RColorBrewer", quietly = TRUE) &&
    requireNamespace("kaleido", quietly = TRUE)

  if (has_3d_deps) {
    hues <- c("Blues", "Greens", "Oranges", "Purples")  # one per Fr scenario, low->high
    fig3d <- plotly::plot_ly()
    for (i in seq_along(fr_scenarios)) {
      fr_val <- fr_scenarios[i]
      sub <- surface_grid %>% filter(Fr == fr_val)
      s_vals <- sort(unique(sub$s))
      a_vals <- sort(unique(sub$alpha))
      z_mat <- matrix(sub$PBR, nrow = length(s_vals), byrow = FALSE)

      fig3d <- fig3d %>%
        plotly::add_surface(
          x = a_vals, y = s_vals, z = z_mat,
          colorscale = list(c(0, 1), c("white", RColorBrewer::brewer.pal(9, hues[i])[7])),
          showscale = FALSE, opacity = 0.85,
          name = paste0("Fr = ", fr_val)
        ) %>%
        plotly::add_trace(
          type = "scatter3d", mode = "markers",
          x = a_vals[1], y = s_vals[1], z = max(sub$PBR),
          marker = list(size = 8, color = RColorBrewer::brewer.pal(9, hues[i])[7]),
          name = paste0("Fr = ", fr_val), showlegend = TRUE
        )
    }
    fig3d <- fig3d %>%
      plotly::layout(
        title = list(text = paste0("PBR response surfaces across s, alpha and Fr (Nmin = ",
                                    format(nmin_assumed, big.mark = ","), ")")),
        scene = list(
          xaxis = list(title = "Age at first breeding (alpha)"),
          yaxis = list(title = "Adult survival (s)"),
          zaxis = list(title = "PBR (bats/year)")
        ),
        legend = list(title = list(text = "Recovery factor scenario"))
      )

    fig_response_surfaces_3d_path <- file.path(fig_dir, "response_surfaces_3d.png")
    export_ok <- tryCatch({
      plotly::save_image(fig3d, file = fig_response_surfaces_3d_path, width = 900, height = 700)
      TRUE
    }, error = function(e) {
      message(
        "Could not export the 3D response-surface figure as a static image (",
        conditionMessage(e), "). The report still includes the 2D faceted version."
      )
      FALSE
    })
    if (export_ok) fig_response_surfaces_3d <- fig_response_surfaces_3d_path
  } else {
    message(
      "Packages 'plotly'/'RColorBrewer'/'kaleido' not all installed -- skipping the static ",
      "3D response-surface print (install.packages(c(\"plotly\", \"RColorBrewer\", \"kaleido\")) to include it). ",
      "The report still includes the 2D faceted version."
    )
  }

  ## ---- 4. Monte Carlo simulation ------------------------------------------
  set.seed(mc_seed)
  sim <- tibble::tibble(
    s     = runif(mc_n_sim, s_range[1], s_range[2]),
    alpha = runif(mc_n_sim, alpha_range[1], alpha_range[2])
  ) %>%
    mutate(
      lambda_max = lambda_max_niel(s, alpha),
      PBR = pbr_from_components(nmin_assumed, fr_corrected, lambda_max)
    )

  pbr_quantiles <- quantile(sim$PBR, c(0.05, 0.25, 0.50, 0.75, 0.95))

  threshold_percentiles <- pbr_thresholds %>%
    mutate(
      percentile_rank = sapply(threshold, function(x) mean(sim$PBR < x) * 100),
      median_ratio = median(sim$PBR) / threshold,
      p95_ratio = quantile(sim$PBR, 0.95) / threshold
    )

  sim_prior_alpha <- tibble::tibble(
    s     = runif(mc_n_sim, s_range[1], s_range[2]),
    alpha = runif(mc_n_sim, alpha_range_prior[1], alpha_range_prior[2])
  ) %>%
    mutate(
      lambda_max = lambda_max_niel(s, alpha),
      PBR = pbr_from_components(nmin_assumed, fr_corrected, lambda_max)
    )

  alpha_label_updated <- sprintf("%s - %s (updated)", alpha_range[1], alpha_range[2])
  alpha_label_prior <- sprintf("%s - %s (prior)", alpha_range_prior[1], alpha_range_prior[2])

  # Built row-by-row (not via group_by()/summarise(), which reorders groups
  # alphabetically) so each row's quantiles and "% below <threshold>"
  # columns are guaranteed to come from the same simulation, in the order
  # the scenarios are listed below (prior range, then the updated range).
  summarise_alpha_scenario <- function(label, pbr_vec) {
    row <- tibble::tibble(
      alpha_range_label = label,
      q05 = quantile(pbr_vec, 0.05), q50 = quantile(pbr_vec, 0.50), q95 = quantile(pbr_vec, 0.95)
    )
    for (i in seq_len(nrow(pbr_thresholds))) {
      col_name <- sprintf("pct_below_%s", pbr_thresholds$facility[i])
      row[[col_name]] <- mean(pbr_vec < pbr_thresholds$threshold[i]) * 100
    }
    row
  }

  alpha_range_comparison <- bind_rows(
    summarise_alpha_scenario(alpha_label_prior, sim_prior_alpha$PBR),
    summarise_alpha_scenario(alpha_label_updated, sim$PBR)
  )

  fig_pbr_density <- file.path(fig_dir, "pbr_density.png")
  ggsave(fig_pbr_density, width = 7, height = 4, dpi = 150, plot = {
    ggplot(sim, aes(x = PBR)) +
      geom_density(fill = "grey80", colour = "grey40") +
      geom_vline(data = pbr_thresholds, aes(xintercept = threshold, colour = facility),
                 linewidth = 0.8, linetype = "dashed") +
      labs(
        title = "Simulated PBR distribution vs. imposed thresholds",
        subtitle = paste0("Fr = ", fr_corrected, "; Nmin = ", format(nmin_assumed, big.mark = ","),
                           "; s ~ U(", s_range[1], ", ", s_range[2], "); alpha ~ U(",
                           alpha_range[1], ", ", alpha_range[2], ")"),
        x = "PBR (bats/year)", y = "Density", colour = "Imposed threshold"
      ) +
      theme_minimal()
  })

  ## ---- 5. Assemble report params -----------------------------------------
  list(
    project_ref = project_ref,
    species_name = species_name,
    species_common_name = species_common_name,
    thresholds = pbr_thresholds,
    fr_original = fr_original,
    fr_corrected = fr_corrected,
    iucn_status_original_guess = iucn_status_original_guess,
    iucn_status_corrected = iucn_status_corrected,
    backsolve_grid = backsolve_grid,
    nmin_assumed = nmin_assumed,
    project_footprint_km2 = project_footprint_km2,
    combined_region_km2 = combined_region_km2,
    interproject_distance_km = interproject_distance_km,
    nmin_plausibility = nmin_plausibility,
    density_proxy_low = density_proxy_low,
    density_proxy_high = density_proxy_high,
    pbr_corrected_fr = pbr_corrected_fr,
    s_range = s_range,
    alpha_range = alpha_range,
    alpha_range_prior = alpha_range_prior,
    elasticity_ranges = elasticity_ranges,
    fig_elasticity_survival = fig_elasticity_survival,
    fig_elasticity_alpha = fig_elasticity_alpha,
    fig_lambda_max = fig_lambda_max,
    fig_response_surfaces = fig_response_surfaces,
    fig_response_surfaces_3d = fig_response_surfaces_3d,
    fig_pbr_density = fig_pbr_density,
    n_sim = mc_n_sim,
    pbr_quantiles = pbr_quantiles,
    threshold_percentiles = threshold_percentiles,
    alpha_range_comparison = alpha_range_comparison
  )
}
