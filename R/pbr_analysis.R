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
source("R/plotly_screenshot.R")

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

  # One colour per benchmark contour (named by its own value) instead of a
  # single colour for both lines -- two same-coloured, unlabelled isolines
  # are not distinguishable from each other. Avoid white: invisible as a
  # legend key swatch against the (white) legend background, even though
  # it reads fine on the dark raster itself. Reused on the two elasticity
  # plots below as well as fig_lambda_max, so a reader can see directly
  # where the paper's lambda_max benchmarks sit within the same (s, alpha)
  # space the elasticity values are plotted over.
  benchmark_labels <- as.character(lambda_max_benchmarks)
  benchmark_colours <- setNames(c("cyan", "chartreuse")[seq_along(lambda_max_benchmarks)], benchmark_labels)
  lambda_contour <- function() {
    list(
      geom_contour(aes(z = lambda_max, colour = after_stat(factor(level))),
                   breaks = lambda_max_benchmarks, linewidth = 0.5),
      scale_colour_manual(name = "Benchmark", values = benchmark_colours)
    )
  }

  fig_elasticity_survival <- file.path(fig_dir, "elasticity_survival.png")
  ggsave(fig_elasticity_survival, width = 7, height = 4.5, dpi = 150, plot = {
    ggplot(sens_grid, aes(x = s, y = alpha)) +
      geom_raster(aes(fill = E_s), interpolate = TRUE) +
      lambda_contour() +
      scale_fill_viridis_c(name = "E[s]", option = "B", direction = -1) +
      labs(x = "Adult survival (s)", y = "Age at first breeding (alpha)",
           title = "Elasticity of lambda_max to adult survival",
           subtitle = lambda_range_txt) +
      theme_minimal() +
      theme(plot.subtitle = element_text(size = 9))
  })

  fig_elasticity_alpha <- file.path(fig_dir, "elasticity_alpha.png")
  ggsave(fig_elasticity_alpha, width = 7, height = 4.5, dpi = 150, plot = {
    ggplot(sens_grid, aes(x = s, y = alpha)) +
      geom_raster(aes(fill = E_alpha), interpolate = TRUE) +
      lambda_contour() +
      scale_fill_viridis_c(name = "E[alpha]", option = "A", direction = -1) +
      labs(x = "Adult survival (s)", y = "Age at first breeding (alpha)",
           title = "Elasticity of lambda_max to age at first breeding",
           subtitle = lambda_range_txt) +
      theme_minimal() +
      theme(plot.subtitle = element_text(size = 9))
  })

  fig_lambda_max <- file.path(fig_dir, "lambda_max.png")
  ggsave(fig_lambda_max, width = 8, height = 4.5, dpi = 150, plot = {
    ggplot(sens_grid, aes(x = s, y = alpha)) +
      geom_raster(aes(fill = lambda_max), interpolate = TRUE) +
      lambda_contour() +
      scale_fill_viridis_c(name = "lambda_max", option = "C") +
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
  ##          report's plot.ly figure) -- a self-contained HTML widget with
  ##          a FIXED camera (so the static export is reproducible across
  ##          runs/machines -- same convention as the 3D turbine coverage
  ##          maps in idf_bsh_dgy/R/coverage_3d_topography.R). The
  ##          self-contained interactive HTML widget is always kept as a
  ##          standalone deliverable (Word/PowerPoint cannot embed a live
  ##          plotly widget, but the HTML opens directly in any browser,
  ##          no server needed); a static PNG for the docx is additionally
  ##          screenshotted via webshot2 (headless Chrome/Edge,
  ##          R/plotly_screenshot.R) if available. Avoids 'kaleido',
  ##          fragile to install across platforms; the 2D faceted version
  ##          above always covers the same content regardless.
  fig_response_surfaces_3d <- NULL
  fig_response_surfaces_3d_html <- NULL
  has_3d_deps <- requireNamespace("plotly", quietly = TRUE) &&
    requireNamespace("RColorBrewer", quietly = TRUE) &&
    requireNamespace("htmlwidgets", quietly = TRUE)

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
          zaxis = list(title = "PBR (bats/year)"),
          aspectmode = "manual",
          aspectratio = list(x = 1, y = 1, z = 0.8),
          # Fixed viewpoint (not plotly's own default, which is not
          # guaranteed stable across versions/machines) -- same idea as
          # .camera_from_bearing() for the 3D turbine coverage maps, though
          # this scene has no compass bearing to anchor to, so the eye
          # position is just a fixed, deliberately chosen 3/4 elevated view
          # showing all three axes clearly.
          camera = list(
            eye = list(x = 1.35, y = -1.35, z = 0.7),
            center = list(x = 0, y = 0, z = -0.15),
            up = list(x = 0, y = 0, z = 1)
          )
        ),
        legend = list(title = list(text = "Recovery factor scenario"), x = 0.82, y = 0.85),
        margin = list(l = 0, r = 0, b = 0, t = 40)
      )

    fig_response_surfaces_3d_path <- file.path(fig_dir, "response_surfaces_3d.png")
    fig_response_surfaces_3d_html_path <- file.path(fig_dir, "response_surfaces_3d.html")
    screenshot_result <- render_plotly_screenshot(
      fig3d, fig_response_surfaces_3d_path, width = 1000, height = 750,
      html_path = fig_response_surfaces_3d_html_path
    )
    fig_response_surfaces_3d <- screenshot_result$png
    fig_response_surfaces_3d_html <- screenshot_result$html
  } else {
    message(
      "Packages 'plotly'/'RColorBrewer'/'htmlwidgets' not all installed -- skipping the 3D response-surface ",
      "figure entirely (install.packages(c(\"plotly\", \"RColorBrewer\", \"htmlwidgets\")) to include it, plus ",
      "'webshot2' for the static image embedded in the docx). ",
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

  ## ---- 5b. Leslie matrix cross-check (Safi 2006; Dekker & Limpens 2024) --
  ## Validated against the published model: female 2-stage matrix
  ## [[0, F], [S_juv, S_adult]] with F = p_breed*litter*sex_ratio*S_juv
  ## reproduces the published lambda = 1.047 exactly (see
  ## R/leslie_dekker_limpens.R, where this was first established, and
  ## references/leslie_matrix_parametrisation.md for the full literature
  ## trail, including the Myotis lucifugus proxy tried before the
  ## species-specific Safi/Dekker & Limpens data were located).
  leslie_F <- leslie_p_breed * leslie_litter * leslie_sex_ratio * leslie_s_juv

  build_dekker_stage_matrix <- function(n_stages, s_juv, s_adult, fecundity) {
    A <- matrix(0, n_stages, n_stages)
    A[1, n_stages] <- fecundity
    for (i in 1:(n_stages - 1)) A[i + 1, i] <- s_juv
    A[n_stages, n_stages] <- s_adult
    A
  }
  leslie_lambda <- function(A) max(Mod(eigen(A, only.values = TRUE)$values))

  A_leslie_validated <- build_dekker_stage_matrix(2, leslie_s_juv, leslie_s_adult_f, leslie_F)
  lambda_leslie_validated <- leslie_lambda(A_leslie_validated)

  leslie_grid <- tidyr::expand_grid(
    s_adult = seq(s_range[1], s_range[2], length.out = 50),
    n_stages = leslie_stage_range
  ) %>%
    mutate(
      lambda_Leslie = mapply(function(s_a, n) leslie_lambda(build_dekker_stage_matrix(n, leslie_s_juv, s_a, leslie_F)),
                              s_adult, n_stages),
      n_stages_label = factor(paste0(n_stages, ifelse(n_stages == 2, " (validated)", "")),
                               levels = paste0(leslie_stage_range, ifelse(leslie_stage_range == 2, " (validated)", "")))
    )

  # Where does the validated structure's lambda sit relative to the PBR
  # benchmarks, across the same adult-survival range used throughout?
  leslie_vs_pbr <- leslie_grid %>%
    filter(n_stages == 2) %>%
    summarise(
      lambda_min = min(lambda_Leslie), lambda_max_val = max(lambda_Leslie),
      reaches_120 = max(lambda_Leslie) >= 1.20, reaches_124 = max(lambda_Leslie) >= 1.24
    )

  fig_leslie_maturation <- file.path(fig_dir, "leslie_maturation_delay.png")
  ggsave(fig_leslie_maturation, width = 7.5, height = 4.5, dpi = 150, plot = {
    ggplot(leslie_grid, aes(x = s_adult, y = lambda_Leslie, colour = n_stages_label)) +
      geom_line(linewidth = 0.9) +
      geom_hline(yintercept = c(1.20, 1.24), linetype = "dashed", colour = c("cyan", "chartreuse")) +
      geom_vline(xintercept = leslie_s_adult_f, linetype = "dotted", colour = "grey40") +
      scale_colour_viridis_d(name = "Pre-reproductive\nstages", option = "C") +
      labs(
        x = "Adult female survival", y = "lambda (Leslie matrix)",
        title = "V. murinus lambda vs adult survival, by maturation delay",
        subtitle = paste0("S_juv = ", leslie_s_juv, " (Safi 2006); dotted: Safi's own S_adult=",
                           leslie_s_adult_f, "; dashed: PBR benchmarks 1.20 (cyan)/1.24 (green)")
      ) +
      theme_minimal() +
      theme(plot.subtitle = element_text(size = 8))
  })

  ## ---- 5c. PVA-lite: 25-year stochastic projection ------------------------
  beta_params <- function(mean, cv) {
    var <- (mean * cv)^2
    common <- mean * (1 - mean) / var - 1
    list(shape1 = mean * common, shape2 = (1 - mean) * common)
  }
  bp_s_adult <- beta_params(leslie_s_adult_f, pva_vital_rate_cv)
  bp_s_juv <- beta_params(leslie_s_juv, pva_vital_rate_cv)
  bp_p_breed <- beta_params(leslie_p_breed, pva_vital_rate_cv)

  simulate_pva_trajectory <- function(n0_juv, n0_adult, n_years, annual_removal) {
    # rbinom()'s size argument silently returns NA for a size that isn't
    # (very nearly) a non-negative integer -- round() alone leaves a
    # representable-double residue (e.g. 1367.000000000000227, from the
    # eigenvector-based stable-stage split below) that trips this; wrapping
    # every size in as.integer(round(.)) clears it for good.
    n_juv <- as.integer(round(n0_juv)); n_adult <- as.integer(round(n0_adult))
    traj <- numeric(n_years + 1)
    traj[1] <- 2 * (n_juv + n_adult)
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
  stable_dist <- Re(eigen(A_leslie_validated)$vectors[, 1]); stable_dist <- stable_dist / sum(stable_dist)
  n0_juv <- round(n0_female * stable_dist[1])
  n0_adult <- n0_female - n0_juv

  pva_scenario_labels <- sprintf("%s threshold (%s/yr)", pbr_thresholds$facility, pbr_thresholds$threshold)
  pva_scenarios <- bind_rows(
    tibble::tibble(scenario = "No additional mortality", threshold = 0),
    tibble::tibble(scenario = pva_scenario_labels, threshold = pbr_thresholds$threshold)
  ) %>%
    mutate(scenario = factor(scenario, levels = c("No additional mortality", pva_scenario_labels)))

  set.seed(pva_seed)
  pva_trajectories <- lapply(seq_len(nrow(pva_scenarios)), function(i) {
    reps <- replicate(pva_n_reps, simulate_pva_trajectory(n0_juv, n0_adult, pva_n_years, pva_scenarios$threshold[i]))
    tibble::tibble(
      scenario = pva_scenarios$scenario[i],
      year = rep(0:pva_n_years, pva_n_reps),
      rep = rep(seq_len(pva_n_reps), each = pva_n_years + 1),
      N = as.vector(reps)
    )
  }) %>% bind_rows()

  pva_summary_by_year <- pva_trajectories %>%
    group_by(scenario, year) %>%
    summarise(median_N = median(N), q05 = quantile(N, 0.05), q25 = quantile(N, 0.25),
              q75 = quantile(N, 0.75), q95 = quantile(N, 0.95), .groups = "drop")

  quasi_ext_threshold <- pva_quasi_extinction_fraction * (2 * n0_female)
  pva_risk_summary <- pva_trajectories %>%
    group_by(scenario, rep) %>%
    summarise(final_N = N[year == pva_n_years], ever_below_threshold = any(N < quasi_ext_threshold), .groups = "drop") %>%
    group_by(scenario) %>%
    summarise(
      median_final_N = median(final_N), q05_final_N = quantile(final_N, 0.05), q95_final_N = quantile(final_N, 0.95),
      p_decline = mean(final_N < 2 * n0_female) * 100, p_quasi_extinction = mean(ever_below_threshold) * 100,
      .groups = "drop"
    )

  fig_pva_projection <- file.path(fig_dir, "pva_25yr_projection.png")
  ggsave(fig_pva_projection, width = 11, height = 4.5, dpi = 150, plot = {
    ggplot(pva_summary_by_year, aes(x = year)) +
      geom_ribbon(aes(ymin = q05, ymax = q95, fill = scenario), alpha = 0.15) +
      geom_ribbon(aes(ymin = q25, ymax = q75, fill = scenario), alpha = 0.3) +
      geom_line(aes(y = median_N, colour = scenario), linewidth = 0.9) +
      geom_hline(yintercept = 2 * n0_female, linetype = "dotted", colour = "grey40") +
      facet_wrap(~scenario, nrow = 1) +
      scale_colour_viridis_d(option = "C", end = 0.8, guide = "none") +
      scale_fill_viridis_d(option = "C", end = 0.8, guide = "none") +
      labs(
        x = "Year", y = "Population (both sexes)",
        title = paste0(pva_n_years, "-year stochastic projection, validated V. murinus demographic model"),
        subtitle = paste0("N0 = ", 2 * n0_female, " (Nmin per facility); shaded: 50%/90% of ",
                           pva_n_reps, " trajectories; dotted: N0")
      ) +
      theme_minimal() +
      theme(strip.text = element_text(face = "bold"), plot.subtitle = element_text(size = 8))
  })

  ## ---- 5d. Three-way comparison summary --------------------------------
  comparison_summary <- tibble::tibble(
    Method = c("PBR / Niel-Lebreton (fixed benchmarks)", "PBR / Niel-Lebreton (Monte Carlo median)",
               "Leslie matrix (validated, Safi's S_adult=0.76)", "PVA-lite (median trend, no removal)"),
    `Growth rate (lambda)` = c(
      sprintf("%.2f - %.2f", 1.20, 1.24),
      sprintf("%.2f", stats::median(sim$lambda_max)),
      sprintf("%.3f", lambda_leslie_validated),
      sprintf("~%.3f (implied)", (pva_summary_by_year %>% filter(scenario == "No additional mortality", year == pva_n_years) %>% pull(median_N) / (2*n0_female))^(1/pva_n_years)
      )
    ),
    `Consistent with imposed thresholds?` = c(
      "By construction (thresholds back-solved from these)",
      sprintf("Thresholds sit at %.0f-%.0fth percentile (conservative)", threshold_percentiles$percentile_rank[1], threshold_percentiles$percentile_rank[2]),
      sprintf("Only reached at S_adult>=%.2f (top of range)", leslie_grid %>% filter(n_stages==2, lambda_Leslie>=1.20) %>% summarise(m=min(s_adult)) %>% pull(m)),
      sprintf("%.0f-%.0f%% probability of decline over %d yrs under imposed thresholds", min(pva_risk_summary$p_decline[pva_risk_summary$scenario!="No additional mortality"]), max(pva_risk_summary$p_decline[pva_risk_summary$scenario!="No additional mortality"]), pva_n_years)
    )
  )

  ## ---- 5e. Inverse demographic boundary analysis ("stress test") --------
  ## Where would the validated Leslie structure's four vital rates need to
  ## sit to reach the PBR benchmarks (1.20/1.24)? Framing (Paulo, 2026-09):
  ## lambda from this matrix is driven by OBSERVED/empirical vital rates
  ## (Safi 2006) -- lambda_observed, not lambda_max, PBR's theoretical
  ## ceiling under favourable, non-resource-limited conditions. Showing
  ## this matrix does not reach 1.20/1.24 at Safi's own values therefore
  ## does NOT by itself show those benchmarks are biologically impossible;
  ## this section instead characterises what demographic combinations would
  ## be needed, and how they compare with the empirical baseline and with
  ## plausible literature ranges (R/leslie_boundary_analysis.R, where this
  ## was first developed and run standalone). No calibration to the
  ## published elasticities (0.530/0.235/0.235) is attempted -- the lambda
  ## match alone (1.0474 vs 1.047) validates the reconstruction.
  boundary_baseline <- list(s_juv = leslie_s_juv, s_adult = leslie_s_adult_f,
                            p_breed = leslie_p_breed, litter = leslie_litter)
  boundary_lambda_of <- function(s_juv, s_adult, p_breed, litter) {
    F_rec <- p_breed * litter * leslie_sex_ratio * s_juv
    leslie_lambda(build_dekker_stage_matrix(2, s_juv, s_adult, F_rec))
  }
  stopifnot(abs(do.call(boundary_lambda_of, boundary_baseline) - lambda_leslie_validated) < 1e-9)

  boundary_ranges <- list(s_adult = leslie_boundary_range_s_adult, s_juv = leslie_boundary_range_s_juv,
                           p_breed = leslie_boundary_range_p_breed, litter = leslie_boundary_range_litter)
  boundary_search_ceiling <- list(s_adult = 0.999, s_juv = 0.99, p_breed = 1.0, litter = 6)
  boundary_param_labels <- c(s_adult = "Adult female survival", s_juv = "Juvenile survival",
                              p_breed = "Breeding fraction", litter = "Litter size")

  boundary_solve_breakeven <- function(param, target, search_upper) {
    f <- function(x) {
      args <- boundary_baseline
      args[[param]] <- x
      do.call(boundary_lambda_of, args) - target
    }
    lo <- boundary_baseline[[param]]
    if (f(lo) >= 0) return(lo)
    if (f(search_upper) < 0) return(NA)
    uniroot(f, c(lo, search_upper))$root
  }

  leslie_breakeven_tbl <- tidyr::expand_grid(param = names(boundary_baseline), target = lambda_max_benchmarks) %>%
    rowwise() %>%
    mutate(
      baseline_value = boundary_baseline[[param]],
      required_value = boundary_solve_breakeven(param, target, boundary_search_ceiling[[param]]),
      plausible_upper = boundary_ranges[[param]][2],
      pct_of_plausible_range = 100 * (required_value - boundary_ranges[[param]][1]) / (plausible_upper - boundary_ranges[[param]][1]),
      within_plausible_range = !is.na(required_value) & required_value <= plausible_upper
    ) %>%
    ungroup() %>%
    mutate(parameter_label = unname(boundary_param_labels[param])) %>%
    select(parameter_label, lambda_target = target, safi_baseline = baseline_value,
           required_value, plausible_upper, pct_of_plausible_range, within_plausible_range)

  # Full 4-way factorial over all vital rates' plausible ranges (n_stages
  # fixed at 2, the validated -- and fastest-maturing, hence most
  # favourable to reaching high lambda -- structure).
  boundary_full_grid <- tidyr::expand_grid(
    s_adult = seq(leslie_boundary_range_s_adult[1], leslie_boundary_range_s_adult[2], length.out = leslie_boundary_grid_resolution),
    s_juv   = seq(leslie_boundary_range_s_juv[1], leslie_boundary_range_s_juv[2], length.out = leslie_boundary_grid_resolution),
    p_breed = seq(leslie_boundary_range_p_breed[1], leslie_boundary_range_p_breed[2], length.out = leslie_boundary_grid_resolution),
    litter  = seq(leslie_boundary_range_litter[1], leslie_boundary_range_litter[2], length.out = leslie_boundary_grid_resolution)
  ) %>%
    mutate(
      lambda = mapply(boundary_lambda_of, s_juv, s_adult, p_breed, litter),
      band = cut(lambda, breaks = c(-Inf, 1.20, 1.24, Inf), labels = c("< 1.20", "1.20-1.24", ">= 1.24"), right = FALSE)
    )
  leslie_boundary_band_shares <- boundary_full_grid %>% count(band) %>% mutate(pct = 100 * n / sum(n))
  leslie_boundary_pct_below_baseline <- 100 * mean(boundary_full_grid$lambda < lambda_leslie_validated)

  surface_grid_boundary <- tidyr::expand_grid(
    s_adult = seq(leslie_boundary_range_s_adult[1], leslie_boundary_range_s_adult[2], length.out = 80),
    litter  = seq(leslie_boundary_range_litter[1], leslie_boundary_range_litter[2], length.out = 80),
    scenario = c("Baseline S_juv & breeding fraction", "Favourable S_juv & breeding fraction")
  ) %>%
    mutate(
      s_juv_use   = ifelse(scenario == "Baseline S_juv & breeding fraction", leslie_s_juv, leslie_boundary_range_s_juv[2]),
      p_breed_use = ifelse(scenario == "Baseline S_juv & breeding fraction", leslie_p_breed, leslie_boundary_range_p_breed[2]),
      lambda = mapply(boundary_lambda_of, s_juv_use, s_adult, p_breed_use, litter)
    )
  boundary_baseline_point <- tibble::tibble(s_adult = leslie_s_adult_f, litter = leslie_litter,
                                             scenario = "Baseline S_juv & breeding fraction")

  fig_leslie_boundary_surface <- file.path(fig_dir, "leslie_boundary_surface.png")
  ggsave(fig_leslie_boundary_surface, width = 10, height = 4.8, dpi = 150, plot = {
    ggplot(surface_grid_boundary, aes(x = s_adult, y = litter, z = lambda)) +
      geom_raster(aes(fill = lambda)) +
      geom_contour(breaks = lambda_max_benchmarks, aes(colour = after_stat(factor(level))), linewidth = 0.8) +
      scale_colour_manual(name = "lambda contour", values = benchmark_colours) +
      scale_fill_viridis_c(option = "D", name = "lambda") +
      geom_point(data = boundary_baseline_point, aes(x = s_adult, y = litter), inherit.aes = FALSE,
                 colour = "red", size = 2.5, shape = 17) +
      facet_wrap(~scenario) +
      labs(
        x = "Adult female survival", y = "Litter size",
        title = "Where would adult survival and litter size need to sit to reach the PBR benchmarks?",
        subtitle = paste0(
          "Red triangle: Safi's empirical baseline (S_adult=", leslie_s_adult_f, ", litter=", leslie_litter,
          "); contours: PBR benchmarks ", paste(lambda_max_benchmarks, collapse = " (cyan) / "), " (green)"
        )
      ) +
      theme_minimal() +
      theme(plot.subtitle = element_text(size = 8), strip.text = element_text(face = "bold"))
  })

  ## ---- 5f. Relative sensitivity of PBR to Nmin vs lambda_max --------------
  ## PBR = 0.5*(lambda_max-1)*Fr*Nmin is linear in BOTH Nmin and
  ## (lambda_max-1), but the two are not equally uncertain in practice: no
  ## independent population estimate exists (Population size (Nmin)
  ## section), while lambda_max is pinned to the paper's own two fixed
  ## benchmarks. A direct, simple comparison: a plausible doubling of Nmin
  ## vs. the full width of the paper's own lambda_max benchmark range,
  ## holding everything else fixed.
  nmin_vs_lambda_sensitivity <- tibble::tibble(
    change = c(
      sprintf("Nmin: %s -> %s (lambda_max = %.2f)", format(nmin_assumed, big.mark = ","),
              format(2 * nmin_assumed, big.mark = ","), min(lambda_max_benchmarks)),
      sprintf("lambda_max: %.2f -> %.2f (Nmin = %s)", min(lambda_max_benchmarks), max(lambda_max_benchmarks),
              format(nmin_assumed, big.mark = ","))
    ),
    pbr_before = c(
      pbr_from_components(nmin_assumed, fr_corrected, min(lambda_max_benchmarks)),
      pbr_from_components(nmin_assumed, fr_corrected, min(lambda_max_benchmarks))
    ),
    pbr_after = c(
      pbr_from_components(2 * nmin_assumed, fr_corrected, min(lambda_max_benchmarks)),
      pbr_from_components(nmin_assumed, fr_corrected, max(lambda_max_benchmarks))
    )
  ) %>%
    mutate(pct_change = 100 * (pbr_after - pbr_before) / pbr_before)

  ## ---- 5g. Nmin x lambda_max PBR surface -----------------------------------
  ## Synthesises the report's two uncertainty axes in one figure: Nmin
  ## (2,000-20,000, spanning the reverse-engineered value and the
  ## migratory-population hypothesis discussed with Paulo) x lambda_max
  ## (1.04-1.24, spanning the validated Leslie matrix's own empirical value
  ## through the paper's fixed benchmarks), Fr held at the corrected value.
  ## Nmin_assumed is marked explicitly as a reverse-engineered scenario, not
  ## an abundance estimate (Reverse-engineering the imposed thresholds).
  nmin_lambda_grid <- tidyr::expand_grid(
    N = seq(2000, 20000, length.out = 200),
    lambda = seq(1.04, 1.24, length.out = 200)
  ) %>%
    mutate(PBR = pbr_from_components(N, fr_corrected, lambda))

  nmin_lambda_breaks <- sort(unique(c(pbr_thresholds$threshold, 200, 240, 500, 1000)))
  nmin_lambda_break_colours <- setNames(
    c("white", "grey75", "cyan", "chartreuse", "yellow", "red")[seq_along(nmin_lambda_breaks)],
    as.character(nmin_lambda_breaks)
  )

  fig_nmin_lambda_surface <- file.path(fig_dir, "nmin_lambda_surface.png")
  ggsave(fig_nmin_lambda_surface, width = 8.5, height = 5.5, dpi = 150, plot = {
    ggplot(nmin_lambda_grid, aes(x = N, y = lambda)) +
      geom_raster(aes(fill = PBR), interpolate = TRUE) +
      geom_contour(aes(z = PBR, colour = after_stat(factor(level))), breaks = nmin_lambda_breaks, linewidth = 0.45) +
      geom_vline(xintercept = nmin_assumed, linetype = "dashed", colour = "white", linewidth = 0.4) +
      annotate(
        "label", x = nmin_assumed, y = 1.235,
        label = paste0("Nmin = ", format(nmin_assumed, big.mark = ","), "\n(reverse-engineered scenario,\nnot an abundance estimate)"),
        colour = "grey20", fill = "white", alpha = 0.85, size = 2.6, hjust = -0.03, vjust = 1, label.size = 0
      ) +
      scale_fill_viridis_c(option = "C", name = "PBR\n(bats/yr)") +
      scale_colour_manual(name = "PBR contour\n(bats/yr)", values = nmin_lambda_break_colours) +
      labs(
        x = "Nmin", y = "lambda_max",
        title = "PBR as a function of Nmin and lambda_max",
        subtitle = paste0(
          "Fr = ", fr_corrected, " throughout. PBR moves the same amount for a given % change in Nmin as for the ",
          "same % change in (lambda_max - 1) -- but Nmin's plausible range is far less constrained."
        )
      ) +
      scale_x_continuous(labels = scales::comma) +
      theme_minimal() +
      theme(plot.subtitle = element_text(size = 7.5))
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
    fig_response_surfaces_3d_html = fig_response_surfaces_3d_html,
    fig_pbr_density = fig_pbr_density,
    n_sim = mc_n_sim,
    pbr_quantiles = pbr_quantiles,
    threshold_percentiles = threshold_percentiles,
    alpha_range_comparison = alpha_range_comparison,
    leslie_s_juv = leslie_s_juv,
    leslie_s_adult_f = leslie_s_adult_f,
    leslie_s_adult_m = leslie_s_adult_m,
    leslie_p_breed = leslie_p_breed,
    leslie_litter = leslie_litter,
    leslie_F = leslie_F,
    lambda_leslie_validated = lambda_leslie_validated,
    leslie_vs_pbr = leslie_vs_pbr,
    fig_leslie_maturation = fig_leslie_maturation,
    n0_total = 2 * n0_female,
    pva_n_years = pva_n_years,
    pva_n_reps = pva_n_reps,
    pva_vital_rate_cv = pva_vital_rate_cv,
    quasi_ext_threshold = quasi_ext_threshold,
    pva_risk_summary = pva_risk_summary,
    fig_pva_projection = fig_pva_projection,
    comparison_summary = comparison_summary,
    leslie_breakeven_tbl = leslie_breakeven_tbl,
    leslie_boundary_band_shares = leslie_boundary_band_shares,
    leslie_boundary_pct_below_baseline = leslie_boundary_pct_below_baseline,
    leslie_boundary_range_s_adult = leslie_boundary_range_s_adult,
    leslie_boundary_range_s_juv = leslie_boundary_range_s_juv,
    leslie_boundary_range_p_breed = leslie_boundary_range_p_breed,
    leslie_boundary_range_litter = leslie_boundary_range_litter,
    fig_leslie_boundary_surface = fig_leslie_boundary_surface,
    nmin_vs_lambda_sensitivity = nmin_vs_lambda_sensitivity,
    fig_nmin_lambda_surface = fig_nmin_lambda_surface
  )
}
