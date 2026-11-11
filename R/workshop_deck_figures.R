##
## Figures for the European-specialist workshop deck (Paulo, 2026-11),
## rebuilt in R from the pptxgenjs v1 draft so the deck can be regenerated
## from this project's own pipeline and restyled from
## inputs/deckSettings_workshop.R alone -- no plotting function below
## should ever need editing just to change a font size, hide a subtitle,
## or resize a point/symbol.
##
## Every function takes `controls` (deck_plot_controls, or a modified
## copy of it) and returns a ggplot object; it does not ggsave itself.
## save_deck_figure() is the one place that writes a PNG, so every figure
## in this deck is saved the same way (bg = "white" always -- see
## R/carrying_capacity_reference.R's own black-background fix this
## project already standardised on).
##

suppressPackageStartupMessages({ library(ggplot2); library(dplyr); library(tidyr); library(scales) })

## ---- Shared theme -- the one place that reads base_size/title_size/
## subtitle_show, so every figure in the deck stays visually consistent
## when Paulo tunes the controls list. ------------------------------------
deck_theme <- function(controls) {
  theme_minimal(base_size = controls$base_size) +
    theme(
      plot.title    = element_text(size = controls$title_size, face = "bold", colour = controls$colour_dark),
      plot.subtitle = element_text(size = controls$base_size * 0.82, colour = controls$colour_muted),
      plot.caption  = element_text(size = controls$base_size * 0.62, colour = controls$colour_muted, hjust = 0, face = "italic"),
      axis.text     = element_text(size = controls$base_size * controls$axis_text_size_mult, colour = "#1A1A1A"),
      axis.title    = element_blank(),
      legend.title  = element_blank(),
      legend.text   = element_text(size = controls$base_size * 0.85),
      legend.position = "bottom",
      panel.grid.minor = element_blank(),
      panel.grid.major.x = element_blank(),
    )
}

## subtitle_show/caption_show are applied here, once, rather than inside
## every plotting function -- flipping either flag in the controls list
## drops it from every figure in the deck in one place.
deck_labs <- function(controls, title, subtitle = NULL, caption = NULL, ...) {
  labs(
    title = title,
    subtitle = if (isTRUE(controls$subtitle_show)) subtitle else NULL,
    caption  = if (isTRUE(controls$caption_show)) caption else NULL,
    ...
  )
}

save_deck_figure <- function(plot, path, controls, width = NULL, height = NULL) {
  ggsave(path, plot = plot,
         width = width %||% controls$fig_width, height = height %||% controls$fig_height,
         dpi = controls$dpi, bg = controls$bg)
  path
}
`%||%` <- function(a, b) if (is.null(a)) b else a

## =========================================================================
## Slide 4 -- turbines curtailed vs. fleet
## =========================================================================
plot_turbine_curtailment <- function(turbines, controls) {
  dt <- turbines %>%
    mutate(not_curtailed = assessed - curtailed) %>%
    select(site, Curtailed = curtailed, `Not curtailed` = not_curtailed) %>%
    pivot_longer(-site, names_to = "status", values_to = "n") %>%
    mutate(status = factor(status, levels = c("Not curtailed", "Curtailed")))

  ggplot(dt, aes(x = site, y = n, fill = status)) +
    geom_col(width = controls$bar_width) +
    geom_text(aes(label = n), position = position_stack(vjust = 0.5),
              colour = "white", size = controls$label_size, fontface = "bold") +
    scale_fill_manual(values = c("Not curtailed" = "#D7E4EC", "Curtailed" = controls$colour_primary)) +
    deck_labs(controls,
      title = "Turbines curtailed vs. fleet",
      subtitle = paste0("Cut-in raised to ", controls$cut_in_speed_ms %||% 6.5, " m/s on the highest-risk subset"),
      caption = "Site labels and turbine counts are anonymised operational facts (not rounded)."
    ) +
    deck_theme(controls) +
    theme(axis.text.y = element_blank(), panel.grid.major.y = element_blank())
}

## =========================================================================
## Slide 5 -- % reduction vs. counterfactual
## =========================================================================
plot_reduction_results <- function(reduction, controls) {
  ggplot(reduction, aes(x = site, y = pct_reduction_vs_counterfactual)) +
    geom_col(width = controls$bar_width, fill = controls$colour_primary) +
    geom_text(aes(label = paste0(pct_reduction_vs_counterfactual, "%")),
              vjust = -0.4, size = controls$label_size * 1.6, fontface = "bold", colour = controls$colour_dark) +
    scale_y_continuous(limits = c(0, 100), expand = expansion(mult = c(0, 0.18))) +
    deck_labs(controls,
      title = "Observed mortality reduction vs. the no-curtailment counterfactual",
      subtitle = "Operational effectiveness, not a randomised causal estimate",
      caption = "Relative figures are more robust to the detection-correction factor than absolute counts."
    ) +
    deck_theme(controls) +
    theme(axis.text.y = element_blank(), panel.grid.major.y = element_blank())
}

## =========================================================================
## Slide 7 -- elasticity of lambda_max (diverging range plot)
## =========================================================================
plot_lambda_elasticity <- function(elasticity, controls) {
  dt <- elasticity %>% mutate(parameter = factor(parameter, levels = rev(parameter)))
  ggplot(dt, aes(y = parameter)) +
    geom_segment(aes(x = low, xend = high, yend = parameter), colour = "grey70", linewidth = controls$line_width * 1.3) +
    geom_point(aes(x = low), size = controls$point_size, colour = controls$colour_secondary) +
    geom_point(aes(x = high), size = controls$point_size, colour = controls$colour_primary) +
    geom_text(aes(x = low, label = sprintf("%.2f", low)), vjust = -1.1, size = controls$label_size, colour = controls$colour_secondary) +
    geom_text(aes(x = high, label = sprintf("%.2f", high)), vjust = -1.1, size = controls$label_size, colour = controls$colour_primary) +
    geom_vline(xintercept = 0, colour = "grey50", linewidth = 0.4) +
    deck_labs(controls,
      title = "Elasticity of lambda_max",
      subtitle = "% change per 1% parameter change -- low (sky) to high (teal) end of the explored range",
      caption = "s in [0.70, 0.95]; alpha in [1.5, 3.5] years (Niel & Lebreton 2005 sensitivity analysis)."
    ) +
    deck_theme(controls)
}

## =========================================================================
## Slide 8 -- implied Nmin range, by lambda_max benchmark
## =========================================================================
plot_nmin_range <- function(nmin_range, controls) {
  dt <- nmin_range %>%
    mutate(lambda_label = paste0("lambda_max = ", lambda_max)) %>%
    group_by(site) %>%
    mutate(site = factor(site, levels = rev(unique(site)))) %>%
    ungroup()

  rng <- dt %>% group_by(site) %>% summarise(lo = min(nmin_required_floor), hi = max(nmin_required_floor), .groups = "drop")

  ggplot() +
    geom_segment(data = rng, aes(x = lo, xend = hi, y = site, yend = site), colour = "grey70", linewidth = controls$line_width * 1.3) +
    geom_point(data = dt, aes(x = nmin_required_floor, y = site, colour = lambda_label), size = controls$point_size + 1) +
    geom_text(data = dt, aes(x = nmin_required_floor, y = site, label = comma(nmin_required_floor)),
              vjust = -1.2, size = controls$label_size) +
    scale_colour_manual(values = setNames(c(controls$colour_secondary, controls$colour_primary), unique(dt$lambda_label))) +
    scale_x_continuous(labels = comma) +
    deck_labs(controls,
      title = "Implied local abundance (Nmin), by lambda_max benchmark",
      subtitle = "No single population size reproduces both sites' thresholds",
      caption = "Under Fr = 0.3 (Near Threatened); values floored to the nearest 100 individuals."
    ) +
    deck_theme(controls)
}

## =========================================================================
## Slide 9 -- Schaefer/logistic sustainable-H curve (conceptual)
## =========================================================================
plot_schaefer_curve <- function(controls, n = 200) {
  dt <- tibble::tibble(alpha = seq(0, 1, length.out = n)) %>%
    mutate(H = alpha * (1 - alpha))

  msy_x <- 0.5
  ggplot(dt, aes(x = alpha, y = H)) +
    geom_line(colour = controls$colour_primary, linewidth = controls$line_width * 1.3) +
    geom_vline(xintercept = msy_x, linetype = "dotted", colour = "grey60") +
    annotate("text", x = msy_x, y = max(dt$H) * 1.04, label = "alpha = 0.5 (MSY)", size = controls$label_size, colour = controls$colour_muted) +
    annotate("text", x = 0.92, y = max(dt$H) * 0.12, label = "alpha near 1\n(near K)", size = controls$label_size, colour = controls$colour_muted) +
    scale_x_continuous(labels = NULL) +
    scale_y_continuous(labels = NULL) +
    deck_labs(controls,
      title = "Sustainable removal H(alpha) = r*K*alpha*(1-alpha)",
      subtitle = "Schematic -- not fitted to this project's data",
      caption = NULL
    ) +
    deck_theme(controls) +
    theme(panel.grid = element_blank())
}

## =========================================================================
## Slide 10 -- exploratory carrying-capacity vs. PBR comparison
## =========================================================================
plot_cc_vs_pbr <- function(cc_comparison, thresholds, controls) {
  dt <- bind_rows(
    cc_comparison %>% transmute(label = scenario, value = H, kind = "Exploratory simulation"),
    thresholds %>% transmute(label = paste0(site, "\nthreshold"), value = threshold, kind = "Current imposed threshold")
  ) %>%
    mutate(label = factor(label, levels = label))

  wrap_labels <- function(x, width = 13) vapply(x, function(s) paste(strwrap(s, width = width), collapse = "\n"), character(1))

  ggplot(dt, aes(x = label, y = value, fill = kind)) +
    geom_col(width = controls$bar_width) +
    geom_text(aes(label = comma(value)), vjust = -0.4, size = controls$label_size * 1.4, fontface = "bold", colour = controls$colour_dark) +
    scale_fill_manual(values = c("Exploratory simulation" = controls$colour_primary, "Current imposed threshold" = controls$colour_dark)) +
    scale_x_discrete(labels = wrap_labels) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.15))) +
    deck_labs(controls,
      title = "Simulated sustainable H vs. imposed thresholds",
      subtitle = "Exploratory -- not part of the reviewed technical note",
      caption = "Age-structured, stochastic simulation at an alpha=0.8 target; values floored."
    ) +
    deck_theme(controls) +
    theme(axis.text.y = element_blank(), panel.grid.major.y = element_blank())
}

## =========================================================================
## Slide 11 -- FCS conceptual spectrum (extinction <-> carrying capacity)
## =========================================================================
plot_fcs_spectrum <- function(controls, pbr_x = 0.18, fcs_x = 0.85) {
  band <- tibble::tibble(x = c(0, 1), y = 0)
  pts <- tibble::tibble(
    x = c(pbr_x, fcs_x), y = 0,
    label = c("PBR / Nmin\n(current practice)", "FCS ideal /\ncarrying-capacity standard"),
    colour = c(controls$colour_dark, controls$colour_accent)
  )
  ggplot() +
    geom_segment(aes(x = 0, xend = 1, y = 0, yend = 0), colour = "#D7E4EC", linewidth = 10, lineend = "round") +
    geom_segment(aes(x = 0, xend = pbr_x, y = 0, yend = 0), colour = controls$colour_primary, linewidth = 10, lineend = "round") +
    geom_point(data = pts, aes(x = x, y = y, colour = colour), size = controls$point_size + 4) +
    scale_colour_identity() +
    geom_text(data = pts, aes(x = x, y = 0.55, label = label, colour = colour), size = controls$label_size, fontface = "bold", lineheight = 0.95) +
    annotate("text", x = 0, y = -0.5, label = "Extinction", hjust = 0, size = controls$label_size, colour = controls$colour_muted) +
    annotate("text", x = 1, y = -0.5, label = "Carrying capacity (K)", hjust = 1, size = controls$label_size, colour = controls$colour_muted) +
    scale_x_continuous(limits = c(-0.05, 1.05)) +
    scale_y_continuous(limits = c(-1, 1.2)) +
    deck_labs(controls, title = "Two reference points, one spectrum", subtitle = NULL, caption = NULL) +
    deck_theme(controls) +
    theme(panel.grid = element_blank(), axis.text = element_blank())
}

## =========================================================================
## Render every figure the deck needs, at the settings in `controls`.
## Returns a named list of file paths (passed straight to the Rmd).
## =========================================================================
build_workshop_deck_figures <- function(fig_dir, controls = deck_plot_controls,
                                          turbines = deck_turbines, reduction = deck_reduction_pct,
                                          elasticity = deck_elasticity, nmin_range = deck_nmin_range,
                                          cc_comparison = deck_cc_comparison, thresholds = deck_thresholds) {
  dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)
  fig_dir <- normalizePath(fig_dir, winslash = "/")
  controls$cut_in_speed_ms <- deck_cut_in_speed_ms

  list(
    fig_turbine_curtailment = save_deck_figure(plot_turbine_curtailment(turbines, controls), file.path(fig_dir, "deck_turbine_curtailment.png"), controls),
    fig_reduction_results   = save_deck_figure(plot_reduction_results(reduction, controls), file.path(fig_dir, "deck_reduction_results.png"), controls),
    fig_lambda_elasticity   = save_deck_figure(plot_lambda_elasticity(elasticity, controls), file.path(fig_dir, "deck_lambda_elasticity.png"), controls),
    fig_nmin_range          = save_deck_figure(plot_nmin_range(nmin_range, controls), file.path(fig_dir, "deck_nmin_range.png"), controls),
    fig_schaefer_curve      = save_deck_figure(plot_schaefer_curve(controls), file.path(fig_dir, "deck_schaefer_curve.png"), controls),
    fig_cc_vs_pbr           = save_deck_figure(plot_cc_vs_pbr(cc_comparison, thresholds, controls), file.path(fig_dir, "deck_cc_vs_pbr.png"), controls),
    fig_fcs_spectrum        = save_deck_figure(plot_fcs_spectrum(controls), file.path(fig_dir, "deck_fcs_spectrum.png"), controls, height = controls$fig_height * 0.75)
  )
}
