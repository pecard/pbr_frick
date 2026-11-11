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

suppressPackageStartupMessages({ library(ggplot2); library(dplyr); library(tidyr); library(scales); library(patchwork) })

## This container's default locale is "C" (not UTF-8), under which the
## graphics device silently renders Greek glyphs (lambda, alpha) as ".."
## in every ggplot-saved PNG, even though the fonts installed (DejaVu
## Sans) have full Greek coverage and the Rmd's own text is unaffected --
## confirmed by direct comparison (Sys.getlocale("LC_CTYPE") "C" vs.
## "C.utf8"). Setting LC_CTYPE once, here, fixes it at the source instead
## of spelling every lambda/alpha out as a word in figure text.
if (Sys.getlocale("LC_CTYPE") %in% c("C", "POSIX")) {
  tryCatch(Sys.setlocale("LC_CTYPE", "C.utf8"), error = function(e) {
    warning("Could not set LC_CTYPE to C.utf8 -- Greek glyphs (λ, α) in figures may render incorrectly: ", conditionMessage(e))
  })
}

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
      title = "Elasticity of λmax",
      subtitle = "% change per 1% parameter change -- low (sky) to high (teal) end of the explored range",
      caption = "s in [0.70, 0.95]; α in [1.5, 3.5] years (Niel & Lebreton 2005 sensitivity analysis)."
    ) +
    deck_theme(controls)
}

## =========================================================================
## Slide 8 -- implied Nmin range, by lambda_max benchmark
## =========================================================================
plot_nmin_range <- function(nmin_range, controls) {
  dt <- nmin_range %>%
    mutate(lambda_label = paste0("λmax = ", lambda_max)) %>%
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
      title = "Implied local abundance (Nmin), by λmax benchmark",
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
    annotate("text", x = msy_x, y = max(dt$H) * 1.04, label = "α = 0.5 (MSY)", size = controls$label_size, colour = controls$colour_muted) +
    annotate("text", x = 0.92, y = max(dt$H) * 0.12, label = "α near 1\n(near K)", size = controls$label_size, colour = controls$colour_muted) +
    scale_x_continuous(labels = NULL) +
    scale_y_continuous(labels = NULL) +
    deck_labs(controls,
      title = "Sustainable removal H(α) = r·K·α·(1−α)",
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
      caption = "Age-structured, stochastic simulation at an α=0.8 target; values floored."
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
## Chapter 7 rationale (Paulo, 2026-11): replaces the old turbine-count
## and reduction-only slides with the fuller operational story --
## Sections 7.1-7.4 of the real technical note.
## =========================================================================

## ---- New slide A -- schematic weekly control chart (Section 7.1's
## mechanism). No real weekly PCFM series is available in this
## environment (gitignored, local-only raw data), so this is explicitly
## a SCHEMATIC of the mechanism -- a generic seasonal cumulative-mortality
## shape with illustrative control limits -- not this project's own real
## weekly trajectory. Flagged as such in the figure itself, consistent
## with this deck's other conceptual (not fitted) diagrams.
plot_control_chart_schematic <- function(controls, n = 200) {
  t <- seq(0, 52, length.out = n)
  seasonal_shape <- function(x) pmax(0, sin(pmax(0, x - 14) / 10 * pi)) * exp(-((x - 24)^2) / 220) +
    pmax(0, sin(pmax(0, x - 34) / 8 * pi)) * exp(-((x - 40)^2) / 90)
  raw <- seasonal_shape(t)
  expected_cum <- cumsum(raw) / sum(raw) * 100
  upper_amber <- pmin(100, expected_cum * 1.15 + 3)
  upper_red   <- pmin(115, expected_cum * 1.32 + 6)
  actual <- pmin(118, expected_cum * 1.30 + ifelse(t > 10, 8, 0) * pmax(0, (t - 10) / 15))

  dt <- tibble::tibble(week = t, expected = expected_cum, amber = upper_amber, red = upper_red, actual = actual)

  ggplot(dt, aes(x = week)) +
    geom_ribbon(aes(ymin = 0, ymax = expected), fill = "#D7E4EC", alpha = 0.6) +
    geom_ribbon(aes(ymin = expected, ymax = amber), fill = controls$colour_accent, alpha = 0.25) +
    geom_ribbon(aes(ymin = amber, ymax = red), fill = "firebrick", alpha = 0.18) +
    geom_line(aes(y = expected), colour = controls$colour_muted, linewidth = 0.8, linetype = "dashed") +
    geom_line(aes(y = actual), colour = controls$colour_dark, linewidth = 1.2) +
    geom_vline(xintercept = 18, linetype = "dotted", colour = controls$colour_primary) +
    annotate("text", x = 18.5, y = 5, label = "Curtailment\nstarts", hjust = 0, size = controls$label_size, colour = controls$colour_primary) +
    annotate("text", x = 46, y = 18, label = "GREEN\n(on pace)", size = controls$label_size * 0.95, colour = controls$colour_muted) +
    annotate("text", x = 46, y = 45, label = "AMBER\n(intensify)", size = controls$label_size * 0.95, colour = "#8a6410") +
    annotate("text", x = 46, y = 80, label = "RED\n(escalate)", size = controls$label_size * 0.95, colour = "firebrick") +
    scale_x_continuous(limits = c(0, 52), expand = c(0, 0)) +
    scale_y_continuous(limits = c(0, 118), expand = c(0, 0), labels = NULL) +
    deck_labs(controls,
      title = "A Poisson-based weekly control chart, not an annual pass/fail line",
      subtitle = "Schematic of the mechanism -- not this project's own real weekly series",
      caption = "Dark line: illustrative observed cumulative mortality. Dashed: expected seasonal pace."
    ) +
    deck_theme(controls) +
    theme(panel.grid = element_blank(), axis.title.y = element_blank())
}

## ---- New slide B -- 2026 mortality position, two questions (Section
## 7.2's two tables): year-to-date vs. the annual reference, and
## year-over-year through the same ISO week. Two-panel patchwork so both
## questions sit on one slide, side by side, exactly as the report poses
## them together to resolve the apparent inconsistency between the two readings.
plot_ytd_vs_yoy <- function(ytd, yoy, controls) {
  p1 <- ytd %>%
    select(site, `Annual reference` = annual_reference, `GenEst-corrected (2026 YTD)` = genest_corrected) %>%
    pivot_longer(-site, names_to = "series", values_to = "value") %>%
    mutate(series = factor(series, levels = c("Annual reference", "GenEst-corrected (2026 YTD)"))) %>%
    ggplot(aes(x = site, y = value, fill = series)) +
    geom_col(position = position_dodge(width = controls$bar_width + 0.08), width = controls$bar_width) +
    geom_text(aes(label = comma(value)), position = position_dodge(width = controls$bar_width + 0.08),
               vjust = -0.4, size = controls$label_size, colour = controls$colour_dark) +
    scale_fill_manual(values = c("Annual reference" = controls$colour_dark, "GenEst-corrected (2026 YTD)" = controls$colour_primary)) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.18))) +
    labs(title = "Year-to-date vs. annual reference", subtitle = "Through ISO week 39", x = NULL, y = NULL, fill = NULL) +
    deck_theme(controls) +
    theme(axis.text.y = element_blank(), panel.grid.major.y = element_blank(), legend.position = "bottom",
          plot.title = element_text(size = controls$base_size + 1, face = "bold", colour = controls$colour_dark),
          plot.subtitle = element_text(size = controls$base_size * 0.75, colour = controls$colour_muted))

  p2 <- yoy %>%
    select(site, `2025` = y2025, `2026` = y2026) %>%
    pivot_longer(-site, names_to = "year", values_to = "value") %>%
    ggplot(aes(x = site, y = value, fill = year)) +
    geom_col(position = position_dodge(width = controls$bar_width + 0.08), width = controls$bar_width) +
    geom_text(aes(label = comma(value)), position = position_dodge(width = controls$bar_width + 0.08),
               vjust = -0.4, size = controls$label_size, colour = controls$colour_dark) +
    geom_text(data = yoy, aes(x = site, y = pmax(y2025, y2026) * 1.18, label = paste0(pct_change, "%")),
               inherit.aes = FALSE, fontface = "bold", size = controls$label_size * 1.2, colour = controls$colour_accent) +
    scale_fill_manual(values = c("2025" = controls$colour_secondary, "2026" = controls$colour_primary)) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.2))) +
    labs(title = "Year-over-year, same weeks", subtitle = "GenEst-corrected mortality", x = NULL, y = NULL, fill = NULL) +
    deck_theme(controls) +
    theme(axis.text.y = element_blank(), panel.grid.major.y = element_blank(), legend.position = "bottom",
          plot.title = element_text(size = controls$base_size + 1, face = "bold", colour = controls$colour_dark),
          plot.subtitle = element_text(size = controls$base_size * 0.75, colour = controls$colour_muted))

  (p1 | p2)
}

## ---- New slide C -- turbine-level retrospective validation (Section
## 7.3): does a simple raw-carcass ranking, corrected with this note's own
## working factor, recover the real/official selection -- and why that
## check is what drove using two project-specific GenEst correction
## factors instead of one shared/flat one.
plot_turbine_validation_summary <- function(validation, controls) {
  validation %>%
    select(site, `% recovered` = pct_official_recovered,
           `% mortality captured` = pct_mortality_captured) %>%
    pivot_longer(-site, names_to = "metric", values_to = "value") %>%
    mutate(metric = factor(metric, levels = c("% recovered", "% mortality captured"))) %>%
    ggplot(aes(x = site, y = value, fill = metric)) +
    geom_col(position = position_dodge(width = controls$bar_width + 0.08), width = controls$bar_width) +
    geom_text(aes(label = paste0(value, "%")), position = position_dodge(width = controls$bar_width + 0.08),
               vjust = -0.4, size = controls$label_size, fontface = "bold", colour = controls$colour_dark) +
    scale_fill_manual(values = c("% recovered" = controls$colour_secondary,
                                  "% mortality captured" = controls$colour_primary)) +
    scale_y_continuous(limits = c(0, 100), expand = expansion(mult = c(0, 0.12))) +
    labs(title = "Turbine selection",
         subtitle = "Independent ranking vs. official selection",
         x = NULL, y = NULL, fill = NULL) +
    deck_theme(controls) +
    theme(axis.text.y = element_blank(), panel.grid.major.y = element_blank(), legend.position = "bottom",
          legend.text = element_text(size = controls$base_size * 0.72))
}

plot_genest_ratio <- function(ratio, rejected_flat, combined_median, controls) {
  dt <- ratio %>% mutate(site = factor(site, levels = rev(site)))
  ggplot(dt) +
    geom_segment(aes(x = min_per_turbine, xend = max_per_turbine, y = site, yend = site),
                 colour = "grey70", linewidth = controls$line_width * 1.3) +
    geom_point(aes(x = min_per_turbine, y = site), size = controls$point_size, colour = controls$colour_secondary) +
    geom_point(aes(x = max_per_turbine, y = site), size = controls$point_size, colour = controls$colour_secondary) +
    geom_point(aes(x = aggregate_ratio, y = site), size = controls$point_size + 1.5, colour = controls$colour_primary) +
    geom_text(aes(x = aggregate_ratio, y = site, label = paste0(aggregate_ratio, "x")),
              vjust = -1.3, size = controls$label_size, fontface = "bold", colour = controls$colour_primary) +
    geom_vline(xintercept = rejected_flat, linetype = "dashed", colour = "firebrick") +
    geom_vline(xintercept = combined_median, linetype = "dotted", colour = controls$colour_muted) +
    scale_x_continuous(limits = c(0, max(dt$max_per_turbine) * 1.12),
                        expand = expansion(mult = c(0, 0))) +
    scale_y_discrete(expand = expansion(add = c(0.6, 0.9))) +
    coord_cartesian(clip = "off") +
    labs(title = "GenEst correction ratio",
         subtitle = paste0("Per-turbine range (grey) vs. project aggregate (teal)\n",
                            "Dashed = rejected flat ", rejected_flat, "x -- Dotted = combined median ", combined_median, "x"),
         x = "Ratio (estimate / raw carcasses)", y = NULL) +
    deck_theme(controls)
}

## Combined view for the deck: validation summary (left) + GenEst ratio
## (right) on one slide, same patchwork convention as plot_ytd_vs_yoy().
plot_turbine_validation_combined <- function(validation, ratio, rejected_flat, combined_median, controls) {
  p1 <- plot_turbine_validation_summary(validation, controls) +
    theme(plot.title = element_text(size = controls$base_size + 1, face = "bold", colour = controls$colour_dark),
          plot.subtitle = element_text(size = controls$base_size * 0.72, colour = controls$colour_muted))
  p2 <- plot_genest_ratio(ratio, rejected_flat, combined_median, controls) +
    theme(plot.title = element_text(size = controls$base_size + 1, face = "bold", colour = controls$colour_dark),
          plot.subtitle = element_text(size = controls$base_size * 0.72, colour = controls$colour_muted))

  (p1 | p2)
}

## =========================================================================
## Render every figure the deck needs, at the settings in `controls`.
## Returns a named list of file paths (passed straight to the Rmd).
## =========================================================================
build_workshop_deck_figures <- function(fig_dir, controls = deck_plot_controls,
                                          turbines = deck_turbines, reduction = deck_reduction_pct,
                                          elasticity = deck_elasticity, nmin_range = deck_nmin_range,
                                          cc_comparison = deck_cc_comparison, thresholds = deck_thresholds,
                                          ytd = deck_ytd_status, yoy = deck_yoy,
                                          turbine_validation = deck_turbine_validation, genest_ratio = deck_genest_ratio) {
  dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)
  fig_dir <- normalizePath(fig_dir, winslash = "/")
  controls$cut_in_speed_ms <- deck_cut_in_speed_ms

  list(
    fig_reduction_results   = save_deck_figure(plot_reduction_results(reduction, controls), file.path(fig_dir, "deck_reduction_results.png"), controls),
    fig_control_chart       = save_deck_figure(plot_control_chart_schematic(controls), file.path(fig_dir, "deck_control_chart.png"), controls),
    fig_ytd_vs_yoy          = save_deck_figure(plot_ytd_vs_yoy(ytd, yoy, controls), file.path(fig_dir, "deck_ytd_vs_yoy.png"), controls,
                                                height = controls$fig_height * 0.9),
    fig_turbine_validation  = save_deck_figure(plot_turbine_validation_combined(turbine_validation, genest_ratio, deck_genest_rejected_flat, deck_genest_combined_median, controls),
                                                file.path(fig_dir, "deck_turbine_validation.png"), controls, height = controls$fig_height * 0.9),
    fig_lambda_elasticity   = save_deck_figure(plot_lambda_elasticity(elasticity, controls), file.path(fig_dir, "deck_lambda_elasticity.png"), controls),
    fig_nmin_range          = save_deck_figure(plot_nmin_range(nmin_range, controls), file.path(fig_dir, "deck_nmin_range.png"), controls),
    fig_schaefer_curve      = save_deck_figure(plot_schaefer_curve(controls), file.path(fig_dir, "deck_schaefer_curve.png"), controls),
    fig_cc_vs_pbr           = save_deck_figure(plot_cc_vs_pbr(cc_comparison, thresholds, controls), file.path(fig_dir, "deck_cc_vs_pbr.png"), controls),
    fig_fcs_spectrum        = save_deck_figure(plot_fcs_spectrum(controls), file.path(fig_dir, "deck_fcs_spectrum.png"), controls, height = controls$fig_height * 0.75)
  )
}
