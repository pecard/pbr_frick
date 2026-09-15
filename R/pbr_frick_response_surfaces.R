# ============================================================
# Rebuilds the "PBR response surfaces across s, alpha and Fr" figure
# (Figure 4 of the earlier Word report) with the 4 Fr scenarios clearly
# identified -- the original plotly version stacked 4 surfaces, each shaded
# with its own PBR-height colour scale, with no legend distinguishing which
# surface belonged to which Fr. This produces two versions:
#  (1) a static, faceted 2D heatmap (one panel per Fr) -- unambiguous by
#      construction, report-ready as a static image;
#  (2) an interactive 3D plotly figure with one solid hue per Fr surface
#      and a proper legend, closer to the original chart type.
# Both use Nmin = 4000 (see pbr_frick_bsh_dgy_thresholds.R for why) and
# lambda_max derived from (s, alpha) via the Niel & Lebreton approximation,
# with the paper's fixed lambda_max = 1.20/1.24 benchmarks and the current
# BSH/DGY thresholds (144/120) overlaid as reference contours. The alpha
# range (1.5-3.5 years) follows the literature review in
# references/alpha_first_breeding_evidence.md.
# ============================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(plotly)
})

lambda_max_niel <- function(s, alpha) {
  disc <- (s - s * alpha - alpha - 1)^2 - 4 * s * alpha^2
  ((s * alpha - s + alpha + 1) + sqrt(pmax(disc, 0))) / (2 * alpha)
}

pbr_from_components <- function(Nmin, Fr, lambda_max) {
  Rmax <- lambda_max - 1
  0.5 * Rmax * Fr * Nmin
}

Nmin_assumed <- 4000
Fr_scenarios <- c(0.1, 0.3, 0.5, 1.0)  # Vulnerable, Near Threatened,
                                        # Least Concern (declining/unknown),
                                        # Least Concern (stable/increasing)

grid <- tidyr::expand_grid(
  s     = seq(0.70, 0.95, length.out = 60),
  alpha = seq(1.5, 3.5, length.out = 60),
  Fr    = Fr_scenarios
) %>%
  mutate(
    lambda_max = lambda_max_niel(s, alpha),
    PBR = pbr_from_components(Nmin_assumed, Fr, lambda_max),
    Fr_label = factor(
      paste0("Fr = ", Fr),
      levels = paste0("Fr = ", Fr_scenarios)
    )
  )

# ---- 1. Static faceted heatmap (one panel per Fr) -------------------------
p_facets <- ggplot(grid, aes(x = s, y = alpha, fill = PBR)) +
  geom_raster(interpolate = TRUE) +
  geom_contour(aes(z = PBR), breaks = c(120, 144), colour = "white", linewidth = 0.35) +
  facet_wrap(~Fr_label, nrow = 1) +
  scale_fill_viridis_c(name = "PBR\n(bats/year)", option = "C") +
  labs(
    title = "PBR response surfaces across adult survival (s), age at first breeding (alpha) and Fr",
    subtitle = paste0(
      "Nmin = ", format(Nmin_assumed, big.mark = ","),
      "; lambda_max derived from s and alpha. White contours: PBR = 120 (DGY) and 144 (BSH), current imposed thresholds"
    ),
    x = "Adult survival (s)", y = "Age at first breeding (alpha, years)"
  ) +
  theme_minimal() +
  theme(legend.position = "right", strip.text = element_text(face = "bold"))

ggsave("figures/pbr_response_surfaces_facets.png", p_facets, width = 13, height = 4, dpi = 150)

# ---- 2. Interactive 3D plotly with a legend distinguishing the 4 surfaces -
# plotly does not give add_surface() traces a legend entry by default, and
# colouring every surface with the same PBR-height colour scale (as in the
# original figure) makes them visually indistinguishable where they
# overlap. Fix: one solid hue per Fr (encoded via a single-colour
# colourscale using surface height only for shading within that hue), plus
# an invisible scatter3d marker per surface to force a legend entry.
hues <- c("Blues", "Greens", "Oranges", "Purples")  # one per Fr scenario, low->high

fig <- plot_ly()
for (i in seq_along(Fr_scenarios)) {
  fr_val <- Fr_scenarios[i]
  sub <- grid %>% filter(Fr == fr_val)
  z_mat <- matrix(sub$PBR, nrow = length(unique(sub$s)), byrow = FALSE)
  s_vals <- sort(unique(sub$s))
  a_vals <- sort(unique(sub$alpha))

  fig <- fig %>%
    add_surface(
      x = a_vals, y = s_vals, z = z_mat,
      colorscale = list(c(0, 1), c("white", RColorBrewer::brewer.pal(9, hues[i])[7])),
      showscale = FALSE,
      opacity = 0.85,
      name = paste0("Fr = ", fr_val),
      hovertemplate = paste0(
        "Fr = ", fr_val,
        "<br>s: %{y:.2f}<br>alpha: %{x:.2f}<br>PBR: %{z:.0f} bats/yr<extra></extra>"
      )
    ) %>%
    add_trace(
      type = "scatter3d", mode = "markers",
      x = a_vals[1], y = s_vals[1], z = max(sub$PBR),
      marker = list(size = 8, color = RColorBrewer::brewer.pal(9, hues[i])[7]),
      name = paste0("Fr = ", fr_val),
      showlegend = TRUE
    )
}

fig <- fig %>%
  layout(
    title = list(text = paste0(
      "PBR response surfaces across s, alpha and Fr (Nmin = ",
      format(Nmin_assumed, big.mark = ","), ")"
    )),
    scene = list(
      xaxis = list(title = "Age at first breeding (alpha)"),
      yaxis = list(title = "Adult survival (s)"),
      zaxis = list(title = "PBR (bats/year)")
    ),
    legend = list(title = list(text = "Recovery factor scenario"))
  )

htmlwidgets::saveWidget(fig, "figures/pbr_response_surfaces_3d.html", selfcontained = TRUE)

cat("Wrote figures/pbr_response_surfaces_facets.png and figures/pbr_response_surfaces_3d.html\n")
