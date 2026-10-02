##
## Illustrative next-period turbine classification for the Paris Expert
## Workshop scenario -- mirrors the real report's evidence hierarchy
## (Core / Expansion candidate / Not flagged / Watch-list), reusing the
## Response period's post-curtailment records as the "since curtailment
## began" evidence, exactly as the real report reuses its post-5-May 2026
## data for the same purpose (no separate future period is simulated).
##

suppressPackageStartupMessages({ library(dplyr); library(ggplot2) })

run_paris_turbine_classification <- function(fig_dir = "outputs/figures", sim,
                                              expansion_mortality_fraction = 0.50,
                                              watchlist_n = 3) {
  dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)
  fig_dir <- normalizePath(fig_dir, winslash = "/")

  post_response <- sim$weekly %>%
    filter(period == "Response", week >= sim$response_start_week) %>%
    group_by(turbine) %>%
    summarise(raw_post = sum(raw_count), .groups = "drop")

  boundary_value <- min(sim$turbine_totals_baseline$total[sim$turbine_totals_baseline$turbine %in% sim$curtailed_turbines])

  tbl <- sim$turbine_totals_baseline %>%
    left_join(post_response, by = "turbine") %>%
    mutate(
      raw_post = tidyr::replace_na(raw_post, 0),
      curtailed = turbine %in% sim$curtailed_turbines,
      group = case_when(
        curtailed ~ "Core (curtailed)",
        !curtailed & total >= expansion_mortality_fraction * boundary_value & raw_post > 0 ~ "Expansion candidate",
        TRUE ~ "Not flagged"
      ),
      core_status = case_when(
        !curtailed ~ NA_character_,
        raw_post > 0 ~ "Residual fatalities since response start -- monitor/consider intensifying",
        TRUE ~ "No residual fatalities observed -- maintain current prescription"
      )
    )

  weakest_core <- tbl %>% filter(group == "Core (curtailed)") %>% arrange(desc(rank)) %>% head(watchlist_n) %>% pull(turbine)
  strongest_candidates <- tbl %>% filter(group == "Expansion candidate") %>% arrange(desc(total)) %>% head(watchlist_n) %>% pull(turbine)

  tbl <- tbl %>% mutate(
    watchlist = turbine %in% c(weakest_core, strongest_candidates),
    group_display = ifelse(watchlist & group == "Core (curtailed)", "Core (curtailed) -- watch-list", group)
  )

  summary_tbl <- tbl %>%
    count(group, name = "n_turbines") %>%
    tidyr::pivot_wider(names_from = group, values_from = n_turbines, values_fill = 0)

  core_residual_n <- sum(tbl$curtailed & tbl$raw_post > 0)
  candidates <- tbl %>% filter(group == "Expansion candidate") %>% arrange(rank) %>% pull(turbine)

  closest_miss <- tbl %>%
    filter(!curtailed, group == "Not flagged", raw_post > 0) %>%
    slice_max(total, n = 1, with_ties = FALSE)
  weak_pct <- if (nrow(closest_miss) == 0) NA_real_ else round(100 * closest_miss$total / boundary_value)

  fig_classification <- file.path(fig_dir, "paris_turbine_classification.png")
  ggsave(fig_classification, width = 8, height = 5.5, dpi = 150, bg = "white", plot = {
    plot_dt <- tbl %>% mutate(group_plot = factor(
      group_display, levels = c("Core (curtailed)", "Core (curtailed) -- watch-list", "Expansion candidate", "Not flagged")
    ))
    ggplot(plot_dt, aes(x = rank, y = raw_post, colour = group_plot)) +
      geom_point(size = 2.5, alpha = 0.85) +
      scale_colour_manual(name = "Next-period group", values = c(
        "Core (curtailed)" = "grey30", "Core (curtailed) -- watch-list" = "firebrick",
        "Expansion candidate" = "forestgreen", "Not flagged" = "grey80"
      )) +
      labs(
        x = "Baseline rank (lower = higher untreated risk)",
        y = "Raw records since response start",
        title = "Next-period turbine classification: Baseline risk vs. recurrence evidence",
        subtitle = paste0(
          "Expansion candidate: Baseline mortality >= ", round(100 * expansion_mortality_fraction),
          "% of boundary value and confirmed recurrence.\n",
          "Watch-list: marginal Core turbines paired with strongest candidates, for monitoring only."
        )
      ) +
      theme_minimal() +
      theme(plot.subtitle = element_text(size = 9.5))
  })

  list(
    turbine_classification = tbl,
    turbine_classification_summary = summary_tbl,
    turbine_classification_core_residual_n = core_residual_n,
    turbine_classification_candidates = candidates,
    turbine_classification_weak_pct = weak_pct,
    turbine_classification_boundary_value = boundary_value,
    fig_turbine_classification = fig_classification
  )
}
