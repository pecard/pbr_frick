##
## Illustrative turbine-selection validation for the Paris Expert Workshop
## scenario -- mirrors the real report's Section 7.3 (does a raw-data
## ranking recover the actual selection?), but reframed for a synthetic
## scenario where the "official" selection has no independent real-world
## source to validate against: here we validate the OBSERVED-mortality
## ranking (what the curtailment decision is actually based on) against
## the TRUE underlying risk ranking (the generative risk_weight, known
## only because this is a simulation) -- the same methodological point
## (does the method recover the right turbines from noisy observed data),
## demonstrated on ground truth instead of an operational history.
##

suppressPackageStartupMessages({ library(dplyr); library(ggplot2) })

run_paris_turbine_validation <- function(fig_dir = "outputs/figures", sim) {
  dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)
  fig_dir <- normalizePath(fig_dir, winslash = "/")

  comparison <- sim$turbine_totals_baseline %>%
    mutate(observed_rank = rank) %>%
    arrange(desc(risk_weight)) %>%
    mutate(true_rank = row_number()) %>%
    select(turbine, true_rank, observed_rank, total, risk_weight)

  n_selected <- sim$n_curtailed
  true_selected <- comparison$turbine[comparison$true_rank <= n_selected]
  observed_selected <- sim$curtailed_turbines
  overlap_n <- length(intersect(true_selected, observed_selected))

  fig_validation <- file.path(fig_dir, "paris_turbine_validation.png")
  ggsave(fig_validation, width = 7.5, height = 6, dpi = 150, bg = "white", plot = {
    ggplot(comparison, aes(x = true_rank, y = observed_rank)) +
      geom_abline(slope = 1, intercept = 0, linetype = "dotted", colour = "grey50") +
      geom_vline(xintercept = n_selected + 0.5, linetype = "dotted", colour = "steelblue") +
      geom_hline(yintercept = n_selected + 0.5, linetype = "dotted", colour = "steelblue") +
      geom_point(aes(colour = turbine %in% observed_selected), size = 2, alpha = 0.8) +
      scale_colour_manual(values = c(`TRUE` = "firebrick", `FALSE` = "grey60"), name = "Curtailed\n(observed)") +
      labs(
        x = "True rank (underlying risk, known only in simulation)",
        y = "Observed rank (Baseline mortality -- what curtailment is based on)",
        title = "Does observed mortality recover the true risk ranking?",
        subtitle = paste0(
          overlap_n, " of ", n_selected, " turbines selected on observed mortality ",
          "are also in the true top ", n_selected, " by underlying risk."
        )
      ) +
      theme_minimal() +
      theme(plot.subtitle = element_text(size = 10))
  })

  list(
    turbine_validation_comparison = comparison,
    turbine_validation_overlap_n = overlap_n,
    turbine_validation_overlap_pct = round(100 * overlap_n / n_selected),
    fig_turbine_validation = fig_validation
  )
}
