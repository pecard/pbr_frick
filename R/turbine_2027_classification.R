##
## 2027 turbine classification (Paulo, 2026-09): implements the evidence
## hierarchy from his post-Version-13 review -- NOT a re-ranked "top N"
## list, but three operational groups built on an explicit rule: 2025
## (the only complete untreated year) stays the risk baseline; 2026 is
## used to CORROBORATE candidates for expansion, never to justify
## removing a currently-curtailed turbine (low mortality at an already-
## curtailed turbine is evidence the intervention may be working, not
## evidence the turbine is now low-risk); any reallocation is a watch-
## list for 2027, not a prescriptive swap, given the short post-
## treatment window (~19 weeks) still available.
##
## Three groups:
##   Core (curtailed): the real, already-selected 17/7 turbines. Split by
##     whether they show ANY post-5-May-2026 V. murinus record: turbines
##     with residual fatalities are flagged for monitoring/possible
##     intensification, not treated as evidence curtailment failed there
##     or that it could be relaxed.
##   Expansion candidates: NOT currently curtailed, but corroborated by
##     BOTH years independently -- 2025 estimated mortality at least half
##     of the lowest currently-curtailed turbine's (i.e. comparable to
##     the turbines that WERE selected, not just "not zero"), AND at
##     least one confirmed V. murinus record after curtailment started in
##     2026 (their mortality is not suppressed by treatment, since they
##     were never curtailed -- genuine untreated recurrence evidence,
##     not noise from a single low-count turbine).
##   Review/watch-list: the weakest-evidence currently-curtailed turbines
##     (closest to the original 50% boundary), paired conceptually with
##     the strongest expansion candidates -- flagged for closer 2027
##     monitoring and joint re-assessment once a full season of evidence
##     accumulates, explicitly NOT a proposed substitution now.
##
## All turbine identifiers are anonymised by official 2025 rank (T-NN),
## consistent with Section 7.3 -- never the real code.
##

suppressPackageStartupMessages({ library(dplyr); library(ggplot2) })

run_turbine_2027_classification <- function(fig_dir = "outputs/figures",
                                             bash_path = "data-raw/BashWPP_Weekly_PCFM_PBR.xlsx",
                                             djangeldy_path = "data-raw/DjangeldyWPP_Weekly_PCFM_PBR.xlsx",
                                             official_bash_path = "data-raw/official_turbine_selection_bash.csv",
                                             official_djangeldy_path = "data-raw/official_turbine_selection_djangeldy.csv",
                                             expansion_mortality_fraction = 0.50,
                                             watchlist_n = 3) {

  dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)
  fig_dir <- normalizePath(fig_dir, winslash = "/")

  d <- load_real_mortality_data(bash_path, djangeldy_path)

  pad_bsh <- function(x) sprintf("BSH%02d", as.integer(sub("^BSH", "", x)))
  pad_dzh <- function(x) x
  official_bash <- read.csv(official_bash_path, stringsAsFactors = FALSE) %>%
    mutate(official_selected = as.logical(official_selected)) %>% rename(turbine = turbine_real)
  official_dzh <- read.csv(official_djangeldy_path, stringsAsFactors = FALSE) %>%
    mutate(official_selected = as.logical(official_selected)) %>% rename(turbine = turbine_real)
  d <- suppressWarnings(d %>% mutate(turbine = ifelse(project == facility_labels[1], pad_bsh(turbine), pad_dzh(turbine))))

  classify_project <- function(official, proj_label) {
    post_may <- d %>% filter(project == proj_label, date >= curtailment_start_date) %>%
      count(turbine, name = "raw_post_may")
    boundary_value <- min(official$official_estimated_mortality[official$official_selected])

    tbl <- official %>%
      left_join(post_may, by = "turbine") %>%
      mutate(
        raw_post_may = tidyr::replace_na(raw_post_may, 0),
        project = proj_label,
        turbine_anon = sprintf("T-%02d", official_rank),
        group = case_when(
          official_selected ~ "Core (curtailed)",
          !official_selected & official_estimated_mortality >= expansion_mortality_fraction * boundary_value &
            raw_post_may > 0 ~ "Expansion candidate",
          TRUE ~ "Not flagged"
        ),
        core_status = case_when(
          !official_selected ~ NA_character_,
          raw_post_may > 0 ~ "Residual fatalities post-5-May -- monitor/consider intensifying",
          TRUE ~ "No residual fatalities observed -- maintain current prescription"
        )
      )

    ## Review/watch-list: the N weakest-evidence Core turbines (closest to
    ## the original boundary, i.e. highest official_rank among curtailed)
    ## paired with the N strongest Expansion candidates (by 2025 estimated
    ## mortality) -- flagged, not swapped.
    weakest_core <- tbl %>% filter(group == "Core (curtailed)") %>% arrange(desc(official_rank)) %>% head(watchlist_n) %>% pull(turbine_anon)
    strongest_candidates <- tbl %>% filter(group == "Expansion candidate") %>%
      arrange(desc(official_estimated_mortality)) %>% head(watchlist_n) %>% pull(turbine_anon)

    tbl %>% mutate(
      watchlist = turbine_anon %in% c(weakest_core, strongest_candidates),
      group = ifelse(watchlist & group == "Core (curtailed)", "Core (curtailed) -- watch-list (marginal)", group)
    )
  }

  tbl_bash <- classify_project(official_bash, facility_labels[1])
  tbl_dzh <- classify_project(official_dzh, facility_labels[2])
  classification_all <- bind_rows(tbl_bash, tbl_dzh)

  summary_by_group <- classification_all %>%
    mutate(group_summary = ifelse(grepl("^Core", group), "Core (curtailed)", group)) %>%
    count(project, group_summary, name = "n_turbines") %>%
    tidyr::pivot_wider(names_from = group_summary, values_from = n_turbines, values_fill = 0)

  flagged <- classification_all %>%
    filter(group == "Expansion candidate" | watchlist) %>%
    arrange(project, desc(group), official_rank) %>%
    select(project, turbine_anon, official_rank, official_estimated_mortality, raw_post_may, group, core_status, watchlist)

  core_residual_summary <- classification_all %>%
    filter(official_selected) %>%
    group_by(project) %>%
    summarise(
      n_core = n(),
      n_with_residual = sum(raw_post_may > 0),
      .groups = "drop"
    )

  fig_2027_classification <- file.path(fig_dir, "turbine_2027_classification.png")
  ggsave(fig_2027_classification, width = 10.5, height = 5.5, dpi = 150, plot = {
    plot_dt <- classification_all %>%
      mutate(group_plot = factor(
        ifelse(grepl("watch-list", group), "Core (curtailed) -- watch-list", gsub(" -- watch-list \\(marginal\\)", "", group)),
        levels = c("Core (curtailed)", "Core (curtailed) -- watch-list", "Expansion candidate", "Not flagged")
      ))
    ggplot(plot_dt, aes(x = official_rank, y = raw_post_may, colour = group_plot)) +
      geom_point(size = 2, alpha = 0.85) +
      facet_wrap(~project, scales = "free_x") +
      scale_colour_manual(name = "2027 group", values = c(
        "Core (curtailed)" = "grey30", "Core (curtailed) -- watch-list" = "firebrick",
        "Expansion candidate" = "forestgreen", "Not flagged" = "grey80"
      )) +
      labs(
        x = "2025 official rank (lower = higher untreated risk)",
        y = "Raw V. murinus records, 5 May-current 2026",
        title = "2027 turbine classification: baseline risk (2025) vs. recurrence evidence (2026)",
        subtitle = paste0(
          "Expansion candidates: comparable 2025 risk to Core AND confirmed 2026 recurrence.\n",
          "Watch-list: marginal Core turbines paired with the strongest candidates, for monitoring only."
        )
      ) +
      theme_minimal() +
      theme(plot.subtitle = element_text(size = 9.5), strip.text = element_text(face = "bold"), legend.position = "bottom")
  })

  list(
    turbine_2027_classification = classification_all,
    turbine_2027_flagged = flagged,
    turbine_2027_summary_by_group = summary_by_group,
    turbine_2027_core_residual_summary = core_residual_summary,
    turbine_2027_expansion_fraction = expansion_mortality_fraction,
    turbine_2027_watchlist_n = watchlist_n,
    fig_turbine_2027_classification = fig_2027_classification
  )
}
