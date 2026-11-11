##
## Interactive PBR simulator (Paulo, 2026-11): sweep one PBR parameter at
## a time -- Nmin, Fr (via IUCN status), adult survival (s), or age at
## first breeding (alpha) -- holding the others fixed at slider values,
## and compare the resulting PBR against Frick et al.'s two FIXED
## precautionary lambda_max benchmarks (1.24 / 1.20) versus the REALISTIC
## lambda_max implied by s and alpha through the Niel & Lebreton (2005)
## demographic-invariant approximation.
##
## Reuses the project's own formulas (R/pbr_functions.R) and real
## threshold VALUES (inputs/pbrSettings_BSH_DGY.R) -- nothing here
## re-derives PBR, it only sweeps and plots what those functions already
## compute. Facility NAMES are this app's own fixed "Project 1"/"Project 2"
## anonymisation, independent of that settings file's facility_labels
## (see the note by project_thresholds below for why).
## Run with: shiny::runApp("shiny_pbr_simulator")
##

library(shiny)
suppressPackageStartupMessages({ library(ggplot2); library(dplyr); library(tidyr); library(scales) })

source("../R/pbr_functions.R")
source("../inputs/pbrSettings_BSH_DGY.R")

## ---- App-level controls -- the knobs for this app's own look/ranges,
## kept separate from pbrSettings_BSH_DGY.R's analysis parameters so
## restyling the simulator never means editing server logic below. -------
app_controls <- list(
  colour_1.24   = "#1C7293",  # teal
  colour_1.20   = "#6CA6C1",  # sky
  colour_realistic = "#E8A33D",  # amber
  colour_threshold = "#0B2545",  # navy
  base_size     = 14,
  s_sweep_range     = c(0.50, 0.99),
  alpha_sweep_range = c(1.0, 5.0),
  fr_sweep_range    = c(0.1, 1.0),
  nmin_sweep_range  = c(200, 8000),
  nmin_default      = 3000,
  surface_grid_res  = 120
)

iucn_choices <- c(
  "Vulnerable" = "Vulnerable",
  "Near Threatened" = "Near Threatened",
  "Least Concern (declining or unknown)" = "Least Concern (declining or unknown)",
  "Least Concern (stable or increasing)" = "Least Concern (stable or increasing)"
)

## Facility NAMES shown here are this app's own fixed anonymisation
## (Project 1 = Bash WPP = 112 bats/yr; Project 2 = Djangeldy WPP = 88
## bats/yr -- the same mapping documented in inputs/pbrSettings_BSH_DGY.R),
## kept independent of that shared file's own facility_labels, which
## Paulo toggles between real and anonymised names depending on which
## deliverable (the real technical note vs. an external/workshop context)
## is being built. Decoupling this way means the simulator never
## accidentally leaks whichever real/anonymised state that file happens
## to be in -- only the real threshold VALUES (not identity-revealing)
## are reused from it. Paulo, 2026-11: anonymise consistently everywhere.
project_thresholds <- tibble::tibble(
  facility  = c("Project 1", "Project 2"),
  threshold = pbr_thresholds$threshold
)

## ======================================================================
## UI
## ======================================================================
ui <- fluidPage(
  titlePanel("PBR simulator -- Frick benchmarks vs. the realistic Niel & Lebreton lambda_max"),
  sidebarLayout(
    sidebarPanel(
      width = 3,
      selectInput("vary_param", "Vary this parameter:",
                  choices = c("Nmin (local abundance)" = "nmin",
                              "Fr (recovery factor, via IUCN status)" = "fr",
                              "Adult survival (s)" = "s",
                              "Age at first breeding (alpha)" = "alpha"),
                  selected = "s"),
      hr(),
      strong("Held fixed (unless swept above):"),
      sliderInput("nmin_fixed", "Nmin", min = app_controls$nmin_sweep_range[1], max = app_controls$nmin_sweep_range[2],
                  value = app_controls$nmin_default, step = 50),
      selectInput("iucn_fixed", "IUCN status (sets Fr)", choices = iucn_choices, selected = "Near Threatened"),
      sliderInput("s_fixed", "Adult survival (s)", min = app_controls$s_sweep_range[1], max = app_controls$s_sweep_range[2],
                  value = leslie_s_adult_f, step = 0.01),
      sliderInput("alpha_fixed", "Age at first breeding (alpha, years)", min = app_controls$alpha_sweep_range[1], max = app_controls$alpha_sweep_range[2],
                  value = 2.0, step = 0.1),
      hr(),
      checkboxGroupInput("show_thresholds", "Show imposed thresholds:",
                          choices = setNames(project_thresholds$facility, project_thresholds$facility),
                          selected = project_thresholds$facility)
    ),
    mainPanel(
      width = 9,
      tabsetPanel(
        tabPanel("Sensitivity sweep",
          br(),
          fluidRow(
            column(4, wellPanel(strong("PBR at lambda_max = 1.24"), br(), textOutput("pbr_124", inline = TRUE))),
            column(4, wellPanel(strong("PBR at lambda_max = 1.20"), br(), textOutput("pbr_120", inline = TRUE))),
            column(4, wellPanel(strong("PBR at realistic lambda_max"), br(), textOutput("pbr_realistic", inline = TRUE)))
          ),
          plotOutput("sweep_plot", height = "500px")
        ),
        tabPanel("lambda_max surface (s x alpha)",
          br(),
          plotOutput("surface_plot", height = "550px")
        )
      )
    )
  )
)

## ======================================================================
## Server
## ======================================================================
server <- function(input, output, session) {

  fr_fixed <- reactive(fr_from_iucn(input$iucn_fixed))
  lambda_realistic_fixed <- reactive(lambda_max_niel(input$s_fixed, input$alpha_fixed))

  ## ---- Sweep data: one row per point along the x-axis of whichever
  ## parameter is being varied, three PBR series (1.24 fixed, 1.20 fixed,
  ## realistic), all other parameters held at their slider values. -------
  sweep_dt <- reactive({
    vp <- input$vary_param
    n <- 80

    if (vp == "nmin") {
      x <- seq(app_controls$nmin_sweep_range[1], app_controls$nmin_sweep_range[2], length.out = n)
      Fr <- fr_fixed(); lam_real <- lambda_realistic_fixed()
      tibble::tibble(
        x = x, x_label = "Nmin (local abundance)",
        pbr_124 = pbr_from_components(x, Fr, 1.24),
        pbr_120 = pbr_from_components(x, Fr, 1.20),
        pbr_realistic = pbr_from_components(x, Fr, lam_real)
      )
    } else if (vp == "fr") {
      x <- seq(app_controls$fr_sweep_range[1], app_controls$fr_sweep_range[2], length.out = n)
      Nmin <- input$nmin_fixed; lam_real <- lambda_realistic_fixed()
      tibble::tibble(
        x = x, x_label = "Fr (recovery factor)",
        pbr_124 = pbr_from_components(Nmin, x, 1.24),
        pbr_120 = pbr_from_components(Nmin, x, 1.20),
        pbr_realistic = pbr_from_components(Nmin, x, lam_real)
      )
    } else if (vp == "s") {
      x <- seq(app_controls$s_sweep_range[1], app_controls$s_sweep_range[2], length.out = n)
      Nmin <- input$nmin_fixed; Fr <- fr_fixed(); alpha_f <- input$alpha_fixed
      tibble::tibble(
        x = x, x_label = "Adult survival (s)",
        pbr_124 = pbr_from_components(Nmin, Fr, 1.24),
        pbr_120 = pbr_from_components(Nmin, Fr, 1.20),
        pbr_realistic = pbr_from_components(Nmin, Fr, lambda_max_niel(x, alpha_f))
      )
    } else {
      x <- seq(app_controls$alpha_sweep_range[1], app_controls$alpha_sweep_range[2], length.out = n)
      Nmin <- input$nmin_fixed; Fr <- fr_fixed(); s_f <- input$s_fixed
      tibble::tibble(
        x = x, x_label = "Age at first breeding (alpha, years)",
        pbr_124 = pbr_from_components(Nmin, Fr, 1.24),
        pbr_120 = pbr_from_components(Nmin, Fr, 1.20),
        pbr_realistic = pbr_from_components(Nmin, Fr, lambda_max_niel(s_f, x))
      )
    }
  })

  current_x <- reactive({
    switch(input$vary_param,
      nmin  = input$nmin_fixed,
      fr    = fr_fixed(),
      s     = input$s_fixed,
      alpha = input$alpha_fixed
    )
  })

  output$pbr_124 <- renderText(sprintf("%.0f bats/yr", pbr_from_components(input$nmin_fixed, fr_fixed(), 1.24)))
  output$pbr_120 <- renderText(sprintf("%.0f bats/yr", pbr_from_components(input$nmin_fixed, fr_fixed(), 1.20)))
  output$pbr_realistic <- renderText(
    sprintf("%.0f bats/yr (lambda_max = %.3f)", pbr_from_components(input$nmin_fixed, fr_fixed(), lambda_realistic_fixed()),
            lambda_realistic_fixed())
  )

  output$sweep_plot <- renderPlot({
    dt <- sweep_dt() %>%
      pivot_longer(c(pbr_124, pbr_120, pbr_realistic), names_to = "series", values_to = "PBR") %>%
      mutate(series = recode(series,
        pbr_124 = "Fixed lambda_max = 1.24 (Frick)",
        pbr_120 = "Fixed lambda_max = 1.20 (Frick)",
        pbr_realistic = "Realistic lambda_max (Niel & Lebreton, current s/alpha)"
      ))

    thr <- project_thresholds %>% filter(facility %in% input$show_thresholds)

    ggplot(dt, aes(x = x, y = PBR, colour = series)) +
      geom_line(linewidth = 1.1) +
      { if (nrow(thr) > 0) geom_hline(data = thr, aes(yintercept = threshold), linetype = "dotted", colour = app_controls$colour_threshold) } +
      { if (nrow(thr) > 0) geom_text(data = thr, aes(x = min(dt$x), y = threshold, label = paste0(facility, " = ", threshold)),
                                      inherit.aes = FALSE, hjust = 0, vjust = -0.5, size = 4, colour = app_controls$colour_threshold) } +
      geom_vline(xintercept = current_x(), linetype = "dashed", colour = "grey50") +
      scale_colour_manual(values = c(
        "Fixed lambda_max = 1.24 (Frick)" = app_controls$colour_1.24,
        "Fixed lambda_max = 1.20 (Frick)" = app_controls$colour_1.20,
        "Realistic lambda_max (Niel & Lebreton, current s/alpha)" = app_controls$colour_realistic
      )) +
      labs(x = unique(dt$x_label), y = "PBR (bats/year)", colour = NULL,
           title = paste0("PBR vs. ", unique(dt$x_label)),
           subtitle = "Dashed grey: current slider value. Dotted navy: imposed thresholds.") +
      theme_minimal(base_size = app_controls$base_size) +
      theme(legend.position = "bottom", plot.subtitle = element_text(colour = "grey40"))
  })

  output$surface_plot <- renderPlot({
    res <- app_controls$surface_grid_res
    grid <- tidyr::expand_grid(
      s = seq(app_controls$s_sweep_range[1], app_controls$s_sweep_range[2], length.out = res),
      alpha = seq(app_controls$alpha_sweep_range[1], app_controls$alpha_sweep_range[2], length.out = res)
    ) %>%
      mutate(lambda_max = lambda_max_niel(s, alpha))

    current_pt <- tibble::tibble(s = input$s_fixed, alpha = input$alpha_fixed,
                                   lambda_max = lambda_realistic_fixed())

    ggplot(grid, aes(x = s, y = alpha, z = lambda_max)) +
      geom_raster(aes(fill = lambda_max)) +
      geom_contour(breaks = c(1.20, 1.24), aes(colour = after_stat(factor(level))), linewidth = 0.9) +
      scale_colour_manual(name = "lambda_max benchmark", values = setNames(c(app_controls$colour_1.20, app_controls$colour_1.24), c("1.2", "1.24"))) +
      scale_fill_viridis_c(option = "D", name = "lambda_max\n(realistic)") +
      geom_point(data = current_pt, aes(x = s, y = alpha), inherit.aes = FALSE,
                 colour = "white", fill = app_controls$colour_realistic, shape = 21, size = 4, stroke = 1.2) +
      labs(x = "Adult survival (s)", y = "Age at first breeding (alpha, years)",
           title = "Where do the current s/alpha sliders sit relative to Frick's benchmarks?",
           subtitle = "White-amber point: current slider position. Lines: lambda_max = 1.20 / 1.24 contours.") +
      theme_minimal(base_size = app_controls$base_size) +
      theme(plot.subtitle = element_text(colour = "grey40"))
  })
}

shinyApp(ui, server)
