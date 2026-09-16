##
## Static PNG export for interactive plotly widgets, via a self-contained
## HTML file + webshot2 (headless Chrome/Edge). Avoids the 'kaleido'
## dependency (fragile to install across platforms/architectures); mirrors
## the same pattern already used for the 3D turbine coverage maps in the
## idf_bsh_dgy project (R/coverage_3d_topography.R, save_coverage_3d_plots()):
## htmlwidgets::saveWidget() with the camera baked into the widget's own
## layout, then a headless-browser screenshot of that static HTML.
##

render_plotly_screenshot <- function(widget, png_path, width = 1000, height = 750, delay = 2) {
  if (!requireNamespace("htmlwidgets", quietly = TRUE) || !requireNamespace("webshot2", quietly = TRUE)) {
    return(NULL)
  }
  html_path <- tempfile(fileext = ".html")
  htmlwidgets::saveWidget(widget, html_path, selfcontained = TRUE)
  tryCatch(
    {
      webshot2::webshot(html_path, png_path, vwidth = width, vheight = height, delay = delay)
      png_path
    },
    error = function(e) {
      message(
        "Could not render a static screenshot of the interactive plot (", conditionMessage(e),
        "). Needs the 'webshot2' package and a Chrome/Edge browser installed."
      )
      NULL
    }
  )
}
