##
## Static PNG export for interactive plotly widgets, via a self-contained
## HTML file + webshot2 (headless Chrome/Edge). Avoids the 'kaleido'
## dependency (fragile to install across platforms/architectures); mirrors
## the same pattern already used for the 3D turbine coverage maps in the
## idf_bsh_dgy project (R/coverage_3d_topography.R, save_coverage_3d_plots()):
## htmlwidgets::saveWidget() with the camera baked into the widget's own
## layout, then a headless-browser screenshot of that static HTML.
##

# html_path: pass a real path (not the default tempfile) to keep the
# self-contained interactive HTML as a deliverable alongside the static
# PNG -- Word/PowerPoint cannot embed a live plotly widget, but the HTML
# file opens directly in any browser (no server needed) and can be shared
# alongside the docx/pptx for anyone who wants to explore the surface
# interactively (rotate, hover for exact values, toggle Fr scenarios).
render_plotly_screenshot <- function(widget, png_path, width = 1000, height = 750, delay = 2,
                                      html_path = tempfile(fileext = ".html")) {
  if (!requireNamespace("htmlwidgets", quietly = TRUE)) {
    return(list(png = NULL, html = NULL))
  }
  htmlwidgets::saveWidget(widget, html_path, selfcontained = TRUE)

  if (!requireNamespace("webshot2", quietly = TRUE)) {
    message("Package 'webshot2' not installed -- the interactive HTML was still saved to ", html_path,
            ", but no static PNG could be rendered for the report.")
    return(list(png = NULL, html = html_path))
  }
  png_result <- tryCatch(
    {
      webshot2::webshot(html_path, png_path, vwidth = width, vheight = height, delay = delay)
      png_path
    },
    error = function(e) {
      message(
        "Could not render a static screenshot of the interactive plot (", conditionMessage(e),
        "). Needs the 'webshot2' package and a Chrome/Edge browser installed. The interactive HTML was still saved to ", html_path, "."
      )
      NULL
    }
  )
  list(png = png_result, html = html_path)
}
