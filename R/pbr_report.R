##
## Builds and renders the PBR technical note (.docx) from the objects
## computed by run_pbr_analysis() (R/pbr_analysis.R), passed in as params.
## The template (report/pbr_technical_note_template.Rmd) does not
## recompute anything -- it only formats params already calculated here.
## Mirrors R/report.R's build_idf_report() in idf_bsh_dgy.
##
## reference_docx (optional): path to a .docx whose styles, header/footer
## (including page numbering, if the .docx has it) and page setup pandoc
## reuses in the generated document -- the CONTENT of that .docx
## (body text/images) is ignored, only styles/sectPr count. NULL (the
## default) keeps rmarkdown::word_document()'s generic styling.
##

build_pbr_report <- function(output_file, report_params,
                              template = "report/pbr_technical_note_template.Rmd",
                              reference_docx = NULL) {

  out_dir <- dirname(output_file)
  dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

  output_options <- if (!is.null(reference_docx)) {
    list(reference_docx = normalizePath(reference_docx))
  } else {
    NULL
  }

  rmarkdown::render(
    input          = template,
    output_dir     = normalizePath(out_dir),
    output_file    = basename(output_file),
    output_options = output_options,
    params         = report_params,
    envir          = new.env()
  )
}
