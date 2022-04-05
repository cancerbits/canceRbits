#' Create a single sample report for scRNA-seq data
#'
#' @param sample_path Path to h5 file
#' @param out_report_path Path of the resulting report
#' @param sample_name Name to use for this sample; default is foo
#' @param out_rds_path Path of the resulting Seurat object (Rds file); default is NULL
#' and no file will. be written
#' @param return_seurat Boolean indicating whether to return the final seurat object;
#' default is FALSE and will return NULL invisibly
#' @param ... parameters passed to rmarkdown::render, e.g. quiet = TRUE
#'
#' @return NULL
#'
#' @section Details:
#' Uses the single_sample_overview rmarkdown template
#'
#' @export
#'
cb_single_sample_report <- function(sample_path,
                                    out_report_path,
                                    sample_name = 'foo',
                                    out_rds_path = NULL,
                                    return_seurat = FALSE, ...) {
  if (!requireNamespace("rmarkdown", quietly = TRUE)) {
    stop(
      "Package \"rmarkdown\" must be installed to use this function.",
      call. = FALSE
    )
  }

  rmd_path <- system.file('rmarkdown', 'templates', 'single_sample_overview',
                          'skeleton', 'skeleton.Rmd', package = 'canceRbits')

  rmarkdown::render(
    input = rmd_path,
    output_dir = dirname(out_report_path),
    knit_root_dir = tempdir(),
    envir = new.env(),
    params = list(sample_path = sample_path, sample_name = sample_name,
                  out_rds_path = out_rds_path),
    output_file = basename(out_report_path),
    ...
  )

  if (return_seurat) {
    return(readRDS(file = out_rds_path))
  } else {
    return(invisible())
  }
}
