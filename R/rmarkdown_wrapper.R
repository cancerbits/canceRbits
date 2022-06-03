#' Create a single sample report for scRNA-seq data
#'
#' @param sample_counts Either a sparse count matrix, a path to h5 file, or a path to an rds file
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
cb_single_sample_report <- function(sample_counts,
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

  # get the path to the rmarkdown template that comes with this package
  rmd_path <- system.file('rmarkdown', 'templates', 'single_sample_overview',
                          'skeleton', 'skeleton.Rmd', package = 'canceRbits')

  # if the input is not a path, write the matrix to a file
  cleanup_in_rds <- FALSE
  if (inherits(x = sample_counts, what = 'dgCMatrix')) {
    tmp_rds_path <- paste0(tempfile(), '.Rds')
    saveRDS(object = sample_counts, file = tmp_rds_path)
    sample_counts <- tmp_rds_path
    cleanup_in_rds <- TRUE
  } else if (!inherits(x = sample_counts, what = 'character')) {
    stop('sample_counts needs to be a dgCMatrix or a character vector of length one')
  }

  # if we want to return a Seurat object, but don't want it written to a
  # permanent file, we need a tmp file
  cleanup_out_rds <- FALSE
  if (return_seurat && is.null(out_rds_path)) {
    out_rds_path <- paste0(tempfile(), '.Rds')
    cleanup_out_rds <- TRUE
  }

  # set up a unique temp knit root dir
  krd <- tempfile()
  dir.create(krd)

  old_option <- options('knitr.duplicate.label')
  options(knitr.duplicate.label = "allow")

  rmarkdown::render(
    input = rmd_path,
    output_dir = dirname(out_report_path),
    knit_root_dir = krd,
    envir = new.env(),
    params = list(sample_path = sample_counts,
                  sample_name = sample_name,
                  out_rds_path = out_rds_path),
    output_file = basename(out_report_path),
    ...
  )

  options(knitr.duplicate.label = old_option)

  # remove the knit root directory
  unlink(x = krd, recursive = TRUE, force = TRUE)

  if (cleanup_in_rds) {
    unlink(sample_counts)
  }

  if (return_seurat) {
    s <- readRDS(file = out_rds_path)
    if (cleanup_out_rds) {
      unlink(out_rds_path)
    }
    return(s)
  } else {
    return(invisible())
  }
}


#' Create an inferCNV report
#'
#' @param counts Count matrix
#' @param conditions The cell conditions, or clusters/group labels; its length must match the number of columns of the count matrix
#' @param out_report_path Path of the resulting report
#' @param out_dir The output directory for inferCNV, will be created if it does not exist; the default is NULL and a temporary directory will be used
#' @param ref_conditions The conditions that make up the reference (e.g. wildtype); default is NULL
#' and the mean of all cells is used as reference
#' @param gene_order_file The file with gene to location mappings; default is NULL and it will be downloaded from
#' https://data.broadinstitute.org/Trinity/CTAT/cnv/gencode_v19_gene_pos.txt
#' @param num_threads The number of CPU threads inferCNV will use; default is 1
#' @param keep_infercnv_scores Boolean, if TRUE infercnv scores are kept in the output directory, else only the
#' infercnv.png file is kept; default is FALSE
#' @param ... parameters passed to rmarkdown::render, e.g. quiet = TRUE
#'
#' @return NULL
#'
#' @section Details:
#' See cb_run_infercnv for details
#'
#' @export
cb_infercnv_report <- function(counts, conditions, out_report_path, out_dir = NULL,
                               ref_conditions = NULL, gene_order_file = NULL,
                               num_threads = 1, keep_infercnv_scores = FALSE, ...) {


  if (!requireNamespace("rmarkdown", quietly = TRUE)) {
    stop(
      "Package \"rmarkdown\" must be installed to use this function.",
      call. = FALSE
    )
  }

  # get the path to the rmarkdown template that comes with this package
  rmd_path <- system.file('rmarkdown', 'templates', 'infercnv',
                          'skeleton', 'skeleton.Rmd', package = 'canceRbits')

  # create temp data file
  data_path <- paste0(tempfile(pattern = 'infercnv_report_data_'), '.Rds')
  saveRDS(object = list(counts = counts, conditions = conditions, ref_conditions = ref_conditions),
          file = data_path)

  # do we need to clean up?
  do_clean_out_dir <- TRUE
  if (!is.null(out_dir) | keep_infercnv_scores) {
    do_clean_out_dir <- FALSE
  }

  # set out_dir if missing
  if (is.null(out_dir)) {
    out_dir <- tempfile(pattern = 'infercnv_out_dir_')
  }

  # make sure all paths are expanded
  out_report_path <- path.expand(out_report_path)
  out_dir <- path.expand(out_dir)
  if (!is.null(gene_order_file)) {
    gene_order_file <- path.expand(gene_order_file)
  }

  # set up a unique temp knit root dir
  krd <- tempfile()
  dir.create(krd)

  old_option <- options('knitr.duplicate.label')
  options(knitr.duplicate.label = "allow")

  # render the report
  rmarkdown::render(
    input = rmd_path,
    output_dir = dirname(out_report_path),
    knit_root_dir = krd,
    envir = new.env(),
    params = list(data_path = data_path,
                  gene_order_file = gene_order_file,
                  out_dir = out_dir,
                  num_threads = num_threads,
                  keep_infercnv_scores = keep_infercnv_scores),
    output_file = basename(out_report_path),
    ...
  )

  options(knitr.duplicate.label = old_option)

  # remove the knit root directory
  unlink(x = krd, recursive = TRUE, force = TRUE)

  # Clean up
  if (do_clean_out_dir) {
    message('Removing the temporary output directory')
    unlink(x = out_dir, recursive = TRUE)
  } else {
    message('inferCNV output files are in ', out_dir)
  }

  unlink(data_path)

  return(invisible())
}

