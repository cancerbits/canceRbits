# functions useful for annotating single-cell data

# write celltypist compatible count matrix and feature and cell files to a specific directory
# output is list with filenames
#' @importFrom Matrix writeMM
#' @importFrom utils write.table
to_mtxgz <- function(path, x) {
  if (!requireNamespace("R.utils", quietly = TRUE)) {
    stop(
      "Package \"R.utils\" must be installed to use this function.",
      call. = FALSE
    )
  }
  if (!inherits(x = x, what = 'dgCMatrix')) {
    stop('Input must be a sparse matrix (dgCMatrix)')
  }

  dir.create(path, showWarnings = FALSE)
  mfile <- file.path(path, "matrix.mtx")
  rfile <- file.path(path, "rownames.tsv")
  cfile <- file.path(path, "colnames.tsv")

  if (any(file.exists(mfile, rfile, cfile))) {
    stop('One or more of the output files already exist')
  }

  writeMM(x, file = mfile)
  R.utils::gzip(mfile)

  write.table(x = rownames(x), file = rfile, quote = FALSE, row.names = FALSE, col.names = FALSE)
  write.table(x = colnames(x), file = cfile, quote = FALSE, row.names = FALSE, col.names = FALSE)

  return(list(matrix = paste0(mfile, '.gz'), rownames = rfile, colnames = cfile))
}

#' Run celltypist
#'
#' @param counts Count matrix with features as rows and cells as columns
#' @param ct_models Name(s) of celltypist model; default is 'Immune_All_High.pkl'
#' @param ct_mode Celltypist label prediction mode, either 'best_match' or 'prob_match'; default is 'best_match'
#' @param verbose Show or hide progress; default is FALSE
#'
#' @return A named list (one entry per model) of data frames of predictions
#'
#' @section Details:
#' This is a simple wrapper that writes the count matrix to a temporary file,
#' calls celltyper, reads the results.
#' Check \url{https://www.celltypist.org} for celltypist specific details.
#' For celltypist to find existing models, use \code{Sys.setenv(CELLTYPIST_FOLDER = path)},
#' where \code{path} must point to the directory that contains the data/models directories.
#'
#' If \code{s} is a Seurat object, \code{ct_out <- celltypist(counts = s$RNA@counts)} would
#' run celltypist to get the annotations, \code{s <- AddMetaData(s, metadata = ct_out[[1]])}
#' would add them to the meta data of \code{s}
#'
#' @export
#'
#' @importFrom utils read.csv
#'
#' @examples
#' \donttest{
#' # Coming soon
#' }
#'
cb_celltypist <- function(counts, ct_models = 'Immune_All_High.pkl',
                         ct_mode = 'best_match', verbose = FALSE) {
  if (!requireNamespace("sys", quietly = TRUE)) {
    stop(
      "Package \"sys\" must be installed to use this function.",
      call. = FALSE
    )
  }
  if (!inherits(x = counts, what = 'dgCMatrix')) {
    stop('Input must be a sparse matrix (dgCMatrix)')
  }

  # check if celltypist is installed
  ret <- system(command = 'which celltypist', intern = FALSE, ignore.stdout = TRUE, ignore.stderr = TRUE)
  if (ret > 0) {
    stop(sprintf('The system call "which celltypist" had exit code %d - please make sure celltypist is installed and paths are set correctly', ret))
  }

  # set up a unique temp dir
  tmp_path <- tempfile(pattern = 'celltypist_')
  dir.create(tmp_path)
  on.exit({
    if (file.exists(tmp_path)) {
      unlink(tmp_path, recursive = TRUE)
    }
  })

  if (verbose) {
    message('Write temporary mtx.gz file')
  }
  mtx <- to_mtxgz(path = tmp_path, x = counts)

  pred_list <- lapply(ct_models, function(ct_model) {
    if (verbose) {
      message('Call celltypist with model ', ct_model)
    }
    res <- sys::exec_wait(cmd = "celltypist",
                          args = c("--indata", mtx$matrix,
                                   "--model", ct_model,
                                   "--transpose-input",
                                   "--gene-file", mtx$rownames,
                                   "--cell-file", mtx$colnames,
                                   "--mode", ct_mode,
                                   "--majority-voting",
                                   "--outdir", tmp_path,
                                   "--quiet"),
                          std_out = verbose, std_err = verbose)
    if (verbose) {
      message('Read results')
    }
    read.csv(file = file.path(tmp_path, 'predicted_labels.csv'), row.names = 1)
  })
  names(pred_list) <- ct_models

  return(pred_list)
}

