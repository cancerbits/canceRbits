#' Run inferCNV
#'
#' @param counts Count matrix
#' @param conditions The cell conditions, or clusters/group labels; its length must match the number of columns of the count matrix
#' @param out_dir The output directory for inferCNV, will be created if it does not exist
#' @param ref_conditions The conditions that make up the reference (e.g. wildtype); default is NULL
#' and the mean of all cells is used as reference
#' @param gene_order_file The file with gene to location mappings; default is NULL and it will be downloaded from
#' https://data.broadinstitute.org/Trinity/CTAT/cnv/gencode_v19_gene_pos.txt
#' @param num_threads The number of CPU threads inferCNV will use; default is 1
#' @param keep_infercnv_scores Boolean, if TRUE infercnv scores are kept in the output directory, else only the
#' infercnv.png file is kept; default is FALSE
#'
#' @return NULL (inferCNV output will be in out_dir)
#'
#' @section Details:
#' Check ?infercnv::CreateInfercnvObject and ?infercnv::run for inferCNV specific help
#'
#' @importFrom utils download.file
#'
#' @export
cb_run_infercnv <- function(counts, conditions, out_dir, ref_conditions = NULL,
                            gene_order_file = NULL, num_threads = 1, keep_infercnv_scores = FALSE) {

  # create output directory if it does not exist
  dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

  # make sure we have a genome position file
  if (is.null(gene_order_file)) {
    gene_order_file <- file.path(out_dir, 'gencode_v19_gene_pos.txt')
  }
  if (!file.exists(gene_order_file)) {
    download.file(url = 'https://data.broadinstitute.org/Trinity/CTAT/cnv/gencode_v19_gene_pos.txt',
                  destfile = gene_order_file)
  }

  # set up meta data
  annot_df <- data.frame(condition = as.character(conditions))
  rownames(annot_df) <- colnames(counts)

  # create inferCNV object and run method
  infercnv_obj = infercnv::CreateInfercnvObject(
    raw_counts_matrix = as.matrix(counts),
    annotations_file = annot_df,
    ref_group_names = as.character(ref_conditions),
    gene_order_file = gene_order_file
  )

  infercnv_obj = infercnv::run(
    infercnv_obj,
    cutoff=0.1,  # use 1 for smart-seq, 0.1 for 10x-genomics
    out_dir=out_dir,
    cluster_by_groups=TRUE,   # cluster
    denoise=FALSE,
    HMM=FALSE,
    save_rds = FALSE,
    save_final_rds = FALSE,
    num_threads = num_threads
  )

  if (keep_infercnv_scores) {
    # convert the infercnv files to Rds
    for (f in list.files(path = out_dir, pattern = '^infercnv\\.(references|observations)\\.txt$', full.names = TRUE)) {
      mat <- read_infercnv_mat(file_path = f)
      saveRDS(object = mat, file = gsub(pattern = '\\.txt$', replacement = '.Rds', x = f))
    }
  }

  # delete everything we don't need
  unlink(x = list.files(path = out_dir, pattern = '\\.(txt|dat|preliminary\\.png)$', full.names = TRUE))

  return(invisible())
}

read_infercnv_mat <- function(file_path) {
  obs_table <- readr::read_delim(file = file_path, delim = ' ', quote = '"', skip = 1,
                                 col_names = FALSE, progress = FALSE, show_col_types = FALSE)
  mat <- t(as.matrix(obs_table[, -1]))
  cell_names <- readr::read_delim(file = file_path, delim = ' ', quote = '"', n_max = 1,
                                  col_names = FALSE, progress = FALSE, show_col_types = FALSE)
  cell_names <- as.character(cell_names[1, ])
  rownames(mat) <- cell_names
  colnames(mat) <- as.character(dplyr::pull(obs_table, 1))
  return(mat)
}
