#' canceRbits: Cancer genomics toolbox
#'
#' Wrapper and helper functions commmonly used by the Developmental
#' Cancer Genomics group at St. Anna Children's Cancer Research Institute.
#'
#' @section GSE functions:
#' These are the main functions used for gene set enrichment:
#' \enumerate{
#'   \item \code{\link{cb_hyper}}
#'   \item \code{\link{cb_fgsea}}
#' }
#'
#' @section QC/filter functions:
#' \enumerate{
#'   \item \code{\link{cb_filter_count_matrix}}
#' }
#'
#' @section Workflows:
#' These functions wrap multiple common commands that are usually used together
#' as part of a workflow
#' \enumerate{
#'   \item \code{\link{cb_seurat_pipeline}}
#' }
#'
#' @section Reports:
#' These functions wrap rmarkdown templates to generate html reports
#' \enumerate{
#'   \item \code{\link{cb_single_sample_report}}
#' }
#'
#' @docType package
#' @name canceRbits
NULL
