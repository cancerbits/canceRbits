#' Find positive cluster markers
#'
#' @param counts Count matrix with features as rows and cells as columns
#' @param grouping Factor or character vector indicating the group membership of the cells
#' @param FDR_cutoff Remove results with FDR below this value; set to NULL or Inf to keep all
#' results; default is 0.05
#'
#' @return A data frame of results
#'
#' @section Details:
#' For each group, this creates pseudobulk data for the cells in that group and
#' for the cells that are not in the group. The two pseudobulk sets are compared
#' against each other using edgeR. Only results with positive fold-change are kept
#' and FDR values are calculated after all groups have been processed.
#'
#' @export
#'
#' @importFrom magrittr %>%
#' @importFrom dplyr filter mutate
#' @importFrom rlang .data
#' @importFrom stats p.adjust
#'
#' @examples
#' \donttest{
#' # Coming soon
#' }
#'
cb_pos_markers <- function(counts, grouping, FDR_cutoff = 0.05) {
  grouping <- droplevels(as.factor(grouping))
  de_lst <- lapply(levels(grouping), function(gl) {
    de_out <- pseudobulk_de(mat = counts,
                            grouping = factor(grouping == gl, levels = c('FALSE', 'TRUE')),
                            G = 3, test_method = 'edgeR', test_type = 'QLF')
    filter(de_out, .data$logFC > 0) %>%
      mutate(group = gl)
  })
  # combine, add FDR, filter
  res <- do.call(rbind, de_lst) %>%
    mutate(FDR = p.adjust(.data$pval, method = 'fdr'))
  if (!is.null(FDR_cutoff)) {
    res <- filter(res, .data$FDR < FDR_cutoff)
  }
  return(as.data.frame(res))
}

pseudobulk_de <- function(mat, grouping, G, test_method, test_type) {
  sel1 <- grouping == levels(grouping)[1]
  sel2 <- grouping == levels(grouping)[2]
  n1 <- sum(sel1)
  n2 <- sum(sel2)
  if (n1 < G | n2 < G) {
    stop('too few cells')
  }
  # create replicate labels per group
  replicate <- numeric(length = n1 + n2)
  replicate[sel1] <- sample(1:n1 %% G)
  replicate[sel2] <- sample(1:n2 %% G) + G
  replicate <- factor(replicate, levels = 0:(2*G-1))
  # create pseudobulk
  pb <- pseudobulk(counts = mat, grouping = replicate)
  # create pseudobulk grouping factor
  pb_group <- factor(1:(2*G) > G)

  if (test_method == 'edgeR') {
    res <- de_edger(mat = pb, grouping = pb_group, test_type = test_type)
  } else {
    stop('test_method unknown')
  }
  return(res)
}

pseudobulk <- function(counts, grouping) {
  mat <- sapply(levels(grouping), function(gr) {
    sparseMatrixStats::rowSums2(x = counts, cols = grouping == gr)
  })
  colnames(mat) <- levels(grouping)
  rownames(mat) <- rownames(counts)
  return(mat)
}

# run edgeR on bulk or pseudobulk data
de_edger <- function(mat, grouping, test_type) {
  design <- stats::model.matrix(~ grouping)

  y <- edgeR::DGEList(counts = mat, group = grouping)
  y <- edgeR::calcNormFactors(y)
  y <- edgeR::estimateDisp(y, design)

  if (test_type == 'exact') {
    test_res <- edgeR::exactTest(y)
  }
  else if (test_type == 'LRT') {
    fit <- edgeR::glmFit(y, design)
    test_res <- edgeR::glmLRT(fit)
  }
  else if (test_type == 'QLF') {
    fit <- edgeR::glmQLFit(y, design)
    test_res <- edgeR::glmQLFTest(fit)
  }
  else {
    stop('test_type must be "exact", "LRT" or "QLF"')
  }

  res <- test_res$table %>% tibble::rownames_to_column(var = 'feature') %>%
    dplyr::rename(pval = .data$PValue) %>%
    dplyr::mutate(FDR = stats::p.adjust(.data$pval, method = 'fdr')) %>%
    dplyr::arrange(.data$FDR, -abs(.data$logFC))
  return(res)
}
