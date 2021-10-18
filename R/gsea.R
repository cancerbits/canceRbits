#' Perform gene set enrichment using hypergeometric test
#'
#' @param signature Character vector of signature genes (the foreground set)
#' @param background Character vector of genes (the background set)
#' @param genesets List of character vectors of pre-defined gene sets (e.g. pathways)
#' @param min_size Minimum size of a gene set to be tested (after intersection with background) (default is 9)
#' @param max_size Maximum size of a gene set to be tested (after intersection with background) (default is 500)
#' @param collapse Boolean indicating whether in the output lower ranking gene sets should be removed if they are 99\% included in a higher ranking one (default is TRUE)
#' @param verbose Boolean indicating whether to show messages (default is TRUE)
#'
#' @return A data frame of results
#'
#' @section  Details:
#' This is a light wrapper around hypeR::hypeR. This functions filters the gene sets
#' given the background and makes sure that there is empty output if there is nothing to be tested.
#' Output is ordered based on p-value. Lower ranking gene sets that are 99\% contained in higher ranking ones can be omitted.
#'
#' This function will always output a data frame. In the case of no results it will have zero rows.
#' This behavior makes it compatible with a group_by, summarize dplyr workflow.
#'
#' @section Notes:
#' The background set should contain all (and only those) genes that were tested/considered to generate the signature gene set (aka foreground).
#'
#' To load local copies of pathways from EnrichR, check out the helper
#' function \code{\link{cb_enrichr_gsets}}.
#' To load them from the web check out \code{\link[hypeR:enrichr_gsets]{hypeR::enrichr_gsets}}.
#'
#' @section TODO:
#' Add documentation and more options to sanitize the geneset labels.
#'
#' @export
#'
#' @examples
#' \donttest{
#' background <- LETTERS
#' genesets <- split(LETTERS, f = round((1:length(LETTERS)) / 4))
#' names(genesets) <- paste0('gene_sets_', 1:length(genesets))
#' signature <- c('N', 'O', 'V', 'W', 'X')
#' cb_hyper(signature = signature,
#'          background = background,
#'          genesets = genesets, min_size = 2,
#'          collapse = FALSE, verbose = FALSE)
#' }
#'
#' @importFrom magrittr %>%
#' @importFrom dplyr mutate arrange
#' @importFrom rlang .data
cb_hyper <- function(signature, background, genesets, min_size = 9,
                      max_size = 500, collapse = TRUE, verbose = TRUE) {
  # make sure signature genes are inlcuded in background
  background <- union(signature, background)
  # intersect gene_sets with background
  genesets <- lapply(genesets, function(x) intersect(x, background))
  genesets_size <- sapply(genesets, length)
  # remove gene sets that are too small or too large
  genesets <- genesets[genesets_size >= min_size & genesets_size <= max_size]

  if (verbose) {
    message('cb_hyper: ', length(signature), ' of ', length(background), ' genes, ',
            length(genesets), ' genesets')
  }

  # return empty data frame if there is nothing to test
  if (length(signature) < 1 || length(genesets) < 1) {
    res <- data.frame(label = as.character(),
                      pval = as.double(),
                      fdr = as.double(),
                      signature = as.integer(),
                      geneset = as.integer(),
                      overlap = as.integer(),
                      background = as.integer(),
                      hits = as.character())
    return(res)
  }

  hyp_res <- hypeR::hypeR(signature = signature,
                          genesets = genesets,
                          test = 'hypergeometric',
                          background = background,
                          quiet = TRUE)
  hyp_res <- hyp_res$as.data.frame() %>% arrange(.data$pval)
  if (collapse & nrow(hyp_res) > 1) {
    n <- nrow(hyp_res)

    # check how much of a gene set is contained in a higher ranking one
    tmp <- matrix(0, n, n)
    for (i in 2:n) {
      pw_i <- hyp_res$label[i]
      size_i <- length(genesets[[pw_i]])
      for (j in 1:(i-1)) {
        pw_j <- hyp_res$label[j]
        n_com <- length(intersect(genesets[[pw_i]], genesets[[pw_j]]))
        tmp[i, j] <- n_com / size_i
      }
    }
    kick_out <- apply(tmp * lower.tri(tmp), 1, max) >= 0.99
    hyp_res <- hyp_res[!kick_out, ]
  }
  mutate(hyp_res, label = gsub(pattern = ' \\(GO:\\d+\\)$', replacement = '', x = .data$label))
}

#' Perform a scored gene set enrichment analysis
#'
#' @param genes Character vector of genes that have been scored
#' @param scores The scores associated with the genes
#' @param genesets List of character vectors of pre-defined gene sets (e.g. pathways)
#' @param min_size Minimum size of a gene set to be tested (after intersection with input genes) (default is 9)
#' @param max_size Maximum size of a gene set to be tested (after intersection with input genes) (default is 500)
#' @param n_leading_edge Number of genes to show in the leadingEdge output column - these are the genes driving the enrichment (default is 10)
#' @param verbose Boolean indicating whether to show messages (default is TRUE)
#' @param ... Otional arguments passed on to \code{\link[fgsea:fgsea]{fgsea::fgsea}}
#'
#' @return A table with GSEA results. Each row corresponds to a tested pathway
#'
#' @section  Details:
#' This is a light wrapper around \code{\link[fgsea:fgsea]{fgsea::fgsea}}. This functions filters the gene sets
#' given the background and makes sure that there is empty output if there is nothing to be tested.
#' Output is ordered based on p-value.
#'
#' This function will always output a data frame. In the case of no results it will have zero rows.
#' This behavior makes it compatible with a group_by, summarize dplyr workflow.
#'
#' To load local copies of pathways from EnrichR, check out the helper
#' function \code{\link{cb_enrichr_gsets}}.
#' To load them from the web check out \code{\link[hypeR:enrichr_gsets]{hypeR::enrichr_gsets}}.
#'
#' @section TODO:
#' Add documentation and more options to sanitize the geneset labels.
#'
#' @export
#'
#' @examples
#' \donttest{
#' background <- LETTERS
#' genesets <- split(LETTERS, f = round((1:length(LETTERS)) / 4))
#' names(genesets) <- paste0('gene_sets_', 1:length(genesets))
#' signature <- c('N', 'O', 'V', 'W', 'X')
#' cb_hyper(signature = signature,
#'          background = background,
#'          genesets = genesets, min_size = 2,
#'          collapse = FALSE, verbose = FALSE)
#' }
#'
#' @importFrom magrittr %>%
#' @importFrom dplyr mutate arrange group_by ungroup
#' @importFrom rlang .data
cb_fgsea <- function(genes, scores, genesets, min_size = 9, max_size = 500,
                     n_leading_edge = 10, verbose = TRUE, ...) {
  # If there is nothing to test, return empty data frame
  if (length(genes) == 0) {
    res <- data.frame(pathway = as.character(),
                      pval = as.double(),
                      padj = as.double(),
                      log2err = as.double(),
                      ES = as.double(),
                      NES = as.double(),
                      size = as.double(),
                      leadingEdge = as.list())
    return(res)
  }

  if (min(scores) < 0) {
    if (max(scores) > 0) {
      score_type <- 'std'
    }
    else {
      if (verbose) {
        message('No positive gene scores given - will set score_type to neg in fgsea')
      }
      score_type <- 'neg'
    }
  } else {
    if (verbose) {
      message('No negative gene scores given - will set score_type to pos in fgsea')
    }
    score_type <- 'pos'
  }

  names(scores) <- genes
  if (verbose) {
    message(sprintf('Run fGSEA with %d genes', length(scores)))
    #message('Top genes', head(names(scores), 10))
    #message('Bottom genes', rev(tail(names(scores), 10)))
  }
  fgsea::fgsea(pathways = genesets,
               stats    = scores,
               scoreType = score_type,
               minSize  = min_size,
               maxSize  = max_size,
               ...) %>%
    mutate(pathway = gsub(pattern = '^GOBP_', replacement = '', x = .data$pathway)) %>%
    group_by(.data$pathway) %>%
    mutate(
      tmp_n = length(unlist(.data$leadingEdge)),
      leadingEdge = paste(unlist(.data$leadingEdge)[0:min(.data$tmp_n, n_leading_edge)], collapse = ','),
      tmp_n = NULL) %>%
    ungroup() %>%
    arrange(.data$pval)
}

#' Load local copy of EnrichR gene sets
#'
#' @param filepath Path to tsv file
#' @param verbose Boolean indicating whether to show messages (default is TRUE)
#'
#' @return A data frame of results
#'
#' @section Details:
#' The input tsv files are meant to be in the format as provided
#' at \url{https://maayanlab.cloud/Enrichr/#libraries}.
#'
#' @export
#'
#' @examples
#' \donttest{
#' filepath <- paste0('/data_synology_rg3/cancerbits/users/',
#'                    'christoph_h/projects/gsea/enrichr_lib/',
#'                    'GO_Biological_Process_2021.tsv')
#' genesets <- cb_enrichr_gsets(filepath = filepath)
#' }
#'
cb_enrichr_gsets <- function(filepath, verbose = TRUE) {
  con = file(filepath, "r")
  lines = readLines(con)
  genesets_raw <- stringr::str_split(string = lines, pattern = '\t\t', simplify = TRUE)
  close(con)

  genesets <- list()
  for (i in 1:nrow(genesets_raw)) {
    x <- genesets_raw[i, ]
    gs_genes <- stringr::str_split(string = x[2], pattern = '\t')[[1]]
    gs_genes <- gsub(pattern = ',[0-9.]*$', replacement = '', x = gs_genes)
    genesets[[i]] <- setdiff(gs_genes, '')
  }
  names(genesets) <- genesets_raw[, 1]

  if (verbose) {
    message('Loaded ', length(genesets), ' genesets')
  }

  return(genesets)
}
