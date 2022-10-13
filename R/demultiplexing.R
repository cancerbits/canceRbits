
#' De-multiplex
#'
#' @param tag_counts Count matrix with cells rows and barcodes (tags) as columns
#' @param min_count Minimum number of total tag counts required per cell
#' @param pseudo_count The pseudo-count that will be added to all observed counts before
#' log-transformation
#'
#' @return A data frame of results
#' The 'tag' column holds the final tag assignment and will match one of the
#' input column names, 'neg', or 'multi'
#'
#' @section Details:
#'
#' @export
#'
#' @examples
#' \donttest{
#' # Coming soon
#' }
#'
cb_demux_gmm <- function(tag_counts, min_count = 30, pseudo_count = 10) {

  # set up return data frame
  res <- data.frame(cell = rownames(tag_counts))
  res$has_NA <- apply(is.na(tag_counts), 1, any)

  # replace any missing values with lowest barcode count
  tag_counts <- apply(tag_counts, 2, function(y) {
    y[is.na(y)] <- min(y, na.rm = TRUE)
    y
  })
  res$counts_sum <- rowSums(tag_counts)

  # filter and log-transform
  #res$tag <- rep(NA_character_, nrow(res))
  res$tag <- rep(NA, nrow(res))
  sel <- res$counts_sum >= min_count
  tcl <- log10(tag_counts[sel, ] + pseudo_count)

  # demux pairs
  dp <- list()
  for (i in 1:(ncol(tcl) - 1)) {
    for (j in (i + 1):ncol(tcl)) {
      tmp <- demux_one_pair(x = tcl[, i], y = tcl[, j],
                            class_names = c('neg', colnames(tcl)[i], colnames(tcl)[j], 'multi'))
      dp[[length(dp) + 1]] <- tmp
    }
  }
  dp <- as.data.frame(dp, col.names = paste0('V', 1:length(dp)))

  cell_tag <- apply(dp, 1, function(x) {
    x <- as.character(x)
    if (any(x == 'multi')) {
      return('multi')
    }
    x <- unique(x[x != 'neg'])
    if (length(x) > 1) {
      return('multi')
    }
    if (length(x) == 1) {
      return(x)
    }
    return('neg')
  })
  res$tag[sel] <- cell_tag
  res$tag <- factor(res$tag, levels = c('neg', colnames(tag_counts), 'multi'))
  res$assigned <- res$tag %in% colnames(tag_counts)

  res$tag_count <- rep(NA_integer_, nrow(res))
  res$tag_count[res$assigned] <- sapply(which(res$assigned), function(i) tag_counts[res$cell[i], as.character(res$tag[i])])
  res$tag_frac <- res$tag_count / res$counts_sum
  return(res)
}

#' Demultiplex one pair of barcode counts
#'
#' @return Character vector or factor of demultiplex classes
#'
#' @param x log-counts of barcode 1
#' @param y log-counts of barcode 2
#' @param N Number of training cells to use for semi-supervised classification
#' @param class_names The names of the output classes
#' @param label_levels If no NULL then the output will be a factor with these levels
#'
#' @keywords internal
#'
demux_one_pair <- function(x, y, N = 50, class_names = c('neg', 'x', 'y', 'multi'), label_levels = NULL) {
  d <- x - y
  m <- (x + y) / 2
  rd <- d / m

  # use semi-supervised approach for initial classification in 1D
  init_labels <- rep(NA, length(y))
  init_labels[order(rd)[1:N]] <- 1
  init_labels[order(-rd)[1:N]] <- 3
  mneg <- mean(rd[!is.na(init_labels)])
  init_labels[order(abs(rd - mneg))[1:N]] <- 2

  mod_rd <- mclust::MclustSSC(data = rd, class = init_labels, G = 3, modelNames = 'E', verbose = FALSE)
  #plot(m, rd, col = mod_rd$classification)

  # generate fake doublets
  dub_n <- min(sum(mod_rd$classification == 1), sum(mod_rd$classification == 3))
  #dub_n <- max(round(dub_n * 0.05), 13)
  dub_x <- sample(x[mod_rd$classification == 3], size = dub_n,
                  prob = x[mod_rd$classification == 3]^2, replace = TRUE)
  dub_y <- sample(y[mod_rd$classification == 1], size = dub_n,
                  prob = y[mod_rd$classification == 1]^2, replace = TRUE)
  dub_d <- dub_x - dub_y
  dub_m <- (dub_x + dub_y) / 2
  dub_rd <- dub_d / dub_m

  # reset the initial labels for the negative group
  init_labels[init_labels == 2] <- NA
  tmp <- m
  tmp[mod_rd$classification != 2] <- Inf
  init_labels[order(tmp)[1:N]] <- 2

  m <- c(m, dub_m)
  rd <- c(rd, dub_rd)
  init_labels <- c(init_labels, rep(4, dub_n))
  #plot(m, rd, col = init_labels)

  mod_2d <- mclust::MclustSSC(data = cbind(m, rd), class = init_labels, G = 4, modelNames = 'VVV', verbose = FALSE)
  #plot(mod_2d, what = 'classification')
  #title(sprintf('%s vs %s', class_names[2], class_names[3]))

  res <- rep(class_names[1], length(rd))
  res[mod_2d$classification == 4] <- class_names[4]
  res[mod_2d$classification == 1] <- class_names[3]
  res[mod_2d$classification == 3] <- class_names[2]
  res <- res[-which(init_labels == 4)]
  if (!is.null(label_levels)) {
    res <- factor(res, levels = label_levels)
  }
  return(res)
}
