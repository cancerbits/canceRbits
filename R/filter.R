#' Filter a scRNA-seq count matrix based on QC metrics
#'
#' @param counts Count matrix with features as rows and cells as columns
#' @param percent_mito_th Max percentage of mitochondrial counts; cell with more will be removed; default is 15
#' @param log_counts_z_th Z-score thresholds for sum of counts and sum of features
#' (these cell statistics will be log-transformed and scaled); default is \code{c(-3, 3)}
#' @param feature_outlier_z_th Z-score thresholds for the residuals of the UMI sum vs number
#' of features relationship; default is \code{c(-5, 5)}
#' @param min_features Min of detected features; cells with less will be removed; default is 300
#' @param mito_pattern Regex pattern used to identify mitochondrial features; default is \code{"^MT[^0-9a-zA-Z]+"}
#' @param return_seurat Boolean indicating whether to return the filtered matrix as Seurat object; default is TRUE
#' @param sample_id Optional string giving the sample name; will be used if return value is a Seurat object; default is NULL
#' @param verbose Boolean indicating whether to show messages; default is TRUE
#'
#' @return A list with components:
#' 1) "filtered" filtered object,
#' 2) "figures" list of figures,
#' 3) "fig_title" figure title
#'
#' @section Details:
#' Four filters are applied in succession:
#' 1) percent mitochondrial reads
#' 2) number of features
#' 3) sum of counts
#' 4) transcripts vs genes
#'
#' Filters 1 to 3 are applied to the cell-level statistics; for 2) and 3) these are log10-transformed before Z-scores are calculated.
#' Filter 4) applies a loess fit to the relationship of the log10-transformed values of number of features as function of sum of counts, the residuals are then z-scored.
#'
#' @export
#'
#' @examples
#' \donttest{
#' filter_out <- canceRbits::cb_filter_count_matrix(sctransform::pbmc)
#' # show figures
#' patchwork::wrap_plots(filter_out$figures) + patchwork::plot_annotation(title = filter_out$fig_title)
#' # show filtered object
#' show(filter_out$filtered)
#' }
#'
#' @importFrom stats loess
#' @importFrom Seurat CreateSeuratObject PercentageFeatureSet
#' @importFrom ggplot2 ggplot geom_histogram geom_vline xlab ylab annotate
#' @importFrom ggplot2 theme_get unit scale_x_continuous theme element_blank
#' @importFrom ggplot2 annotation_logticks geom_point scale_color_manual scale_y_continuous
#' @importFrom ggplot2 aes
#' @importFrom rlang .data
#'

cb_filter_count_matrix <- function(
  counts,
  percent_mito_th = 15,
  log_counts_z_th = c(-3, 3),
  feature_outlier_z_th = c(-5, 5),
  min_features = 300,
  mito_pattern = "^MT[^0-9a-zA-Z]+",
  return_seurat = TRUE,
  sample_id = NULL,
  verbose = TRUE) {
  
  s <- CreateSeuratObject(counts = counts)
  
  s[['orig.ident']] <- sample_id
  s[['percent.mito']] <- PercentageFeatureSet(s, pattern = mito_pattern)
  
  keep <- 1:ncol(s)

  ## --- 1) percent mitochondrial reads ---
  keep_this <- s[['percent.mito']][,1] <= percent_mito_th
  txt <- sprintf("Remove %d of %d\n(%1.1f%%) cells", sum(!keep_this), length(keep), sum(!keep_this) / length(keep) * 100)
  p1 <- ggplot(s@meta.data, aes(x = .data$percent.mito)) +
    geom_histogram(binwidth = 1) +
    geom_vline(xintercept = percent_mito_th, color = 'red') +
    xlab('Mitochondrial reads in %') +
    ylab('Cells') +
    annotate('label', x = Inf, y = Inf, label = txt, size = theme_get()$text$size * 0.75 / ggplot2::.pt, 
             vjust = 1, hjust = 1, label.r = unit(0, 'cm'))
  keep <- keep[keep_this]

  ## --- 2) number of features ---
  nFeature_RNA <- s@meta.data$nFeature_RNA[keep]
  nFeature_RNA_logscaled <- scale(log10(nFeature_RNA))
  keep_this <- nFeature_RNA >= min_features &
    nFeature_RNA_logscaled >= log_counts_z_th[1] &
    nFeature_RNA_logscaled <= log_counts_z_th[2]
  th <- range(nFeature_RNA[keep_this]) + c(-1e-10, 1e-10)
  txt <- sprintf("Remove %d of %d\n(%1.1f%%) cells", sum(!keep_this), length(keep), sum(!keep_this) / length(keep) * 100)
  p2 <- ggplot(s@meta.data[keep, ], aes(x = .data$nFeature_RNA)) +
    geom_histogram(binwidth = 33) +
    geom_vline(xintercept = th, color = 'red') +
    xlab('Number of genes detected') +
    ylab('Cells') +
    annotate('label', x = Inf, y = Inf, label = txt, size = theme_get()$text$size * 0.75 / ggplot2::.pt, 
             vjust = 1, hjust = 1, label.r = unit(0, 'cm'))
  keep <- keep[keep_this]

  ## --- 3) number of molecules ---
  nCount_RNA <- s@meta.data$nCount_RNA[keep]
  nCount_RNA_logscaled <- scale(log10(nCount_RNA))
  keep_this <- nCount_RNA_logscaled >= log_counts_z_th[1] &
               nCount_RNA_logscaled <= log_counts_z_th[2]
  th <- range(nCount_RNA[keep_this]) + c(-1e-10, 1e-10)
  txt <- sprintf("Remove %d of %d\n(%1.1f%%) cells", sum(!keep_this), length(keep), sum(!keep_this) / length(keep) * 100)
  p3 <- ggplot(s@meta.data[keep, ], aes(x = log10(.data$nCount_RNA))) +
    geom_histogram(binwidth = 0.04) +
    geom_vline(xintercept = log10(th), color = 'red') +
    xlab('Number of transcripts') +
    ylab('Cells') +
    scale_x_continuous(breaks = log10(10^(0:7)), labels = as.character(10^(0:7))) +
    theme(panel.grid.minor = element_blank()) +
    annotation_logticks(sides = 'b') +
    annotate('label', x = Inf, y = Inf, label = txt, size = theme_get()$text$size * 0.75 / ggplot2::.pt, 
             vjust = 1, hjust = 1, label.r = unit(0, 'cm'))
  keep <- keep[keep_this]

  ## --- 4) transcripts vs genes ---
  md <- s@meta.data[keep, ]
  mod <- loess(log10(nFeature_RNA) ~ log10(nCount_RNA), data = md, span = 1)
  md$nFeature_outlier <- scale(mod$residuals) < feature_outlier_z_th[1] | 
                         scale(mod$residuals) > feature_outlier_z_th[2]
  keep_this <- !md$nFeature_outlier
  txt <- sprintf("Remove %d of %d\n(%1.1f%%) cells", sum(!keep_this), length(keep), sum(!keep_this) / length(keep) * 100)
  p4 <- ggplot(md, aes(x = log10(.data$nCount_RNA), y = log10(.data$nFeature_RNA), color = .data$nFeature_outlier)) +
    geom_point(size = 0.5) + 
    scale_color_manual(values = c('grey35', 'red'), guide = 'none') +
    xlab('Transcripts') + 
    ylab('Genes') +
    scale_y_continuous(breaks = log10((1:55)*1000), labels = (1:55)*1000) +
    scale_x_continuous(breaks = log10(10^(0:7)), labels = as.character(10^(0:7))) +
    theme(panel.grid.minor = element_blank()) +
    annotation_logticks() +
    annotate('label', x = Inf, y = -Inf, label = txt, size = theme_get()$text$size * 0.75 / ggplot2::.pt, 
             vjust = 0, hjust = 1, label.r = unit(0, 'cm'))
  keep <- keep[keep_this]

  keep_cells <- colnames(s)[keep]
  txt <- sprintf('Sample %s; Keeping %d of %d cells (%1.2f percent)\n',
                 sample_id, length(keep_cells), ncol(s),
                 length(keep_cells)/ncol(s)*100)
  if (verbose) message(txt)

  s <- s[, keep_cells, drop = FALSE]
  
  # ✅ Seurat v5 fix: use GetAssayData instead of slot access
  if (!return_seurat) {
    s <- GetAssayData(s, slot = "counts", assay = "RNA")
  }

  return(list(filtered = s, figures = list(p1, p2, p3, p4), fig_title = txt))
}
