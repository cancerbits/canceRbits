
#' Run a Seurat workflow
#'
#' @param x Sparse matrix or Seurat object
#' @param max_pc The number of principal components to use; default is 15
#' @param cluster_res Resolution parameter used with FindClusters; default is 0.3
#' @param metric Distance metric to use for FindNeighbors and RunUMAP; default is manhattan
#' @param k_param k.param parameter for FindNeighbors; default is 20
#' @param n_neighbors n.neighbors parameter for RunUMAP; default is 40
#' @param verbose Boolean indicating whether to print progress; is passed on
#' to most functions used internally; default is FALSE
#'
#' @return A list with items
#' * s - the final Seurat object
#' * figures - a list of figures
#' * fig_title - a figure title
#'
#' @section Details:
#' Performs a fairly standard Seurat analysis workflow including
#' * SCTransform
#' * RunPCA
#' * RunUMAP
#' * FindNeighbors
#' * FindClusters
#'
#' It will also find positive markers of the clusters and save the results in
#' object$RNA@misc$markers
#'
#' @section TODO:
#' Expose more of the parameters that are currently hard coded. Add examples.
#'
#' @export
#'
#' @importFrom Seurat CreateSeuratObject SCTransform RunPCA RunUMAP FindNeighbors
#' @importFrom Seurat FindClusters GetResidual DimPlot WhichCells DoHeatmap VlnPlot
#' @importFrom Seurat NoLegend
#' @importFrom rlang .data
#' @importFrom magrittr %>%
#' @importFrom dplyr mutate filter
#' @importFrom ggplot2 ggplot geom_bar
#'
#' @examples
#' \donttest{
#' # Coming soon
#' }
#' @md
cb_seurat_pipeline <- function(x, max_pc = 15, metric = 'manhattan',
                               k_param = 20, n_neighbors = 40,
                               cluster_res = 0.3, verbose = FALSE) {
  if (inherits(x = x, what = 'dgCMatrix')) {
    s <- CreateSeuratObject(counts = x)
  } else if (inherits(x = x, what = 'Seurat')) {
    s <- x
  } else {
    stop('x needs to be sparse matrix (dgCMatrix) or Seurat object')
  }

  s <- SCTransform(s, verbose = verbose, method = "qpoisson")
  s <- RunPCA(s, npcs = max_pc, verbose = verbose)
  dims <- 1:max_pc

  s <- RunUMAP(s, dims = dims, n.neighbors = n_neighbors, metric = metric,
               verbose = verbose)
  s <- FindNeighbors(s, reduction = "pca", dims = dims, verbose = verbose,
                     k.param = k_param, nn.method = "annoy", annoy.metric = metric)
  s <- FindClusters(s, resolution = cluster_res, verbose = verbose)

  # change to Seurat v5 Assay
  counts_mat <- GetAssayData(s, assay = "RNA", slot = "counts")
  cluster_labels <- s[["seurat_clusters", drop = TRUE]]

  s[["RNA"]]@misc$markers <- cb_pos_markers(counts = counts_mat, grouping = cluster_labels)

  s[["RNA"]]@misc$top_markers <- s[["RNA"]]@misc$markers %>%
    dplyr::group_by(.data$group) %>%
    dplyr::filter(.data$logFC > 0.5, .data$logCPM > 3) %>%
    dplyr::mutate(grp_rank = 1:dplyr::n()) %>%
    dplyr::filter(.data$grp_rank <= 5)

  s <- GetResidual(s, features = s[["RNA"]]@misc$top_markers$feature, verbose = verbose)
  s$log10_nCount_RNA <- log10(s[["nCount_RNA"]])

  # Plotting
  p1 <- DimPlot(s, label = TRUE, repel = TRUE)
  cells <- WhichCells(s, downsample = 100)
  p2 <- DoHeatmap(s, features = s[["RNA"]]@misc$top_markers$feature, 
                  slot = "scale.data", cells = cells) + NoLegend()
  p3 <- ggplot(s@meta.data, aes(.data$seurat_clusters, fill = .data$seurat_clusters)) + 
    geom_bar() + NoLegend()
  p4 <- VlnPlot(s, features = c('percent.mito', 'log10_nCount_RNA', 'nFeature_RNA'), pt.size = 0)

  common_title <- sprintf("Unsupervised clustering %s, %d cells", s@meta.data$orig.ident[1], ncol(s))

  return(list(s = s, figures = list(p1 = p1, p2 = p2, p3 = p3, p4 = p4), fig_title = common_title))
}

#' Load count matrix given file
#'
#' @param path File path (h5 or rds)
#' @param modality Slot to use in case file has multiple count matrices
#'
#' @return Sparse count matrix
#'
#' @section Details:
#' Looks at file extension to use the appropriate function to load the data.
#' Supported at the moment: H5, HDF5, RDS (capitalization does not matter)
#'
#' @export
#'
#' @importFrom tools file_ext
#' @importFrom Seurat Read10X_h5
cb_load_counts <- function(path, modality = 'Gene Expression') {
  ext <- toupper(file_ext(path))
  if (any(ext %in% c('H5', 'HDF5'))) {
    counts <- Read10X_h5(filename = path)
    if (is.list(counts)) {
      counts <- counts[[modality]]
    }
  } else if (ext == 'RDS') {
    object <- readRDS(file = path)
    if (inherits(x = object, what = 'dgCMatrix')) {
      counts <- object
    } else {
      stop('Expected dgCMatrix in Rds file')
    }
  } else {
    stop('No or unknown extension in path')
  }
  return(counts)
}
