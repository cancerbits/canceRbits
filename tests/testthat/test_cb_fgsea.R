
test_that("cb_fgsea returns expected results on simple input", {
  background <- LETTERS
  genesets <- split(LETTERS, f = round((1:length(LETTERS)) / 4))
  names(genesets) <- paste0('gene_set_', 1:length(genesets))
  genes <- LETTERS
  scores <- 1:(length(LETTERS))

  set.seed(75820)
  fgsea_res <- cb_fgsea(genes = genes,
                        scores = scores,
                        genesets = genesets,
                        min_size = 2,
                        verbose = FALSE)

  expect_equal(dim(fgsea_res), c(7, 8))
  expect_equal(fgsea_res$pathway[1], 'gene_set_7')
  expect_equal(fgsea_res$NES[1], 2.319275, tolerance = 1e-6)
})
