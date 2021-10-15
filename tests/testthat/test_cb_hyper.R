
test_that("cb_hyper returns expected results on simple input", {
  background <- LETTERS
  genesets <- split(LETTERS, f = round((1:length(LETTERS)) / 4))
  names(genesets) <- paste0('gene_sets_', 1:length(genesets))
  signature <- c('N', 'O', 'V', 'W', 'X')
  hyper_res <- cb_hyper(signature = signature,
                        background = background,
                        genesets = genesets, min_size = 2,
                        collapse = FALSE, verbose = FALSE)
  expect_equal(dim(hyper_res), c(7, 8))
  expect_equal(hyper_res$label[1], 'gene_sets_7')
  expect_equal(hyper_res[1, 2], 0.034)

  genesets$gene_sets_8 <- setdiff(x = genesets$gene_sets_7, y = 'V')
  hyper_res <- cb_hyper(signature = signature,
                        background = background,
                        genesets = genesets, min_size = 2,
                        collapse = TRUE, verbose = FALSE)
  expect_equal(dim(hyper_res), c(7, 8))
  expect_equal(hyper_res$label[1], 'gene_sets_7')
  expect_equal(hyper_res[1, 2], 0.034)
})
