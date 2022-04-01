test_that("cb_filter_count_matrix returns expected results on simple input", {
  filter_out <- canceRbits::cb_filter_count_matrix(sctransform::pbmc)
  expect_equal(dim(filter_out$filtered), c(914, 120))
})
