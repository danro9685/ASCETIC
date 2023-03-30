context("ASCETIC")

data("resExampleSingleSamplesResampling")
test_that("ASCETIC produces correct output", {
    expect_equal(names(resExampleSingleSamplesResampling),c("dataset","ccf_dataset","ranking_estimate","poset","inference"))
})

data("resExamplePhylogeniesDatasetBootstrap")
test_that("ASCETIC produces correct output", {
    expect_equal(names(resExamplePhylogeniesDatasetBootstrap),c("dataset","models","ranking_estimate","poset","inference"))
})
