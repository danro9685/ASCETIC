context("ASCETIC")

data("res_example_single_samples_resampling")
test_that("ASCETIC produces correct output", {
    expect_equal(names(res_example_single_samples_resampling),c("dataset","ccf_dataset","ranking_estimate","poset","inference"))
})

data("res_example_phylogenies_dataset_bootstrap")
test_that("ASCETIC produces correct output", {
    expect_equal(names(res_example_phylogenies_dataset_bootstrap),c("dataset","ccf_dataset","ranking_estimate","poset","inference"))
})
