#' @name dataset_example_single_samples
#' @title Binary matrix for 10 hepatocellular carcinoma sequencing samples.
#' @description Binary matrix for 10 hepatocellular carcinoma sequencing samples from Nguyen, Bastien, et al., 2022.
#' @docType data
#' @usage data(dataset_example_single_samples)
#' @format Binary matrix where rows are samples and columns are mutations.
#' @source Nguyen, Bastien, et al. "Genomic characterization of metastatic patterns from prospective clinical sequencing of 25,000 patients." Cell 185.3 (2022): 563-575.
#' @return Binary matrix where rows are samples and columns are mutations.
NULL

#' @name ccf_dataset_example_single_samples
#' @title Cancer cell fractions matrix for 10 hepatocellular carcinoma sequencing samples.
#' @description Cancer cell fractions matrix for 10 hepatocellular carcinoma sequencing samples from Nguyen, Bastien, et al., 2022.
#' @docType data
#' @usage data(ccf_dataset_example_single_samples)
#' @format Matrix of cancer cell fractions where rows are samples and columns are mutations.
#' @source Nguyen, Bastien, et al. "Genomic characterization of metastatic patterns from prospective clinical sequencing of 25,000 patients." Cell 185.3 (2022): 563-575.
#' @return Matrix of cancer cell fractions where rows are samples and columns are mutations.
NULL

#' @name vaf_dataset_example_single_samples
#' @title R data.frame for 10 hepatocellular carcinoma sequencing samples.
#' @description R data.frame for 10 hepatocellular carcinoma sequencing samples from Nguyen, Bastien, et al., 2022.
#' @docType data
#' @usage data(vaf_dataset_example_single_samples)
#' @format R data.frame with 8 columns: 1) SAMPLE_ID; 2) GENE_ID; 3) REF_COUNT; 4) ALT_COUNT; 5) COPY_NUMBER; 6) NORMAL_PLOIDY; 7) VAF_ESTIMATE; 8) CCF_ESTIMATE.
#' @source Nguyen, Bastien, et al. "Genomic characterization of metastatic patterns from prospective clinical sequencing of 25,000 patients." Cell 185.3 (2022): 563-575.
#' @return R data.frame with 8 columns: 1) SAMPLE_ID; 2) GENE_ID; 3) REF_COUNT; 4) ALT_COUNT; 5) COPY_NUMBER; 6) NORMAL_PLOIDY; 7) VAF_ESTIMATE; 8) CCF_ESTIMATE.
NULL

#' @name dataset_example_phylogenies
#' @title Binary matrix for 10 myeloid malignancies samples sequenced with Tapestri.
#' @description Binary matrix for 10 myeloid malignancies samples sequenced with Tapestri from Miles, Linde A., et al., 2020.
#' @docType data
#' @usage data(dataset_example_phylogenies)
#' @format Binary matrix where rows are samples and columns are mutations.
#' @source Miles, Linde A., et al. "Single-cell mutation analysis of clonal evolution in myeloid malignancies." Nature 587.7834 (2020): 477-482.
#' @return Binary matrix where rows are samples and columns are mutations.
NULL

#' @name models_phylogenies
#' @title Phylogenetic models estimated for each sample included in dataset.
#' @description Phylogenetic models estimated for each sample included in dataset from Miles, Linde A., et al., 2020.
#' @docType data
#' @usage data(models_phylogenies)
#' @format R list of 10 elements providing the phylogenetic models estimated for the samples in the dataset.
#' @source Miles, Linde A., et al. "Single-cell mutation analysis of clonal evolution in myeloid malignancies." Nature 587.7834 (2020): 477-482.
#' @return R list of 10 elements providing the phylogenetic models estimated for the samples in the dataset.
NULL

#' @name res_example_single_samples_resampling
#' @title Example of results obtained by running the function perform.ascetic.ccf.dataset.resampling on dataset_example_single_samples.
#' @description Example of results by the function perform.ascetic.ccf.dataset.resampling on dataset_example_single_samples from Nguyen, Bastien, et al., 2022.
#' @docType data
#' @usage data(res_example_single_samples_resampling)
#' @format Results obtained by running the function perform.ascetic.ccf.dataset.resampling on dataset_example_single_samples.
#' @return Results obtained by running the function perform.ascetic.ccf.dataset.resampling on dataset_example_single_samples.
NULL

#' @name res_example_phylogenies_dataset_bootstrap
#' @title Example of results obtained by running the function perform.ascetic.phylogenies.dataset.bootstrap on dataset_example_phylogenies
#' @description Example of results obtained by running the function perform.ascetic.phylogenies.dataset.bootstrap on dataset_example_phylogenies from Miles, Linde A., et al., 2020.
#' @docType data
#' @usage data(res_example_phylogenies_dataset_bootstrap)
#' @format Results by perform.ascetic.phylogenies.dataset.bootstrap on dataset_example_phylogenies.
#' @return Results by perform.ascetic.phylogenies.dataset.bootstrap on dataset_example_phylogenies.
NULL
