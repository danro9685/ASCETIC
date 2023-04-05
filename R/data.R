#' @name datasetExampleSingleSamples
#' @title Binary matrix for 10 hepatocellular carcinoma sequencing samples.
#' @description Binary matrix for 10 hepatocellular carcinoma sequencing samples from Nguyen, Bastien, et al., 2022.
#' @docType data
#' @usage data(datasetExampleSingleSamples)
#' @format Binary matrix where rows are samples and columns are mutations.
#' @source Nguyen, Bastien, et al. "Genomic characterization of metastatic patterns from prospective clinical sequencing of 25,000 patients." Cell 185.3 (2022): 563-575.
#' @return Binary matrix where rows are samples and columns are mutations.
NULL

#' @name ccfDatasetExampleSingleSamples
#' @title Cancer cell fractions matrix for 10 hepatocellular carcinoma sequencing samples.
#' @description Cancer cell fractions matrix for 10 hepatocellular carcinoma sequencing samples from Nguyen, Bastien, et al., 2022.
#' @docType data
#' @usage data(ccfDatasetExampleSingleSamples)
#' @format Matrix of cancer cell fractions where rows are samples and columns are mutations.
#' @source Nguyen, Bastien, et al. "Genomic characterization of metastatic patterns from prospective clinical sequencing of 25,000 patients." Cell 185.3 (2022): 563-575.
#' @return Matrix of cancer cell fractions where rows are samples and columns are mutations.
NULL

#' @name vafDatasetExampleSingleSamples
#' @title R data.frame for 10 hepatocellular carcinoma sequencing samples.
#' @description R data.frame for 10 hepatocellular carcinoma sequencing samples from Nguyen, Bastien, et al., 2022.
#' @docType data
#' @usage data(vafDatasetExampleSingleSamples)
#' @format R data.frame with 8 columns: 1) SAMPLE_ID; 2) GENE_ID; 3) REF_COUNT; 4) ALT_COUNT; 5) COPY_NUMBER; 6) NORMAL_PLOIDY; 7) VAF_ESTIMATE; 8) CCF_ESTIMATE.
#' @source Nguyen, Bastien, et al. "Genomic characterization of metastatic patterns from prospective clinical sequencing of 25,000 patients." Cell 185.3 (2022): 563-575.
#' @return R data.frame with 8 columns: 1) SAMPLE_ID; 2) GENE_ID; 3) REF_COUNT; 4) ALT_COUNT; 5) COPY_NUMBER; 6) NORMAL_PLOIDY; 7) VAF_ESTIMATE; 8) CCF_ESTIMATE.
NULL

#' @name datasetExamplePhylogenies
#' @title Binary matrix for 10 myeloid malignancies samples sequenced with Tapestri.
#' @description Binary matrix for 10 myeloid malignancies samples sequenced with Tapestri from Miles, Linde A., et al., 2020.
#' @docType data
#' @usage data(datasetExamplePhylogenies)
#' @format Binary matrix where rows are samples and columns are mutations.
#' @source Miles, Linde A., et al. "Single-cell mutation analysis of clonal evolution in myeloid malignancies." Nature 587.7834 (2020): 477-482.
#' @return Binary matrix where rows are samples and columns are mutations.
NULL

#' @name modelsPhylogenies
#' @title Phylogenetic models estimated for each sample included in dataset.
#' @description Phylogenetic models estimated for each sample included in dataset from Miles, Linde A., et al., 2020.
#' @docType data
#' @usage data(modelsPhylogenies)
#' @format R list of 10 elements providing the phylogenetic models estimated for the samples in the dataset.
#' @source Miles, Linde A., et al. "Single-cell mutation analysis of clonal evolution in myeloid malignancies." Nature 587.7834 (2020): 477-482.
#' @return R list of 10 elements providing the phylogenetic models estimated for the samples in the dataset.
NULL

#' @name resExampleSingleSamplesResampling
#' @title Example of results obtained by running the function perform.ascetic.ccf.dataset.resampling on datasetExampleSingleSamples.
#' @description Example of results by the function perform.ascetic.ccf.dataset.resampling on datasetExampleSingleSamples from Nguyen, Bastien, et al., 2022.
#' @docType data
#' @usage data(resExampleSingleSamplesResampling)
#' @format Results obtained by running the function perform.ascetic.ccf.dataset.resampling on datasetExampleSingleSamples.
#' @return Results obtained by running the function perform.ascetic.ccf.dataset.resampling on datasetExampleSingleSamples.
NULL

#' @name resExamplePhylogeniesDatasetBootstrap
#' @title Example of results obtained by running the function perform.ascetic.phylogenies.dataset.bootstrap on datasetExamplePhylogenies.
#' @description Example of results obtained by running the function perform.ascetic.phylogenies.dataset.bootstrap on datasetExamplePhylogenies from Miles, Linde A., et al., 2020.
#' @docType data
#' @usage data(resExamplePhylogeniesDatasetBootstrap)
#' @format Results by perform.ascetic.phylogenies.dataset.bootstrap on datasetExamplePhylogenies.
#' @return Results by perform.ascetic.phylogenies.dataset.bootstrap on datasetExamplePhylogenies.
NULL

#' @name dataExamples
#' @title List data structure to run examples
#' @description List data structure to run examples.
#' @docType data
#' @usage data(dataExamples)
#' @format List data structure to run examples.
#' @source List data structure to run examples.
#' @return List data structure to run examples.
NULL

#' @name amlExample
#' @title Example dataset for 366 patients with Acute Myeloid Leukaemia.
#' @description 366 patients with Acute Myeloid Leukaemia obtained from Tyner, Jeffrey W., et al. 2018.
#' @docType data
#' @usage data(amlExample)
#' @format List of two elements, i.e., survival data and evolutionary steps extracted by ASCETIC.
#' @source Tyner, Jeffrey W., et al. "Functional genomic landscape of acute myeloid leukaemia." Nature 562.7728 (2018): 526-531.
#' @return List of two elements, i.e., survival data and evolutionary steps extracted by ASCETIC.
NULL
