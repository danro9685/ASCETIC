#' Perform the ASCETIC inference framework on single samples (using CCF) datasets
#' with re-sampling for a robust estimation of the agony ranking.
#'
#' @examples
#' set.seed(12345)
#' data(datasetExampleSingleSamples)
#' data(ccfDatasetExampleSingleSamples)
#' data(vafDatasetExampleSingleSamples)
#' resExampleSingleSamplesResampling <- asceticCCFResampling(
#'                                                dataset = datasetExampleSingleSamples,
#'                                                ccfDataset = ccfDatasetExampleSingleSamples,
#'                                                vafDataset = vafDatasetExampleSingleSamples,
#'                                                nsampling = 5,
#'                                                regularization = "aic",
#'                                                command = "hc",
#'                                                restarts = 0 )
#'
#' @title asceticCCFResampling
#' @param dataset Binary matrix where rows are samples and columns are mutations.
#' Each cell of the matrix is 1 if the related mutation was observed in the sample; 0 otherwise.
#' Values reported in the variable named dataset must be consistent with the ones reported in ccfDataset and vafDataset.
#' @param ccfDataset Matrix where rows are samples and columns are mutations.
#' Each cell of the matrix is the cancer cell fraction (CCF) estimated for the related mutation when observed in the sample.
#' The CCF value is 0 if the mutation was not observed in the sample.
#' Values reported in ccfDataset must be consistent with the ones reported in dataset and vafDataset.
#' @param vafDataset R data.frame with 8 columns: 1) SAMPLE_ID, sample name.
#'                                                2) GENE_ID, gene name.
#'                                                3) REF_COUNT, total counts for reference allele.
#'                                                4) ALT_COUNT, total counts for alternate allele.
#'                                                5) COPY_NUMBER, total copy number estimate.
#'                                                6) NORMAL_PLOIDY, ploidy for normal sample; this is either 1 for mutations on sex chromosomes or 2.
#'                                                7) VAF_ESTIMATE, variant allele frequency (VAF) estimate.
#'                                                8) CCF_ESTIMATE, cancer cell fraction (CCF) estimate.
#' Values reported in vafDataset must be consistent with the ones reported in dataset and ccfDataset.
#' @param nsampling Number of re-sampling to be performed for a robust estimation of the agony ranking.
#' Higher values lead to improved estimates, but require higher computational burden; default value is 100.
#' @param regularization Regularization to be used for the maximum likelihood estimation.
#' Possible values are aic for the Akaike information criterion and bic for the Bayesian information criterion.
#' For the complete list of options, we refer to the manual of the bnlearn package.
#' @param command Optimization technique to be used for maximum likelihood estimation.
#' Valid values are either hc for Hill Climbing or tabu for Tabu Search.
#' @param restarts Number of restarts to be performed during the maximum likelihood estimation when Hill Climbing optimization technique is used.
#' Higher values lead to improved estimates, but require higher computational burden; default value is 10.
#' This parameter is ignored if tabu search is selected.
#' @return A list of 5 elements: 1) dataset, input dataset.
#'                               2) ccfDataset, input ccfDataset.
#'                               3) rankingEstimate, ranking among mutations estimated by agony.
#'                                  Lower rankings correspond to early mutations.
#'                               4) poset, partially order set among mutations estimated by ASCETIC from the agony ranking.
#'                               5) inference, inferred ASCETIC evolutionary model for each selected regularization.
#' @export asceticCCFResampling
#' @import bnlearn
#' @import Rcpp
#' @importFrom stats rbeta
#' @importFrom utils read.table write.table
#' @useDynLib ASCETIC, .registration = TRUE
#' @exportPattern "^[[:alpha:]]+"
#'
asceticCCFResampling <- function(dataset,
           ccfDataset,
           vafDataset,
           nsampling = 100,
           regularization = c("aic", "bic"),
           command = "hc",
           restarts = 10) {

    # estimate the probability raising models
    prModel <- .estimatePrModelSingleSamples(ccfDataset)
    prNull <- .estimatePrNull(dataset)
    
    # estimate agony ranking nsampling times by resampling to obtain a robust and stable estimation
    cat(0, "\n")
    resamplingAgonyRankingEstimate <-
      matrix(0, nrow = ncol(ccfDataset), ncol = 2)
    rownames(resamplingAgonyRankingEstimate) <- colnames(ccfDataset)
    colnames(resamplingAgonyRankingEstimate) <- c("variable", "rank")
    resamplingAgonyRankingEstimate[, "variable"] <- 1:ncol(ccfDataset)
    cont <- 0
    for (i in seq_len(nsampling)) {
      # assess temporal priority with minimum agony for nsampling sampling iterations
      currCcfDataset <-
        matrix(0, nrow = length(unique(vafDataset$SAMPLE_ID)), ncol = length(unique(vafDataset$GENE_ID)))
      rownames(currCcfDataset) <- sort(unique(vafDataset$SAMPLE_ID))
      colnames(currCcfDataset) <- sort(unique(vafDataset$GENE_ID))
      for (j in seq_len(nrow(vafDataset))) {
        vafDatasetSamples <-
          rbeta(
            n = 1,
            shape1 = vafDataset$ALT_COUNT[j],
            shape2 = vafDataset$REF_COUNT[j],
            ncp = 0
          )
        currCcfEstimate <- as.numeric(((
          vafDataset$NORMAL_PLOIDY[j]
          + (vafDataset$COPY_NUMBER[j] - vafDataset$NORMAL_PLOIDY[j])
        ) * vafDatasetSamples))
        currCcfDataset[vafDataset$SAMPLE_ID[j], vafDataset$GENE_ID[j]] <-
          currCcfEstimate
      }
      currCcfDataset[which(currCcfDataset < 0)] <- 0
      currCcfDataset[which(currCcfDataset > 1)] <- 1
      # compute a set of total orderings, one per patient
      totalOrderings <- .inferOrderInCcfDataset(currCcfDataset)
      # compute a set of time orderings among each event, given the total orderings
      agonyArcs <- .buildAgonyInput(totalOrderings)
      # estimate a best agony poset given the time orderings
      if (!is.null(agonyArcs)) {
        # estimate a best agony ranking
        agonyRanking <- agony(agonyArcs)
        cont <- cont + 1
        for (j in 1:nrow(agonyRanking)) {
          resamplingAgonyRankingEstimate[agonyRanking[j, 1], "rank"] <-
            resamplingAgonyRankingEstimate[agonyRanking[j, 1], "rank"]
          + agonyRanking[j, 2]
        }
      }
      cat(i / nsampling, "\n")
    }
    
    # the final estimation is the mean (approximated to integer) of the nsampling rankings
    if (cont > 0) {
      resamplingAgonyRankingEstimate[, "rank"] <-
        round(resamplingAgonyRankingEstimate[, "rank"] / cont)
    }
    
    # FIRST: assess temporal priority with minimum agony
    agonyPoset <-
      .buildRankingAdjMatrix(resamplingAgonyRankingEstimate, ncol(ccfDataset))
    # SECOND: verify probability raising
    agonyPoset <- .applyPr(agonyPoset, prModel, prNull)
    
    # perform the inference
    agonyInference <- list()
    for (reg in regularization) {
      agonyInference[[reg]] <-
        .performLikelihoodFit(dataset, agonyPoset, reg, command, restarts)
    }
    
    # create the data structures
    rownames(agonyPoset) <- as.character(1:ncol(dataset))
    colnames(agonyPoset) <- as.character(1:ncol(dataset))
    poset <- agonyPoset
    for (reg in regularization) {
      rownames(agonyInference[[reg]]) <- as.character(1:ncol(dataset))
      colnames(agonyInference[[reg]]) <- as.character(1:ncol(dataset))
    }
    inference <- agonyInference
    
    # create the data structures with the results
    results <- list(
      dataset = dataset,
      ccfDataset = ccfDataset,
      rankingEstimate = resamplingAgonyRankingEstimate,
      poset = poset,
      inference = inference
    )
    
    return(results)
    
}

#' Perform the ASCETIC inference framework on multiple samples (using phylogenies as inputs) datasets
#' with bootstrap for a robust estimation of the agony ranking.
#'
#' @examples
#' set.seed(12345)
#' data(datasetExamplePhylogenies)
#' data(modelsPhylogenies)
#' resExamplePhylogeniesDatasetBootstrap <- asceticPhylogeniesBootstrap(
#'                                                     dataset = datasetExamplePhylogenies,
#'                                                     models = modelsPhylogenies,
#'                                                     nsampling = 5,
#'                                                     regularization = "aic",
#'                                                     command = "hc",
#'                                                     restarts = 0 )
#'
#' @title asceticPhylogeniesBootstrap
#' @param dataset Binary matrix where rows are samples and columns are mutations.
#' Each cell of the matrix is 1 if the related mutation was observed in the sample; 0 otherwise.
#' @param models Phylogenetic models estimated for each sample included in dataset.
#' @param nsampling Number of re-sampling to be performed for a robust estimation of the agony ranking.
#' Higher values lead to improved estimates, but require higher computational burden; default value is 100.
#' @param regularization Regularization to be used for the maximum likelihood estimation.
#' Possible values are aic for the Akaike information criterion and bic for the Bayesian information criterion.
#' For the complete list of options, we refer to the manual of the bnlearn package.
#' @param command Optimization technique to be used for maximum likelihood estimation.
#' Valid values are either hc for Hill Climbing or tabu for Tabu Search.
#' @param restarts Number of restarts to be performed during the maximum likelihood estimation when Hill Climbing optimization technique is used.
#' Higher values lead to improved estimates, but require higher computational burden; default value is 10.
#' This parameter is ignored if tabu search is selected.
#' @return A list of 5 elements: 1) dataset, input dataset.
#'                               2) models, input phylogenetic models
#'                               3) rankingEstimate, ranking among mutations estimated by agony.
#'                                  Lower rankings correspond to early mutations.
#'                               4) poset, partially order set among mutations estimated by ASCETIC from the agony ranking.
#'                               5) inference, inferred ASCETIC evolutionary model for each selected regularization.
#' @export asceticPhylogeniesBootstrap
#' @import bnlearn
#' @import Rcpp
#' @importFrom utils read.table write.table
#' @useDynLib ASCETIC, .registration = TRUE
#' @exportPattern "^[[:alpha:]]+"
#'
asceticPhylogeniesBootstrap <- function(dataset,
           models,
           nsampling = 100,
           regularization = c("aic", "bic"),
           command = "hc",
           restarts = 10) {

    # estimate the probability raising models
    prModel <- .estimatePrModelMultipleSamples(models, colnames(dataset))
    prNull <- .estimatePrNull(dataset)
    
    # estimate agony ranking nsampling times by bootstrap to obtain a robust and stable estimation
    cat(0, "\n")
    bootstrapAgonyRankingEstimate <-
      matrix(0, nrow = ncol(dataset), ncol = 2)
    rownames(bootstrapAgonyRankingEstimate) <- colnames(dataset)
    colnames(bootstrapAgonyRankingEstimate) <- c("variable", "rank")
    bootstrapAgonyRankingEstimate[, "variable"] <- 1:ncol(dataset)
    cont <- 0
    for (i in 1:nsampling) {
      # assess temporal priority with minimum agony for nsampling bootstrap iterations
      currModels <- sample(1:length(models), replace = TRUE)
      currModels <- models[currModels]
      for (j in 1:length(currModels)) {
        currm <- currModels[[j]]
        currModels[[j]] <-
          matrix(0,
                 nrow = length(colnames(dataset)),
                 ncol = length(colnames(dataset)))
        rownames(currModels[[j]]) <- colnames(dataset)
        colnames(currModels[[j]]) <- colnames(dataset)
        currModels[[j]][rownames(currm), colnames(currm)] <- currm
      }
      # compute a set of time orderings among each event, given a set of models inferred from multiple samples
      agonyArcs <-
        .buildAgonyInputMultipleSamples(currModels, colnames(dataset))
      # estimate a best agony poset given the time orderings
      if (!is.null(agonyArcs)) {
        # estimate a best agony ranking and save the results to file
        agonyRanking <- agony(agonyArcs)
        cont <- cont + 1
        for (j in 1:nrow(agonyRanking)) {
          bootstrapAgonyRankingEstimate[agonyRanking[j, 1], "rank"] <-
            bootstrapAgonyRankingEstimate[agonyRanking[j, 1], "rank"]
          + agonyRanking[j, 2]
        }
      }
      cat(i / nsampling, "\n")
    }
    
    # the final estimation is the mean (approximated to integer) of the nsampling rankings
    if (cont > 0) {
      bootstrapAgonyRankingEstimate[, "rank"] <-
        round(bootstrapAgonyRankingEstimate[, "rank"] / cont)
    }
    
    # FIRST: assess temporal priority with minimum agony
    agonyPoset <-
      .buildRankingAdjMatrix(bootstrapAgonyRankingEstimate, ncol(dataset))
    # SECOND: verify probability raising
    agonyPoset <- .applyPr(agonyPoset, prModel, prNull)
    
    # perform the inference
    agonyInference <- list()
    for (reg in regularization) {
      agonyInference[[reg]] <-
        .performLikelihoodFit(dataset, agonyPoset, reg, command, restarts)
    }
    
    # create the data structures with the posets and the inference results
    rownames(agonyPoset) <- as.character(1:ncol(dataset))
    colnames(agonyPoset) <- as.character(1:ncol(dataset))
    poset <- agonyPoset
    for (reg in regularization) {
      rownames(agonyInference[[reg]]) <- as.character(1:ncol(dataset))
      colnames(agonyInference[[reg]]) <- as.character(1:ncol(dataset))
    }
    inference <- agonyInference
    
    # Create the data structures with the results
    results <- list(
      dataset = dataset,
      models = models,
      rankingEstimate = bootstrapAgonyRankingEstimate,
      poset = poset,
      inference = inference
    )
    
    return(results)
    
}

#' Perform the ASCETIC inference framework on single samples (using CCF) datasets.
#'
#' @examples
#' set.seed(12345)
#' data(datasetExampleSingleSamples)
#' data(ccfDatasetExampleSingleSamples)
#' resExampleSingleSamples = asceticCCF(
#'                               dataset = datasetExampleSingleSamples,
#'                               ccfDataset = ccfDatasetExampleSingleSamples,
#'                               regularization = "aic",
#'                               command = "hc",
#'                               restarts = 0 )
#'
#' @title asceticCCF
#' @param dataset Binary matrix where rows are samples and columns are mutations.
#' Each cell of the matrix is 1 if the related mutation was observed in the sample; 0 otherwise.
#' Values reported in the variable named dataset must be consistent with the ones reported in ccfDataset.
#' @param ccfDataset Matrix where rows are samples and columns are mutations.
#' Each cell of the matrix is the cancer cell fraction (CCF) estimated for the related mutation when observed in the sample.
#' The CCF value is 0 if the mutation was not observed in the sample.
#' Values reported in ccfDataset must be consistent with the ones reported in dataset.
#' @param regularization Regularization to be used for the maximum likelihood estimation.
#' Possible values are aic for the Akaike information criterion and bic for the Bayesian information criterion.
#' For the complete list of options, we refer to the manual of the bnlearn package.
#' @param command Optimization technique to be used for maximum likelihood estimation.
#' Valid values are either hc for Hill Climbing or tabu for Tabu Search.
#' @param restarts Number of restarts to be performed during the maximum likelihood estimation when Hill Climbing optimization technique is used.
#' Higher values lead to improved estimates, but require higher computational burden; default value is 10.
#' This parameter is ignored if tabu search is selected.
#' @return A list of 4 elements: 1) dataset, input dataset.
#'                               2) ccfDataset, input ccfDataset.
#'                               3) poset, partially order set among mutations estimated by ASCETIC from the agony ranking.
#'                               4) inference, inferred ASCETIC evolutionary model for each selected regularization.
#' @export asceticCCF
#' @import bnlearn
#' @import Rcpp
#' @importFrom utils read.table write.table
#' @useDynLib ASCETIC, .registration = TRUE
#' @exportPattern "^[[:alpha:]]+"
#'
asceticCCF <- function(dataset,
           ccfDataset,
           regularization = c("aic", "bic"),
           command = "hc",
           restarts = 10) {

    # estimate the probability raising models
    prModel <- .estimatePrModelSingleSamples(ccfDataset)
    prNull <- .estimatePrNull(dataset)
    
    # compute the agony based poset
    # First: assess temporal priority with minimum agony
    agonyPoset <- .estimateAgonyPosetSingleSamples(ccfDataset)
    # Second: verify probability raising
    agonyPoset <- .applyPr(agonyPoset, prModel, prNull)
    
    # perform the inference
    agonyInference = list()
    for (reg in regularization) {
      agonyInference[[reg]] = .performLikelihoodFit(dataset, agonyPoset, reg, command, restarts)
    }
    
    # create the data structures with the posets and the inference results
    rownames(agonyPoset) <- as.character(1:ncol(dataset))
    colnames(agonyPoset) <- as.character(1:ncol(dataset))
    poset <- agonyPoset
    for (i in regularization) {
      rownames(agonyInference[[i]]) <- as.character(1:ncol(dataset))
      colnames(agonyInference[[i]]) <- as.character(1:ncol(dataset))
    }
    
    # create the data structures with the results
    results <-
      list(
        dataset = dataset,
        ccfDataset = ccfDataset,
        poset = poset,
        inference = agonyInference
      )
    
    return(results)
    
}

#' Perform the ASCETIC inference framework on multiple samples (using phylogenies as inputs) datasets.
#'
#' @examples
#' set.seed(12345)
#' data(datasetExamplePhylogenies)
#' data(modelsPhylogenies)
#' resExamplePhylogeniesDataset = asceticPhylogenies(
#'                                          dataset = datasetExamplePhylogenies,
#'                                          models = modelsPhylogenies,
#'                                          regularization = "aic",
#'                                          command = "hc",
#'                                          restarts = 0)
#'
#' @title asceticPhylogenies
#' @param dataset Binary matrix where rows are samples and columns are mutations.
#' Each cell of the matrix is 1 if the related mutation was observed in the sample; 0 otherwise.
#' @param models Phylogenetic models estimated for each sample included in dataset.
#' @param regularization Regularization to be used for the maximum likelihood estimation.
#' Possible values are aic for the Akaike information criterion and bic for the Bayesian information criterion.
#' For the complete list of options, we refer to the manual of the bnlearn package.
#' @param command Optimization technique to be used for maximum likelihood estimation.
#' Valid values are either hc for Hill Climbing or tabu for Tabu Search.
#' @param restarts Number of restarts to be performed during the maximum likelihood estimation when Hill Climbing optimization technique is used.
#' Higher values lead to improved estimates, but require higher computational burden; default value is 10.
#' This parameter is ignored if tabu search is selected.
#' @return A list of 4 elements: 1) dataset, input dataset.
#'                               2) models, input phylogenetic models
#'                               3) poset, partially order set among mutations estimated by ASCETIC from the agony ranking.
#'                               4) inference, inferred ASCETIC evolutionary model for each selected regularization.
#' @export asceticPhylogenies
#' @import bnlearn
#' @import Rcpp
#' @importFrom utils read.table write.table
#' @useDynLib ASCETIC, .registration = TRUE
#' @exportPattern "^[[:alpha:]]+"
#'
asceticPhylogenies <- function(dataset,
           models,
           regularization = c("aic", "bic"),
           command = "hc",
           restarts = 10) {

    # estimate the probability raising models
    prModel <-
      .estimatePrModelMultipleSamples(models, colnames(dataset))
    prNull <- .estimatePrNull(dataset)
    
    # compute the agony-based poset
    # First: assess temporal priority with minimum agony
    agonyPoset <-
      .estimateAgonyPosetMultipleSamples(models, colnames(dataset))
    # Second: verify probability raising
    agonyPoset <- .applyPr(agonyPoset, prModel, prNull)
    
    # perform the inference
    agonyInference = list()
    for (reg in regularization) {
      agonyInference[[reg]] = .performLikelihoodFit(dataset, agonyPoset, reg, command, restarts)
    }
    
    # create the data structures with the posets and the inference results
    rownames(agonyPoset) <- as.character(1:ncol(dataset))
    colnames(agonyPoset) <- as.character(1:ncol(dataset))
    poset <- agonyPoset
    for (i in regularization) {
      rownames(agonyInference[[i]]) <- as.character(1:ncol(dataset))
      colnames(agonyInference[[i]]) <- as.character(1:ncol(dataset))
    }

    # create the data structures with the results
    results <-
      list(
        dataset = dataset,
        models = models,
        poset = poset,
        inference = agonyInference
      )
    
    return(results)
    
}
