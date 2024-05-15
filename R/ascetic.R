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
    resamplingAgonyRankingEstimate[, "variable"] <- seq_len(ncol(ccfDataset))
    cont <- 0
    
    # warning: side-effects of for loops
    # no possible refactoring with apply-style
    for (i in seq_len(nsampling)) {
      # assess temporal priority with minimum agony for nsampling sampling iterations
      currCcfDataset <-
        matrix(0, nrow = length(unique(vafDataset$SAMPLE_ID)), ncol = length(unique(vafDataset$GENE_ID)))
      rownames(currCcfDataset) <- sort(unique(vafDataset$SAMPLE_ID))
      colnames(currCcfDataset) <- sort(unique(vafDataset$GENE_ID))
      for (j in seq_len(nrow(vafDataset))) {
        vafDatasetSamples <-
          rbeta(
            n <- 1,
            shape1 <- vafDataset$ALT_COUNT[j],
            shape2 <- vafDataset$REF_COUNT[j],
            ncp <- 0
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
        agonyRanking <- agony(agonyArcs, .get_seed())
        cont <- cont + 1
        for (j in seq_len(nrow(agonyRanking))) {
          resamplingAgonyRankingEstimate[agonyRanking[j, 1], "rank"] <-
            resamplingAgonyRankingEstimate[agonyRanking[j, 1], "rank"] + agonyRanking[j, 2]
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
    agonyInference <- setNames(
      lapply(regularization, function(reg) {
        result <- .performLikelihoodFit(dataset, agonyPoset, reg, command, restarts)
        rownames(result) <- colnames(result) <- as.character(seq_len(ncol(dataset)))
        return(result)
      }),
      regularization
    )
    
    # create the data structures
    rownames(agonyPoset) <- as.character(seq_len(ncol(dataset)))
    colnames(agonyPoset) <- as.character(seq_len(ncol(dataset)))
    poset <- agonyPoset
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

#' Perform the ASCETIC inference framework on multiple samples (using mutational trees as inputs) datasets
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
#' @param models Phylogenetic models estimated for each sample included in dataset given as input 
#' to ASCETIC as a adjacency matrices representing mutational trees.
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
#'                               2) models, input mutational trees
#'                               3) rankingEstimate, ranking among mutations estimated by agony.
#'                                  Lower rankings correspond to early mutations.
#'                               4) poset, partially order set among mutations estimated by ASCETIC from the agony ranking.
#'                               5) inference, inferred ASCETIC evolutionary model for each selected regularization.
#' @export asceticPhylogeniesBootstrap
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
  bootstrapAgonyRankingEstimate[, "variable"] <- seq_len(ncol(dataset))
  cont <- 0

  # use lapply for bootstrap iterations
  bootstrap_results <- lapply(seq_len(nsampling), function(i) {
    # assess temporal priority with minimum agony for nsampling bootstrap iterations
    currModels <- sample(seq_len(length(models)), replace = TRUE)
    currModels <- models[currModels]
    for (j in seq_len(length(currModels))) {
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
    agonyArcs <- .buildAgonyInputMultipleSamples(currModels, colnames(dataset))
    # estimate a best agony poset given the time orderings
    if (!is.null(agonyArcs)) {
      # estimate a best agony ranking and save the results to file
      agonyRanking <- agony(agonyArcs, .get_seed())
      cont <- cont + 1
      for (j in seq_len(nrow(agonyRanking))) {
        bootstrapAgonyRankingEstimate[agonyRanking[j, 1], "rank"] <-
          bootstrapAgonyRankingEstimate[agonyRanking[j, 1], "rank"] + agonyRanking[j, 2]
      }
    }
    cat(i / nsampling, "\n")
    # return the results for this iteration
    return(list(
      currModels = currModels,
      agonyArcs = agonyArcs
    ))
  })

  # combine the results from the lapply output
  for (i in seq_len(nsampling)) {
    currModels <- bootstrap_results[[i]]$currModels
    agonyArcs <- bootstrap_results[[i]]$agonyArcs

    # the final estimation is the mean (approximated to integer) of the nsampling rankings
    if (!is.null(agonyArcs)) {
      # estimate a best agony ranking and save the results to file
      agonyRanking <- agony(agonyArcs, .get_seed())
      cont <- cont + 1
      for (j in seq_len(nrow(agonyRanking))) {
        bootstrapAgonyRankingEstimate[agonyRanking[j, 1], "rank"] <-
          bootstrapAgonyRankingEstimate[agonyRanking[j, 1], "rank"] + agonyRanking[j, 2]
      }
    }
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
  agonyInference <- setNames(
    lapply(regularization, function(reg) {
      result <- .performLikelihoodFit(dataset, agonyPoset, reg, command, restarts)
      rownames(result) <- colnames(result) <- as.character(seq_len(ncol(dataset)))
      return(result)
    }),
    regularization
  )

  # create the data structures with the posets and the inference results
  rownames(agonyPoset) <- as.character(seq_len(ncol(dataset)))
  colnames(agonyPoset) <- as.character(seq_len(ncol(dataset)))
  poset <- agonyPoset
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
#' resExampleSingleSamples <- asceticCCF(
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
    agonyInference <- setNames(
      lapply(regularization, function(reg) {
        result <- .performLikelihoodFit(dataset, agonyPoset, reg, command, restarts)
        rownames(result) <- colnames(result) <- as.character(seq_len(ncol(dataset)))
        return(result)
      }),
      regularization
    )
    
    # create the data structures with the posets and the inference results
    rownames(agonyPoset) <- as.character(seq_len(ncol(dataset)))
    colnames(agonyPoset) <- as.character(seq_len(ncol(dataset)))
    poset <- agonyPoset
    
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

#' Perform the ASCETIC inference framework on multiple samples (using mutational trees as inputs) datasets.
#'
#' @examples
#' set.seed(12345)
#' data(datasetExamplePhylogenies)
#' data(modelsPhylogenies)
#' resExamplePhylogeniesDataset <- asceticPhylogenies(
#'                                          dataset = datasetExamplePhylogenies,
#'                                          models = modelsPhylogenies,
#'                                          regularization = "aic",
#'                                          command = "hc",
#'                                          restarts = 0)
#'
#' @title asceticPhylogenies
#' @param dataset Binary matrix where rows are samples and columns are mutations.
#' Each cell of the matrix is 1 if the related mutation was observed in the sample; 0 otherwise.
#' @param models Phylogenetic models estimated for each sample included in dataset given as input 
#' to ASCETIC as a adjacency matrices representing mutational trees.
#' @param regularization Regularization to be used for the maximum likelihood estimation.
#' Possible values are aic for the Akaike information criterion and bic for the Bayesian information criterion.
#' For the complete list of options, we refer to the manual of the bnlearn package.
#' @param command Optimization technique to be used for maximum likelihood estimation.
#' Valid values are either hc for Hill Climbing or tabu for Tabu Search.
#' @param restarts Number of restarts to be performed during the maximum likelihood estimation when Hill Climbing optimization technique is used.
#' Higher values lead to improved estimates, but require higher computational burden; default value is 10.
#' This parameter is ignored if tabu search is selected.
#' @return A list of 4 elements: 1) dataset, input dataset.
#'                               2) models, input mutational trees
#'                               3) poset, partially order set among mutations estimated by ASCETIC from the agony ranking.
#'                               4) inference, inferred ASCETIC evolutionary model for each selected regularization.
#' @export asceticPhylogenies
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
    agonyInference <- setNames(
      lapply(regularization, function(reg) {
        result <- .performLikelihoodFit(dataset, agonyPoset, reg, command, restarts)
        rownames(result) <- colnames(result) <- as.character(seq_len(ncol(dataset)))
        return(result)
      }),
      regularization
    )
    
    # create the data structures with the posets and the inference results
    rownames(agonyPoset) <- as.character(seq_len(ncol(dataset)))
    colnames(agonyPoset) <- as.character(seq_len(ncol(dataset)))
    poset <- agonyPoset

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

#' Perform the assessment via cross-validation of the model inferred by the ASCETIC framework on 
#' single samples (using CCF) datasets.
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
#'                                                nsampling = 3,
#'                                                regularization = c("aic","bic"),
#'                                                command = "hc",
#'                                                restarts = 0 )
#' resExampleCCFAssessment <- asceticCCFAssessment(
#'                                        inference = resExampleSingleSamplesResampling,
#'                                        niterations = 3,
#'                                        vafDataset = vafDatasetExampleSingleSamples,
#'                                        nsampling = 3)
#'
#' @title asceticCCFAssessment
#' @param inference Model inferred with ASCETIC using either the function asceticCCFResampling or the function asceticCCF
#' @param niterations Number of cross-validation iterations to be performed for a robust assessment of ASCETIC model.
#' Higher values lead to improved estimates, but require higher computational burden; default value is 100.
#' @param vafDataset R data.frame with 8 columns: 1) SAMPLE_ID, sample name.
#'                                                2) GENE_ID, gene name.
#'                                                3) REF_COUNT, total counts for reference allele.
#'                                                4) ALT_COUNT, total counts for alternate allele.
#'                                                5) COPY_NUMBER, total copy number estimate.
#'                                                6) NORMAL_PLOIDY, ploidy for normal sample; this is either 1 for mutations on sex chromosomes or 2.
#'                                                7) VAF_ESTIMATE, variant allele frequency (VAF) estimate.
#'                                                8) CCF_ESTIMATE, cancer cell fraction (CCF) estimate.
#' Values reported in vafDataset must be consistent with the ones reported in dataset and ccfDataset.
#' If vafDataset is not provided, sampling cannot be performed.
#' @param nsampling Number of re-sampling to be performed for a robust estimation of the agony ranking.
#' Higher values lead to improved estimates, but require higher computational burden; default value is 100.
#' @return A list of 3 elements for which the estimate by cross-validation is performed: 
#'                               1) rankingEstimate, ranking among mutations estimated by agony.
#'                                  Lower rankings correspond to early mutations. This is returned only if nsampling > 0.
#'                               2) poset, partially order set among mutations estimated by ASCETIC from the agony ranking.
#'                               3) inference, inferred ASCETIC evolutionary model for each selected regularization.
#' @export asceticCCFAssessment
#'
asceticCCFAssessment <- function(inference, niterations = 100, vafDataset = NULL, nsampling = 100) {

    # structures to save the results
    niterations <- max(1,niterations) # at least 1 iteration needs to be performed
    ranking_estimate <- matrix(NA, nrow = ncol(inference$dataset), ncol = niterations)
    rownames(ranking_estimate) <- colnames(inference$dataset)
    colnames(ranking_estimate) <- paste0("Iteration ",1:ncol(ranking_estimate))
    poset <- matrix(0, nrow = ncol(inference$dataset), ncol = ncol(inference$dataset))
    rownames(poset) <- colnames(inference$dataset)
    colnames(poset) <- colnames(inference$dataset)
    reg <- list()
    for(i in names(inference$inference)) {
      reg[[i]] <- poset
    }

    # perform the assessment of the inferred ASCETIC model via cross-validation
    valid_samples <- 1:nrow(inference$dataset)
    cv_size <- round(length(valid_samples) * 0.80) # consider 80% samples at each cross-validation step
    cat("Starting cross-validation...","\n")
    if(is.null(vafDataset)) { # if vafDataset is not provided, I cannot perform sampling
      nsampling <- 0
    }
    for(i in 1:niterations) {

      # build the inputs for the current iteration
      selected_samples <- sample(valid_samples, size = cv_size, replace = FALSE)
      curr_dataset <- inference$dataset[selected_samples,]
      curr_ccf_dataset <- inference$ccfDataset[selected_samples,]
      curr_vaf_dataset <- vafDataset[which(vafDataset$SAMPLE_ID%in%rownames(curr_dataset)),]
      
      # perform the inference for the current iteration
      if(nsampling > 0) {
        cv_ascetic_results <- asceticCCFResampling(dataset = curr_dataset, 
                                            ccfDataset = curr_ccf_dataset, 
                                            vafDataset = curr_vaf_dataset, 
                                            nsampling = nsampling, 
                                            regularization = names(inference$inference), 
                                            command = "hc", 
                                            restarts = 10)
      }
      else {
        cv_ascetic_results <- asceticCCF(dataset = curr_dataset, 
                                      ccfDataset = curr_ccf_dataset, 
                                      regularization = names(inference$inference), 
                                      command = "hc", 
                                      restarts = 10)
      }

      # save the results for the current iteration
      if(nsampling > 0) {
        ranking_estimate[,i] <- as.numeric(cv_ascetic_results$rankingEstimate[,"rank"])
      }
      poset <- (poset + cv_ascetic_results$poset)
      for(j in names(inference$inference)) {
        reg[[j]] <- (reg[[j]] + cv_ascetic_results$inference[[j]])
      }
      cat("Cross-validation progress: ",(i/niterations),"\n")
        
    }

    # save the results
    results <- list()
    if(nsampling > 0) {
      results[["ranking"]] <- (rowSums(ranking_estimate) / niterations)
    }
    else {
      results[["ranking"]] <- NA
    }
    results[["poset"]] <- (poset / niterations)
    for(i in names(inference$inference)) {
      results[["inference"]][[i]] <- (reg[[i]] / niterations)
    }
      
    return(results)

}

#' Perform the assessment via cross-validation of the model inferred by the ASCETIC framework on 
#' multiple samples (using mutational trees as inputs) datasets.
#'
#' @examples
#' set.seed(12345)
#' data(datasetExamplePhylogenies)
#' data(modelsPhylogenies)
#' resExamplePhylogeniesDatasetBootstrap <- asceticPhylogeniesBootstrap(
#'                                                     dataset = datasetExamplePhylogenies,
#'                                                     models = modelsPhylogenies,
#'                                                     nsampling = 3,
#'                                                     regularization = c("aic","bic"),
#'                                                     command = "hc",
#'                                                     restarts = 0 )
#' resExamplePhylogeniesAssessment <- asceticPhylogeniesAssessment(
#'                                           inference = resExamplePhylogeniesDatasetBootstrap,
#'                                           niterations = 3,
#'                                           nsampling = 3)
#'
#' @title asceticPhylogeniesAssessment
#' @param inference Model inferred with ASCETIC using either the function asceticPhylogeniesBootstrap or the function asceticPhylogenies.
#' @param niterations Number of cross-validation iterations to be performed for a robust assessment of ASCETIC model.
#' Higher values lead to improved estimates, but require higher computational burden; default value is 100.
#' @param nsampling Number of re-sampling to be performed for a robust estimation of the agony ranking.
#' Higher values lead to improved estimates, but require higher computational burden; default value is 100.
#' @return A list of 3 elements for which the estimate by cross-validation is performed: 
#'                               1) rankingEstimate, ranking among mutations estimated by agony.
#'                                  Lower rankings correspond to early mutations. This is returned only if nsampling > 0.
#'                               2) poset, partially order set among mutations estimated by ASCETIC from the agony ranking.
#'                               3) inference, inferred ASCETIC evolutionary model for each selected regularization.
#' @export asceticPhylogeniesAssessment
#'
asceticPhylogeniesAssessment <- function(inference, niterations = 100, nsampling = 100) {

    # structures to save the results
    niterations <- max(1,niterations) # at least 1 iteration needs to be performed
    ranking_estimate <- matrix(NA, nrow = ncol(inference$dataset), ncol = niterations)
    rownames(ranking_estimate) <- colnames(inference$dataset)
    colnames(ranking_estimate) <- paste0("Iteration ",1:ncol(ranking_estimate))
    poset <- matrix(0, nrow = ncol(inference$dataset), ncol = ncol(inference$dataset))
    rownames(poset) <- colnames(inference$dataset)
    colnames(poset) <- colnames(inference$dataset)
    reg <- list()
    for(i in names(inference$inference)) {
      reg[[i]] <- poset
    }

    # perform the assessment of the inferred ASCETIC model via cross-validation
    valid_samples <- 1:nrow(inference$dataset)
    cv_size <- round(length(valid_samples) * 0.80) # consider 80% samples at each cross-validation step
    cat("Starting cross-validation...","\n")
    for(i in 1:niterations) {

      # build the inputs for the current iteration
      selected_samples <- sample(valid_samples, size = cv_size, replace = FALSE)
      curr_dataset <- inference$dataset[selected_samples,]
      curr_models <- inference$models[names(inference$models)[which(names(inference$models)%in%rownames(curr_dataset))]]
      
      # perform the inference for the current iteration
      if(nsampling > 0) {
        cv_ascetic_results <- asceticPhylogeniesBootstrap(dataset = curr_dataset, 
                                               models = curr_models, 
                                               nsampling = nsampling, 
                                               regularization = names(inference$inference), 
                                               command = "hc", 
                                               restarts = 10)
      }
      else {
        cv_ascetic_results <- asceticPhylogenies(dataset = curr_dataset, 
                                        models = curr_models, 
                                        regularization = names(inference$inference), 
                                        command = "hc", 
                                        restarts = 10)
      }

      # save the results for the current iteration
      if(nsampling > 0) {
        ranking_estimate[,i] <- as.numeric(cv_ascetic_results$rankingEstimate[,"rank"])
      }
      poset <- (poset + cv_ascetic_results$poset)
      for(j in names(inference$inference)) {
        reg[[j]] <- (reg[[j]] + cv_ascetic_results$inference[[j]])
      }
      cat("Cross-validation progress: ",(i/niterations),"\n")
        
    }

    # save the results
    results <- list()
    if(nsampling > 0) {
      results[["ranking"]] <- (rowSums(ranking_estimate) / niterations)
    }
    else {
      results[["ranking"]] <- NA
    }
    results[["poset"]] <- (poset / niterations)
    for(i in names(inference$inference)) {
      results[["inference"]][[i]] <- (reg[[i]] / niterations)
    }
      
    return(results)

}
