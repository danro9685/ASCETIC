# Estimate the model for the probability raising, i.e., j --> i, P(t_j<t_i | i,j),
# from a dataset of cancer cell fractions.
#
# @examples
# data(dataExamples)
# dataExample <- dataExamples[["estimate.pr.model.single.samples"]]
# res <- .estimatePrModelSingleSamples(ccfDataset = dataExample)
#
# @title .estimatePrModelSingleSamples
# @param ccfDataset Matrix where rows are samples and columns are mutations.
#
.estimatePrModelSingleSamples <- function(ccfDataset) {

  # data structure to save results
  prModelEstimate <-
    matrix(0,
           nrow = ncol(ccfDataset),
           ncol = ncol(ccfDataset))
  
  # consider any pair of events in the dataset
  for (i in seq_len(nrow(prModelEstimate))) {
    for (j in i:nrow(prModelEstimate)) {
      if (i != j) {
        # if i and j are not the same event
        
        # get the subset of the dataset where events i and j are co-occurring
        coOccurringSamples <-
          which(ccfDataset[, i] > 0 & ccfDataset[, j] > 0)
        
        # if there is at least 1 sample where the 2 genes are co-occurring
        if (length(coOccurringSamples) > 0) {
          currI <- as.vector(ccfDataset[coOccurringSamples, i])
          currJ <- as.vector(ccfDataset[coOccurringSamples, j])
          prModelEstimate[i, j] <-
            (length(which(currI > currJ)) / length(coOccurringSamples))
          prModelEstimate[j, i] <-
            (length(which(currJ > currI)) / length(coOccurringSamples))
        }
        else {
          prModelEstimate[i, j] <- 0
          prModelEstimate[j, i] <- 0
        }
        
      }
    }
  }
  
  return(prModelEstimate)
  
}

# Estimate the model for the probability raising, i.e., j --> i, P(t_j<t_i | i,j),
# from a phylogenetic models available for multiple distinct patients.
#
# @examples
# data(dataExamples)
# dataExample <- dataExamples[["estimate.pr.model.multiple.samples"]]
# res <- .estimatePrModelMultipleSamples(models = dataExample$models, events = dataExample$events)
#
# @title .estimatePrModelMultipleSamples
# @param models Phylogenetic models estimated for each sample included in dataset.
# @param events Driver genes to be considered in the inference.
#
.estimatePrModelMultipleSamples <- function(models, events) {

  # data structure to save results
  prModelEstimate <-
    matrix(0, nrow = length(events), ncol = length(events))
  
  # data structures to save the counts for each event
  countsTimeOrderings <-
    matrix(0, nrow = length(events), ncol = length(events))
  countsCoOccurrence <-
    matrix(0, nrow = length(events), ncol = length(events))
  
  # consider any pair of events in the dataset
  for (i in seq_len(length(models))) {
    # compute the time ordering counts, i.e., P(t_j<t_i,i,j)
    currModel <- models[[i]]
    timeOrderings <- which(currModel == 1, arr.ind = TRUE)
    if (nrow(timeOrderings) > 0) {
      timeOrderingsParents <- rownames(currModel)[timeOrderings[, "row"]]
      timeOrderingsChildren <-
        colnames(currModel)[timeOrderings[, "col"]]
      for (j in seq_len(length(timeOrderingsParents))) {
        countsTimeOrderings[which(events == timeOrderingsParents[j]), which(events == timeOrderingsChildren[j])] <-
          countsTimeOrderings[which(events == timeOrderingsParents[j]), which(events == timeOrderingsChildren[j])] + 1
      }
    }
    
    # increment the co-occurrence count, i.e., P(i,j)
    countsCoOccurrence[which(events %in% colnames(currModel)), which(events %in% colnames(currModel))] <-
      countsCoOccurrence[which(events %in% colnames(currModel)), which(events %in% colnames(currModel))] + 1
    
  }
  
  # save the results
  prModelEstimate <- countsTimeOrderings / countsCoOccurrence
  prModelEstimate[which(is.nan(prModelEstimate), arr.ind = TRUE)] <-
    0
  
  return(prModelEstimate)
  
}

# Estimate the null model for the probability raising, i.e., j --> i, P(i, not j).
#
# @examples
# data(dataExamples)
# dataExample <- dataExamples[["estimate.pr.null"]]
# res <- .estimatePrNull(dataExample)
#
# @title .estimatePrNull
# @param dataset Binary matrix where rows are samples and columns are mutations.
#
.estimatePrNull <- function(dataset) {

  # data structure to save results
  prNullEstimate <-
    matrix(0, nrow = ncol(dataset), ncol = ncol(dataset))
  
  # compute marginal and joint probabilities
  probsEstimate <- .estimateProbs(dataset)
  
  # consider any pair of events in the dataset
  for (i in seq_len(nrow(prNullEstimate))) {
    for (j in i:nrow(prNullEstimate)) {
      if (i != j) {
        prNullEstimate[i, j] <-
          probsEstimate$marginalProbs[i, 1] - probsEstimate$jointProbs[i, j]
        prNullEstimate[j, i] <-
          probsEstimate$marginalProbs[j, 1] - probsEstimate$jointProbs[j, i]
      }
    }
  }
  
  return(prNullEstimate)
  
}

# Estimate marginal and joint probabilities from the dataset.
#
# @examples
# data(dataExamples)
# dataExample <- dataExamples[["estimate.probs"]]
# res <- .estimateProbs(dataExample)
#
# @title .estimateProbs
# @param dataset Binary matrix where rows are samples and columns are mutations.
#
.estimateProbs <- function(dataset) {

  # data structure to save results
  pairCount <- matrix(0, nrow = ncol(dataset), ncol = ncol(dataset))
  
  # compute marginal and joint probabilities
  for (i in seq_len(ncol(dataset))) {
    for (j in seq_len(ncol(dataset))) {
      val1 <- dataset[, i]
      val2 <- dataset[, j]
      pairCount[i, j] <- (t(val1) %*% val2)
    }
  }
  marginalProbs <-
    matrix(as.matrix(diag(pairCount) / nrow(dataset)),
           nrow = ncol(dataset),
           ncol = 1)
  jointProbs <- as.matrix(pairCount / nrow(dataset))
  
  # save results
  results <-
    list(marginalProbs = marginalProbs, jointProbs = jointProbs)
  return(results)
  
}

# Estimate the poset based on the best agony ranking inferred from a dataset of cancer cell fractions.
#
# @title .estimateAgonyPosetSingleSamples
# @param ccfDataset Matrix where rows are samples and columns are mutations.
#
.estimateAgonyPosetSingleSamples <- function(ccfDataset) {

  # compute a set of total orderings, one per patient
  totalOrderings <- .inferOrderInCcfDataset(ccfDataset)
  
  # compute a set of time orderings among events, given the total orderings
  agonyArcs <- .buildAgonyInput(totalOrderings)
  
  # estimate the best agony poset given the time orderings
  if (!is.null(agonyArcs)) {
    agonyPoset <- .computeAgonyPoset(agonyArcs, ncol(ccfDataset))
  } else {
    agonyPoset <-
      matrix(1,
             nrow = ncol(ccfDataset),
             ncol = ncol(ccfDataset))
    diag(agonyPoset) <- 0
  }
  
  return(agonyPoset)
  
}

# Estimate the poset based on a best agony ranking from a set of models inferred from multiple samples.
#
# @title .estimateAgonyPosetMultipleSamples
# @param models Phylogenetic models estimated for each sample included in dataset.
# @param events Driver genes to be considered in the inference.
#
.estimateAgonyPosetMultipleSamples <- function(models, events) {

  # compute a set of time orderings among events, given a set of models inferred from multiple samples
  agonyArcs <- .buildAgonyInputMultipleSamples(models, events)
  
  # estimate the best agony poset given the time orderings
  if (!is.null(agonyArcs)) {
    agonyPoset <- .computeAgonyPoset(agonyArcs, length(events))
  }
  else {
    agonyPoset <-
      matrix(1, nrow = length(models), ncol = length(events))
    diag(agonyPoset) <- 0
  }
  
  return(agonyPoset)
  
}

# Estimate a total ordering for a dataset given the observed cancer cell fractions in each sample.
#
# @examples
# data(dataExamples)
# dataExample <- dataExamples[["infer.order.in.ccf.dataset"]]
# res <- .inferOrderInCcfDataset(dataExample)
#
# @title .inferOrderInCcfDataset
# @param ccfDataset Matrix where rows are samples and columns are mutations.
#
.inferOrderInCcfDataset <- function(ccfDataset) {

  totalOrder <- list()
  
  for (i in seq_len(nrow(ccfDataset))) {
    currRes <- .inferTotalOrder(ccfDataset[i,])
    if (length(currRes) > 1) {
      totalOrder[[(length(totalOrder) + 1)]] <- currRes
    }
  }
  
  return(totalOrder)
  
}

# Estimate a total ordering for a sample given the observed cancer cell fractions.
#
# @examples
# data(dataExamples)
# dataExample <- dataExamples[["infer.total.order"]]
# res <- .inferTotalOrder(dataExample)
#
# @title .inferTotalOrder
# @param ccfEntry Cancer cell fractions for each sample.
#
.inferTotalOrder <- function(ccfEntry) {

  totalOrderSort <-
    sort(ccfEntry, decreasing = TRUE, index.return = TRUE)
  totalOrder <- list()
  totalOrder[["value"]] <-
    totalOrderSort$x[which(totalOrderSort$x > 0)]
  totalOrder[["position"]] <-
    totalOrderSort$ix[which(totalOrderSort$x > 0)]
  
  return(totalOrder)
  
}

# Create the input for computing the best agony ranking given a total ordering.
#
# @examples
# data(dataExamples)
# dataExample <- dataExamples[["build.agony.input"]]
# res <- .buildAgonyInput(dataExample)
#
# @title .buildAgonyInput
# @param agonyOrders Agony-based rankings.
#
.buildAgonyInput <- function(agonyOrders) {

  myArcs <- NULL
  
  if (length(agonyOrders) > 0) {
    for (i in seq_len(length(agonyOrders))) {
      currOrder <- agonyOrders[[i]]
      if (length(unique(currOrder$value)) > 1) {
        currOrderVal <- unique(currOrder$value)
        for (j in seq_len(length(currOrderVal[-length(currOrderVal)]))) {
          currParents <- which(currOrder$value == currOrderVal[j])
          currChildren <-
            which(currOrder$value == currOrderVal[(j + 1)])
          for (a in names(currParents)) {
            for (b in names(currChildren)) {
              currEdge <- c(currOrder$position[j], currOrder$position[(j + 1)])
              myArcs <- rbind(myArcs, currEdge)
            }
          }
        }
      }
    }
    if (!is.null(myArcs)) {
      colnames(myArcs) <- c("parent", "child")
      rownames(myArcs) <- as.character(1:nrow(myArcs))
    }
  }
  
  return(myArcs)
  
}

# Compute a set of time orderings among events, given a set of models inferred from multiple samples.
#
# @examples
# data(dataExamples)
# dataExample <- dataExamples[["build.agony.input.multiple.samples"]]
# res <- .buildAgonyInputMultipleSamples(dataExample$models,dataExample$events)
#
# @title .buildAgonyInputMultipleSamples
# @param models Phylogenetic models estimated for each sample included in dataset.
# @param events Driver genes to be considered in the inference.
#
.buildAgonyInputMultipleSamples <- function(models, events) {

  myArcs <- NULL
  
  if (length(models) > 0) {
    for (i in seq_len(length(models))) {
      currModel <- models[[i]]
      timeOrderings <- which(currModel == 1, arr.ind = TRUE)
      timeOrderingsParents <-
        rownames(currModel)[timeOrderings[, "row"]]
      timeOrderingsChildren <-
        colnames(currModel)[timeOrderings[, "col"]]
      if (length(timeOrderingsParents) > 0) {
        for (j in seq_len(length(timeOrderingsParents))) {
          # consider each arch in the input model
          currParent <- which(events == timeOrderingsParents[j])
          currChild <- which(events == timeOrderingsChildren[j])
          # add this ordering to the arcs to be considered
          if (length(currParent) > 0 && length(currChild) > 0) {
            currEdge <- c(currParent, currChild)
            myArcs <- rbind(myArcs, currEdge)
          }
        }
      }
    }
  }
  if (!is.null(myArcs)) {
    colnames(myArcs) <- c("parent", "child")
    rownames(myArcs) <- as.character(1:nrow(myArcs))
  }
  
  return(myArcs)
  
}

# Estimate the best agony ranking poset given a set of temporal observations.
#
# @title .computeAgonyPoset
# @param agonyArcs Agony-based rankings.
# @param numEvents Number of driver genes considered for the inference.
#
.computeAgonyPoset <- function(agonyArcs, numEvents) {

  # data structure to save results
  poset <- matrix(0, nrow = numEvents, ncol = numEvents)
  
  if (!is.null(agonyArcs)) {
    # estimate the best agony ranking and save the results to file
    agonyRanking <- agony(agonyArcs, .get_seed())
    
    # compute the poset based on the best agony ranking
    poset <- .buildRankingAdjMatrix(agonyRanking, numEvents)
    
  }
  
  return(poset)
  
}

# Compute the poset based on the best agony ranking.
#
# @examples
# data(dataExamples)
# dataExample <- dataExamples[["build.ranking.adj.matrix"]]
# res <- .buildRankingAdjMatrix(dataExample$agony_ranking,dataExample$num_events)
#
# @title .buildRankingAdjMatrix
# @param agonyRanking Agony-based rankings.
# @param numEvents Number of driver genes considered for the inference.
#
.buildRankingAdjMatrix <- function(agonyRanking, numEvents) {

  # data structure to save results
  adjMatrix <- matrix(0, nrow = numEvents, ncol = numEvents)
  
  if (length(agonyRanking) > 0) {
    currRank <- (max(agonyRanking[, 2]) - 1)
    while (currRank >= 0) {
      # consider each ranking from max to min + 1 (nodes at minimum ranking have no predecessors)
      currChildren <-
        agonyRanking[which(agonyRanking[, 2] == (currRank + 1)), 1]
      for (myRank in 0:currRank) {
        # the predecessors of nodes at ranking curr_rank + 1 are all the nodes at ranking from 0 to curr_rank included
        currParents <-
          agonyRanking[which(agonyRanking[, 2] == myRank), 1]
        adjMatrix[currParents, currChildren] <- 1
      }
      currRank <- currRank - 1
    }
  } else {
    adjMatrix <- matrix(1, nrow = numEvents, ncol = numEvents)
    diag(adjMatrix) <- 0
  }
  
  return(adjMatrix)
  
}

# Enforce probability raising constraints.
#
# @examples
# data(dataExamples)
# dataExample <- dataExamples[["apply.pr"]]
# res <- .applyPr( poset = dataExample$poset,
#                  prModel = dataExample$pr_model,
#                  prNull = dataExample$pr_null)
#
# @title .applyPr
# @param poset Partially order set to be considered during the inference.
# @param prModel ASCETIC statistical model based on Suppes' theory of probabilistic causation.
# @param prNull ASCETIC null model based on Suppes' theory of probabilistic causation.
#
.applyPr <- function(poset, prModel, prNull) {

  for (i in seq_len(nrow(poset))) {
    for (j in seq_len(ncol(poset))) {
      # Consider arc i --> j if it is in the poset.
      if (poset[i, j] == 1) {
        if (prModel[i, j] <= prNull[i, j]) {
          poset[i, j] <- 0
        }
      }
    }
  }
  
  return(poset)
  
}

# Perform structure learning by maximum likelihood given a poset.
#
# @examples
# data(dataExamples)
# dataExample <- dataExamples[["perform.likelihood.fit"]]
# res <- .performLikelihoodFit( dataset = dataExample$dataset,
#                               poset = dataExample$poset,
#                               regularization = "aic",
#                               command = "hc",
#                               restarts = 0)
#
# @title .performLikelihoodFit
# @param dataset Binary matrix where rows are samples and columns are mutations.
# @param poset Partially order set to be considered during the inference.
# @param regularization Regularization to be used for the maximum likelihood estimation.
# @param command Optimization technique to be used for maximum likelihood estimation.
# @param restarts Number of restarts to be performed during the maximum likelihood estimation when
# Hill Climbing optimization technique is used.
#
.performLikelihoodFit <- function(dataset,
           poset,
           regularization,
           command,
           restarts) {

    # initialization
    adjMatrix <- poset
    rownames(adjMatrix) <- as.character(1:nrow(adjMatrix))
    colnames(adjMatrix) <- as.character(1:ncol(adjMatrix))
    adjMatrixFit <-
      matrix(0, nrow = nrow(adjMatrix), ncol = ncol(adjMatrix))
    rownames(adjMatrixFit) <- as.character(1:ncol(adjMatrixFit))
    colnames(adjMatrixFit) <- as.character(1:nrow(adjMatrixFit))
    
    # create a categorical data frame from the dataset
    data <- .asCategoricalDataset(dataset)
    
    # create the blacklist based on the poset
    cont <- 0
    parent <- -1
    child <- -1
    for (i in rownames(adjMatrix)) {
      for (j in colnames(adjMatrix)) {
        if (i != j) {
          if (adjMatrix[i, j] == 0) {
            # [i,j] refers to causation i --> j
            cont <- cont + 1
            if (cont == 1) {
              parent <- i
              child <- j
            } else {
              parent <- c(parent, i)
              child <- c(child, j)
            }
          }
        }
      }
    }
    
    # perform the reconstruction by likelihood fit with regularization
    if (cont > 0) {
      blacklist <- data.frame(from = parent, to = child)
      if (command == "hc") {
        myNet <-
          hc(data,
             score = regularization,
             blacklist = blacklist,
             restart = restarts)
      } else if (command == "tabu") {
        myNet <- tabu(data, score = regularization, blacklist = blacklist)
      }
    } else {
      if (command == "hc") {
        myNet <- hc(data, score = regularization, restart = restarts)
      } else if (command == "tabu") {
        myNet <- tabu(data, score = regularization)
      }
    }
    myArcs <- myNet$arcs
    
    # build the adjacency matrix of the reconstructed topology
    if (length(nrow(myArcs)) > 0 && nrow(myArcs) > 0) {
      for (i in seq_len(nrow(myArcs))) {
        # [i,j] refers to causation i --> j
        adjMatrixFit[myArcs[i, 1], myArcs[i, 2]] <- 1
      }
    }
    
    return(adjMatrixFit)
    
}

# Create a categorical data structure.
#
# @examples
# data(dataExamples)
# dataExample <- dataExamples[["as.categorical.dataset"]]
# res <- as.categorical.dataset(dataExample)
#
# @title as.categorical.dataset
# @param dataset Binary matrix where rows are samples and columns are mutations.
#
.asCategoricalDataset <- function(dataset) {

  # create a categorical data frame from the dataset
  data <- matrix("missing", nrow = nrow(dataset), ncol = ncol(dataset))
  data[dataset == 1] <- "observed"
  data <- data.frame(data, stringsAsFactors = TRUE)

  data <- as.data.frame(lapply(data, function(x) {
    levels(x) <- c("missing", "observed")
    x
  }))
  
  # renaming
  colnames(data) <- as.character(seq_len(ncol(data)))
  rownames(data) <- as.character(seq_len(nrow(data)))
  
  return(data)
  
}

# Get a number to be used as seed
#
#
# @title get.seed
#
.get_seed <- function() {
  sample.int(.Machine$integer.max, 1)
}
