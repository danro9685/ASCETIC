#' Select the evolutionary steps from ASCETIC analysis that are significantly associated with survival data.
#'
#' @examples
#' set.seed(12345)
#' data(amlExample)
#' resExampleEvosigs = evoSigs( survivalData = amlExample$survival_data,
#'                              evolutionarySteps = amlExample$evolutionary_steps )
#'
#' @title evoSigs
#' @param survivalData Survival data. Notice that no missing values are allowed.
#' @param evolutionarySteps Evolutionary steps per patient as inferred by ASCETIC.
#' @return A list of 3 elements: 1) survivalAnalysis, results from survival analysis.
#'                               2) evolutionarySteps, evolutionary steps significantly associated with survival data.
#'                               3) clustersPrevalence, prevalence per cluster of the selected evolutionary steps.
#' @export evoSigs
#'
evoSigs <- function(survivalData, evolutionarySteps) {

  # create analysis data structure
  analysisData <- data.frame(
    cbind(as.matrix(survivalData),
          as.matrix(evolutionarySteps)),
    check.rows = FALSE,
    stringsAsFactors = FALSE
  )
  colnames(analysisData) <-
    c("status", "times", colnames(evolutionarySteps))
  
  # select the significant evolutionary steps based on association to survival data
  survivalAnalysis <- list()
  isValid <- FALSE
  tryCatch({
    lassoCov <-
      analysisData[, colnames(analysisData)[3:ncol(analysisData)]]
    selFeatures <- paste0("lassoCov$", colnames(lassoCov))
    x <-
      as.matrix(cbind(rep(1, nrow(analysisData)), analysisData[,-c(1, 2)]))
    rownames(x) <- 1:nrow(x)
    colnames(x) <- c("(Intercept)", selFeatures)
    y <- Surv(as.numeric(analysisData$times),
              as.numeric(analysisData$status))
    cvFit <-
      cv.glmnet(x,
                y,
                family = "cox",
                maxit = 1000000,
                alpha = 1)
    survivalAnalysis[["data"]] <- analysisData
    survivalAnalysis[["fit"]] <- cvFit
    isValid <- TRUE
  }, error = function(e) {
  })
  
  # evaluate the selected features evolutionary steps
  if (isValid == TRUE) {
    currData <- survivalAnalysis[["data"]]
    currFit <- survivalAnalysis[["fit"]]
    betaTemp <- matrix(NA, nrow = nrow(currData), ncol = 1)
    Coefficients <-
      as.numeric(coef(currFit, s = currFit$lambda.min)[selFeatures,])
    coeffNames <-
      gsub("lassoCov\\$", "", names(coef(currFit, s = currFit$lambda.min)[selFeatures,]))
    coeffNames <- coeffNames[which(Coefficients != 0)]
    Coefficients <- Coefficients[which(Coefficients != 0)]
    x <- Coefficients
    for (k in 1:nrow(currData)) {
      y <- currData[k, coeffNames, drop = FALSE]
      betaTemp[k] <- as.numeric(x %*% t(y))
    }
    beta <- as.numeric(betaTemp)
    survivalAnalysis$beta <- beta
    significantCoefficients <- list()
    coeffPos <- as.numeric(Coefficients[which(Coefficients > 0)])
    names(coeffPos) <- coeffNames[which(Coefficients > 0)]
    coeffNeg <- as.numeric(Coefficients[which(Coefficients < 0)])
    names(coeffNeg) <- coeffNames[which(Coefficients < 0)]
    significantCoefficients[["positive"]] <- coeffPos
    significantCoefficients[["negative"]] <- coeffNeg
    survivalAnalysis$significantCoefficients <-
      significantCoefficients
    range <- sort(unique(fivenum(survivalAnalysis$beta)[2:4]))
    survivalAnalysis$clusters <- rep("C1", length(beta))
    if (length(range) == 1) {
      survivalAnalysis$clusters[survivalAnalysis$beta >= range] <- "C2"
    }
    if (length(range) == 2) {
      survivalAnalysis$clusters[survivalAnalysis$beta >= range[1]] <- "C2"
      survivalAnalysis$clusters[survivalAnalysis$beta >= range[2]] <-
        "C3"
    }
    if (length(range) == 3) {
      survivalAnalysis$clusters[survivalAnalysis$beta >= range[1]] <- "C2"
      survivalAnalysis$clusters[survivalAnalysis$beta >= range[2]] <-
        "C3"
      survivalAnalysis$clusters[survivalAnalysis$beta >= range[3]] <-
        "C4"
    }
    clusters <- rep(NA, length(beta))
    clustPos <- 0
    for (currClustPos in sort(unique(survivalAnalysis$clusters))) {
      clustPos <- clustPos + 1
      clusters[which(survivalAnalysis$clusters == currClustPos)] <-
        paste0("C", clustPos)
    }
    survivalAnalysis$clusters <- NULL
    survivalAnalysis$clusters <- clusters
  }
  if (length(unique(survivalAnalysis$clusters)) < 2) {
    isValid <- FALSE
  }
  
  # save the extracted evolutionary signatures
  evolutionarySignatures <- list()
  if (isValid == FALSE) {
    message(
      "No evolutionary step could be significantly associated to differences in prognosis...",
      "\n"
    )
  }
  if (isValid == TRUE) {
    survivalAnalysis <- survivalAnalysis
    selectedEvolutionarySteps <-
      c(
        survivalAnalysis$significantCoefficients$positive,
        survivalAnalysis$significantCoefficients$negative
      )
    selectedEvolutionarySteps <-
      selectedEvolutionarySteps[sort(unique(names(selectedEvolutionarySteps)))]
    clusters <- survivalAnalysis$clusters
    freqFeatures <-
      array(NA, c(length(unique(clusters)), length(evolutionarySteps)))
    rownames(freqFeatures) <- sort(unique(clusters))
    colnames(freqFeatures) <- sort(names(evolutionarySteps))
    for (i in rownames(freqFeatures)) {
      currSamplesData <-
        evolutionarySteps[which(clusters == i), , drop = FALSE]
      currEstimate <-
        colSums(currSamplesData) / nrow(currSamplesData)
      freqFeatures[i, names(currEstimate)] <-
        as.numeric(currEstimate)
    }
    evolutionarySignatures <-
      list(
        survivalAnalysis = survivalAnalysis,
        evolutionarySteps = selectedEvolutionarySteps,
        clustersPrevalence = freqFeatures
      )
  }
  
  return(evolutionarySignatures)
  
}
