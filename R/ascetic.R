#' Perform the ASCETIC inference framework on single samples (using CCF) datasets 
#' with re-sampling for a robust estimation of the agony ranking.
#'
#' @examples
#' set.seed(12345)
#' data(dataset_example_single_samples)
#' data(ccf_dataset_example_single_samples)
#' data(vaf_dataset_example_single_samples)
#' res_example_single_samples_resampling = perform.ascetic.ccf.dataset.resampling( 
#'                                                 dataset = dataset_example_single_samples, 
#'                                                 ccf.dataset = ccf_dataset_example_single_samples, 
#'                                                 vaf.dataset = vaf_dataset_example_single_samples, 
#'                                                 nsampling = 5, 
#'                                                 regularization = "aic", 
#'                                                 command = "hc", 
#'                                                 restarts = 0 )
#'
#' @title perform.ascetic.ccf.dataset.resampling
#' @param dataset Binary matrix where rows are samples and columns are mutations.
#' Each cell of the matrix is 1 if the related mutation was observed in the sample; 0 otherwise.
#' Values reported in the variable named dataset must be consistent with the ones reported in ccf.dataset and vaf.dataset.
#' @param ccf.dataset Matrix where rows are samples and columns are mutations.
#' Each cell of the matrix is the cancer cell fraction (CCF) estimated for the related mutation when observed in the sample.
#' The CCF value is 0 if the mutation was not observed in the sample.
#' Values reported in ccf.dataset must be consistent with the ones reported in dataset and vaf.dataset.
#' @param vaf.dataset R data.frame with 8 columns: 1) SAMPLE_ID, sample name.
#'                                                 2) GENE_ID, gene name.
#'                                                 3) REF_COUNT, total counts for reference allele.
#'                                                 4) ALT_COUNT, total counts for alternate allele.
#'                                                 5) COPY_NUMBER, total copy number estimate.
#'                                                 6) NORMAL_PLOIDY, ploidy for normal sample; this is either 1 for mutations on sex chromosomes or 2.
#'                                                 7) VAF_ESTIMATE, variant allele frequency (VAF) estimate.
#'                                                 8) CCF_ESTIMATE, cancer cell fraction (CCF) estimate.
#' Values reported in vaf.dataset must be consistent with the ones reported in dataset and ccf.dataset.
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
#'                               2) ccf_dataset, input ccf_dataset.
#'                               3) ranking_estimate, ranking among mutations estimated by agony.
#'                                  Lower rankings correspond to early mutations.
#'                               4) poset, partially order set among mutations estimated by ASCETIC from the agony ranking.
#'                               5) inference, inferred ASCETIC evolutionary model for each selected regularization.
#' @export perform.ascetic.ccf.dataset.resampling
#' @import bnlearn
#' @import Rcpp
#' @importFrom stats rbeta
#' @importFrom utils read.table write.table
#' @useDynLib ASCETIC, .registration = TRUE
#' @exportPattern "^[[:alpha:]]+"
#'
perform.ascetic.ccf.dataset.resampling <- function( dataset, ccf.dataset, vaf.dataset, nsampling = 100, regularization = c("aic","bic"), command = "hc", restarts = 10 ) {

    # estimate the probability raising models
    pr_model = estimate.pr.model.single.samples(ccf.dataset)
    pr_null = estimate.pr.null(dataset)

    # estimate agony ranking nsampling times by resampling to obtain a robust and stable estimation
    cat(0,"\n")
    resampling_agony_ranking_estimate = array(0,c(ncol(ccf.dataset),2))
    rownames(resampling_agony_ranking_estimate) = colnames(ccf.dataset)
    colnames(resampling_agony_ranking_estimate) = c("variable","rank")
    resampling_agony_ranking_estimate[,"variable"] = 1:ncol(ccf.dataset)
    cont = 0
    for(i in 1:nsampling) {
        # assess temporal priority with minimum agony for nsampling sampling iterations
        curr_ccf.dataset = array(0,c(length(unique(vaf.dataset$SAMPLE_ID)),length(unique(vaf.dataset$GENE_ID))))
        rownames(curr_ccf.dataset) = sort(unique(vaf.dataset$SAMPLE_ID))
        colnames(curr_ccf.dataset) = sort(unique(vaf.dataset$GENE_ID))
        for(j in 1:nrow(vaf.dataset)) {
            vaf_dataset_samples = rbeta(n=1,shape1=vaf.dataset$ALT_COUNT[j],shape2=vaf.dataset$REF_COUNT[j],ncp=0)
            curr_ccf_estimate = as.numeric(((vaf.dataset$NORMAL_PLOIDY[j]+(vaf.dataset$COPY_NUMBER[j]-vaf.dataset$NORMAL_PLOIDY[j]))*vaf_dataset_samples))
            curr_ccf.dataset[vaf.dataset$SAMPLE_ID[j],vaf.dataset$GENE_ID[j]] = curr_ccf_estimate
        }
        curr_ccf.dataset[which(curr_ccf.dataset>1)] = 1
        # compute a set of total orderings, one per patient
        total_orderings = infer.order.in.ccf.dataset(curr_ccf.dataset)
        # compute a set of time orderings among each event, given the total orderings
        agony_arcs = build.agony.input(total_orderings)
        # estimate a best agony poset given the time orderings
        if(!is.null(agony_arcs)) {
            # safely create the directoy to save the files to compute the best agony ranking
            # agony_files = paste0(getwd(),"/agony_files")
            # unlink(agony_files,recursive=TRUE,force=TRUE)
            # dir.create(agony_files,showWarnings=FALSE)
            # save the orderings to file
            # write.table(agony_arcs,file=paste0(agony_files,"/inputs.txt"),quote=FALSE,row.names=FALSE,col.names=FALSE)
            
            
            # estimate a best agony ranking and save the results to file
            
            #print("agony_arcs")
            #print(agony_arcs)
            
            
            agony_ranking = agony(agony_arcs)
            # read the estimated best agony ranking
            # agony_ranking = read.table(file=paste0(agony_files,"/outputs.txt"),check.names=FALSE,stringsAsFactors=FALSE)
            
            # remove the created files and directories
            # unlink(agony_files,recursive=TRUE,force=TRUE)
            cont = cont + 1
            for(j in 1:nrow(agony_ranking)) {
                resampling_agony_ranking_estimate[agony_ranking[j,1],"rank"] = resampling_agony_ranking_estimate[agony_ranking[j,1],"rank"] + agony_ranking[j,2]
            }
        }
        cat(i/nsampling,"\n")
    }

    # the final estimation is the mean (approximated to integer) of the nsampling rankings
    if(cont>0) {
        resampling_agony_ranking_estimate[,"rank"] = round(resampling_agony_ranking_estimate[,"rank"]/cont)
    }

    # FIRST: assess temporal priority with minimum agony
    agony_poset = build.ranking.adj.matrix(resampling_agony_ranking_estimate,ncol(ccf.dataset))
    # SECOND: verify probability raising
    agony_poset = apply.pr(agony_poset,pr_model,pr_null)

    # perform the inference
    agony_inference = list()
    for(reg in regularization) {
        agony_inference[[reg]] = perform.likelihood.fit(dataset,agony_poset,reg,command,restarts)
    }
    
    # create the data structures with the posets and the inference results
    rownames(agony_poset) = as.character(1:ncol(dataset))
    colnames(agony_poset) = as.character(1:ncol(dataset))
    poset = agony_poset
    for(reg in regularization) {
        rownames(agony_inference[[reg]]) = as.character(1:ncol(dataset))
        colnames(agony_inference[[reg]]) = as.character(1:ncol(dataset))
    }
    inference = agony_inference
        
    # create the data structures with the results
    rownames(dataset) = as.character(1:nrow(dataset))
    colnames(dataset) = as.character(1:ncol(dataset))
    rownames(ccf.dataset) = as.character(1:nrow(ccf.dataset))
    colnames(ccf.dataset) = as.character(1:ncol(ccf.dataset))
    rownames(pr_model) = as.character(1:ncol(dataset))
    colnames(pr_model) = as.character(1:ncol(dataset))
    rownames(pr_null) = as.character(1:ncol(dataset))
    colnames(pr_null) = as.character(1:ncol(dataset))
    results = list(dataset=dataset,ccf_dataset=ccf.dataset,ranking_estimate=resampling_agony_ranking_estimate,poset=poset,inference=inference)

    return(results)

}

#' Perform the ASCETIC inference framework on multiple samples (using phylogenies as inputs) datasets 
#' with bootstrap for a robust estimation of the agony ranking.
#'
#' @examples
#' set.seed(12345)
#' data(dataset_example_phylogenies)
#' data(models_phylogenies)
#' res_example_phylogenies_dataset_bootstrap = perform.ascetic.phylogenies.dataset.bootstrap( 
#'                                                     dataset = dataset_example_phylogenies, 
#'                                                     models = models_phylogenies, 
#'                                                     nsampling = 5, 
#'                                                     regularization = "aic", 
#'                                                     command = "hc", 
#'                                                     restarts = 0 )
#'
#' @title perform.ascetic.phylogenies.dataset.bootstrap
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
#'                               3) ranking_estimate, ranking among mutations estimated by agony.
#'                                  Lower rankings correspond to early mutations.
#'                               4) poset, partially order set among mutations estimated by ASCETIC from the agony ranking.
#'                               5) inference, inferred ASCETIC evolutionary model for each selected regularization.
#' @export perform.ascetic.phylogenies.dataset.bootstrap
#' @import bnlearn
#' @import Rcpp
#' @importFrom utils read.table write.table
#' @useDynLib ASCETIC, .registration = TRUE
#' @exportPattern "^[[:alpha:]]+"
#'
perform.ascetic.phylogenies.dataset.bootstrap <- function( dataset, models, nsampling = 100, regularization = c("aic","bic"), command = "hc", restarts = 10 ) {

    # estimate the probability raising models
    pr_model = estimate.pr.model.multiple.samples(models,colnames(dataset))
    pr_null = estimate.pr.null(dataset)

    # estimate agony ranking nsampling times by bootstrap to obtain a robust and stable estimation
    cat(0,"\n")
    bootstrap_agony_ranking_estimate = array(0,c(ncol(dataset),2))
    rownames(bootstrap_agony_ranking_estimate) = colnames(dataset)
    colnames(bootstrap_agony_ranking_estimate) = c("variable","rank")
    bootstrap_agony_ranking_estimate[,"variable"] = 1:ncol(dataset)
    cont = 0
    for(i in 1:nsampling) {
        # assess temporal priority with minimum agony for nsampling bootstrap iterations
        curr_models = sample(1:length(models),replace=TRUE)
        curr_models = models[curr_models]
        for(j in 1:length(curr_models)) {
            currm = curr_models[[j]]
            curr_models[[j]] = array(0,c(length(colnames(dataset)),length(colnames(dataset))))
            rownames(curr_models[[j]]) = colnames(dataset)
            colnames(curr_models[[j]]) = colnames(dataset)
            curr_models[[j]][rownames(currm),colnames(currm)] = currm
        }
        # compute a set of time orderings among each event, given a set of models inferred from multiple samples
        agony_arcs = build.agony.input.multiple.samples(curr_models,colnames(dataset))
        # estimate a best agony poset given the time orderings
        if(!is.null(agony_arcs)) {
            # safely create the directoy to save the files to compute the best agony ranking
            #agony_files = paste0(getwd(),"/agony_files")
            #unlink(agony_files,recursive=TRUE,force=TRUE)
            #dir.create(agony_files,showWarnings=FALSE)
            # save the orderings to file
            # write.table(agony_arcs,file=paste0(agony_files,"/inputs.txt"),quote=FALSE,row.names=FALSE,col.names=FALSE)
            

            # estimate a best agony ranking and save the results to file
            agony_ranking = agony(agony_arcs)
            # read the estimated best agony ranking
            # agony_ranking = read.table(file=paste0(agony_files,"/outputs.txt"),check.names=FALSE,stringsAsFactors=FALSE)
            # remove the created files and directories
            #unlink(agony_files,recursive=TRUE,force=TRUE)
            cont = cont + 1
            for(j in 1:nrow(agony_ranking)) {
                bootstrap_agony_ranking_estimate[agony_ranking[j,1],"rank"] = bootstrap_agony_ranking_estimate[agony_ranking[j,1],"rank"] + agony_ranking[j,2]
            }
        }
        cat(i/nsampling,"\n")
    }

    # the final estimation is the mean (approximated to integer) of the nsampling rankings
    if(cont>0) {
        bootstrap_agony_ranking_estimate[,"rank"] = round(bootstrap_agony_ranking_estimate[,"rank"]/cont)
    }

    # FIRST: assess temporal priority with minimum agony
    agony_poset = build.ranking.adj.matrix(bootstrap_agony_ranking_estimate,ncol(dataset))
    # SECOND: verify probability raising
    agony_poset = apply.pr(agony_poset,pr_model,pr_null)

    # perform the inference
    agony_inference = list()
    for(reg in regularization) {
        agony_inference[[reg]] = perform.likelihood.fit(dataset,agony_poset,reg,command,restarts)
    }
    
    # create the data structures with the posets and the inference results
    rownames(agony_poset) = as.character(1:ncol(dataset))
    colnames(agony_poset) = as.character(1:ncol(dataset))
    poset = agony_poset
    for(reg in regularization) {
        rownames(agony_inference[[reg]]) = as.character(1:ncol(dataset))
        colnames(agony_inference[[reg]]) = as.character(1:ncol(dataset))
    }
    inference = agony_inference
        
    # create the data structures with the results
    rownames(dataset) = as.character(1:nrow(dataset))
    colnames(dataset) = as.character(1:ncol(dataset))
    names(models) = as.character(1:length(models))
    rownames(pr_model) = as.character(1:ncol(dataset))
    colnames(pr_model) = as.character(1:ncol(dataset))
    rownames(pr_null) = as.character(1:ncol(dataset))
    colnames(pr_null) = as.character(1:ncol(dataset))
    results = list(dataset=dataset,models=models,ranking_estimate=bootstrap_agony_ranking_estimate,poset=poset,inference=inference)

    return(results)

}

#' Perform the ASCETIC inference framework on single samples (using CCF) datasets.
#'
#' @examples
#' set.seed(12345)
#' data(dataset_example_single_samples)
#' data(ccf_dataset_example_single_samples)
#' res_example_single_samples = perform.ascetic.ccf.dataset( 
#'                                      dataset = dataset_example_single_samples, 
#'                                      ccf.dataset = ccf_dataset_example_single_samples, 
#'                                      regularization = "aic", 
#'                                      command = "hc", 
#'                                      restarts = 0 )
#'
#' @title perform.ascetic.ccf.dataset
#' @param dataset Binary matrix where rows are samples and columns are mutations.
#' Each cell of the matrix is 1 if the related mutation was observed in the sample; 0 otherwise.
#' Values reported in the variable named dataset must be consistent with the ones reported in ccf.dataset.
#' @param ccf.dataset Matrix where rows are samples and columns are mutations.
#' Each cell of the matrix is the cancer cell fraction (CCF) estimated for the related mutation when observed in the sample.
#' The CCF value is 0 if the mutation was not observed in the sample.
#' Values reported in ccf.dataset must be consistent with the ones reported in dataset.
#' @param regularization Regularization to be used for the maximum likelihood estimation.
#' Possible values are aic for the Akaike information criterion and bic for the Bayesian information criterion.
#' For the complete list of options, we refer to the manual of the bnlearn package.
#' @param command Optimization technique to be used for maximum likelihood estimation.
#' Valid values are either hc for Hill Climbing or tabu for Tabu Search.
#' @param restarts Number of restarts to be performed during the maximum likelihood estimation when Hill Climbing optimization technique is used.
#' Higher values lead to improved estimates, but require higher computational burden; default value is 10.
#' This parameter is ignored if tabu search is selected.
#' @return A list of 4 elements: 1) dataset, input dataset.
#'                               2) ccf_dataset, input ccf_dataset.
#'                               3) poset, partially order set among mutations estimated by ASCETIC from the agony ranking.
#'                               4) inference, inferred ASCETIC evolutionary model for each selected regularization.
#' @export perform.ascetic.ccf.dataset
#' @import bnlearn
#' @import Rcpp
#' @importFrom utils read.table write.table
#' @useDynLib ASCETIC, .registration = TRUE
#' @exportPattern "^[[:alpha:]]+"
#'
perform.ascetic.ccf.dataset <- function( dataset, ccf.dataset, regularization = c("aic","bic"), command = "hc", restarts = 10 ) {

    # estimate the probability raising models
    pr_model = estimate.pr.model.single.samples(ccf.dataset)
    pr_null = estimate.pr.null(dataset)

    # compute the agony based poset
    # FIRST: assess temporal priority with minimum agony
    agony_poset = estimate.agony.poset.single.samples(ccf.dataset)
    # SECOND: verify probability raising
    agony_poset = apply.pr(agony_poset,pr_model,pr_null)

    # perform the inference
    agony_inference = list()
    for(reg in regularization) {
        agony_inference[[reg]] = perform.likelihood.fit(dataset,agony_poset,reg,command,restarts)
    }
    
    # create the data structures with the posets and the inference results
    rownames(agony_poset) = as.character(1:ncol(dataset))
    colnames(agony_poset) = as.character(1:ncol(dataset))
    poset = agony_poset
    for(reg in regularization) {
        rownames(agony_inference[[reg]]) = as.character(1:ncol(dataset))
        colnames(agony_inference[[reg]]) = as.character(1:ncol(dataset))
    }
    inference = agony_inference
        
    # create the data structures with the results
    rownames(dataset) = as.character(1:nrow(dataset))
    colnames(dataset) = as.character(1:ncol(dataset))
    rownames(ccf.dataset) = as.character(1:nrow(ccf.dataset))
    colnames(ccf.dataset) = as.character(1:ncol(ccf.dataset))
    rownames(pr_model) = as.character(1:ncol(dataset))
    colnames(pr_model) = as.character(1:ncol(dataset))
    rownames(pr_null) = as.character(1:ncol(dataset))
    colnames(pr_null) = as.character(1:ncol(dataset))
    results = list(dataset=dataset,ccf_dataset=ccf.dataset,poset=poset,inference=inference)

    return(results)

}

#' Perform the ASCETIC inference framework on multiple samples (using phylogenies as inputs) datasets.
#'
#' @examples
#' set.seed(12345)
#' data(dataset_example_phylogenies)
#' data(models_phylogenies)
#' res_example_phylogenies_dataset = perform.ascetic.phylogenies.dataset( 
#'                                                   dataset = dataset_example_phylogenies, 
#'                                                   models = models_phylogenies, 
#'                                                   regularization = "aic", 
#'                                                   command = "hc", 
#'                                                   restarts = 0 )
#'
#' @title perform.ascetic.phylogenies.dataset
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
#' @export perform.ascetic.phylogenies.dataset
#' @import bnlearn
#' @import Rcpp
#' @importFrom utils read.table write.table
#' @useDynLib ASCETIC, .registration = TRUE
#' @exportPattern "^[[:alpha:]]+"
#'
perform.ascetic.phylogenies.dataset <- function( dataset, models, regularization = c("aic","bic"), command = "hc", restarts = 10 ) {

    # estimate the probability raising models
    pr_model = estimate.pr.model.multiple.samples(models,colnames(dataset))
    pr_null = estimate.pr.null(dataset)

    # compute the agony based poset
    # FIRST: assess temporal priority with minimum agony
    agony_poset = estimate.agony.poset.multiple.samples(models,colnames(dataset))
    # SECOND: verify probability raising
    agony_poset = apply.pr(agony_poset,pr_model,pr_null)

    # perform the inference
    agony_inference = list()
    for(reg in regularization) {
        agony_inference[[reg]] = perform.likelihood.fit(dataset,agony_poset,reg,command,restarts)
    }
    
    # create the data structures with the posets and the inference results
    rownames(agony_poset) = as.character(1:ncol(dataset))
    colnames(agony_poset) = as.character(1:ncol(dataset))
    poset = agony_poset
    for(reg in regularization) {
        rownames(agony_inference[[reg]]) = as.character(1:ncol(dataset))
        colnames(agony_inference[[reg]]) = as.character(1:ncol(dataset))
    }
    inference = agony_inference

    # create the data structures with the results
    rownames(dataset) = as.character(1:nrow(dataset))
    colnames(dataset) = as.character(1:ncol(dataset))
    names(models) = as.character(1:length(models))
    rownames(pr_model) = as.character(1:ncol(dataset))
    colnames(pr_model) = as.character(1:ncol(dataset))
    rownames(pr_null) = as.character(1:ncol(dataset))
    colnames(pr_null) = as.character(1:ncol(dataset))
    results = list(dataset=dataset,models=models,poset=poset,inference=inference)

    return(results)

}
