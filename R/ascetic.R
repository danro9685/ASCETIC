#' Perform the ASCETIC inference framework on single samples (using CCF) datasets 
#' with re-sampling for a robust estimation of the agony ranking.
#'
#' @examples
#' set.seed(12345)
#' data(dataset_example_single_samples)
#' data(ccf_dataset_example_single_samples)
#' data(vaf_dataset_example_single_samples)
#' res_example_single_samples_resampling = perform.ascetic.ccf.dataset.resampling( dataset = dataset_example_single_samples, 
#'                                                                                 ccf.dataset = ccf_dataset_example_single_samples, 
#'                                                                                 vaf.dataset = vaf_dataset_example_single_samples, 
#'                                                                                 nsampling = 5, 
#'                                                                                 regularization = "aic", 
#'                                                                                 command = "hc", 
#'                                                                                 restarts = 0 )
#'
#' @title perform.ascetic.ccf.dataset.resampling
#' @param dataset Binary matrix where rows are samples and columns are mutations.
#' Each cell of the matrix is 1 if the related mutation was observed in the sample; 0 otherwise.
#' Values reported in dataset must be consistent with the ones reported in ccf.dataset and vaf.dataset.
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
            agony_files = paste0(getwd(),"/agony_files")
            unlink(agony_files,recursive=TRUE,force=TRUE)
            dir.create(agony_files,showWarnings=FALSE)
            # save the orderings to file
            write.table(agony_arcs,file=paste0(agony_files,"/inputs.txt"),quote=FALSE,row.names=FALSE,col.names=FALSE)
            # estimate a best agony ranking and save the results to file
            agony(agony_files)
            # read the estimated best agony ranking
            agony_ranking = read.table(file=paste0(agony_files,"/outputs.txt"),check.names=FALSE,stringsAsFactors=FALSE)
            # remove the created files and directories
            unlink(agony_files,recursive=TRUE,force=TRUE)
            cont = cont + 1
            for(j in 1:nrow(agony_ranking)) {
                resampling_agony_ranking_estimate[agony_ranking[j,"V1"],"rank"] = resampling_agony_ranking_estimate[agony_ranking[j,"V1"],"rank"] + agony_ranking[j,"V2"]
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
#' res_example_phylogenies_dataset_bootstrap = perform.ascetic.phylogenies.dataset.bootstrap( dataset = dataset_example_phylogenies, 
#'                                                                                            models = models_phylogenies, 
#'                                                                                            nsampling = 5, 
#'                                                                                            regularization = "aic", 
#'                                                                                            command = "hc", 
#'                                                                                            restarts = 0 )
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
            agony_files = paste0(getwd(),"/agony_files")
            unlink(agony_files,recursive=TRUE,force=TRUE)
            dir.create(agony_files,showWarnings=FALSE)
            # save the orderings to file
            write.table(agony_arcs,file=paste0(agony_files,"/inputs.txt"),quote=FALSE,row.names=FALSE,col.names=FALSE)
            # estimate a best agony ranking and save the results to file
            agony(agony_files)
            # read the estimated best agony ranking
            agony_ranking = read.table(file=paste0(agony_files,"/outputs.txt"),check.names=FALSE,stringsAsFactors=FALSE)
            # remove the created files and directories
            unlink(agony_files,recursive=TRUE,force=TRUE)
            cont = cont + 1
            for(j in 1:nrow(agony_ranking)) {
                bootstrap_agony_ranking_estimate[agony_ranking[j,"V1"],"rank"] = bootstrap_agony_ranking_estimate[agony_ranking[j,"V1"],"rank"] + agony_ranking[j,"V2"]
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
#' res_example_single_samples = perform.ascetic.ccf.dataset( dataset = dataset_example_single_samples, 
#'                                                           ccf.dataset = ccf_dataset_example_single_samples, 
#'                                                           regularization = "aic", 
#'                                                           command = "hc", 
#'                                                           restarts = 0 )
#'
#' @title perform.ascetic.ccf.dataset
#' @param dataset Binary matrix where rows are samples and columns are mutations.
#' Each cell of the matrix is 1 if the related mutation was observed in the sample; 0 otherwise.
#' Values reported in dataset must be consistent with the ones reported in ccf.dataset.
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
#' res_example_phylogenies_dataset = perform.ascetic.phylogenies.dataset( dataset = dataset_example_phylogenies, 
#'                                                                        models = models_phylogenies, 
#'                                                                        regularization = "aic", 
#'                                                                        command = "hc", 
#'                                                                        restarts = 0 )
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

# estimate the model for the probability raising, i.e., j --> i, P(t_j<t_i|i,j), from a dataset of cancer cell fractions
estimate.pr.model.single.samples <- function( ccf.dataset ) {
    
    pr.model.estimate = array(0,c(ncol(ccf.dataset),ncol(ccf.dataset)))
    
    # consider any pair of events in the dataset
    for (i in 1:nrow(pr.model.estimate)) {
        for (j in i:nrow(pr.model.estimate)) {
            if(i!=j) {
                # get the subset of the dataset where events i and j are co-occurring
                co_occurring_samples = which(ccf.dataset[,i]>0&ccf.dataset[,j]>0)
                # if there is at least 1 sample where the 2 genes are co-occurring
                if(length(co_occurring_samples)>0) {
                    curr_i = as.vector(ccf.dataset[co_occurring_samples,i])
                    curr_j = as.vector(ccf.dataset[co_occurring_samples,j])
                    pr.model.estimate[i,j] = length(which(curr_i>curr_j))/length(co_occurring_samples)
                    pr.model.estimate[j,i] = length(which(curr_j>curr_i))/length(co_occurring_samples)
                }
                else {
                    pr.model.estimate[i,j] = 0
                    pr.model.estimate[j,i] = 0
                }
            }
        }
    }
    
    return(pr.model.estimate)
    
}

# estimate the model for the probability raising, i.e., j --> i, P(t_j<t_i|i,j), from a inferences on multiple patients
estimate.pr.model.multiple.samples <- function( models, events ) {
    
    pr.model.estimate = array(0,c(length(events),length(events)))
    
    # get the counts for each event
    counts.time.orderings = array(0,dim(pr.model.estimate))
    counts.co_occurance = array(0,dim(pr.model.estimate))
    
    for(i in 1:length(models)) {
        
        curr_model = models[[i]]
        # compute the time ordering counts, i.e., P(t_j<t_i,i,j)
        time_orderings = which(curr_model==1,arr.ind=TRUE)

        if(length(time_orderings)>0) {
            time_orderings_parents = rownames(curr_model)[time_orderings[,"row"]]
            time_orderings_children = colnames(curr_model)[time_orderings[,"col"]]
            for(j in 1:length(time_orderings_parents)) {
                counts.time.orderings[which(events==time_orderings_parents[j]),which(events==time_orderings_children[j])] = counts.time.orderings[which(events==time_orderings_parents[j]),which(events==time_orderings_children[j])] + 1
            }
        }
        
        # compute the co_occurance count increment, i.e., P(i,j)
        counts.co_occurance[which(events%in%colnames(curr_model)),which(events%in%colnames(curr_model))] = counts.co_occurance[which(events%in%colnames(curr_model)),which(events%in%colnames(curr_model))] + 1
        
    }
    
    pr.model.estimate = counts.time.orderings / counts.co_occurance
    pr.model.estimate[which(is.nan(pr.model.estimate),arr.ind=TRUE)] = 0
    
    return(pr.model.estimate)
    
}

# estimate the null model for the probability raising, i.e., j --> i, P(i,not j)
estimate.pr.null <- function( dataset ) {
    
    pr.null.estimate = array(0,c(ncol(dataset),ncol(dataset)))

    # compute marginal and joint probabilities
    probs.estimate = estimate.probs(dataset)
    
    # consider any pair of events in the dataset
    for (i in 1:nrow(pr.null.estimate)) {
        for (j in i:nrow(pr.null.estimate)) {
            if(i!=j) {
                pr.null.estimate[i,j] = probs.estimate$marginal.probs[i,1] - probs.estimate$joint.probs[i,j]
                pr.null.estimate[j,i] = probs.estimate$marginal.probs[j,1] - probs.estimate$joint.probs[j,i]
            }
        }
    }
    
    return(pr.null.estimate)
    
}

# estimate marginal and joint probabilities from the dataset
estimate.probs <- function( dataset ) {

    pair.count = array(0,dim=c(ncol(dataset),ncol(dataset)))

    for (i in 1:ncol(dataset)) {
        for (j in 1:ncol(dataset)) {
            val1 = dataset[,i]
            val2 = dataset[,j]
            pair.count[i,j] = (t(val1)%*%val2)
        }
    }

    marginal.probs = array(as.matrix(diag(pair.count)/nrow(dataset)),dim=c(ncol(dataset),1))
    joint.probs = as.matrix(pair.count/nrow(dataset))

    return(list(marginal.probs=marginal.probs,joint.probs=joint.probs))

}

# estimate the poset based on a best agony ranking from a dataset of cancer cell fractions
estimate.agony.poset.single.samples <- function( ccf.dataset ) {

    # compute a set of total orderings, one per patient
    total_orderings = infer.order.in.ccf.dataset(ccf.dataset)

    # compute a set of time orderings among each event, given the total orderings
    agony_arcs = build.agony.input(total_orderings)

    # estimate a best agony poset given the time orderings
    if(!is.null(agony_arcs)) {
        agony_poset = compute.agony.poset(agony_arcs,ncol(ccf.dataset))
    }
    else {
        agony_poset = array(1,c(ncol(ccf.dataset),ncol(ccf.dataset)))
        diag(agony_poset) = 0
    }

    return(agony_poset)

}

# estimate the poset based on a best agony ranking from a set of models inferred from multiple samples
estimate.agony.poset.multiple.samples <- function( models, events ) {

    # compute a set of time orderings among each event, given a set of models inferred from multiple samples
    agony_arcs = build.agony.input.multiple.samples(models,events)

    # estimate a best agony poset given the time orderings
    if(!is.null(agony_arcs)) {
        agony_poset = compute.agony.poset(agony_arcs,length(events))
    }
    else {
        agony_poset = array(1,c(ncol(ccf.dataset),length(events)))
        diag(agony_poset) = 0
    }

    return(agony_poset)

}

# estimate a total ordering for a dataset given the observed cancer cell fractions in each sample
infer.order.in.ccf.dataset <- function( ccf.dataset  ) {

    total.order = list()

    for(i in 1:nrow(ccf.dataset)) {
        curr_res = infer.total.order(ccf.dataset[i,])
        if(length(curr_res)>1) {
            total.order[[(length(total.order)+1)]] = curr_res
        }
    }

    return(total.order)

}

# estimate a total ordering for a sample given the observed cancer cell fractions
infer.total.order <- function( ccf.entry ) {

    total.order_sort = sort(ccf.entry,decreasing=TRUE,index.return=TRUE)
    total.order = list()
    total.order[["value"]] = total.order_sort$x[which(total.order_sort$x>0)]
    total.order[["position"]] = total.order_sort$ix[which(total.order_sort$x>0)]

    return(total.order)

}

# create the input file for computing a best agony ranking given a total ordering
build.agony.input <- function( agony_orders ) {

    my_arcs = NULL

    if(length(agony_orders)>0) {
        for(i in 1:length(agony_orders)) {
            curr_order = agony_orders[[i]]
            if(length(unique(curr_order$value))>1) {
                curr_order_val = unique(curr_order$value)
                for (j in 1:length(curr_order_val[-length(curr_order_val)])) {
                    curr_parents = which(curr_order$value==curr_order_val[j])
                    curr_children = which(curr_order$value==curr_order_val[(j+1)])
                    for(a in names(curr_parents)) {
                        for(b in names(curr_children)) {
                            curr_edge = c(curr_order$position[j],curr_order$position[(j+1)])
                            my_arcs = rbind(my_arcs,curr_edge)
                        }
                    }
                }
            }
        }
        if(!is.null(my_arcs)) {
            colnames(my_arcs) = c("parent","child")
            rownames(my_arcs) = as.character(1:nrow(my_arcs))
        }
    }

    return(my_arcs)

}

# compute a set of time orderings among each event, given a set of models inferred from multiple samples
build.agony.input.multiple.samples <- function( models, events ) {

    my_arcs = NULL

    if(length(models)>0) {
        
        for(i in 1:length(models)) {
            
            curr_model = models[[i]]
            time_orderings = which(curr_model==1,arr.ind=TRUE)
            time_orderings_parents = rownames(curr_model)[time_orderings[,"row"]]
            time_orderings_children = colnames(curr_model)[time_orderings[,"col"]]
            
            if(length(time_orderings_parents)>0) {
                for(j in 1:length(time_orderings_parents)) {
                    # consider each arch in the input model
                    curr_parent = which(events==time_orderings_parents[j])
                    curr_child = which(events==time_orderings_children[j])
                    # add this ordering to the arcs to be considered
                    if(length(curr_parent)>0 && length(curr_child)>0) {
                        curr_edge = c(curr_parent,curr_child)
                        my_arcs = rbind(my_arcs,curr_edge)
                    }
                    
                }
            }
            
        }
        
    }
    
    if(!is.null(my_arcs)) {
        colnames(my_arcs) = c("parent","child")
        rownames(my_arcs) = as.character(1:nrow(my_arcs))
    }

    return(my_arcs)

}

# estimate an optimal agony ranking poset given a set of temporal observations
compute.agony.poset <- function( agony_arcs, num_events, agony_files = paste0(getwd(),"/agony_files") ) {

    poset = array(0,c(num_events,num_events))

    # safely create the directoy to save the files to compute the best agony ranking
    unlink(agony_files,recursive=TRUE,force=TRUE)
    dir.create(agony_files,showWarnings=FALSE)

    if(!is.null(agony_arcs)) {
        # save the orderings to file
        write.table(agony_arcs,file=paste0(agony_files,"/inputs.txt"),quote=FALSE,row.names=FALSE,col.names=FALSE)
        # estimate a best agony ranking and save the results to file
        agony(agony_files)
        # read the estimated best agony ranking
        agony_ranking = read.table(file=paste0(agony_files,"/outputs.txt"),check.names=FALSE,stringsAsFactors=FALSE)
        # compute the poset based on the best agony ranking
        poset = build.ranking.adj.matrix(agony_ranking,num_events)
    }

    # remove the created files and directories
    unlink(agony_files,recursive=TRUE,force=TRUE)

    return(poset)

}

# compute the poset based on the best agony ranking
build.ranking.adj.matrix <- function( agony_ranking, num_events ) {
    
    adj.matrix = array(0,c(num_events,num_events))

    if(length(agony_ranking)>0) {
        curr_rank = (max(agony_ranking[,2])-1)
        while(curr_rank>=0) {
            # consider each ranking from max to min + 1 (nodes at minimum ranking have no predecessors)
            curr_children = agony_ranking[which(agony_ranking[,2]==(curr_rank+1)),1]
            for (my_rank in 0:curr_rank) {
                # the predecessors of nodes at ranking curr_rank + 1 are all the nodes at ranking from 0 to curr_rank included
                curr_parents = agony_ranking[which(agony_ranking[,2]==my_rank),1]
                adj.matrix[curr_parents,curr_children] = 1
            }
            curr_rank = curr_rank - 1
        }
    }
    else {
        adj.matrix = array(1,c(ncol(num_events),ncol(num_events)))
        diag(adj.matrix) = 0
    }
    
    return(adj.matrix)
    
}

# enforce probability raising constraints
apply.pr <- function( poset, pr_model, pr_null ) {

    for(i in 1:nrow(poset)) {
        for(j in 1:ncol(poset)) {
            # consider arc i --> j, if it is in the poset
            if(poset[i,j]==1) {
                if(pr_model[i,j]<=pr_null[i,j]) {
                    poset[i,j] = 0
                }
            }
        }
    }
    
    return(poset)

}

# perform structure learning by maximum likelihood given a poset
perform.likelihood.fit <- function( dataset, poset, regularization, command, restarts ) {

    # initialization
    adj.matrix = poset
    rownames(adj.matrix) = as.character(1:nrow(adj.matrix))
    colnames(adj.matrix) = as.character(1:ncol(adj.matrix))
    adj.matrix.fit = array(0, c(nrow(adj.matrix),ncol(adj.matrix)))
    rownames(adj.matrix.fit) = as.character(1:ncol(adj.matrix.fit))
    colnames(adj.matrix.fit) = as.character(1:nrow(adj.matrix.fit))

    # create a categorical data frame from the dataset
    data = as.categorical.dataset(dataset)

    # create the blacklist based on the poset
    cont = 0
    parent = -1
    child = -1
    for (i in rownames(adj.matrix)) {
        for (j in colnames(adj.matrix)) {
            if(i!=j) {
                if (adj.matrix[i,j] == 0) {
                    # [i,j] refers to causation i --> j
                    cont = cont + 1
                    if (cont == 1) {
                        parent = i
                        child = j
                    }
                    else {
                        parent = c(parent,i)
                        child = c(child,j)
                    }
                }
            }
        }
    }

    # perform the reconstruction by likelihood fit with regularization
    if(cont>0) {
        blacklist = data.frame(from=parent,to=child)
        if(command == "hc") {
            my.net = hc(data,score=regularization,blacklist=blacklist,restart=restarts)
        }
        else if (command == "tabu") {
            my.net = tabu(data,score=regularization,blacklist=blacklist)
        }
    }
    else {
        if(command == "hc") {
            my.net = hc(data,score=regularization,restart=restarts)
        }
        else if (command == "tabu") {
            my.net = tabu(data,score=regularization)
        }
    }
    my.arcs = my.net$arcs
    
    # build the adjacency matrix of the reconstructed topology
    if(length(nrow(my.arcs))>0 && nrow(my.arcs)>0) {
        for (i in 1:nrow(my.arcs)) {
            # [i,j] refers to causation i --> j
            adj.matrix.fit[my.arcs[i,1],my.arcs[i,2]] = 1
        }
    }
    return(adj.matrix.fit)

}

# create a categorical data structure
as.categorical.dataset <- function( dataset ) {

    # create a categorical data frame from the dataset
    data = array("missing",c(nrow(dataset),ncol(dataset)))
    for (i in 1:nrow(dataset)) {
        for (j in 1:ncol(dataset)) {
            if (dataset[i,j]==1) {
                data[i,j] = "observed"
            }
        }
    }
    data = data.frame(data,stringsAsFactors=TRUE)
    for (n in names(data)) {
        levels(data[[n]]) = c("missing","observed")
    }

    # renaming
    colnames(data) = as.character(1:ncol(data))
    rownames(data) = as.character(1:nrow(data))
    return(data)

}
