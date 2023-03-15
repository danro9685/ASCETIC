#' Estimate the model for the probability raising, i.e., j --> i, P(t_j<t_i|i,j), from a dataset of cancer cell fractions.
#'
#' @examples
#' data(data_examples)
#' data_example = data_examples[["estimate.pr.model.single.samples"]]
#' res = estimate.pr.model.single.samples(ccf.dataset=data_example)
#'
#' @title estimate.pr.model.single.samples
#' @param ccf.dataset Matrix where rows are samples and columns are mutations.
#'
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

#' Estimate the model for the probability raising, i.e., j --> i, P(t_j<t_i|i,j), from a inferences on multiple patients.
#'
#' @examples
#' data(data_examples)
#' data_example = data_examples[["estimate.pr.model.multiple.samples"]]
#' res = estimate.pr.model.multiple.samples(data_example$models,data_example$events)
#'
#' @title estimate.pr.model.multiple.samples
#' @param models Phylogenetic models estimated for each sample included in dataset.
#' @param events Driver genes to be considered in the inference.
#'
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

#' Estimate the null model for the probability raising, i.e., j --> i, P(i,not j).
#'
#' @examples
#' data(data_examples)
#' data_example = data_examples[["estimate.pr.null"]]
#' res = estimate.pr.null(data_example)
#'
#' @title estimate.pr.null
#' @param dataset Binary matrix where rows are samples and columns are mutations.
#'
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

#' Estimate marginal and joint probabilities from the dataset.
#'
#' @examples
#' data(data_examples)
#' data_example = data_examples[["estimate.probs"]]
#' res = estimate.probs(data_example)
#'
#' @title estimate.probs
#' @param dataset Binary matrix where rows are samples and columns are mutations.
#'
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

#' Estimate the poset based on a best agony ranking from a dataset of cancer cell fractions.
#'
#' @title estimate.agony.poset.single.samples
#' @param ccf.dataset Matrix where rows are samples and columns are mutations.
#'
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

#' Estimate the poset based on a best agony ranking from a set of models inferred from multiple samples.
#'
#' @title estimate.agony.poset.multiple.samples
#' @param models Phylogenetic models estimated for each sample included in dataset.
#' @param events Driver genes to be considered in the inference.
#'
estimate.agony.poset.multiple.samples <- function( models, events ) {

    # compute a set of time orderings among each event, given a set of models inferred from multiple samples
    agony_arcs = build.agony.input.multiple.samples(models,events)

    # estimate a best agony poset given the time orderings
    if(!is.null(agony_arcs)) {
        agony_poset = compute.agony.poset(agony_arcs,length(events))
    }
    else {
        agony_poset = array(1,c(length(models),length(events)))
        diag(agony_poset) = 0
    }

    return(agony_poset)

}

#' Estimate a total ordering for a dataset given the observed cancer cell fractions in each sample.
#'
#' @examples
#' data(data_examples)
#' data_example = data_examples[["infer.order.in.ccf.dataset"]]
#' res = infer.order.in.ccf.dataset(data_example)
#'
#' @title infer.order.in.ccf.dataset
#' @param ccf.dataset Matrix where rows are samples and columns are mutations.
#'
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

#' Estimate a total ordering for a sample given the observed cancer cell fractions.
#'
#' @examples
#' data(data_examples)
#' data_example = data_examples[["infer.total.order"]]
#' res = infer.total.order(data_example)
#'
#' @title infer.total.order
#' @param ccf.entry Cancer cell fractions for each sample.
#'
infer.total.order <- function( ccf.entry ) {

    total.order_sort = sort(ccf.entry,decreasing=TRUE,index.return=TRUE)
    total.order = list()
    total.order[["value"]] = total.order_sort$x[which(total.order_sort$x>0)]
    total.order[["position"]] = total.order_sort$ix[which(total.order_sort$x>0)]

    return(total.order)

}

#' Create the input file for computing a best agony ranking given a total ordering.
#'
#' @examples
#' data(data_examples)
#' data_example = data_examples[["build.agony.input"]]
#' res = build.agony.input(data_example)
#'
#' @title build.agony.input
#' @param agony_orders Agony-based rankings.
#'
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

#' Compute a set of time orderings among each event, given a set of models inferred from multiple samples.
#'
#' @examples
#' data(data_examples)
#' data_example = data_examples[["build.agony.input.multiple.samples"]]
#' res = build.agony.input.multiple.samples(data_example$models,data_example$events)
#'
#' @title build.agony.input.multiple.samples
#' @param models Phylogenetic models estimated for each sample included in dataset.
#' @param events Driver genes to be considered in the inference.
#'
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

#' Estimate an optimal agony ranking poset given a set of temporal observations.
#'
#' @title compute.agony.poset
#' @param agony_arcs Agony-based rankings.
#' @param num_events Number of driver genes considered for the inference.
#' @param agony_files Location of agony files.
#'
compute.agony.poset <- function( agony_arcs, num_events, agony_files = paste0(getwd(),"/agony_files") ) {

    poset = array(0,c(num_events,num_events))

    # safely create the directoy to save the files to compute the best agony ranking
    unlink(agony_files,recursive=TRUE,force=TRUE)
    dir.create(agony_files,showWarnings=FALSE)

    if(!is.null(agony_arcs)) {
        # save the orderings to file
        write.table(agony_arcs,file=paste0(agony_files,"/inputs.txt"),quote=FALSE,row.names=FALSE,col.names=FALSE)
        

        # estimate a best agony ranking and save the results to file
        agony(agony_arcs, paste0(agony_files,"/outputs.txt"))
        # read the estimated best agony ranking
        agony_ranking = read.table(file=paste0(agony_files,"/outputs.txt"),check.names=FALSE,stringsAsFactors=FALSE)
        # compute the poset based on the best agony ranking
        poset = build.ranking.adj.matrix(agony_ranking,num_events)
    }

    # remove the created files and directories
    unlink(agony_files,recursive=TRUE,force=TRUE)

    return(poset)

}

#' Compute the poset based on the best agony ranking.
#'
#' @examples
#' data(data_examples)
#' data_example = data_examples[["build.ranking.adj.matrix"]]
#' res = build.ranking.adj.matrix(agony_ranking=data_example$agony_ranking,num_events=data_example$num_events)
#'
#' @title build.ranking.adj.matrix
#' @param agony_ranking Agony-based rankings.
#' @param num_events Number of driver genes considered for the inference.
#'
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

#' Enforce probability raising constraints.
#'
#' @examples
#' data(data_examples)
#' data_example = data_examples[["apply.pr"]]
#' res = apply.pr(poset=data_example$poset,pr_model=data_example$pr_model,pr_null=data_example$pr_null)
#'
#' @title apply.pr
#' @param poset Partially order set to be considered during the inference.
#' @param pr_model ASCETIC statistical model based on Suppes' theory of probabilistic causation.
#' @param pr_null ASCETIC null model based on Suppes' theory of probabilistic causation.
#'
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

#' Perform structure learning by maximum likelihood given a poset.
#'
#' @examples
#' data(data_examples)
#' data_example = data_examples[["perform.likelihood.fit"]]
#' res = perform.likelihood.fit( dataset = data_example$dataset, 
#'                               poset = data_example$poset, 
#'                               regularization = "aic", 
#'                               command = "hc", 
#'                               restarts = 0)
#'
#' @title perform.likelihood.fit
#' @param dataset Binary matrix where rows are samples and columns are mutations.
#' @param poset Partially order set to be considered during the inference.
#' @param regularization Regularization to be used for the maximum likelihood estimation.
#' @param command Optimization technique to be used for maximum likelihood estimation.
#' @param restarts Number of restarts to be performed during the maximum likelihood estimation when Hill Climbing optimization technique is used.
#'
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

#' Create a categorical data structure.
#'
#' @examples
#' data(data_examples)
#' data_example = data_examples[["as.categorical.dataset"]]
#' res = as.categorical.dataset(data_example)
#'
#' @title as.categorical.dataset
#' @param dataset Binary matrix where rows are samples and columns are mutations.
#'
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
