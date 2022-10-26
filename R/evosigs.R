#' Select evolutionary steps from ASCETIC analysis which are significantly associated with survival data.
#'
#' @examples
#' set.seed(12345)
#' data(aml_example)
#' res_example_evosigs = evosigs( 
#'                               survival_data = aml_example$survival_data, 
#'                               evolutionary_steps = aml_example$evolutionary_steps )
#'
#' @title evosigs
#' @param survival_data Survival data. Notice that no missing values are allowed.
#' @param evolutionary_steps Evolutionary steps per patient as inferred by ASCETIC.
#' @return A list of 3 elements: 1) survival_analysis, results from survival analysis.
#'                               2) evolutionary_steps, evolutionary steps significantly associated with survival data.
#'                               3) clusters_prevalence, prevalence per cluster of the selected evolutionary steps.
#' @export evosigs
#' @importFrom glmnet cv.glmnet
#' @importFrom survival Surv survfit
#' @importFrom survminer ggsurvplot
#'
evosigs <- function( survival_data, evolutionary_steps ) {

    # create analyis data structure
    analysis_data = data.frame(cbind(as.matrix(survival_data),as.matrix(evolutionary_steps)),check.rows=FALSE,stringsAsFactors=FALSE)
    colnames(analysis_data) = c("status","times",colnames(evolutionary_steps))

    # select significant evolutionary steps based on association to survival data
    survival_analysis = list()
    is.valid = FALSE
    tryCatch(expr = {
        y = Surv(as.numeric(analysis_data$times),as.numeric(analysis_data$status))
        lasso_cov = analysis_data[,colnames(analysis_data)[3:ncol(analysis_data)]]
        x = as.matrix(cbind(rep(1,nrow(analysis_data)),analysis_data[,-c(1,2)]))
        rownames(x) = 1:nrow(x)
        colnames(x) = c("(Intercept)",paste0("lasso_cov$",colnames(lasso_cov)))
        cv.fit = cv.glmnet(x,y,family="cox",maxit=1000000,alpha=1)
        survival_analysis[["data"]] = analysis_data
        survival_analysis[["fit"]] = cv.fit
        is.valid = TRUE
    }, error = function(e) {})

    # evaluate selected features evolutionary steps
    if(is.valid==TRUE) {
        curr_data = survival_analysis[["data"]]
        curr_fit = survival_analysis[["fit"]]
        beta_temp = matrix(NA,nrow=dim(curr_data)[1],ncol=1)
        Coefficients = as.numeric(coef(curr_fit,s=curr_fit$lambda.min)[-1,])
        coeff_names = gsub("lasso_cov\\$","",names(coef(curr_fit,s=curr_fit$lambda.min)[-1,]))
        coeff_names = coeff_names[which(Coefficients!=0)]
        Coefficients = Coefficients[which(Coefficients!=0)]
        for(k in 1:nrow(curr_data)) {
            x = Coefficients
            y = curr_data[k,coeff_names,drop=FALSE]
            beta_temp[k] = as.numeric(x%*%t(y))
        }
        beta = as.numeric(beta_temp)
        survival_analysis$beta = beta
        significant_coefficients = list()
        coeff_pos = as.numeric(Coefficients[which(Coefficients>0)])
        names(coeff_pos) = coeff_names[which(Coefficients>0)]
        coeff_neg = as.numeric(Coefficients[which(Coefficients<0)])
        names(coeff_neg) = coeff_names[which(Coefficients<0)]
        significant_coefficients[["positive"]] = coeff_pos
        significant_coefficients[["negative"]] = coeff_neg
        survival_analysis$significant_coefficients = significant_coefficients
        range = sort(unique(fivenum(survival_analysis$beta)[2:4]))
        survival_analysis$clusters = rep("C1",length(beta))
        if(length(range)==1) {
            survival_analysis$clusters[survival_analysis$beta>=range] = "C2"
        }
        if(length(range)==2) {
            survival_analysis$clusters[survival_analysis$beta>=range[1]] = "C2"
            survival_analysis$clusters[survival_analysis$beta>=range[2]] = "C3"
        }
        if(length(range)==3) {
            survival_analysis$clusters[survival_analysis$beta>=range[1]] = "C2"
            survival_analysis$clusters[survival_analysis$beta>=range[2]] = "C3"
            survival_analysis$clusters[survival_analysis$beta>=range[3]] = "C4"
        }
        clusters = rep(NA,length(beta))
        clust_pos = 0
        for(curr_clust_pos in sort(unique(survival_analysis$clusters))) {
            clust_pos = clust_pos + 1
            clusters[which(survival_analysis$clusters==curr_clust_pos)] = paste0("C",clust_pos)
        }
        survival_analysis$clusters = NULL
        survival_analysis$clusters = clusters
    }
    if(length(unique(survival_analysis$clusters))<2) {
        is.valid = FALSE
    }

    # perform Kaplan-Meier analysis
    if(is.valid==TRUE) {
        survival_data_km = survival_analysis$data
        K = survival_analysis$clusters
        print(ggsurvplot(survfit(Surv(as.numeric(survival_data_km$times),as.numeric(survival_data_km$status))~K,data=survival_data_km),
            xlab = "Time",
            ylab = "Survival", mark.time = TRUE, pval = TRUE, ggtheme = theme_bw(),
            title = "Kaplan-Meier analysis",
            palette = "Dark2",
            font.main = 18,
            font.x = 18,
            font.y = 18,
            font.caption = 18,
            font.legend = 18,
            font.tickslab = 18,
            risk.table = TRUE))
    }

    # save the extracted evolutionary signatures
    evolutionary_signatures = list()
    if(is.valid==FALSE) {
        message("No evolutionary step could be significantly associated to differences in prognosis...","\n")
    }
    if(is.valid==TRUE) {
        survival_analysis = survival_analysis
        selected_evolutionary_steps = c(survival_analysis$significant_coefficients$positive,survival_analysis$significant_coefficients$negative)
        selected_evolutionary_steps = selected_evolutionary_steps[sort(unique(names(selected_evolutionary_steps)))]
        clusters = survival_analysis$clusters
        freq_features = array(NA,c(length(unique(clusters)),length(evolutionary_steps)))
        rownames(freq_features) = sort(unique(clusters))
        colnames(freq_features) = sort(names(evolutionary_steps))
        for(i in rownames(freq_features)) {
            curr_samples_data = evolutionary_steps[which(clusters==i),,drop=FALSE]
            curr_estimate = colSums(curr_samples_data)/nrow(curr_samples_data)
            freq_features[i,names(curr_estimate)] = as.numeric(curr_estimate)
        }
        evolutionary_signatures = list(survival_analysis=survival_analysis,evolutionary_steps=selected_evolutionary_steps,clusters_prevalence=freq_features)
    }

    return(evolutionary_signatures)

}
