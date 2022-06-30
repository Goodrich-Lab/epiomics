#' Perform 'omics wide association study
#' @description 
#' Implements a meet in the middle analysis for identifying omics associated 
#' with both exposures and outcomes, as described by 
#' 
#' @import data.table
#' @export 
#' @param df Dataframe
#' @param exposure Name of the exposure of interest. Can be either continuous or 
#' dichotomous. 
#' @param outcome Name of the outcome of interest. Can be either continuous or 
#' dichotomous. For dichotomous variables, must set \code{outcome_family} to
#'  "logistic", and values must be either 0/1 or a factor with the first level 
#'  representing the reference group.
#' @param omics Names of all omics features in the dataset 
#' @param covars Names of covariates (can be NULL)
#' @param outcome_family "gaussian" for linear models (via lm) or "binomial" for 
#' logistic (via glm) 
#' @param confidence_level Confidence level for marginal significance 
#' (defaults to 0.95)
#' 
#' @returns 
#' A list of three dataframes, containing: 
#' 1) Results from the Exposure-Omics Wide Association Study  
#' 2) Results from the Omics-Outcome Wide Association Study
#' 3) Overlapping significant features from 1 and 2.     
#' For each omics wide association, results are provided in a data frame with 6 columns:  
#' feature_name: name of the omics feature   
#' estimate: the model estimate for the feature. For linear models, this is the 
#' beta: for logistic models, this is the log odds. 
#' se: Standard error of the estimate
#' p_value: p-value for the estimate
#' adjusted_pval: FDR adjusted p-value
#' threshold: Marginal significance, based on unadjusted p-values 
#' 
#' @examples 
#'  # Simulate dataset
#' set.seed(4656)
#' n_omic_ftrs = 100
#' n_ids = 400
#' # Simulate omics
#' omics_df <- matrix(nrow = n_ids,
#'                    ncol = n_omic_ftrs)
#' omics_df <- apply(omics_df, MARGIN = 2, FUN = function(x){rnorm(n_omic_ftrs)})
#' omics_df <- as.data.frame(omics_df)
#' colnames(omics_df) <- paste0("feature_", colnames(omics_df))
#' # Simulate covariates and outcomes
#' cov_out <- data.frame(id = c(1:n_ids),
#'                       sex = sample(c("male", "female"),
#'                                    n_ids, replace=TRUE,prob=c(.5,.5)),
#'                       age = rnorm(10, 10, 2),
#'                       exposure = rlnorm(n_ids, meanlog = 2.3, sdlog = 1),
#'                       disease = sample(0:1, n_ids, replace=TRUE,prob=c(.9,.1)),
#'                       weight =  rlnorm(n_ids, meanlog = 3, sdlog = 0.2))
#' 
#' # Create Test Data
#' test_data <- cbind(cov_out, omics_df)
#' 
#' # Get names of omics
#' colnames_omic_fts <- colnames(test_data)[grep("feature_",
#'                                               colnames(test_data))]
#' 
#' # Meet in the middle with a dichotomous outcome
#' res <- meet_in_middle(df = test_data,
#'                       exposure = "exposure", 
#'                       outcome = "disease", 
#'                       omics = colnames_omic_fts,
#'                       covars = c("age", "sex"), 
#'                       outcome_family = "binomial")
#' 
#' # Meet in the middle with a continuous outcome 
#' res <- meet_in_middle(df = test_data,
#'                       exposure = "exposure", 
#'                       outcome = "weight", 
#'                       omics = colnames_omic_fts,
#'                       covars = c("age", "sex"), 
#'                       outcome_family = "gaussian")
#' 
#' # Meet in the middle with a continuous outcome and no covariates
#' res <- meet_in_middle(df = test_data,
#'                       exposure = "exposure", 
#'                       outcome = "weight", 
#'                       omics = colnames_omic_fts,
#'                       outcome_family = "gaussian")
#' 
meet_in_middle <- compiler::cmpfun(
  function(df, 
           exposure,
           outcome,
           omics,
           covars = NULL,
           outcome_family = "gaussian",
           confidence_level = 0.95){
    alpha = 1-confidence_level
    
    df <- data.table::as.data.table(df)
    # exposure_omics_owas
    exposure_omics_owas <- epiomics::owas(df = df,
                                          var = exposure, 
                                          omics = omics, 
                                          covars = covars,
                                          var_exposure_or_outcome = "exposure", 
                                          family = "gaussian",
                                          confidence_level = confidence_level)
    
    
    # omics_outcome_owas
    omics_outcome_owas <- epiomics::owas(df = df,
                                         var = outcome, 
                                         omics = omics,
                                         covars = covars,
                                         var_exposure_or_outcome = "outcome", 
                                         family = outcome_family,
                                         confidence_level = confidence_level)
    
    # Find overlap
    x_o_fts <- exposure_omics_owas[exposure_omics_owas$p_value<alpha]$feature_name
    o_y_fts <- omics_outcome_owas[omics_outcome_owas$p_value<alpha]$feature_name
    overlap_fts <- x_o_fts[x_o_fts %in% o_y_fts]
    
    if(length(overlap_fts)>0){
      # Subset only overlapping sig
      x_o <- exposure_omics_owas[exposure_omics_owas$feature_name %in% overlap_fts]
      names(x_o) <- paste0(colnames(x_o), c("", rep("_exp_omic", ncol(x_o)-1 )))
      o_y <- omics_outcome_owas[omics_outcome_owas$feature_name %in% overlap_fts]
      names(o_y) <- paste0(colnames(o_y), c("", rep("_omic_out", ncol(x_o)-1 )))
      
      # Merge
      overlap <- merge(x_o, o_y, by = "feature_name")
    } else {overlap = NULL}
    
    final_results <- list(exposure_omics_owas = exposure_omics_owas, 
                          omics_outcome_owas = omics_outcome_owas, 
                          overlap = overlap)
    
    return(final_results)
    
  }
)