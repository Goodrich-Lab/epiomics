#' Perform 'omics wide association study
#' @description 
#' Implements a meet in the middle analysis for identifying omics associated 
#' with both exposures and outcomes
#' 
#' @import data.table
#' @export 
#' @param df Dataframe
#' @param exposure Name of the exposure of interest. Can be either continuous or 
#' dichotomous. 
#' @param outcome Name of the outcome of interest. Can be either continuous or 
#' dichotomous. For dichotomous variables, must set \code{outcome_model} to
#'  "logistic", and values must be either 0/1 or a factor with the first level 
#'  representing the reference group.
#' @param omics Names of all omics features in the dataset 
#' @param covars Names of covariates (can be NULL)
#' @param outcome_model_type "linear" for linear models (via lm) or "logistic" for 
#' logistic (via glm) 
#' @param confidence_level Confidence level for marginal significance 
#' (defaults to 0.95)
#' 
#' @returns 
#' A list of two, for the exposure-omics and omics-outcome MWAS. Each element
#' contains a data frame with 6 columns:  
#' feature_name: name of the omics feature   
#' estimate: the model estimate for the feature. For linear models, this is the 
#' beta; for logistic models, this is the log odds. 
#' se: Standard error of the estimate
#' p_value: p-value for the estimate
#' adjusted_pval: FDR adjusted p-value
#' threshold: Marginal significance, based on unadjusted p-values 
#' 
meet_in_middle <- compiler::cmpfun(
  function(df, 
           exposure,
           outcome,
           omics,
           covars = NULL,
           outcome_model_type = "linear",
           confidence_level = 0.95){
    alpha = 1-confidence_level
    # exposure_omics_owas
    exposure_omics_owas <- epiomics::owas(df = df,
                                          var = exposure, 
                                          omics = omics, 
                                          covars = covars,
                                          var_exposure_or_outcome = "exposure", 
                                          model_type = "linear",
                                          confidence_level = confidence_level)
    
    
    # omics_outcome_owas
    omics_outcome_owas <- epiomics::owas(df = df,
                                         var = "outcome", 
                                         omics = omics,
                                         covars = covars,
                                         var_exposure_or_outcome = "outcome", 
                                         model_type = outcome_model_type,
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
    }
    
    final_results <- list(exposure_omics_owas = exposure_omics_owas, 
                          omics_outcome_owas = omics_outcome_owas, 
                          overlap = overlap)
    
    
    
    return(final_results)
    
  }
)