#' Perform 'omics wide association study
#' @description 
#' Implements an omics wide association study with the option of using the 
#' 'omics data as either the dependent variable (i.e., for performing an 
#' exposure --> 'omics analysis) or using the 'omics as the independent variable 
#' (i.e., for performing an 'omics --> outcome analysis). Allows for either 
#' continuous or dichotomous outcomes, and provides the option to adjust for 
#' covariates. 
#' 
#' @import data.table
#' @export 
#' @param df Dataset
#' @param var Name of the variable of interest- this is usually either an 
#' exposure variable or an outcome variable. Can be either continuous or 
#' dichotomous. For dichotomous variables, must set \code{model_type} to
#'  "logistic", and values must be either 0/1 or a factor with the first level 
#'  representing the reference group.
#' @param omics Names of all omics features in the dataset 
#' @param covars Names of covariates (can be NULL)
#' @param var_exposure_or_outcome Is the variable of interest an exposure 
#' (independent variable) or outcome (dependent variable)? Must be either
#' "exposure" or "outcome"
#' @param model_type "linear" for linear models (via lm) or "logistic" for 
#' logistic (via glm) 
#' @param confidence_level Confidence level for marginal significance 
#' (defaults to 0.95, or an alpha of 0.05)
#' 
#' @returns 
#' A data frame with 6 columns:  
#' feature_name: name of the omics feature   
#' estimate: the model estimate for the feature. For linear models, this is the 
#' beta; for logistic models, this is the log odds. 
#' se: Standard error of the estimate
#' p_value: p-value for the estimate
#' adjusted_pval: FDR adjusted p-value
#' threshold: Marginal significance, based on unadjusted p-values 
#' 
#' @examples 
#' # Simulate dataset
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
#'                                    n_ids, replace=T,prob=c(.5,.5)),
#'                       age = rnorm(10, 10, 2),
#'                       pfos = rlnorm(n_ids, meanlog = 2.3, sdlog = 1),
#'                       disease = sample(0:1, n_ids, replace=T,prob=c(.9,.1)),
#'                       weight =  rlnorm(n_ids, meanlog = 3, sdlog = 0.2))
#' 
#' # Create Test Data
#' test_data <- cbind(cov_out, omics_df)
#' 
#' # Get names of omics
#' colnames_omic_fts <- colnames(test_data)[grep("feature_",
#'                                               colnames(test_data))]
#' 
#' # Run function with continuous exposure as the variable of interest
#' owas(df = test_data, 
#'      var = "pfos", 
#'      omics = colnames_omic_fts, 
#'      covars = c("age", "sex"), 
#'      var_exposure_or_outcome = "exposure", 
#'      model_type = "linear")
#' 
#' 
#' # Run function with dichotomous outcome as the variable of interest
#' owas(df = test_data, 
#'      var = "disease", 
#'      omics = colnames_omic_fts, 
#'      covars = c("age", "sex"), 
#'      var_exposure_or_outcome = "outcome", 
#'      model_type = "logistic")
#'  
#' 
owas <- compiler::cmpfun(
  function(df, 
           var,
           omics, 
           covars,
           var_exposure_or_outcome, 
           model_type = "linear", 
           confidence_level = 0.95){
    alpha = 1-confidence_level
    # Change data frame to data table for speed
    df <- data.table(df)
    
    # Pivot longer  
    dt_l = melt.data.table(
      df,
      id.vars = c(var, covars),
      measure.vars = omics,
      variable.name = "feature_name",
      value.name = "feature_value"
    )
    
    # Set formula for model ------------------
    # depending on whether variable of interest is the exposure or the outcome 
    if(var_exposure_or_outcome == "exposure"){
      # If variable is exposure: 
      if(is.null(covars)){
        mod_formula <- paste0("feature_value~", var)
      } else {
        mod_formula <- paste0("feature_value~", 
                              paste0(covars, collapse = "+"), "+",
                              var)
      } 
      
    } else if(var_exposure_or_outcome == "outcome"){
      # If variable is outcome: 
      if(is.null(covars)){
        mod_formula <- paste0(var, "~ feature_value")
      } else{
        mod_formula <- paste0(var, "~", 
                              paste0(covars, collapse = "+"),
                              "+feature_value")
      }
      
    } else {
      stop("var_exposure_or_outcome must be either \"exposure\" or \"outcome\" ")
    }
    
    
    # Run models -------------------------
    if(model_type == "linear"){
      # Linear models:
      res <- dt_l[, 
                  {fit <- lm(mod_formula, data = .SD) 
                  coef(summary(fit))[nrow(coef(summary(fit))), # Select last row
                                     c(1, 2, 4)] # Select Estimate, Std Error, and p_val
                  }, 
                  by = feature_name]
      
      # Add column for estimate 
      res <- cbind(res, c("estimate", "se", "p_value"))
      
      # Pivot wider
      final_results <- dcast(data = res, 
                             feature_name ~ V2, 
                             value.var = "V1")[,c(1, 2, 4, 3)]
      
    } else if(model_type == "logistic"){
      # Logistic models
      res <- dt_l[, 
                  {fit <- glm(mod_formula, data = .SD, family=binomial(link='logit')) 
                  coef(summary(fit))[nrow(coef(summary(fit))), # Select last row
                                     c(1, 2, 4)] # Select Estimate, Std Error, and p_val
                  }, 
                  by = feature_name]
      
      # Add column for estimate 
      res <- cbind(res, c("estimate", "se", "p_value"))
      
      # Pivot wider
      final_results <- dcast(data = res, 
                             feature_name ~ V2, 
                             value.var = "V1")[,c(1, 2, 4, 3)]
      
    } else {
      stop("model_type must be either \"linear\" or \"logistic\" ")
    }
    
    # Calculate adjusted p value
    final_results$adjusted_pval = p.adjust(final_results$p_value, method = "fdr")
    
    final_results$threshold = ifelse(final_results$adjusted_pval < alpha, 
                                     "Significant",
                                     "Non-significant")
    
    return(final_results)
    
  }
)
