#' Perform 'omics wide association study using qgcomp
#' @description 
#' Implements an omics wide association study using QGComp to model associations of 
#' exposure mixtures with each individual 'omics feature as an outcome
#' 'omics data as either the dependent variable. Allows for either 
#' continuous or dichotomous outcomes, and provides the option to adjust for 
#' covariates. 
#' 
#' @import data.table qgcomp
#' @importFrom stats binomial coef glm lm p.adjust  
#' @export 
#' @param df Dataset
#' @param expnms Name of the exposures. Can be either continuous or 
#' dichotomous. For dichotomous variables, must set \code{q} to
#'  "NULL", and values must be either 0/1.
#' @param omics Names of all omics features in the dataset 
#' @param covars Names of covariates (can be NULL)
#' @param q NULL or number of quantiles used to create quantile indicator variables representing the exposure variables. Defaults to 4If NULL, then gcomp proceeds with un-transformed version of exposures in the input datasets (useful if data are already transformed, or for performing standard g-computation). 
#' @param confidence_level Confidence level for marginal significance 
#' (defaults to 0.95, or an alpha of 0.05)
#' 
#' @returns 
#' A data frame with 6 columns:  
#' feature_name: name of the omics feature   
#' psi: the model estimate for the feature. For linear models, this is the 
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
#' omics_df <- apply(omics_df, MARGIN = 2, FUN = function(x){rnorm(n_ids)})
#' omics_df <- as.data.frame(omics_df)
#' colnames(omics_df) <- paste0("feature_", colnames(omics_df))
#' # Simulate covariates and outcomes
#' cov_out <- data.frame(id = c(1:n_ids), 
#'                       sex = sample(c("male", "female"), 
#'                                   n_ids, replace=TRUE,prob=c(.5,.5)),
#'                      age = rnorm(10, 10, 2),
#'                      exposure1 = rlnorm(n_ids, meanlog = 2.3, sdlog = 1),
#'                      exposure2 = rlnorm(n_ids, meanlog = 2.3, sdlog = 1),
#'                      exposure3 = rlnorm(n_ids, meanlog = 2.3, sdlog = 1),
#'                      disease = sample(0:1, n_ids, replace=TRUE,prob=c(.9,.1)),
#'                      weight =  rlnorm(n_ids, meanlog = 3, sdlog = 0.2))
#' 
#' # Create Test Data
#' test_data <- cbind(cov_out, omics_df)
#' 
#' # Get names of omics
#' colnames_omic_fts <- colnames(test_data)[grep("feature_",
#'                                               colnames(test_data))]
#' 
#' # Names of exposures in mixture
#'  exposure_names = c("exposure1", "exposure2", "exposure3")
#' 
#' # Run function without covariates
#' out <- qgcomp_owas(df = test_data,
#'                    expnms = exposure_names,
#'                    omics = colnames_omic_fts,
#'                    q = 4, 
#'                    confidence_level = 0.95) 
#' 
#' 
#' # Run analysis with covariates
#' out <- qgcomp_owas(df = test_data,
#'                    expnms = c("exposure1", "exposure2", "exposure3"),
#'                    covars = c("weight", "age", "sex"),
#'                    omics = colnames_omic_fts,
#'                    q = 4, 
#'                    confidence_level = 0.95) 
#'  
#' 
qgcomp_owas <- compiler::cmpfun(
  function(df, 
           expnms,
           omics, 
           covars = NULL,
           q = 4, 
           confidence_level = 0.95){
    
    alpha = 1-confidence_level
    
    # Check for issues in data 
    # Check if variable of interest is in data
    if(FALSE %in% (expnms %in% colnames(df))){ 
      stop(paste0("Not all exposure variables are found in the data. Check data.") ) 
    }    
    # Check if all omics features are in the data
    if(FALSE %in% (omics %in% colnames(df))){ 
      stop("Not all omics variables are found in the data. Check omics column names.")  
    }    
    # Check if covars are in data
    if(FALSE %in% (covars %in% colnames(df))){ 
      stop("Not all covariates are found in the data. Check covariate column names.") 
    }    
    
    
    # Change data frame to data table for speed
    df <- data.table(df)
    
    # Pivot longer  
    dt_l = melt.data.table(
      df,
      id.vars = c(expnms, covars),
      measure.vars = omics,
      variable.name = "feature_name",
      value.name = "feature_value"
    )
    
    # Run models -------------------------
    # Split to list
    dt_list <- split(dt_l, f = dt_l$feature_name)
    # Apply qgcomp function
    
      res <- lapply(dt_list, 
                    FUN = function(x){ 
                      dt_l <- x[,c(colnames(x) %in% c(expnms, "feature_value", covars)),
                                with = FALSE]
                      fit <- qgcomp::qgcomp(feature_value~.,
                                            expnms = expnms, 
                                            q = q,
                                            alpha = 0.05,
                                            data = dt_l)
                      c(psi = fit$coef[2],
                                  # se = fit$
                                  lcl_psi = fit$ci[1],
                                  ucl_psi = fit$ci[2],
                                  p_value = fit$pval[2])
                    })
    
    
    # Add column for estimate 
    res_2 <- t(as.data.frame(res)) 
    res_2 <- as.data.frame(res_2)
    res_2$feature <- rownames(res_2)
    rownames(res_2) <- NULL
    
    # Reorder (select feature then all other cols)
    final_results <- res_2[colnames(res_2)[c(5, 1:4)]]
    
    # Add column for estimate 
    colnames(final_results) <- c("feature", "psi", #"se",
                                 "lcl_psi", "ucl_psi", "p_value") 
    
    # Calculate adjusted p value
    final_results$adjusted_pval = p.adjust(final_results$p_value, method = "fdr")
    
    final_results$threshold = ifelse(final_results$adjusted_pval < alpha, 
                                     "Significant",
                                     "Non-significant")
    
    return(final_results)
    
  }
  
)
