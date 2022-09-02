
# Test that owas works with multiple outcomes ----------------
test_that("owas works with multiple traits", {
  set.seed(4656)
  n_omic_ftrs = 100
  n_ids = 400
  # Simulate omics
  omics_df <- matrix(nrow = n_ids,
                     ncol = n_omic_ftrs)
  omics_df <- apply(omics_df, 
                    MARGIN = 2, 
                    FUN = function(x){rnorm(n_ids)})
  omics_df <- as.data.frame(omics_df)
  colnames(omics_df) <- paste0("feature_", colnames(omics_df))
  
  # Simulate covariates and outcomes
  cov_out <- data.frame(id = c(1:n_ids),
                        sex = sample(c("male", "female"),
                                     n_ids, replace=TRUE,prob=c(.5,.5)),
                        age = rnorm(10, 10, 2),
                        weight =  rlnorm(n_ids, meanlog = 3, sdlog = 0.2),
                        exposure1 = rlnorm(n_ids, meanlog = 2.3, sdlog = 1),
                        exposure2 = rlnorm(n_ids, meanlog = 2.3, sdlog = 1),
                        exposure3 = rlnorm(n_ids, meanlog = 2.3, sdlog = 1),
                        disease1 = sample(0:1, n_ids, replace=TRUE, prob=c(.9,.1)),
                        disease2 = sample(0:1, n_ids, replace=TRUE,prob=c(.9,.1)))
  
  # Create Test Data
  test_data <- cbind(cov_out, omics_df)
  
  # Get names of omics
  colnames_omic_fts <- colnames(test_data)[grep("feature_",
                                                colnames(test_data))]
  
  # Run function with multiple exposures as the variable of interest
  expnms = c("exposure1", "exposure2", "exposure3")
  
  # Continuous exposure, covars
  cont_owas_out <- owas(
    df = test_data,
    var = expnms,
    omics = colnames_omic_fts,
    covars = c("age", "sex"),
    var_exposure_or_outcome = "exposure",
    family = "gaussian")
  
  
  # Test that function returns expected dimensions
  testthat::expect_equal(object = ncol(cont_owas_out), 
                         expected = 8)
  testthat::expect_equal(object = nrow(cont_owas_out), 
                         expected = length(colnames_omic_fts)*length(expnms))
  
  # Continuous exposure, no covars
  owas(df = test_data,
       var = expnms,
       omics = colnames_omic_fts,
       # covars = c("age", "sex"),
       var_exposure_or_outcome = "exposure",
       family = "gaussian")
  
  # Continuous outcome, with covars
  owas(df = test_data,
       var = expnms,
       omics = colnames_omic_fts,
       covars = c("age", "sex"),
       var_exposure_or_outcome = "outcome",
       family = "gaussian")
  
  # Continuous outcome, no covars
  owas(df = test_data,
       var = expnms,
       omics = colnames_omic_fts,
       # covars = c("age", "sex"),
       var_exposure_or_outcome = "outcome",
       family = "gaussian")
  
  
})




# Test that owas works with single outcome ----------------
test_that("owas works with single continuous outcome", {
  set.seed(4656)
  n_omic_ftrs = 100
  n_ids = 400
  
  # Simulate omics
  omics_df <- matrix(nrow = n_ids,
                     ncol = n_omic_ftrs)
  omics_df <- apply(omics_df, MARGIN = 2, FUN = function(x){rnorm(n_ids)})
  omics_df <- as.data.frame(omics_df)
  colnames(omics_df) <- paste0("feature_", colnames(omics_df))
  # Simulate covariates and outcomes
  cov_out <- data.frame(id = c(1:n_ids),
                        sex = sample(c("male", "female"),
                                     n_ids, replace=TRUE,prob=c(.5,.5)),
                        age = rnorm(10, 10, 2),
                        weight =  rlnorm(n_ids, meanlog = 3, sdlog = 0.2),
                        exposure1 = rlnorm(n_ids, meanlog = 2.3, sdlog = 1),
                        exposure2 = rlnorm(n_ids, meanlog = 2.3, sdlog = 1),
                        exposure3 = rlnorm(n_ids, meanlog = 2.3, sdlog = 1),
                        disease1 = sample(0:1, n_ids, replace=TRUE, prob=c(.9,.1)),
                        disease2 = sample(0:1, n_ids, replace=TRUE,prob=c(.9,.1)))
  
  # Create Test Data
  test_data <- cbind(cov_out, omics_df)
  
  # Get names of omics
  colnames_omic_fts <- colnames(test_data)[grep("feature_",
                                                colnames(test_data))]
  
  # Run function with continuous exposure as the variable of interest
  cont_owas_out <- owas(df = test_data,
                        var = "exposure1",
                        omics = colnames_omic_fts,
                        covars = c("age", "sex"),
                        var_exposure_or_outcome = "exposure",
                        family = "gaussian")
  
  # Test that function returns expected dimensions
  testthat::expect_equal(object = ncol(cont_owas_out), 
                         expected = 8)
  testthat::expect_equal(object = nrow(cont_owas_out), 
                         expected = 100)
  
  # Run function with dichotomous outcome as the variable of interest
  # owas(df = test_data,
  #                var = "disease1",
  #                omics = colnames_omic_fts,
  #                covars = c("age", "sex"),
  #                var_exposure_or_outcome = "outcome",
  #                family = "binomial")
})


# Test that meet in middle works ----------------
test_that("meet_in_middle works", {
  set.seed(4656)
  n_omic_ftrs = 100
  n_ids = 400
  
  # Simulate omics
  omics_df <- matrix(nrow = n_ids,
                     ncol = n_omic_ftrs)
  omics_df <- apply(omics_df, MARGIN = 2, FUN = function(x){rlnorm(n_ids)})
  omics_df <- as.data.frame(omics_df)
  colnames(omics_df) <- paste0("feature_", colnames(omics_df))
  # Simulate covariates and outcomes
  cov_out <- data.frame(id = c(1:n_ids),
                        sex = sample(c("male", "female"),
                                     n_ids, replace=TRUE,prob=c(.5,.5)),
                        age = rnorm(10, 10, 2),
                        weight =  rlnorm(n_ids, meanlog = 3, sdlog = 0.2),
                        exposure1 = rlnorm(n_ids, meanlog = 2.3, sdlog = 1),
                        exposure2 = rlnorm(n_ids, meanlog = 2.3, sdlog = 1),
                        exposure3 = rlnorm(n_ids, meanlog = 2.3, sdlog = 1),
                        disease1 = sample(0:1, n_ids, replace=TRUE, prob=c(.9,.1)),
                        disease2 = sample(0:1, n_ids, replace=TRUE,prob=c(.9,.1)))
  
  
  # Create Test Data
  test_data <- cbind(cov_out, omics_df)
  
  # Get names of omics
  colnames_omic_fts <- colnames(test_data)[grep("feature_",
                                                colnames(test_data))]
  
  # Run function with continuous exposure as the variable of interest
  cont_meet_in_middle <- meet_in_middle(df = test_data,
                                        exposure = "exposure1", 
                                        outcome = "disease1", 
                                        omics = colnames_omic_fts,
                                        # covars = c("age", "sex"), 
                                        outcome_family = "binomial", 
                                        confidence_level = 0.95)
  #############################################################
  # Test that function returns expected dimensions
  testthat::expect_equal(object = ncol(cont_meet_in_middle$overlap), 
                         expected = 15)
  testthat::expect_equal(object = nrow(cont_meet_in_middle$overlap), 
                         expected = 1)
  
  #Test that this grabs the correct estimates: 
  testthat::expect_equal(
    object = cont_meet_in_middle$overlap$estimate_exp_omic[1],
    expected = as.numeric(coef(lm(feature_V63 ~ exposure1, 
                                  data = test_data))[2]))  
  
  #Test that this grabs the correct estimates: 
  testthat::expect_equal(
    object = cont_meet_in_middle$overlap$estimate_omic_out[1],
    expected = as.numeric(coef(glm(disease1 ~ feature_V63, 
                                  data = test_data, 
                                  family = binomial(link = "logit")))[2]))  
  
})
