
# Test that owas works with multiple outcomes ----------------
test_that("owas works with multiple traits", {
  # Load Example Data
  data("example_data")
  
  # Get names of omics
  colnames_omic_fts <- colnames(example_data)[grep("feature_",
                                                   colnames(example_data))]
  
  # Run function with multiple exposures as the variable of interest
  expnms = c("exposure1", "exposure2", "exposure3")
  
  # Continuous exposure, covars
  cont_owas_out <- owas(
    df = example_data,
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
  owas(df = example_data,
       var = expnms,
       omics = colnames_omic_fts,
       # covars = c("age", "sex"),
       var_exposure_or_outcome = "exposure",
       family = "gaussian")
  
  # Continuous outcome, with covars
  owas(df = example_data,
       var = expnms,
       omics = colnames_omic_fts,
       covars = c("age", "sex"),
       var_exposure_or_outcome = "outcome",
       family = "gaussian")
  
  # Continuous outcome, no covars
  owas(df = example_data,
       var = expnms,
       omics = colnames_omic_fts,
       # covars = c("age", "sex"),
       var_exposure_or_outcome = "outcome",
       family = "gaussian")
  
  
})


# Test that owas works with single outcome ----------------
test_that("owas works with single continuous outcome", {
  # Load Example Data
  data("example_data")
  
  # Get names of omics
  colnames_omic_fts <- colnames(example_data)[grep("feature_",
                                                   colnames(example_data))]
  
  # Run function with continuous exposure as the variable of interest
  cont_owas_out <- owas(df = example_data,
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
  # owas(df = example_data,
  #                var = "disease1",
  #                omics = colnames_omic_fts,
  #                covars = c("age", "sex"),
  #                var_exposure_or_outcome = "outcome",
  #                family = "binomial")
})


# Test that meet in middle works ----------------
test_that("meet_in_middle works", {
  
  # Load Example Data
  data("example_data")
  
  # Get names of omics
  colnames_omic_fts <- colnames(example_data)[grep("feature_",
                                                   colnames(example_data))]
  
  # Run function with continuous exposure as the variable of interest
  cont_meet_in_middle <- meet_in_middle(df = example_data,
                                        exposure = "exposure1", 
                                        outcome = "disease1", 
                                        omics = colnames_omic_fts,
                                        # covars = c("age", "sex"), 
                                        outcome_family = "binomial", 
                                        confidence_level = 0.8)
  #############################################################
  # Test that function returns expected dimensions
  testthat::expect_equal(object = ncol(cont_meet_in_middle$overlap), 
                         expected = 15)
  testthat::expect_equal(object = nrow(cont_meet_in_middle$overlap), 
                         expected = 4)
  
  #Test that this grabs the correct estimates: 
  testthat::expect_equal(
    object = cont_meet_in_middle$overlap$estimate_exp_omic[1],
    expected = as.numeric(coef(lm(feature_V47 ~ exposure1, 
                                  data = example_data))[2]))  
  
  #Test that this grabs the correct estimates: 
  testthat::expect_equal(
    object = cont_meet_in_middle$overlap$estimate_omic_out[1],
    expected = as.numeric(
      coef(glm(disease1 ~ feature_V47, 
               data = example_data, 
               family = binomial(link = "logit")))[2]))  
  
})



# Test that owas gives correct errors ----------------
test_that("owas gives correct errors", {
  # Load Example Data
  data("example_data")
  
  # Get names of omics
  colnames_omic_fts <- colnames(example_data)[grep("feature_",
                                                   colnames(example_data))]
  
  # Run function with multiple exposures as the variable of interest
  disease_names = c("diseasecat1", "diseasecat2", "diseasecat3") 
  
  ## Error that ref_group must be specified ----
  error_message <- testthat::capture_error(
    cont_owas_out <- owas(df = example_data,
                          var = disease_names,
                          omics = colnames_omic_fts,
                          covars = c("age", "sex"),
                          var_exposure_or_outcome = "exposure",
                          family = "gaussian")
  )
  
  # Test that function returns expected dimensions
  testthat::expect_equal(object = error_message$message, 
                         expected = 'If var is character or factor, ref_group must be specified')
  
  
  # Erro: two categories
  error_message <- testthat::capture_error(
    cont_owas_out <- owas(df = example_data,
                          var = disease_names,
                          omics = colnames_omic_fts,
                          covars = c("age", "sex"),
                          var_exposure_or_outcome = "exposure",
                          family = "gaussian", 
                          ref_group = "healthy")
  )
  
  # Test that function returns expected dimensions
  testthat::expect_equal(object = error_message$message, 
                         expected = 'Currently var can only contain a maximum of two unique categories')
  
  # Error: same type
  error_message <- testthat::capture_error(
    cont_owas_out <- owas(df = example_data,
                          var = c(disease_names, "disease1"),
                          omics = colnames_omic_fts,
                          covars = c("age", "sex"),
                          var_exposure_or_outcome = "exposure",
                          family = "gaussian", 
                          ref_group = "healthy")
  )
  
  # Test that function returns expected dimensions
  testthat::expect_equal(object = error_message$message, 
                         expected = 'All variables in \'var\' must be the same type')
  
  
  
  
})
