
# Test that owas works with multiple outcomes ----------------
test_that("meet_in_middle works", {
  # Load Example Data
  data("example_data")
  
  # Get names of omics
  colnames_omic_fts <- colnames(example_data)[grep("feature_",
                                                   colnames(example_data))][1:10]
  
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