
# Test that owas works with multiple outcomes ----------------
test_that("owas works with multiple traits", {
  # Load Example Data
  data("example_data")

  # Get names of omics
  colnames_omic_fts <- colnames(example_data)[grep("feature_",
                                                colnames(example_data))][1:10]

  # Names of exposures in mixture
   exposure_names = c("exposure1", "exposure2", "exposure3")

  # Run function without covariates
  out <- owas_qgcomp(df = example_data,
                     expnms = exposure_names,
                     omics = colnames_omic_fts,
                     q = 4,
                     confidence_level = 0.95)


  # Run analysis with covariates
  out <- owas_qgcomp(df = example_data,
                     expnms = c("exposure1", "exposure2", "exposure3"),
                     covars = c("weight", "age", "sex"),
                     omics = colnames_omic_fts,
                     q = 4,
                     confidence_level = 0.95)
  
  # Test that function returns expected dimensions
  testthat::expect_equal(object = ncol(out), 
                         expected = 11)
  testthat::expect_equal(object = nrow(out), 
                         expected = length(colnames_omic_fts))
  
})
