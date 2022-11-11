
# Test that owas_clogit works ----------------
test_that("owas_clogit works", {
  
  # Load Example Data
  data("example_data")
  
  # Get names of omics
  colnames_omic_fts <- colnames(example_data)[grep("feature_",
                                                colnames(example_data))]
  
  
  # No covars
  owas_clogit_out <- owas_clogit(
    df = example_data,
    cc_status = "cc_status",
    cc_set = "case_control_set", 
    omics = colnames_omic_fts, 
    covars = NULL,
    confidence_level = 0.95, 
    conf_int = FALSE)
  
  # Test that function returns expected dimensions
  testthat::expect_equal(object = ncol(owas_clogit_out), 
                         expected = 8)
  testthat::expect_equal(object = nrow(owas_clogit_out), 
                         expected = length(colnames_omic_fts))
})

