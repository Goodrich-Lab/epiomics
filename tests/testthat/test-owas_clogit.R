
# Test that owas_clogit works ----------------
test_that("owas_clogit works", {
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
                        sex = rep(c("male", "female"),n_ids/2),
                        case_control_set = rep(c(1:(n_ids/2)), 2),
                        cc_status = c(rep(0,n_ids/2), rep(1,n_ids/2)),
                        age = rep(rnorm(n = n_ids/2, mean = 65, sd = 5) |> 
                                    round(digits = 0), 2),
                        weight =  rlnorm(n_ids, meanlog = 5.1, sdlog = 0.15),
                        exposure1 = rlnorm(n_ids, meanlog = 2.3, sdlog = 1),
                        exposure2 = rlnorm(n_ids, meanlog = 2.3, sdlog = 1),
                        exposure3 = rlnorm(n_ids, meanlog = 2.3, sdlog = 1))
  
  
  # Create Test Data
  test_data <- cbind(cov_out, omics_df)
  
  # Get names of omics
  colnames_omic_fts <- colnames(test_data)[grep("feature_",
                                                colnames(test_data))]
  
  
  # No covars
  owas_clogit_out <- owas_clogit(
    df = test_data,
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

