
# Test that owas_mixed_effects works ----------------
test_that("owas_mixed_effects works", {
  set.seed(4656)
  n_omic_ftrs = 10
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
                        exposure3 = rlnorm(n_ids, meanlog = 2.3, sdlog = 1), 
                        disease1 = sample(0:1, n_ids, replace=TRUE, prob=c(.9,.1)))
  
  
  # Create Test Data
  test_data <- cbind(cov_out, omics_df)
  
  # Get names of omics
  colnames_omic_fts <- colnames(test_data)[grep("feature_",
                                                colnames(test_data))]
  
  
  # Continuous exp -------------
  ## with covars and conf int------------
  owas_mixed_effects_out <- owas_mixed_effects(
    df = test_data,
    var = "cc_status",
    omics = colnames_omic_fts,
    random_effects = "1|case_control_set",
    covars = c("age", "sex"),
    var_exposure_or_outcome = "exposure",
    family = "gaussian",
    conf_int = TRUE)
  
  ## without covars and with conf int------------
  owas_mixed_effects_out <- owas_mixed_effects(
    df = test_data,
    var = "cc_status",
    omics = colnames_omic_fts,
    random_effects = "1|case_control_set",
    var_exposure_or_outcome = "exposure",
    family = "gaussian",
    conf_int = TRUE)

  ## with covars, without conf int------------
  owas_mixed_effects_out <- owas_mixed_effects(
    df = test_data,
    var = "cc_status",
    omics = colnames_omic_fts,
    random_effects = "1|case_control_set",
    covars = c("age", "sex"),
    var_exposure_or_outcome = "exposure",
    family = "gaussian",
    conf_int = FALSE)
  
  ## without covars, without conf int------------
  owas_mixed_effects_out <- owas_mixed_effects(
    df = test_data,
    var = "cc_status",
    omics = colnames_omic_fts,
    random_effects = "1|case_control_set",
    var_exposure_or_outcome = "exposure",
    family = "gaussian",
    conf_int = FALSE, 
    REML = TRUE)
  
  
  # Test continuous analysis
  ### Run single model ---------------
  fit <- lmerTest::lmer(feature_V1 ~ cc_status + (1|case_control_set), 
                     data = test_data, REML = TRUE) |>
    summary() |> 
    coef()
  
  res_manual <- fit[nrow(fit),] |> t() |> as.data.frame()
  
  est_manual <- res_manual$Estimate
  
  ### Test that single model and mixed model give same estimates
  est_fxn <- owas_mixed_effects_out[
    owas_mixed_effects_out$feature_name == "feature_V1", 
    ]$estimate
  
  testthat::expect_equal(object = est_fxn, 
                         expected = est_manual)
  
  # Categorical exp -------------------------
  ## with covars and conf int ---------------
  owas_mixed_effects_out <- owas_mixed_effects(
    df = test_data,
    var = "disease1",
    omics = colnames_omic_fts,
    random_effects = "1|case_control_set",
    # covars = c("age", "sex"),
    var_exposure_or_outcome = "outcome", 
    family = "binomial",
    conf_int = TRUE)
  
  ### Run single model ---------------
  fit <- lme4::glmer(disease1 ~ feature_V1 + (1| case_control_set), 
                     data = test_data, 
                     family = binomial(link = "logit")) |>
    summary() |> 
    coef()
  
  res_manual <- fit[nrow(fit),] |> t() |> as.data.frame()
  
  est_manual <- res_manual$Estimate
  
  ### Test that single model and mixed model give same estimates
  est_fxn <- owas_mixed_effects_out[owas_mixed_effects_out$feature_name == "feature_V1", ]$estimate
  
  testthat::expect_equal(object = est_fxn, 
                         expected = est_manual)
  
})