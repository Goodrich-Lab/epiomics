# epiomics: Functions for omics data analysis in observational studies

**epiomics** provides a collection of fast and flexible functions for the analysis of 'omics data in observational studies.

## Installation

You can install epiomics from [GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("JAGoodrich/epiomics")

library(epiomics)
```

## Omics wide association study (owas)

The basis of many omics analysis in epidemiology begin with an omics wide association study. The function `owas()` implements an omics wide association study with the option of using the 'omics data as either the dependent variable (i.e., for performing an exposure --\> 'omics analysis) or using the 'omics as the independent variable (i.e., for performing an 'omics --\> outcome analysis). `owas()` provides the option to adjust for covariates, and allows for either continuous or dichotomous outcomes. `owas()` can also handle multiple variables of interest (ie, multiple exposures or multiple traits). 

Start with simulating data:

``` r
# Simulate dataset
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
                                               

# Get names of traits
trait_nms = c("disease1", "disease2")
```

### Run owas with continuous exposure as the variable of interest

``` r
owas(df = test_data, 
     var = "exposure1", 
     omics = colnames_omic_fts, 
     covars = c("age", "sex"), 
     var_exposure_or_outcome = "exposure", 
     family = "gaussian")
     
# Equivalent: 
owas(df = test_data, 
     var = "exposure1", 
     omics = colnames_omic_fts, 
     covars = c("age", "sex"), 
     var_exposure_or_outcome = "exposure")  
```

### Run owas with dichotomous outcome as the variable of interest
``` r
owas(df = test_data, 
     var = "disease1", 
     omics = colnames_omic_fts, 
     covars = c("age", "sex"), 
     var_exposure_or_outcome = "outcome", 
     family = "binomial")
```

### Run owas with multiple continuous exposures as the variable of interest
``` r
# Get names of exposures
expnms = c("exposure1", "exposure2", "exposure3")

owas(df = test_data, 
     var = expnms, 
     omics = colnames_omic_fts, 
     covars = c("age", "sex"), 
     var_exposure_or_outcome = "exposure", 
     family = "gaussian")
```


## Meet in the Middle

The function `meet_in_middle()` conducts meet in the middle screening between an exposure, omics, and an outcome, as described by Cadiou et al., 2021. This function provides the option to adjust for covariates, and allows for either continuous or dichotomous outcomes. Examples are based on the simulated data created above. 

### Meet in the middle with a dichotomous outcome

``` r
res <- meet_in_middle(df = test_data,
                      exposure = "exposure1", 
                      outcome = "disease1", 
                      omics = colnames_omic_fts,
                      covars = c("age", "sex"), 
                      outcome_family = "binomial")
res
``` 


### Meet in the middle with a continuous outcome 

``` r
res <- meet_in_middle(df = test_data,
                      exposure = "exposure1", 
                      outcome = "weight", 
                      omics = colnames_omic_fts,
                      covars = c("age", "sex"), 
                      outcome_family = "gaussian")

res
``` 


### Meet in the middle with a continuous outcome and no covariates

``` r 
res <- meet_in_middle(df = test_data,
                      exposure = "exposure1", 
                      outcome = "weight", 
                      omics = colnames_omic_fts,
                      outcome_family = "gaussian")
res
```
