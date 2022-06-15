# epiomics: Functions for omics data analysis in observational studies

**epiomics** provides a collection of fast and flexible functions for the analysis of 'omics data in observational studies. 

## Installation

You can install epiomics from [GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("JAGoodrich/epiomics")
```

## 'Omics wide association study (owas) 

The basis of many omics analysis in epidemiology begin with an 'omics wide association study. The function `owas()` implements an omics wide association study with the option of using the 'omics data as either the dependent variable (i.e., for performing an exposure --> 'omics analysis) or using the 'omics as the independent variable (i.e., for performing an 'omics --> outcome analysis). `owas()` provides the option to adjust for covariates, and allows for either continuous or dichotomous outcomes.  


To illustrate this function, start with simulating data: 

