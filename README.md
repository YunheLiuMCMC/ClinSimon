# ClinSimon

## Intended Usage
Implements methods for managing enrollment challenges in Simon's 
    Two-Stage Design for phase II clinical trials. Provides functions for:
    (1) adaptive threshold adjustment to handle under-enrollment,
    (2) sample size recalibration for over-enrollment scenarios, and
    (3) comprehensive post-inference analysis tools including power calculation,
    confidence intervals, and hypothesis testing. These methods help maintain 
    statistical rigor while adapting to real-world enrollment variations in
    clinical trials.

## Installation Instruction

My repository is public, so you can use following R code to download:   
```r
devtools::install_github("YunheLiuMCMC/ClinSimon")
```
After making this package Vignette, you can use following R code to download:
```r
devtools::install_github("YunheLiuMCMC/ClinSimon"", build_vignettes = TRUE)
```

## What is Left

- Write the function ATSS_Design( ) to provide an Adaptive Threshold and Sample Size Simon Design

- Write the function SimonAnalysis( ) to calculate the Uniformly minimum-variance
unbiased estimator (UMVUE), Confidence Intervals (Clopper-Pearson, Jung exact, 
and Mid-$p$) and $p$-Value given the design parameters obtained from the Adaptive Threshold Simon Design (ATS Simon)

- Make the Vignette
