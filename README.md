# ClinSimon

## Overview

**ClinSimon** provides methods to address enrollment challenges in Simon's Two-Stage Design for phase II clinical trials. The package offers tools to:

1. **Adaptive threshold adjustment**: Manage under-enrollment scenarios.
2. **Sample size recalibration**: Handle over-enrollment cases.
3. **Post-inference analysis**: Perform power calculations, confidence interval estimation, and hypothesis testing.

These methods help maintain statistical rigor while accommodating real-world enrollment variations in clinical trials.

## Installation Instructions

You can install the **ClinSimon** package directly from GitHub using the following R code:
```r
# Uncomment and run the following line if you don't already have devtools installed:
# install.packages("devtools")

# Install ClinSimon from GitHub:
devtools::install_github("YunheLiuMCMC/ClinSimon")
```

If the package includes a vignette, you can install it with the vignette built:
```r
# Install ClinSimon with vignette support:
devtools::install_github("YunheLiuMCMC/ClinSimon", build_vignettes = TRUE)
```

## What is Left

- Write the function ATSS_Design( ) to provide an Adaptive Threshold and Sample Size Simon Design

- Write the function SimonAnalysis( ) to calculate the Uniformly minimum-variance
unbiased estimator (UMVUE), Confidence Intervals (Clopper-Pearson, Jung exact, 
and Mid-$p$) and $p$-Value given the design parameters obtained from the Adaptive Threshold Simon Design (ATS Simon)

- Make the Vignette
