---
output: html_document
---

# Final Project - R package

## Intended Usage
Providing an Adaptive Threshold Simon Design (ATS Simon) method for Simon's 
two-stage design in oncology trials when the realized sample sizes in the $1^{st}$ 
and/or $2^{nd}$ stage(s) are different from the planned sample sizes in the $1^{st}$ 
and/or $2^{nd}$ stage(s). The proposed ATS Simon design tries to follow 
sample sizes of the original design, to that end, this design updates the original
thresholds of $(r_1, r)$ in the $1^{st}$ and/or the $2^{nd}$ stage(s) to satisfy the 
type I error rate as the original planned design (note: power will decrease if the
realized sample size is smaller than the original one).

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
