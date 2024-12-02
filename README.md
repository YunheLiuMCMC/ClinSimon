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

If you want to install the package with a vignette, you should install the `clinfun`
package first. Then you can install it with the vignette built:
```r
#Install the package clinfun first:
install.packages("clinfun")
# Install ClinSimon with vignette support:
devtools::install_github("YunheLiuMCMC/ClinSimon", build_vignettes = TRUE)
```

## Example
```r
# Load the ClinSimon package
library(ClinSimon)

# Generate a Simon two stage design from package Clinfun
library(Clinfun)
# Specify the parameters and constraints
trial = ph2simon(0.25, 0.45, 0.1, 0.1)
# Print
trial

# Simon 2-stage Phase II design 

# Unacceptable response rate:  0.25 
# Desirable response rate:  0.45 
# Error rates: alpha =  0.1 ; beta =  0.1 

#            r1 n1  r  n EN(p0) PET(p0)   qLo   qHi
# Minimax     5 23 13 39  31.50  0.4685 0.752 1.000
# Admissible  3 15 13 40  28.47  0.4613 0.026 0.752
# Optimal     3 14 14 44  28.36  0.5213 0.000 0.026
```

Right now, one clinician decided to use Optimal deisgn to recruit patients.
Suppose we currently have outcome data for only 11 patients in the $1^{st}$ stage
and assume the $2^{nd}$ stage sample size for evaluable patients remains as 
planned, We can use `ATS_Design()` to provide updated thresholds for the 
interim analysis, without	needing	to wait for the	number of evaluable	patients 
to reach 14. This means the input *`n1_star`* is 11 and *`n_star`* is 41 (=11+30 
instead of the original 14+30 as the total sample size), as the original optimal
design specifies a sample size of 30 for the second stage.  

```r
ATS_Design(n1=14,n=44,n1_star=11,n_star=41,r1=3,r=14,p0=0.25,p1=0.45,alpha=0.1)

#                                 r1* r* n1* n* alpha(n*) Type I Power EN(p0) PET(p0)
# Adaptive Threshold Simon Design   2 14  11 41     0.088   0.06 0.854 27.344   0.455
```

The updated design parameters by ATS Simon design method are (2, 14). We also 
output actual sample sizes of stage 1 and the total sample size (11, 41). The 
type I error rate of this updated design is 0.06, which is below the original 
type I error constraint of 0.1. As we know, if alpha decreases, power will also 
decrease. Therefore, the current power of 85.4% is lower than the original 
design's 90% as expected, but still close to the original power. Additionally, 
the probability of early termination is 0.455, which is close to the original 
design's 0.521.  


## Vignette
For a detailed overview and examples, check out the [ClinSimon vignette](https://github.com/YunheLiuMCMC/ClinSimon/blob/main/vignettes/ClinSimon.pdf).

## License
This package is free and open source software, licensed under GPL
(&gt;=3). See the [LICENSE](https://github.com/YunheLiuMCMC/ClinSimon/blob/main/LICENSE.md) 
file for more details.

