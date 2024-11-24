# Introduction
The `ClinSimon` package serves three purposes, including:

* Providing an Adaptive Threshold Simon Design (ATS Simon) method for Simon's 
two-stage design in oncology trials when the realized sample sizes in the $1^{st}$ 
and/or $2^{nd}$ stage(s) are different from the planned sample sizes in the $1^{st}$ 
and/or $2^{nd}$ stage(s). The proposed ATS Simon design tries to follow 
sample sizes of the original design, to that end, this design updates the original
thresholds of $(r_1, r)$ in the $1^{st}$ and/or the $2^{nd}$ stage(s) to satisfy the 
type I error rate as the original planned design (note: power will decrease if the
realized sample size is smaller than the original one).

* Providing an Adaptive Threshold and Sample Size Simon Design (ATSS Simon) method 
for Simon's two-stage design in oncology trials when the realized sample sizes 
in the $1^{st}$ and/or $2^{nd}$ stage(s) are different from the planned 
sample sizes in the $1^{st}$ and/or $2^{nd}$ stage(s). The proposed ATSS 
Simon method updates not only the original threshold of $(r_{1}^{*}, r^*)$ 
but the original sample sizes of $(n_{1}^{*}, n^*)$ to satisfy the type I error 
rate and power requirements as the original planned design (note:	unlike the 
ATS Simon design, the sample size here will also be updated to satisfy the 
original power). In addition, the ATSS Simon design also satisfies the other 
criteria as in the originally planned design, such as minimizing the average 
sample size under the null hypothesis *$H_0$*.

* Providing	comprehensive	post-trial inference tools at the end of the trial 
for Simon's two-stage design when under- or over-enrollment occurs. This 
includes the computation of point estimates, confidence intervals, and 
p-values for the proposed ATS and ATSS Simon design methods.

This vignette introduces functionalities in {ClinSimon} tailored for	
designing	single-arm clinical trials using the proposed method. Examples featuring
Simon's	two-stage designs from {clinfun} demonstrate the application of these 
functions within our package.

To begin, install and load {ClinSimon} and {clinfun}:

```r
library(clinfun)
library(ClinSimon)
```






