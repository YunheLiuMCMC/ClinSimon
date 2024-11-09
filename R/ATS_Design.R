##' ATS_Design( ) provides an Adaptive Threshold Simon Design (ATS Simon) method
##' for Simon’s two-stage design in oncology trials when the realized sample sizes
##' in the first stage and/or the second stage(s) are different from the planned
##' sample sizes in the first stage and/or the second stage(s). The Proposed ATS
##' Simon design aims to adhere to sample sizes of the original design, to that end,
##' this design updates the original thresholds of (r1, r) in the first and/or
##' the second stages to satisfy the type I error rate as the original planned
##' design (note: power will decrease if the realized sample size is smaller than
##' the original one).
##'
##' @title Adaptive Threshold Simon	Design
##'
##' @param n1 The planned number of patients in stage 1
##' @param n  The planned number of patients in stages 1 and 2
##' @param n1_star The actual number of patients in stage 1
##' @param n_star The actual total number of patients in stages 1 and 2
##' @param r1 Original design threshold in stage 1
##' @param r Original design threshold in stage 2
##' @param p0 Unacceptable efficacy rate
##' @param p1 Desirable efficacy rate
##' @param alpha Original type-I error rate
##'
##' @return a data frame includes the Adaptive Threshold Simon Design (ATS Simon)
##' first stage threshold r1*, second stage threshold r*, actual first stage
##' patients n1*, actual total sample sizes of the two stages patients n*, updated
##' type I error constraint alpha(n*), attained type-I error and Power, Average
##' sample size under null hypothesis EN(p0) and Probability of early termination
##' under null hypothesis PET(p0).
##'
##' @references Yunhe Liu, & Haitao Pan. (2024). \emph{Clinical Trial Design Methods
##' for Managing Under- and Over-Enrollment in Simon's Two-Stage Design, Submitted.}
##'
##' @examples
##' # Adaptive Threshold Simon Design Case 1
##' ATS_Design(19, 36, 17, 34, 3, 10, 0.20, 0.40, 0.1)
##' #                                 r1* r* n1* n* alpha(n*) Type I Power EN(p0) PET(p0)
##' # Adaptive Threshold Simon Design   3 10  17 34     0.091  0.059 0.847 24.669   0.549
##'
##' # Adaptive Threshold Simon Design Case 2
##' ATS_Design(14, 44, 11, 41, 3, 14, 0.25, 0.45, 0.1)
##' #                                 r1* r* n1* n* alpha(n*) Type I Power EN(p0) PET(p0)
##' # Adaptive Threshold Simon Design   2 14  11 41     0.088   0.06 0.854 27.344   0.455
##'
##' @export


ATS_Design <- function(n1,n,n1_star,n_star,r1,r,p0,p1,alpha) {
  # Compute P(Y<=x|n,p), where Y~Bin(n,p)
  cd <- function(x, p, n) {
    epsilon <- 0.00000001
    if (p < epsilon) {
      p <- epsilon
    }
    if (p > 1 - epsilon) {
      p <- 1 - epsilon
    }
    if (x < 0) {
      y <- 0
    } else if (x > n) {
      y <- 1
    } else {
      u <- floor(x)
      y <- pbinom(u, n, p)
    }
    return(y)
  }

  # Compute P(Y=x|n,p), where Y~Bin(n,p)
  pd <- function(x, p, n) {
    epsilon <- 0.00000001
    if (p < epsilon) {
      p <- epsilon
    }
    if (p > 1 - epsilon) {
      p <- 1 - epsilon
    }
    if (x < 0) {
      y <- 0
    } else if (x > n) {
      y <- 0
    } else if (x != floor(x)) {
      y <- 0
    } else {
      y <- dbinom(x, n, p)
    }
    return(y)
  }

  # Compute the power for given design (r1, r, n1, n) under response rate p
  power <- function(r1, r, n1, n, p) {
    y <- 0
    n2 <- n - n1
    for (i in (r1 + 1):n1) {
      y <- y + pd(i, p, n1) * (1 - cd(r - i, p, n2))
    }
    return(y)
  }

  #Given the planned sample sizes of n1 and n, the planned PET,
  #and the actual sample size of n1_star at the first stage,
  #this subroutine finds the largest integer r1 such that P(Y1<=r1 | p0,n1_star) is closest to planned PET.
  first_cutoff<-function(PET, n1, n, n1_star, p0){
    m <- -1
    flag <- 0
    while ((flag == 0) & (m < n1_star + 1)) {
      y <- cd(m, p0, n1_star)
      if (y < PET) {
        m <- m + 1
      } else {
        flag <- 1
      }
    }
    uu <- c(0, 0)
    PET_l <- cd(m-1, p0, n1_star)
    PET_u <- cd(m, p0, n1_star)
    if (abs(PET_l-PET) <  abs(PET_u-PET)){
      uu[1]<-m-1
    } else{
      uu[1]<-m
    }

    PET_modified<-cd(uu[1], p0, n1_star)
    uu[2] <- PET_modified
    return(uu)
  }


  new_design<-function(n1, n, n1_star, n_star, r1, r, p0, p1, alpha){
    PET <- pbinom(r1, n1, p0)
    uu <- first_cutoff(PET, n1, n, n1_star, p0)
    a <- uu[1]
    n2_star <- n_star - n1_star
    flag <- 0
    c <- a - 1

    critical_v<-qnorm(1-alpha/2)
    prop<-sqrt(n_star/(n))
    alpha_spending<-2-2*pnorm(critical_v/prop)
    alpha_spending<-min(alpha_spending,alpha)

    while ((flag == 0) & (c < n_star + 1)) {
      c <- c + 1
      al <- power(a, c, n1_star, n_star, p0)
      if (al <= alpha_spending) {
        ##Change flag, so we stop searching when achieving the target Type I error
        ##We can ensure the highest power under this situation
        flag <- 1
        aaa <- a
        aal <- al
        ppw <- power(a, c, n1_star, n_star, p1)
        aave <- n2_star * (1 - cd(a, p0, n1_star)) + n1_star
      }
    }
    vv <- c(0, 0, 0, 0, 0, 0, 0, 0, 0)
    if (flag == 1) {
      vv[1] <- aaa
      vv[2] <- c
      vv[3] <- n1_star
      vv[4] <- n_star
      vv[5] <- alpha_spending
      vv[6] <- aal
      vv[7] <- ppw
      vv[8] <- aave
      vv[9] <- uu[2]
    }
    Redesign <- round(as.data.frame(t(vv)),3)
    colnames(Redesign) <- c("r1*", "r*", "n1*", "n*", "alpha(n*)",
                            "Type I", "Power", "EN(p0)", "PET(p0)")
    rownames(Redesign) <- c("Adaptive Threshold Simon Design")
    return (Redesign)
  }
  ########################## Final Results #####################################
  AlterDesign <- new_design(n1, n, n1_star, n_star, r1, r, p0, p1, alpha)
  return(AlterDesign)
}








