##' ATSS_Design_Stage2( ) provides an Adaptive Threshold and Sample Size Simon Design
##' (ATSS Simon) method for Simon's two stage design in oncology trials when the
##' realized sample sizes in the second stage is different from the planned
##' sample sizes in the second stage from interim analysis new design. Further
##' adjustment of the threshold at the second stage is needed. So, we update again
##' the second stage threshold r* to satisfy the type I error rate given
##' the interim analysis design first stage threshold r1* and actual two stages
##' sample sizes (n1*, n**).
##'
##' @title Adaptive Threshold and Sample Size Simon Design Two Stages
##'
##' @param p0 Unacceptable efficacy rate
##' @param p1 Desirable efficacy rate
##' @param r1_star Interim analysis design threshold in stage 1
##' @param n1_star The actual number of patients in stage 1
##' @param n_double_star The actual total number of patients in stages 1 and 2
##' @param alpha Original Type-I error rate
##'
##' @return a data frame includes the Adaptive Threshold and Sample Size Simon Design
##' interim analysis design adjusted first stage threshold r1*,
##' Adaptive Threshold and Sample Simon Design stage 2 new design adjusted second
##' stage threshold r*, actual number of patients in the first stage n1*,
##' actual total number of patients in stages 1 and 2 n**, attained Type-I error
##' and Power, Average sample size under null hypothesis EN(p0) and Probability
##' of early termination under null hypothesis PET(p0).
##'
##' @examples
##' # Adaptive Threshold and Sample Size Simon Design two stages analysis case 1
##' ATSS_Design_Stage2(0.05, 0.20, 1, 20, 33, 0.10)
##' #                     r1* r* n1* n** Type I Power EN(p0) PET(p0)
##' # ATSS_Design_Stage2   1  3  20  33   0.07 0.888 23.434   0.736
##'
##' # Adaptive Threshold and Sample Size Simon Design two stages analysis case 2
##' ATSS_Design_Stage2(0.10, 0.30, 2, 18, 24, 0.10)
##' #                   r1* r* n1* n** Type I Power EN(p0) PET(p0)
##' #ATSS_Design_Stage2   2  4  18  24   0.08 0.876 19.597   0.734
##'
##' @export


ATSS_Design_Stage2 <- function(p0,p1,r1_star,n1_star,n_double_star, alpha) {
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

  # Compute the power for given design (r1, r2, n1, n) under response rate p
  power <- function(r1, r2, n1, n, p) {
    y <- 0
    n2 <- n - n1
    for (i in (r1 + 1):n1) {
      y <- y + pd(i, p, n1) * (1 - cd(r2 - i, p, n2))
    }
    return(y)
  }

  ##Find c_star_star given r1_star, n1_star, and n_double_star such that
  ##P(y1>r1_star, Y>c_star_star|p0,n1_star,n_double_star)<=alpha.
  c_double_star <- function(p0,p1,r1_star,n1_star,n_double_star,alpha){
    n1 <- n1_star
    n <- n_double_star
    a <- r1_star
    flag <- 0
    ff <- 0
    cc <- a-1
    while ((ff == 0) & (cc < n + 1)) {
      cc <- cc+1
      al <- power(a,cc,n1,n,p0)
      pw <- power(a,cc,n1,n,p1)
      if (al<=alpha) {
        ave <- n1+(1-cd(a,p0,n1_star))*(n-n1)
        ff <- 1
        flag <- 1
        aaa <- a
        ccc <- cc
        nnn1 <- n1
        nnn <- n
        aal <- al
        ppw <- pw
        aave <- ave
        PET <- pbinom(a, n1, p0)
      }
    }
    vv <- array(0, dim = c(9,1))
    vv[1] <- flag
    if (flag == 1) {
      vv[2] <- aaa
      vv[3] <- ccc
      vv[4] <- nnn1
      vv[5] <- nnn
      vv[6] <- aal
      vv[7] <- ppw
      vv[8] <- aave
      vv[9] <- PET
    }
    Redesign <- round(as.data.frame(t(vv[-1,])),3)
    colnames(Redesign) <- c("r1*", "r*", "n1*", "n**",
                            "Type I", "Power", "EN(p0)", "PET(p0)")
    rownames(Redesign) <- c("ATSS_Design_Stage2")
    return (Redesign)
  }
  ########################## Final Results #####################################
  ATSS_Design_Stage2 <- c_double_star(p0,p1,r1_star,n1_star,n_double_star,alpha)
  return(ATSS_Design_Stage2)
}





