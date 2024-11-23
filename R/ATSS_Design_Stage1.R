##' @title Adaptive Threshold and Sample Size Simon Design Interim Analysis
##'
##' @description
##' ATSS_Design_Stage1( ) provides an Adaptive Threshold and Sample Size Simon Design
##' (ATSS Simon) method for Simon's two stage design in oncology trials when the
##' realized sample sizes in the first stage is different from the planned sample
##' sizes in the first stage. When under-enrollment or over-enrollment occurs at
##' the first stage, we identify the design parameters (r1*, r*, n*) based
##' on the actual sample size n1* ar the first stage to satisfy the type I error
##' rate and power. In addition, the ATSS Simon design also satisfies the other
##' criteria as in the originally planned design, such as minimizing the average
##' sample size under the null hypothesis H0.
##'
##' @param p0 Unacceptable efficacy rate
##' @param p1 Desirable efficacy rate
##' @param n1_star The actual number of patients in stage 1
##' @param alpha Original Type-I error rate
##' @param beta Original Type-II error rate
##'
##' @return a data frame includes the Adaptive Threshold and Sample Simon Design
##' interim analysis' adjusted first stage threshold r1*, second stage threshold
##' r*, actual number of patients in the first stage n1*, new design planned two
##' stages' patients n*, attained Type-I error rate and Power, Average sample
##' size under null hypothesis EN(p0) and Probability of early termination
##' under null hypothesis PET(p0).
##'
##' @examples
##' # Adaptive Threshold and Sample Size Simon Design interim analysis case 1
##' ATSS_Design_Stage1(0.05, 0.20, 20, 0.10, 0.10)
##' #                    r1* r* n1* n* Type I Power EN(p0) PET(p0)
##' # ATSS_Design_Stage1   1  3  20 35   0.08 0.901 23.962   0.736
##'
##' # Adaptive Threshold and Sample Size Simon Design interim analysis case 2
##' ATSS_Design_Stage1(0.10, 0.30, 18, 0.10, 0.10)
##' #                    r1* r* n1* n* Type I Power EN(p0) PET(p0)
##' # ATSS_Design_Stage1   2  4  18 26  0.099 0.904  20.13   0.734
##'
##' @export


ATSS_Design_Stage1 <- function(p0, p1, n1_star, alpha, beta) {

  # Compute P(Y<=x|n,p), where Y~Bin(n,p)
  cd <- function(x, p, n) {
    epsilon <- 1e-8

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
    epsilon <- 1e-8

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

  # Compute power for given design under response rate p
  power <- function(r1, r2, n1, n, p) {
    y <- 0
    n2 <- n - n1

    for (i in (r1 + 1):n1) {
      y <- y + pd(i, p, n1) * (1 - cd(r2 - i, p, n2))
    }

    return(y)
  }

  # Find largest integer satisfying beta constraint
  a_up <- function(p1, n1, beta) {
    flag <- 0
    m <- -1

    while ((flag == 0) & (m < n1 + 1)) {
      y <- cd(m, p1, n1)
      if (y <= beta) {
        m <- m + 1
      } else {
        flag <- 1
        uu <- m - 1
      }
    }

    return(uu)
  }

  # Find optimal design parameters:
  # Finds design (a_star,c_star,n1_star,n_star) for given n1_star at stage 1
  # such that type I and II error constraints are satisfied while minimizing
  # average sample size
  modify <- function(p0, p1, n1_star, alpha, beta) {
    flag <- 0
    ave <- 100
    n1 <- n1_star
    m <- 100

    for (n in (n1 + 1):m) {
      n2 <- n - n1
      u <- a_up(p1, n1, beta)

      for (i in 0:u) {
        aa <- (1 - cd(i, p0, n1)) * n2 + n1

        if (aa < ave) {
          ff <- 0
          cc <- i - 1

          while ((ff == 0) & (cc < n + 1)) {
            cc <- cc + 1
            al <- power(i, cc, n1, n, p0)
            pw <- power(i, cc, n1, n, p1)

            if ((pw >= 1 - beta) & (al <= alpha)) {
              ave <- aa
              ff <- 1
              flag <- 1
              aaa <- i
              ccc <- cc
              nnn1 <- n1
              nnn <- n
              aal <- al
              ppw <- pw
              aave <- ave
              PET <- pbinom(i, n1, p0)
            } else if (pw < 1 - beta) {
              cc <- n + 1
            }
          }
        }
      }
    }

    vv <- array(0, dim = c(9, 1))
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

    Redesign <- round(as.data.frame(t(vv[-1, ])), 3)
    colnames(Redesign) <- c(
      "r1*", "r*", "n1*", "n*",
      "Type I", "Power", "EN(p0)", "PET(p0)"
    )
    rownames(Redesign) <- c("ATSS_Design_Stage1")

    return(Redesign)
  }

  # Calculate final results
  Redesign_Stage_one <- modify(p0, p1, n1_star, alpha, beta)

  return(Redesign_Stage_one)
}
