

SimonAnalysis <- function(m,s,n1,n2,r1,r,alpha,quantile,CI_option,p0) {
  ############################# UMVUE Estimate #################################
  sum_numer<-function(lower,upper,n1,n2,s){
    sum<-0
    for (x1 in lower:upper){
      sum <- sum + choose(n1-1, x1-1)*choose(n2, s-x1)
    }

    return(sum)
  }

  sum_denom<-function(lower,upper,n1,n2,s){
    sum<-0
    for (x1 in lower:upper){
      sum <- sum + choose(n1, x1)*choose(n2, s-x1)
    }

    return(sum)
  }

  point_est<-function(m,s,n1,n2,r1,r){
    if(m==1&s>r1){
      response_est <- paste("The number of response people should be
                            less than stage 1 cut-off if stopping stage is 1.")
    }
    else if(m==2&s<r1){
      response_est <- paste("The number of response people should be
                            larger than stage 1 cut-off if stopping stage is 2.")
    }
    else if(m==1&s<=r1){
      response_est <- s/n1
      #response_est <- round(response_est,3)
    }
    else if(m==2&s>r1){
      upper <- min(s, n1)
      lower <- max(r1+1, s-n2)
      response_est <- sum_numer(lower,upper,n1,n2,s)/sum_denom(lower,upper,n1,n2,s)
      #response_est <- round(response_est,3)
    }
    else{
      response_est <- paste("The input number of stage is over 2")
    }
    return(response_est)
  }

  ############################# Clopper_Pearson CI #############################
  sum_pearson<-function(lower_bound,upper_bound,n1,n2,alpha,p){
    sum<-0
    for (k in lower_bound:upper_bound){
      sum <- sum + choose(n1+n2, k)*p^k*(1-p)^(n1+n2-k)
    }
    sum <- sum-alpha/2
    return(sum)
  }

  CP_CI_est<-function(s,n1,n2,alpha){
    lower_left  <- s
    upper_left  <- n1+n2
    lower_root  <- stats::uniroot(sum_pearson, interval = c(0, 1), lower_bound=lower_left,
                                  upper_bound=upper_left, n1=n1, n2=n2, alpha=alpha)
    lower_bound <- lower_root$root
    lower_bound <- round(lower_bound,3)

    lower_right <- 0
    upper_right <- s
    upper_root  <- stats::uniroot(sum_pearson, interval = c(0, 1), lower_bound=lower_right,
                                  upper_bound=upper_right, n1=n1, n2=n2, alpha=alpha)
    upper_bound <- upper_root$root
    upper_bound <- round(upper_bound,3)

    interval_est <- c(lower_bound, upper_bound)
    return(interval_est)
  }

  ############################# Jung's exact CI ################################
  sum_joint<-function(lower,upper,n1,n2,s,p){
    sum<-0
    for (x1 in lower:upper){
      sum <- sum + choose(n1, x1)*choose(n2, s-x1)
    }
    sum <- sum*p^s*(1-p)^(n1+n2-s)
    return(sum)
  }

  joint_prob<-function(m,s,n1,n2,r1,r,p){
    if(m==1){
      prob <- choose(n1,s)*p^s*(1-p)^(n1-s)
    }
    else if(m==2){
      upper <- min(n1, s)
      lower <- r1+1
      prob <- sum_joint(lower,upper,n1,n2,s,p)
    }
    else{
      prpb <- paste("The input number of stage is over 2")
    }

    return(prob)
  }

  Jung_conf_interval<-function(m,s,n1,n2,r1,r,limit,p){
    matrix_ms <- matrix(nrow = n1+n2+1, ncol = 4)
    matrix_ms[,1] <- c(rep(1, times=r1+1), rep(2, times=n1+n2-r1))
    matrix_ms[,2] <- c(0:r1, (r1+1):(n1+n2))

    for (i in 1:(n1+n2+1)){
      stop_stage <- matrix_ms[i,1]
      resp_total <- matrix_ms[i,2]
      matrix_ms[i,3] <- point_est(stop_stage,resp_total,n1,n2,r1,r)
    }

    for (j in 1:(n1+n2+1)){
      stop_stage <- matrix_ms[j,1]
      resp_total <- matrix_ms[j,2]
      matrix_ms[j,4] <- joint_prob(stop_stage,resp_total,n1,n2,r1,r,p)
    }

    umvue <- point_est(m,s,n1,n2,r1,r)
    subset_matrix_ms <- matrix_ms[matrix_ms[,3] >=  umvue, ]
    estimate <- sum(subset_matrix_ms[,4])-limit
    return(estimate)
  }

  Jung_CI_est<-function(m,s,n1,n2,r1,r,quantile){
    lower_root  <- stats::uniroot(Jung_conf_interval, interval = c(0, 1), m=m, s=s, n1=n1,
                                  n2=n2, r1=r1, r=r, limit=quantile[1])
    lower_bound <- lower_root$root
    lower_bound <- round(lower_bound,3)

    upper_root  <- stats::uniroot(Jung_conf_interval, interval = c(0, 1), m=m, s=s, n1=n1,
                                  n2=n2, r1=r1, r=r, limit=quantile[2])
    upper_bound <- upper_root$root
    upper_bound <- round(upper_bound,3)

    interval_est <- c(lower_bound, upper_bound)
    return(interval_est)
  }

  ############################# Mid-p CI #######################################
  MIDp_conf_interval<-function(m,s,n1,n2,r1,r,limit,p){
    matrix_ms <- matrix(nrow = n1+n2+1, ncol = 4)
    matrix_ms[,1] <- c(rep(1, times=r1+1), rep(2, times=n1+n2-r1))
    matrix_ms[,2] <- c(0:r1, (r1+1):(n1+n2))

    for (i in 1:(n1+n2+1)){
      stop_stage <- matrix_ms[i,1]
      resp_total <- matrix_ms[i,2]
      matrix_ms[i,3] <- point_est(stop_stage,resp_total,n1,n2,r1,r)
    }

    for (j in 1:(n1+n2+1)){
      stop_stage <- matrix_ms[j,1]
      resp_total <- matrix_ms[j,2]
      matrix_ms[j,4] <- joint_prob(stop_stage,resp_total,n1,n2,r1,r,p)
    }

    umvue <- point_est(m,s,n1,n2,r1,r)
    subset_matrix_ms <- matrix_ms[matrix_ms[,3] >  umvue, ]
    real_ms          <- 1/2*matrix_ms[matrix_ms[,1]==m & matrix_ms[,2]==s, ]
    subset_matrix_ms <- rbind(subset_matrix_ms,real_ms)
    estimate <- sum(subset_matrix_ms[,4])-limit
    return(estimate)
  }

  MIDp_CI_est<-function(m,s,n1,n2,r1,r,quantile){
    lower_root  <- stats::uniroot(MIDp_conf_interval, interval = c(0, 1), m=m, s=s, n1=n1,
                                  n2=n2, r1=r1, r=r, limit=quantile[1])
    lower_bound <- lower_root$root
    lower_bound <- round(lower_bound,3)

    upper_root  <- stats::uniroot(MIDp_conf_interval, interval = c(0, 1), m=m, s=s, n1=n1,
                                  n2=n2, r1=r1, r=r, limit=quantile[2])
    upper_bound <- upper_root$root
    upper_bound <- round(upper_bound,3)

    interval_est <- c(lower_bound, upper_bound)
    return(interval_est)
  }

  ############################# p VALUE ########################################
  p_cal<-function(m,s,n1,n2,r1,r,p0){
    umvue <- point_est(m,s,n1,n2,r1,r)
    p <- p0
    matrix_ms <- matrix(nrow = n1+n2+1, ncol = 4)
    matrix_ms[,1] <- c(rep(1, times=r1+1), rep(2, times=n1+n2-r1))
    matrix_ms[,2] <- c(0:r1, (r1+1):(n1+n2))

    for (i in 1:(n1+n2+1)){
      m <- matrix_ms[i,1]
      s <- matrix_ms[i,2]
      matrix_ms[i,3] <- point_est(m,s,n1,n2,r1,r)
    }

    for (j in 1:(n1+n2+1)){
      m <- matrix_ms[j,1]
      s <- matrix_ms[j,2]
      matrix_ms[j,4] <- joint_prob(m,s,n1,n2,r1,r,p)
    }

    subset_matrix_ms <- matrix_ms[matrix_ms[,3] >= umvue, ]
    p_value <- sum(subset_matrix_ms[,4])
    p_value <- round(p_value,4)
    return(p_value)
  }

  ############################# Final Results ##################################
  umvue <- point_est(m,s,n1,n2,r1,r)
  CP_CI <- CP_CI_est(s,n1,n2,alpha)
  Jung_CI <- Jung_CI_est(m,s,n1,n2,r1,r,quantile)
  MIDp_CI <- MIDp_CI_est(m,s,n1,n2,r1,r,quantile)
  p_VALUE <- p_cal(m,s,n1,n2,r1,r,p0)

  if(CI_option == "CP"){
    CI <- CP_CI
  } else if(CI_option == "Jung"){
    CI <- Jung_CI
  } else if(CI_option == "MIDp"){
    CI <- MIDp_CI
  } else{
    CI <- NA
  }

  Inference <- c(umvue, CI, p_VALUE)
  Inference <- round(as.data.frame(t(Inference)),3)
  colnames(Inference) <- c("UMVUE", "CI(lower)", "CI(upper)", "p_Val")
  rownames(Inference) <- c("Post-Trial Inference")
  return(Inference)
}


