
### dgm with endogenous availability process ###


source("functions/utillity.R")

trun <- function(x){
  if(x < 1 & x> -1){
    ans <- x
  } else if (x >= 1){
    ans <- 1
  } else if(x <= -1){
    ans <- -1
  }
  return(ans)
}

dgm_endoAvail <- function(sample_size, total_T, ft, beta_1, beta_2, gt, alpha,
                          gamma3, gamma4,
                          tau, pt, mean_var = 1, j_t = NULL, h_At = NULL) {
  # variable A represent the treatment level  (0 , 1 , 2)
  # variable A1 and A2 is indicator function whether A is equal 1 or 2 respectively
  # variable I is a binary indicating availability of the participant at timepoint t 
  # prob_t indicates the probability of someone available at time point t, assume constant for now
  # tau is E(I_t) which we assume is equal 1 (always available at any given time) unless specifiedd otherwise
  # if there is no exp_var provided, assume error is Normal 0,1
  # AA = 0.3
  
  
  if (is.null(j_t)) {
    j_t <- rep(mean_var, total_T)
  }
  
  if (is.null(h_At)) {
    h_At <- matrix(rep(c(1, 0, 0), total_T), ncol = dim(pt)[2], byrow = TRUE)
  }
  
  
  # treatment name(only for 2 level of treatment for now)
  a_val <- c(0, 1, 2)
  
  df_names <- c("userid", "time", "Y", "A", "A1", "A2", "I", "prob_A","prob_I", "eps", "m_At", "n_et")
  
  dta <- data.frame(matrix(NA, nrow = sample_size * total_T, ncol = length(df_names)))
  names(dta) <- df_names
  
  dta$userid <- rep(1:sample_size, each = total_T)
  dta$time <- rep(1:total_T, times = sample_size)
  
  for (t in 1:total_T) {
    prob_a <- pt[t,]
    # row index for the rows corresponding to day t for every subject
    row_index <- seq(from = t, by = total_T, length = sample_size)
    # generate Availability 
    
    if (t ==1) {
      # generate Availability 
      dta$I[row_index] <- rbinom(sample_size, 1, prob = tau[1])
      dta$prob_I[row_index] <- tau[1]
      dta$m_At[row_index] <- 0
      dta$n_et[row_index] <- 0
    }else{
      m_At <- n_et <- I_temp <- rep(NA, sample_size)
      prob_I_temp <- rep(NA, sample_size)
      A1_temp <- dta$A1[row_index_lag1]
      A2_temp <- dta$A2[row_index_lag1]
      eps_temp <- dta$eps[row_index_lag1]
      # for generate availability for each individual depends on previous trt
      for (i in 1:sample_size) {
        m_At[i] <-( (A1_temp[i] - prob_a[2] * tau[t]) + (A2_temp[i] - prob_a[3] * tau[t])) 
        n_et[i] <- trun(eps_temp[i])
        prob_I_temp[i] <- tau[t] + 
          gamma3 * m_At[i] +
          gamma4 * (n_et[i])
        stopifnot(prob_I_temp[i] <= 1)
        stopifnot(prob_I_temp[i] >= 0)
        I_temp[i] <- rbinom(1, 1, prob_I_temp[i])
      }
      dta$m_At[row_index] <- m_At
      dta$n_et[row_index] <- n_et
      dta$I[row_index] <- I_temp
      dta$prob_I[row_index] <- prob_I_temp
    }
    # generate treatment 
    dta$A[row_index] <- sample(x = a_val, size = sample_size, replace = TRUE, prob = prob_a)
    # set treatment to 0 when availability is 0 
    dta$A[row_index] <- dta$A[row_index] * dta$I[row_index]
    # record probability for each treatment
    dta$prob_A[row_index] <- ifelse(dta$A[row_index] == 0, prob_a[1], 
                                    ifelse(dta$A[row_index] == 1, prob_a[2], prob_a[3]))
    # create pseudo variable for indicator function
    dta$A1[row_index] <- ifelse(dta$A[row_index] == 1, 1, 0)
    dta$A2[row_index] <- ifelse(dta$A[row_index] == 2, 1, 0)
    # generate random error 
    h_Ai <- h_At[t, 1] + h_At[t, 2] * dta$A1[row_index]+h_At[t, 3] * dta$A2[row_index]
    dta$eps[row_index] <- rnorm(sample_size, mean = 0, sd = 1) * sqrt(j_t[t]) * sqrt(h_Ai)
    dta$Y[row_index] <- (dta$A1[row_index] * as.numeric(ft[t, ] %*% beta_1)  +
                           dta$A2[row_index] * as.numeric(ft[t, ] %*% beta_2) +
                           as.numeric(gt[t, ] %*% alpha) +   
                           dta$eps[row_index])
    
    row_index_lag1 <- row_index
  }
  
  return(dta)
}
