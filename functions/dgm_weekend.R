### dgm with weekend effect on gt  ###
source("functions/utillity.R")


is_weekend <- function(t, n_obs = 1){
  # this function tells whether time points is weekend or weekday
  # return 1 if its weekend, 0 otherwise
  # t is the time point input 
  # n_obs is the number of observation perday (see Liao et all)
  
  # convert t into week 1 (since its periodical)
  t_mod <- t %% (n_obs * 7)
  
  # start on monday
  # if t_mod is greater than 1 or less or equal to n_obs * 5 its weekdays, or else weekend
  ans <- ifelse(t_mod <= (n_obs * 5) & t_mod >= 1, 0, 1)
  
  return(ans)
}

dgm_weekend<- function(sample_size, total_T, ft, beta_1, beta_2, gt, alpha,
                          eta1, eta2, gamma2,
                          tau, pt, mean_var = 1, n_obs = 1) {
  # variable A represent the treatment level  (0 , 1 , 2)
  # variable A1 and A2 is indicator function whether A is equal 1 or 2 respectively
  # variable I is a binary indicating availability of the participant at timepoint t 
  # prob_t indicates the probability of someone available at time point t, assume constant for now
  # tau is E(I_t) which we assume is equal 1 (always available at any given time) unless specifiedd otherwise
  # if there is no exp_var provided, assume error is Normal 0,1
  # Y_noEffect for debugging purposes

  # treatment name(only for 2 level of treatment for now)
  a_val <- c(0, 1, 2)
  
  df_names <- c("userid", "time", "Y", "A", "A1", "A2", "I", "prob_A","prob_I", "eps", "weekend", "Y_noEffect")
  
  dta <- data.frame(matrix(NA, nrow = sample_size * total_T, ncol = length(df_names)))
  names(dta) <- df_names
  
  dta$userid <- rep(1:sample_size, each = total_T)
  dta$time <- rep(1:total_T, times = sample_size)
  
  for (t in 1:total_T) {
    prob_a <- pt[t,]
    # row index for the rows corresponding to day t for every subject
    row_index <- seq(from = t, by = total_T, length = sample_size)
    
    # check if its weekdays or weekend 
    dta$weekend[row_index] <- is_weekend(t, n_obs = 1)
    
    # generate Availability 
    dta$I[row_index] <- rbinom(sample_size, 1, prob = tau[t])
    dta$prob_I[row_index] <- ifelse(dta$I[row_index] == 0, 1 - tau[t], tau[t])
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
    dta$eps[row_index] <- rnorm(sample_size, mean = 0, sd = 1) 
    
    
    dta$Y_noEffect[row_index] <- dta$A1[row_index] * as.numeric(ft[t, ] %*% beta_1)  +
                           dta$A2[row_index] * as.numeric(ft[t, ] %*% beta_2) +
                           as.numeric(gt[t, ] %*% alpha) +   
                           dta$eps[row_index]
    
    dta$Y[row_index] <- dta$A1[row_index] * as.numeric(ft[t, ] %*% beta_1)  +
      dta$A2[row_index] * as.numeric(ft[t, ] %*% beta_2) +
      as.numeric(gt[t, ] %*% alpha)+ gamma2 * dta$weekend[row_index] +   
      dta$eps[row_index]
  }
  
  return(dta)
}