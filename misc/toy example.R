rm(list = ls(all = TRUE))
source("wcls_categorical_treatment.R")

library(rootSolve)

# --------------------- Make synthetic data -------------------------------------- 

#' Data generative model for a simple demonstration of wcls estimator with multilevel treatment
#'
#' @param sample_size number of individual in data generative model 
#' @param total_T number of observation for each individual 
#'
#' @return long dataset generated
#' @export
#'
#' @examples
dgm_1 <- function(sample_size, total_T) {
  # variable A represent the treatment level  (0 , 1 , 2)
  # variable A1 and A2 is indicator function whether A is equal 1 or 2 respectively
  # variable I is a binary indicating availability of the participant at timepoint t 
  # prob_t indicates the probability of someone available at time point t, assume constant for now
  # S_t is moderator variable, that takes three values with equal probability

  baseline_Y_S0 <- 0.2
  baseline_Y_S1 <- 0.5
  baseline_Y_S2 <- 0.4
  
  
  # intercept
  beta_1 <- 0.1 # when At = 1 
  beta_3 <- 0.45  # when At = 2
  
  # Moderator effect 
  beta_2 <- 0.3 # when At = 1
  beta_4 <- 0.1 # when At = 2
  
  
  prob_a <- c(0.2, 0.5, 0.3)
  prob_I <- 1
  a_val <- c(0, 1, 2)
  
  
  df_names <- c("userid", "time", "Y", "A", "A1", "A2", "I", "S", "S2", "prob_A","prob_I","Y_A0", "eps")
  
  dta <- data.frame(matrix(NA, nrow = sample_size * total_T, ncol = length(df_names)))
  names(dta) <- df_names
  
  dta$userid <- rep(1:sample_size, each = total_T)
  dta$time <- rep(1:total_T, times = sample_size)
  
  for (t in 1:total_T) {
    # row index for the rows corresponding to day t for every subject
    row_index <- seq(from = t, by = total_T, length = sample_size)
    dta$A[row_index] <- sample(x = a_val, size = sample_size, replace = TRUE, prob = prob_a)
    dta$prob_A[row_index] <- ifelse(dta$A[row_index] == 0, prob_a[1], 
                                    ifelse(dta$A[row_index] == 1, prob_a[2], prob_a[3]))
    dta$A1[row_index] <- ifelse(dta$A[row_index] == 1, 1, 0)
    dta$A2[row_index] <- ifelse(dta$A[row_index] == 2, 1, 0)
    dta$I[row_index] <- rbinom(sample_size, 1, prob = prob_I)
    dta$prob_I[row_index] <- ifelse(dta$I[row_index] == 0, 1 - prob_I, prob_I)
    dta$S[row_index] <- sample(c(0,1,2), sample_size, replace = TRUE)
    dta$S2[row_index] <- ifelse(dta$S[row_index] == 2, 1, 0)
    dta$Y_A0[row_index] <- ifelse(dta$S[row_index] == 0, baseline_Y_S0, 
                                       ifelse(dta$S[row_index] == 1, baseline_Y_S1, baseline_Y_S2))
    dta$eps[row_index] <- rnorm(sample_size, mean = 0, sd = 1)
    dta$Y[row_index] <- (dta$A1[row_index] * (beta_1 + beta_2 * dta$S[row_index])  +
                           dta$A2[row_index] * (beta_3 + beta_4 * dta$S[row_index]) +
                           dta$Y_A0[row_index]+
                           dta$time[row_index]+
                           dta$eps[row_index])
  }
  
  return(dta)
}


#---------------------- Try on single data set-------------------------------
set.seed(2024)
control_vars <- c("S", "time")
moderator_vars <- c("S")
avail_varname = "I"

dta <- dgm_1(sample_size = 500, total_T = 20)

fit  <- wcls_categorical_treatment(
  dta = dta,
  id_varname = "userid",
  decision_time_varname = "time",
  treatment_varname = "A",
  outcome_varname = "Y",
  control_varname = c("S", "time"),
  moderator_varname = c("S"),
  rand_prob_varname = "prob_A",
  estimator_initial_value = NULL,
  trt_level = 2,
  pmatrix_tilde = NULL,
  avail_varname = avail_varname
)

fit

# case where pmatrix_tilde different for some time point 
# in this scenario we assume probability of no treatment is greater for the first 2 decision point (DP)
probs_list <- list(
  c(0.60, 0.20, 0.20),  # DP 1
  c(0.40, 0.3, 0.3),   # DP 2
  c(0.33, 0.33, 0.33)
)

pmatrix_tilde <- matrix(rep(NA, sample_size * total_T * 3), ncol = 3)

for (t in 1: total_T) {
  row_index <- seq(from = t, by = total_T, length = sample_size)
  
  if (t <=3) {
    pmatrix_tilde[row_index, ] <- matrix(rep(probs_list[[t]], each = sample_size), nrow = sample_size)
  }else {
    pmatrix_tilde[row_index, ] <- matrix(rep(probs_list[[3]], each = sample_size), nrow = sample_size)
  }
}
pmatrix_tilde

fit  <- wcls_categorical_treatment(
  dta = dta,
  id_varname = "userid",
  decision_time_varname = "time",
  treatment_varname = "A",
  outcome_varname = "Y",
  control_varname = c("S", "time"),
  moderator_varname = c("S"),
  rand_prob_varname = "prob_A",
  estimator_initial_value = NULL,
  trt_level = 2,
  pmatrix_tilde = pmatrix_tilde,
  avail_varname = avail_varname
)
