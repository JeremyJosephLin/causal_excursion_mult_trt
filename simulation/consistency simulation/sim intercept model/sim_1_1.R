# This simulation is to check consistency of our estimator under control misspecification
# Zt is a moderator variable, here the f(H_t) is correctly specified,
rm(list = ls(all = TRUE))
library(doParallel)
library(doSNOW)
library(utils)
library(foreach)
library(doRNG)
library(rootSolve) # for solver function multiroot()
max_cores <- 16
CL <- makeCluster(min(detectCores() - 1, max_cores))
registerDoParallel(CL)


# --------------------- Build data frame -------------------------------------- 
dgm_1 <- function(sample_size, total_T) {
  # same DGM as dgm_binary above, but faster
  # variable A represent the treatment level  (0 , 1 , 2)
  # variable A1 and A2 is indicator function whether A is equal 1 or 2 respectively
  # variable I is a binary indicating availability of the participant at timepoint t 
  # prob_t indicates the probability of someone available at time point t, assume constant for now
  # Z_t is moderator variable, that takes three values with equal probability

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
  prob_I <- 0.6
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
                           dta$eps[row_index])
  }
  
  return(dta)
}


get_alpha_beta_from_multiroot_result <- function(root, p, q)
{
  if (p == 1) {
    beta_root <- root$root[q+1]
  } else {
    beta_root <- as.matrix(root$root[(q+1) : (q+p)])
  }
  if (q == 1) {
    alpha_root <- root$root[1]
  } else {
    alpha_root <- as.matrix(root$root[1:q])
  }
  return(list(alpha = alpha_root, beta = beta_root))
}

interact <- function(ind_matrix, ind_names, cov, cov_names){
  # This is a function to create dummy matrix of treatment and 
  ncol_interaction <- ncol(cov) * ncol(ind_matrix)
  df_names <- rep(NA, ncol_interaction)
  Interaction_matrix <- 
    matrix(rep(NA, ncol_interaction * nrow(cov)), ncol = ncol_interaction)
  
  count = 0
  
  for (i in 1:ncol(ind_matrix)) {
    for (j in 1:ncol(cov)) {
      count = count + 1 
      df_names[count] <- paste0(cov_names[j], " : ", ind_names[i])
      Interaction_matrix [,count] <- cov[,j] * ind_matrix[,i]
    }
  }
  colnames(Interaction_matrix) <- df_names
  return(Interaction_matrix)
}

find_change_location <- function(v){
  #' Find the locations at which the value changes from the previous one
  # Used when total observation per individual is different, copied from TQ github
  
  n <- length(v)
  if (n <= 1) {
    stop("The vector need to have length > 1.")
  }
  return(c(1, 1 + which(v[1:(n-1)] != v[2:n])))
}

#---------------------------- TQ Estimator -------------------------------------
# TQ estimator can be found in R_Code/ estimator_MEE.R
wcls_categorical_treatment <- function(
    dta,
    id_varname,
    decision_time_varname,
    treatment_varname,
    outcome_varname,
    control_varname,
    moderator_varname,
    rand_prob_varname,
    avail_varname = NULL,
    p_tilde = NULL,
    estimator_initial_value = NULL
)
{
  ### 1. preparation ###
  
  sample_size <- length(unique(dta[, id_varname]))
  total_person_decisionpoint <- nrow(dta)
  total_T <- total_person_decisionpoint / sample_size
  A <- dta[, treatment_varname]

  p_t <- dta[, rand_prob_varname]
  
  # assign p_tilde equal for each treatment if p_tilde is not specified
  if (is.null(p_tilde)) {
    nA <- length(unique(A)) # number of treatment
    p_tilde <- rep(1/nA, nA)
  } else {
    nA <- length(p_tilde)
  }
  
  p_t_tilde <- ifelse(A == 0, p_tilde[1], 
                      ifelse(A ==1, p_tilde[2], p_tilde[3]))
  
  # indicator matrix, each column i represent  I(At = i) 
  ind_matrix <- matrix(rep(NA, (length(A) * (nA - 1) )), nrow = length(A))
  # this is the W = I - p_tilde function in our estimator
  ind_center <- matrix(rep(NA, (length(A) * (nA - 1) )), nrow = length(A))
  # column names of indicator matrix
  ind_names <- paste0("trt = ", 1:(nA-1))
  for (trt in (1:(nA - 1))) {
    ind_matrix[,trt] <- ifelse(A == trt, 1, 0)
    ind_center[,trt] <- ind_matrix[,trt] - p_tilde[trt + 1]
  }
  
  Y <- dta[, outcome_varname]
  
  # St is the covariates of moderator 
  St <- as.matrix(cbind(1, dta[, moderator_varname]))
  St_names <- c("Intercept", moderator_varname)
  
  # Ht is the covariate of control 
  Ht <- as.matrix(cbind(1, dta[, control_varname[1]]))
  Ht_names <- c("Intercept", control_varname)
  
  
  
  # Xdm what we are interested, Zdm is the one we are not interested 
  Xdm <- interact(ind_matrix, ind_names, St, St_names)  # X (moderator) design matrix, intercept added
  Zdm <- Ht  # Z (control) design matrix, intercept added
  # Wdm is the right most vector of our estimating equation 
  Wdm <- interact(ind_center, ind_names, St, St_names)
  Wnames <- colnames(Wdm)
  if (is.null(avail_varname)) {
    avail <- rep(1, total_person_decisionpoint)
  } else {
    avail <- dta[, avail_varname]
  }
  
  weight <- p_t_tilde / p_t
  
  p <-  (nA - 1) * (length(moderator_varname) + 1) # dimension of beta, need to generalize it later
  q <- length(control_varname) + 1 # dimension of alpha
  
  # for now manually add the name for moderator variable  
  Znames <- c("Intercept", control_varname)
  colnames(Zdm) <- Znames
  
  ### 2. estimation ###
  
  estimating_equation <- function(theta) {
    alpha <- as.matrix(theta[1:q])
    beta <- as.matrix(theta[(q+1):(q+p)])
    
    Zdm_alpha <- Zdm %*% alpha
    Wdm_beta <-  Wdm %*% beta
    
    residual <- Y - Zdm_alpha - Wdm_beta
    
    ef <- rep(NA, length(theta)) # value of estimating function
    for (i in 1:q) {
      ef[i] <- sum( weight * residual * avail * Zdm[, i])
    }
    for (i in 1:p) {
      ef[q + i] <- sum( weight * residual * avail * Wdm[,i])
    }
    
    ef <- ef / sample_size
    return(ef)
  }
  
  if (is.null(estimator_initial_value)) {
    estimator_initial_value <- rep(0, length = p + q)
  }
  
  # browser()
  solution <- tryCatch(
    {
      multiroot(estimating_equation, estimator_initial_value)
    },
    error = function(cond) {
      message("\nCatched error in multiroot inside efficient_ee():")
      message(cond)
      return(list(root = rep(NaN, p + q), msg = cond,
                  f.root = rep(NaN, p + q)))
    })
  
  estimator <- get_alpha_beta_from_multiroot_result(solution, p, q)
  
  alpha_hat <- as.vector(estimator$alpha)
  names(alpha_hat) <- Znames   # give alpha variable names
  beta_hat <- as.vector(estimator$beta)
  names(beta_hat) <- Wnames
  
  ### 3. asymptotic variance ###
  
  ### 3.1 Compute M_n matrix (M_n is the empirical expectation of the derivative of the estimating function) ###
  
  Mn_summand <- array(NA, dim = c(total_person_decisionpoint, p+q, p+q))
  # Mn_summand is  D^{(t),T} \frac{\partial r^(t)}{\partial \theta^T}
  # see May2023week2/variance derivation
  
  D_term_collected <- matrix(NA, nrow = p+q, ncol = total_person_decisionpoint)
  partialr_partialtheta_collected <- matrix(NA, nrow = total_person_decisionpoint, ncol = p+q)
  
  for (it in 1:total_person_decisionpoint) {
    # this is to make R code consistent whether X_it, Z_it contains more entries or is just the intercept.
    if (p == 1) {
      Xbeta <- Xdm[it, ] * beta_hat
    } else {
      Xbeta <- as.numeric(Xdm[it, ] %*% beta_hat)
    }
    if (q == 1) {
      Zalpha <- Zdm[it, ] * alpha_hat
    } else {
      Zalpha <- as.numeric(Zdm[it, ] %*% alpha_hat)
    }
    
    pre_multiplier <-  weight[it]
    
    # D_term = D^{(t),T} (dim = (nA * p+q) * 1)
    D_term <- pre_multiplier * avail[it]  * c(Zdm[it, ], Wdm[it, ])
    D_term_collected[, it] <- D_term
    
    # partialr_partialtheta = \frac{\partial r^(t)}{\partial \theta^T}
    partialr_partialtheta <- - c(Zdm[it, ], Wdm[it, ]) 
    partialr_partialtheta_collected[it, ] <- partialr_partialtheta
    
    Mn_summand[it, , ] <-  D_term %o% partialr_partialtheta
  }
  Mn <- apply(Mn_summand, c(2,3), sum) / sample_size
  Mn_inv <- solve(Mn)
  
  ### 3.2 Compute \Sigma_n matrix (\Sigma_n is the empirical variance of the estimating function) ###
  Zdm_alpha <- Zdm %*% alpha_hat
  Wdm_beta <-  Wdm %*% beta_hat
  residual <-  as.matrix(Y - Zdm_alpha - Wdm_beta)
  
  Sigman_summand <- array(NA, dim = c(sample_size, p+q, p+q))
  # Sigman_summand is  \sum_{t=1}^T ( D^{(t),T} r^(t) )^{\otimes 2}
  
  person_first_index <- c(find_change_location(dta[, id_varname]), total_person_decisionpoint + 1)
  
  
  for (i in 1:sample_size) {
    D_term_i <- D_term_collected[, person_first_index[i] : (person_first_index[i+1] - 1)]
    r_term_i <- matrix(residual[person_first_index[i] : (person_first_index[i+1] - 1)], ncol = 1)
    
    Sigman_summand[i, , ] <- D_term_i %*% r_term_i %*% t(r_term_i) %*% t(D_term_i)
  }
  Sigman <- apply(Sigman_summand, c(2,3), sum) / sample_size
  
  varcov <- Mn_inv %*% Sigman %*% t(Mn_inv) / sample_size
  varcovnames <- c(Znames,Wnames)
  colnames(varcov) <- varcovnames
  rownames(varcov) <- varcovnames
  
  alpha_se <- sqrt(diag(varcov)[1:q]) 
  names(alpha_se) <- Znames
  
  beta_se <- sqrt(diag(varcov)[(q+1):(q+p)])
  names(beta_se) <- Wnames
  
  ### 4. small sample correction ###
  
  Sigman_tilde <- 0
  for (i in 1:sample_size) {
    D_term_i <- D_term_collected[, person_first_index[i]:(person_first_index[i + 1] - 1)]
    r_term_i <- matrix(residual[person_first_index[i] : (person_first_index[i+1] - 1)], ncol = 1)
    partialr_partialtheta_i <- partialr_partialtheta_collected[person_first_index[i]:(person_first_index[i + 1] - 1), ]
    H_ii <- partialr_partialtheta_i %*% Mn_inv %*% D_term_i / sample_size
    Ii_minus_Hii_inv <- solve(diag(nrow(H_ii)) - H_ii)
    
    Sigman_tilde <- Sigman_tilde + D_term_i %*% Ii_minus_Hii_inv %*% r_term_i %*% t(r_term_i) %*% t(Ii_minus_Hii_inv) %*% t(D_term_i)
  }
  Sigman_tilde <- Sigman_tilde / sample_size
  
  varcov_adjusted <- Mn_inv %*% Sigman_tilde %*% t(Mn_inv) / sample_size
  colnames(varcov_adjusted) <- varcovnames
  rownames(varcov_adjusted) <- varcovnames
  
  alpha_se_adjusted <- sqrt(diag(varcov_adjusted)[1:q])
  names(alpha_se_adjusted) <- Znames
  
  beta_se_adjusted <- sqrt(diag(varcov_adjusted)[(q + 1):(q + p)])
  names(beta_se_adjusted) <- Wnames
  
  ### 5. calculate confidence interval
  
  conf_int <- cbind(beta_hat - 1.96 * beta_se, beta_hat + 1.96 * beta_se)
  conf_int_adjusted_z <- cbind(beta_hat - 1.96 * beta_se_adjusted, beta_hat + 1.96 * beta_se_adjusted)
  c <- qt(1 - 0.05/2, df = sample_size - p - q)
  conf_int_adjusted_t <- cbind(beta_hat - c * beta_se_adjusted,
                             beta_hat + c * beta_se_adjusted)
  colnames(conf_int) <- colnames(conf_int_adjusted_t)<- colnames(conf_int_adjusted_z) <- c("2.5 %", "97.5 %")
  row.names(conf_int) <- rownames(conf_int_adjusted_t) <- rownames(conf_int_adjusted_z) <- Wnames
  
  
  return(list(beta_hat = beta_hat, alpha_hat = alpha_hat, 
              beta_se = beta_se, alpha_se = alpha_se, 
              beta_se_adjusted = beta_se_adjusted, alpha_se_adjusted = alpha_se_adjusted,
              varcov= varcov,
              varcov_adjusted = varcov_adjusted, 
              conf_int = conf_int, 
              conf_int_adjusted_z = conf_int_adjusted_z,
              conf_int_adjusted_t = conf_int_adjusted_t,
              dims = list(p = p, q = q), 
              p_tilde = p_tilde)
         )
  
}

#---------------------- Try on single data set-------------------------------
control_vars <- c("S")
moderator_vars <- c("S")
avail_varname = "I"

dta <- dgm_1(sample_size = 50, total_T = 20)

fit  <- wcls_categorical_treatment(
  dta = dta,
  id_varname = "userid",
  decision_time_varname = "time",
  treatment_varname = "A",
  outcome_varname = "Y",
  control_varname = control_vars,
  moderator_varname = moderator_vars,
  rand_prob_varname = "prob_A",
  estimator_initial_value = NULL,
  avail_varname = avail_varname
)

fit
# ------------------------- run simulations ------------------------------------ 

control_vars <- c("S")
moderator_vars <- c("S")
avail_varname = "I"

#sample_size <- 31
sample_sizes <- c(15, 20, 25, 30, 40, 50)
total_T <- 15
nsim <- 1000


# result 
result_list <- c()

# Counter 
count = 1 





# for each variation of sample sizes and total T do n_sim number of simulation 
for (i in 1:length(sample_sizes)) {

    print(Sys.time())
    print(paste0("Sample size ", sample_sizes[i], " total T ", total_T))
    
    # changing the number of seed for every simulation
    set.seed(count)
    
    ### important ###
    
    pb <- txtProgressBar(max = nsim, style = 3)
    progress <- function(n) setTxtProgressBar(pb, n)
    opts <- list(progress = progress)
    
    ### end important ###
    
    
    result <-
      foreach(
        isim = 1:nsim,
        .combine = "c",
        .errorhandling = 'remove',
        .packages = c('lme4', 'pracma', 'rootSolve')
      ) %dorng% {
        if (isim %% 10 == 0) {
          cat(paste("Starting iteration", isim, "\n"))
        }  
        
        
        # Generate data
        dta <- dgm_1(sample_size = sample_sizes[i], total_T = total_T)
        
        # Fit TQ estimator
        fit <- wcls_categorical_treatment(
          dta = dta,
          id_varname = "userid",
          decision_time_varname = "time",
          treatment_varname = "A",
          outcome_varname = "Y",
          control_varname = control_vars,
          moderator_varname = moderator_vars,
          rand_prob_varname = "prob_A",
          avail_varname = avail_varname,
          p_tilde = c(0.2, 0.5,  0.3),
          estimator_initial_value = NULL
        )
        

        # Store list of list(2 result)
        # we dont append , so output get updated every simulation
        #?? why make it list of list then
        output <- list(
          list(beta_hat = fit$beta_hat, 
               beta_se = fit$beta_se, 
               varcov = fit$varcov,
               beta_se_adjusted = fit$beta_se_adjusted, 
               varcov_adjusted = fit$varcov_adjusted,
               ci_unadj = fit$conf_int, 
               ci_adj_z = fit$conf_int_adjusted_z,
               ci_adj_t = fit$conf_int_adjusted_t,
               p_tilde = fit$p_tilde,
               data = dta
          ))
      } # end of foreach loop 
    
    
    ### important ###
    close(pb)
    ### end important ###
    
    # update result list
    result_list <- c(result_list, list(list(sample_size = sample_sizes[i], 
                                            total_T = total_T,
                                            result = result)))
    
    # update iteration counter
    count <- count + 1
} # end of loop

# save result list 
saveRDS(object = result_list, file = "sim.rds")

