library(rootSolve)
# update 11/14/2025 pmatrix_tilde dimension is updated to (sample_size * total_T * nA) x nA to acommodate ptilde for each individual 

#------------------ wcls categorical treatment ---------------------------------
# TQ estimator can be found in R_Code/ estimator_MEE.R
#' Estimates the marginal excursion effect for cntinuous outcome MRT 
#' with categorical treatment
#'
#' This estimator returns the estimates for the marginal excursion effect estimator
#' and provides the estimated variance and
#' standard error for the estimators, with small sample correction for
#' an MRT with binary outcome using the "Hat" matrix in the variance estimate
#' and t-distribution or F-distribution critical value with corrected degrees of freedom.
#' 
#' @param dta a data set in a long format
#' @param id_varname the variable name that specifies the subject IDs
#' @param decision_time_varname the variable name that specifies the decision points
#' @param treatment_varname the variable name that specifies the assigned treatments for subjects
#' @param outcome_varname the variable name that specifies the outcomes for subjects
#' @param control_varname a vector of variable names used to reduce noise, this could be NULL
#' @param moderator_varname a vector of variable names of the effect modifiers, this could be NULL
#' @param rand_prob_varname the variable name that specifies the treatment randomizing probability in the data set
#' @param avail_varname the variable name that specifies the availability of the subjects,
#'                      default to be always available at any decision points using NULL
#' @param trt_level variable that specifies the number of categorical treatment 
#' @param pmatrix_tilde a matrix that pre-specifies the probability at each decision point for each individual
#'                      default to be always equal for each treatments using NULL
#' @param estimator_initial_value a numeric vector of the initial value for the estimator,
#'                                its length should be the sum of the length of control and moderator variables plus 2
#'                                default to be all 0's using NULL
#'
#' @returnReturns the estimated beta with its intercept and alpha with its intercept,
#'         the standard error of the estimated beta with its intercept and alpha with its intercept,
#'         the adjusted standard error of the estimated beta with its intercept and alpha with its intercept for small sample,
#'         the estimated variance-covariance matrix for the estimated beta with its intercept and alpha with its intercept,
#'         the estimated variance-covariance matrix for the estimated beta with its intercept and alpha with its intercept for small sample,
#'         the 95 percent confidence interval for beta_hat, and the adjusted 95 percent confidence interval for beta_hat,
#'         the dimension of the moderated variables, and the dimension of the control variables,
#'         the value of the probability weight at each decision point
#' @export 
#'
#' @exampleswcls_categorical_treatment(dta = dta,
#'                                          id_varname = "userid",
#'                                          decision_time_varname = "time",
#'                                          treatment_varname = "A",
#'                                          outcome_varname = "Y",
#'                                          control_varname = c("S", "time"),
#'                                          moderator_varname = c("S"),
#'                                          rand_prob_varname = "prob_A",
#'                                          estimator_initial_value = NULL,
#'                                          trt_level = 2,
#'                                          pmatrix_tilde = NULL,
#'                                          avail_varname = avail_varname)
wcls_categorical_treatment <- function(
    dta,
    id_varname,
    decision_time_varname,
    treatment_varname,
    outcome_varname,
    control_varname,   #intercept will be automatically added
    moderator_varname, #intercept will be automatically added
    rand_prob_varname,
    avail_varname = NULL,
    trt_level, # exclude control group
    pmatrix_tilde = NULL,
    estimator_initial_value = NULL
)
{
  ### 1. preparation ###
  
  sample_size <- length(unique(dta[, id_varname]))
  total_person_decisionpoint <- nrow(dta)
  total_T <- total_person_decisionpoint / sample_size
  A <- dta[, treatment_varname]
  
  p_t <- dta[, rand_prob_varname]
  nA <- trt_level + 1
  # assign pmatrix_tilde equal for each treatment if pmatrix_tilde is not specified
  if (is.null(pmatrix_tilde)) {
    pmatrix_tilde <- matrix(rep(1/nA, nA* total_T * sample_size), ncol = nA)
  } 
  
  # indicator matrix, each column i represent  I(At = i) 
  ind_matrix <- matrix(rep(NA, (length(A) * (nA - 1) )), nrow = length(A))
  # this is the W = I - p_tilde function in our estimator
  ind_center <- matrix(rep(NA, (length(A) * (nA - 1) )), nrow = length(A))
  # column names of indicator matrix
  ind_names <- paste0("trt = ", 1:(nA-1))
  #initialize p_t_tilde
  p_t_tilde = rep(NA, nrow(dta))
  
  for (i in 1:total_T) {
    # extract row index for all individual at time i 
    row_index <- which(dta$time == i)
    # extract p_t for at time point i
    temp <- dta[row_index, treatment_varname]
    p_t_tilde[row_index] = match_p_tilde(temp, pmatrix_tilde[row_index,])
    ind <- build_W_mat(temp,pmatrix_tilde[i,] )
    ind_matrix[which(dta$time == i),] = ind$ind_matrix
    ind_center[which(dta$time == i),] = ind$ind_center
  }
  
  
  Y <- dta[, outcome_varname]
  
  # St is the covariates of moderator 
  St <- as.matrix(cbind(1, dta[, moderator_varname]))
  St_names <- c("Intercept", moderator_varname)
  
  # Ht is the covariate of control 
  Ht <- as.matrix(cbind(1, dta[, control_varname]))
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
              dims = list(p = p, q = q, n = sample_size, nA = nA), 
              pmatrix_tilde = pmatrix_tilde)
  )
  
}



# utility function to run 

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

match_p_tilde <- function(u, pt_tilde){
  # this function match treatment assignment with p_tilde at a given time point
  # u is trt assignment for all individual at a given time point
  # p_tilde is the assignment at given time point 
  temp_p_tilde = rep(NA, length(u))
  for (row in 1:length(u)) {
    # assign p_tilde for each individual
    # grab p_tilde for each individual(row)
    p_tilde <- pt_tilde[row, ]
    for (j in 1:length(p_tilde)) {
      # match with corresponding p_tilde
       if (u[row] == (j-1)) {
        temp_p_tilde[row] = p_tilde[j]
      }
    }
  }
  return(temp_p_tilde)
}

build_W_mat <- function(u, p_tilde){
  # this function match treatment assignment with p_tilde at a given time point
  # u is trt assignment for all individual at a given time point
  # p_tilde is the assignment at given time point 
  n_ind = length(u)
  n_trt = length(p_tilde) - 1
  temp_ind_center <- temp_ind_matrix <- matrix(rep(NA, n_ind*n_trt), ncol = n_trt)
  
  for (ind in 1:n_ind) {
    # assign p_tilde for each individual
    for (trt in 1:n_trt) {
      temp_ind_matrix[,trt] <- ifelse(u == trt, 1, 0)
    }
    p_tilde_mat = matrix(rep(p_tilde[-1], n_ind), ncol = n_trt, byrow = TRUE)# exclude ctrl
    temp_ind_center = temp_ind_matrix - p_tilde_mat
  }
  return(list(ind_matrix = temp_ind_matrix, ind_center = temp_ind_center))
}

