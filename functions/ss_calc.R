## Sample size calculator ##

#' This function is to calculate the sample size needed to guarantee 1-b power 
#' For the derivation of the samplesize calculator see Goodnotes-Power Calculation
#'
#'
#' @param f_t Parametric form of Moderator variable (S_t in the paper, must be a subset of g_t)
#' @param g_t Parametric form of Control Variable (H_t in the paper). 
#' @param beta p by 1 vector Coefficient for ATE (ft^T* beta)
#' @param alpha q by 1 vecotro coefftient, for gt^T*alpha
#' @param p_t (k + 1) by 1 vector of probability treatment 
#' @param mean_var sigma_square bar in the ss derivation
#' @param gamma  desired type 1 error 
#' @param b desired type 2 error (desired power equal 1 - b)
#' @param L Contrast matrix for Hypothesis testing 
#' @param tau Availability probabability for each decision point (nT by 1 vector), 
#'            set to 1 by default (Ind always available)
#' @param exact Boolean variable, Set to False if rounding is desired
#'
#' @return Number of sample size (Individual) needed to guarantee 1-b power
#' @export
#'
#' @examples  mrt_mult_trt_ss(ft_w, gt_w, beta_w, alpha_w, pt,
#'            mean_var = 1, gamma, b,L, tau = taut_w)
#'            
mrt_mult_trt_ss <- function(f_t,
                            g_t, 
                            beta,
                            alpha,
                            p_t, 
                            mean_var, 
                            gamma, 
                            b, 
                            L,
                            tau = rep(1,dim(f_t)[1]), 
                            exact = FALSE
){
  nA <- dim(pt)[2]
  t <- dim(f_t)[1]
  
  # degree of freedom of L*beta
  p <- dim(L)[1]
  # number of variable in g_t
  q <- dim(g_t)[2]
  
  M <- compute_sigma(f_t, beta, p_t, L, nA, tau = tau)
  M_matrix <- M$M_matrix
  
  # Set up the function to solve using uniroot
  power_f <- function(n){
    
    right_hand_side <- pf(q=qf(p=(1-gamma), df1=p, df2=n-q-p), 
                          df1=p, 
                          df2=n-q-p, 
                          ncp=compute_ncp(n, beta,  M_matrix, mean_var, L, tau = tau))
    left_hand_side <- b
    return(right_hand_side - left_hand_side)
  }
  
  # find min n to achieve desired power
  sample_size <- uniroot(power_f, lower=p+q+1, upper=1000000)$root
  
  # round up if non-exact size is requested
  if(exact == FALSE){
    sample_size <- ceiling(sample_size)
  }
  
  return(sample_size)
}



#' Function to compute Variance Covariance Matrix for S.S calculator
#'
#' @param f_t Parametric form of Moderator variable (S_t in the paper, must be a subset of g_t)
#' @param beta p by 1 vector Coefficient for ATE (ft^T* beta)
#' @param p_t (k + 1) by 1 vector of probability treatment 
#' @param L Contrast matrix for Hypothesis testing 
#' @param nA number of Treatment + control (K+1)
#' @param tau Availability probabability for each decision point (nT by 1 vector), 
#'            set to 1 by default (Ind always available)
#'
#' @return Var-Covariance MAtrix
#' @export
#'
#' @examples compute_sigma(f_t, beta, p_t, L, nA, tau = tau)
compute_sigma <- function(f_t,
                          beta,            
                          p_t, 
                          L,
                          nA,
                          tau
) {
  p <- length(beta)
  ## The M and Sigma matrices (needed to compute lambda)
  M_matrix <- matrix(data=0, nrow=p, ncol=p) 
  
  t <- dim(f_t)[1]
  # For each decision point 
  # (T = total # of decision points is taken as the length of p_t, for now)
  i = 1
  for (i in 1:t){
    
    # breaking down steps to identify bug and improve robustness of code
    this_f_t <- as.numeric(f_t[i, ])
    this_p_t <- as.numeric(p_t[i,])
    tau_t <- as.numeric(tau[i])
    if(length(this_f_t) != p/(nA-1)){
      stop("Incorrect dimensions for f_t.")
    }
    this_pt <- this_p_t[2:(p+1)]
    # create p_matrix
    p_matrix <- matrix(data = 0, nrow = (dim(p_t)[2]-1), ncol = (dim(p_t)[2]-1))
    
    for (row in 1:(dim(p_t)[2]-1)) {
      for (col in 1:(dim(p_t)[2]-1)) {
        if (row == col) {
          # diagonal element of p_matrix
          p_matrix[row, col] = this_pt[row] * (1 - this_pt[row])
        }else{
          # off diagonal of p_matrix
          p_matrix[row, col] = -1 * this_pt[row] * this_pt[col]
        }
      }
    } 
    # calculate f_t %*% f_t^T
    this_f_t_f_t <- outer(this_f_t, this_f_t)
    
    # calculate the sigma matrix at time t 
    this_sigma <- kronecker(p_matrix, this_f_t_f_t) * tau_t # kroenecker product(see asynmptotic notes)
    # running sum 
    M_matrix <- M_matrix + this_sigma
  }
  
  
  return(list(M_matrix = M_matrix))
  
}


#' function to calculate the non centralize parameter for Chi-Square
#'
#' @param x proposed sample size (n)
#' @param beta p by 1 vector Coefficient for ATE (ft^T* beta)
#' @param M_matrix Variance Covariance Matrix
#' @param mean_var sigma_square bar in the ss derivation
#' @param L Contrast matrix for Hypothesis testing 
#' @param tau Availability probabability for each decision point (nT by 1 vector), 
#'            set to 1 by default (Ind always available)
#'
#' @return non centrality parameter
#' @export
#'
#' @examples
compute_ncp <- function(x, beta, M_matrix, mean_var, L, tau){
  if(det(M_matrix) == 0){
    stop("sigma_matrix must be nonsingular")
  }
  
  return(as.numeric(x /(mean_var)  * t(L %*% beta) %*%
                      solve(L%*% solve(M_matrix) %*% t(L)) %*% 
                      L %*% beta))
}

