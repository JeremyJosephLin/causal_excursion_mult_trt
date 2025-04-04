rm(list = ls(all = TRUE))
# call functions
source("functions/utillity.R")
source("functions/dgm.R")
source("functions/wcls_categorical_treatment.R")

# grab factorial design
SD <- readRDS("simulation/2. WA-a violated/sim 2.1a/Simulation_design.RDS")

# get command parameters
args <- commandArgs(trailingOnly = TRUE)
isetting <- as.integer(args[1]) # settings 
nsim <- as.integer(args[2])
nsetting <- as.integer(args[3])

# # Debugging 
# isetting = 1
# nsim= 1
# nsetting = 1

setting_start <- (isetting - 1)* nsetting + 1
setting_end <- isetting * nsetting

print(paste0("setting value ", isetting, " nsim ", nsim, " nsetting ", nsetting, " setting start ",setting_start, " setting end ",setting_end))

set.seed(isetting)
# 2. Simulation(parallel) -------------------------------------------------------
options(warn = 0)


library(rootSolve) # for solver function multiroot()
result_list_collected <- list()

for (i in setting_start:setting_end) {
  result <- list()
  
  start_time <- Sys.time()
  current_time <- Sys.time()
  print(Sys.time())
  print(paste0(round(
    difftime(current_time, start_time, units = "hours"), 2
  ),
  " hours has lapsed."))
  cat("i =", i, "/", nrow(SD), "\n")
  
  #print(SD[i, ])
  
  ## gamma, b, m
  gamma <- SD$gamma[i]
  b <- SD$b[i]
  m <- SD$m[i]
  
  ## pt
  pt <- construct_pt(m, SD$pt_shape[i])
  
  ## taut_t
  taut_t <-
    construct_taut_theta(m , SD$AvgTau_t[i], SD$tau_t_shape[i])
  
  ## gt_t, alpha_t
  gt_t <- construct_ftgt(m, SD$gt_t_shape[i])
  alpha_t <-
    as.numeric(SD[i, c("alpha_t0", "alpha_t1", "alpha_t2")[1:SD$q[i]]])
  
  ## ft_t, beta1_t, beta2_t
  ft_t <- construct_ftgt(m, SD$ft_t_shape[i])
  beta1_t <-
    as.numeric(SD[i, c("beta1_t0", "beta1_t1", "beta1_t2")[1:SD$p[i]]])
  beta2_t <-
    as.numeric(SD[i, c("beta2_t0", "beta2_t1", "beta2_t2")[1:SD$p[i]]])
  beta_t <- c(beta1_t, beta2_t)
  L <- construct_L(SD$ft_w_shape[i])
  
  ## taut_w
  # (not used in the following)
  
  ## gt_w, alpha_w
  gt_w <- construct_ftgt(m, SD$gt_w_shape[i])
  alpha_w <-
    as.numeric(SD[i, c("alpha_w0", "alpha_w1", "alpha_w2")[1:SD$q_w[i]]])
  
  ## ft_w
  ft_w <- construct_ftgt(m, SD$ft_w_shape[i])
  beta1_w <-
    as.numeric(SD[i, c("beta1_w0", "beta1_w1", "beta1_w2")[1:SD$p_w[i]]])
  beta2_w <-
    as.numeric(SD[i, c("beta2_w0", "beta2_w1", "beta2_w2")[1:SD$p_w[i]]])
  beta_w <- c(beta1_w, beta2_w)
  
  ## n, p_w, q_w
  n <- SD$n[i]
  p_w <- SD$p_w[i]
  q_w <- SD$q_w[i]
  p_star <- SD$p_star[i] # This is the degree of freedom of T test
  p <- 2 * p_w # this is the dimension of beta
  q <- q_w
  
  ##### simulation begins #####
  for (i_sim in 1:nsim) {
    # Generate data
    dta <-
      dgm(
        SD$n[i],
        total_T = m,
        ft = ft_t,
        beta_1 = beta1_t,
        beta_2 = beta2_t,
        gt = gt_t,
        alpha = alpha_t,
        tau = taut_t,
        pt = pt
      )
    
    dta <- add_analysis_vars(dta, n, ft_w, gt_w)
    control_varname <- get_control_varname(gt_w)
    moderator_varname <- get_moderator_varname(ft_w)
    avail_varname = "I"
    
    # Fit TQ estimator
    fit_TQ <- wcls_categorical_treatment(
      dta = dta,
      id_varname = "userid",
      decision_time_varname = "time",
      treatment_varname = "A",
      outcome_varname = "Y",
      control_varname = control_varname,
      moderator_varname = moderator_varname,
      rand_prob_varname = "prob_A",
      avail_varname = avail_varname,
      trt_level = 2,
      estimator_initial_value = NULL
    )
    
    
    beta_hat <- fit_TQ$beta_hat
    var_beta <- fit_TQ$varcov[(q + 1):(p + q), (q + 1):(p + q)]
    var_beta_adj <-
      fit_TQ$varcov_adjusted[(q + 1):(p + q), (q + 1):(p + q)]
    
    # without small sample correction
    test_stat_unadj <-
      t(L %*% beta_hat) %*% solve(L %*% var_beta %*% t(L)) %*% (L %*% beta_hat)
    
    # with small sample correction
    test_stat_adj <-
      t(L %*% beta_hat) %*% solve(L %*% var_beta_adj %*% t(L)) %*% (L %*% beta_hat)
    
    output <- list(
      list(
        fit = fit_TQ,
        test_stat_unadj = test_stat_unadj,
        test_stat_adj = test_stat_adj
      )
    )
    result <- c(result,output)
  }
  result_list_collected <- c(result_list_collected, list(list(
    setting = i,
    result = result
  )))
}

outName <- paste("results_sim_2_1_setting_",i,"_",nSeed,".RDS",sep = "")
saveRDS(result,file=outName)
