# sim 4.5a : misspecification of gt, constant true eo, quadratic working eo 
rm(list = ls())
source("functions/utillity.R")
source("functions/ss_calc.R")
source("functions/dgm.R")

# 1. Simulation design setup ----------------------------------------------

# SD: simulation design
SD <- expand.grid(gamma = 0.05,  # type I error
                  b = 0.2, # type II error
                  # m = c(30, 60),
                  m = c(30),
                  pt_shape = c("constant equal pk"),
                  # pt_shape = c("constant 0.5"),
                  ft_t_shape = c("constant"),
                  ft_t_theta = c(0),
                  ft_t_theta2 = 0,
                  gt_t_shape = c("constant"),
                  gt_t_theta = 0, 
                  # for specifying gt when gt_t_shape is linear or quadratic
                  # ATE_t = c(1.2, 1.3, 1.4, 1.5),
                  ATE1_t = c(1.2, 1.3),
                  ATE1_t_empirical = NA,
                  ATE2_t = c(1.4),
                  ATE2_t_empirical = NA,
                  ASPNC_t = c(0.2, 0.4, 0.8, 1, 2),
                  ASPNC_t_empirical = NA,
                  #AvgTau_t = c(0.5, 0.8, 1),
                  AvgTau_t = c(1),
                  tau_t_shape = c("constant"),
                  beta1_t0 = NA,
                  beta1_t1 = NA,
                  beta1_t2 = NA,
                  beta2_t0 = NA,
                  beta2_t1 = NA,
                  beta2_t2 = NA,
                  alpha_t0 = NA,
                  alpha_t1 = NA,
                  alpha_t2 = NA,
                  p = NA,
                  q = NA,
                  ft_w_shape = NA,
                  ft_w_theta = NA,
                  ft_w_theta2 = NA,
                  gt_w_shape = c("quadratic_theta"),
                  gt_w_theta = seq(from = -1, to = 1, by = 0.05),
                  beta1_w0 = NA,
                  beta1_w1 = NA,
                  beta1_w2 = NA,
                  beta2_w0 = NA,
                  beta2_w1 = NA,
                  beta2_w2 = NA,
                  alpha_w0 = NA,
                  alpha_w1 = NA,
                  alpha_w2 = NA,
                  ATE1_w = NA,
                  ATE1_w_empirical = NA,
                  ATE2_w = NA,
                  ATE2_w_empirical = NA,
                  ASPNC_w = NA,
                  ASPNC_w_empirical = NA,
                  AvgTau_w = NA,
                  tau_w_shape = NA,
                  p_w = NA,
                  q_w = NA,
                  n = NA,
                  p_star = NA,
                  
                  stringsAsFactors = FALSE)

# remove simulation settings where pt*ft is not in the linear span of gt
SD <- SD[!(SD$ft_t_shape == "linear_theta" &
             SD$gt_t_shape == "constant"), ]
SD <- SD[!(SD$ft_t_shape == "quadratic_theta" &
             SD$gt_t_shape != "quadratic_theta"), ]

SD <- SD[!(SD$ft_t_shape == "constant" &
             SD$ft_t_theta != 0), ]
SD <- SD[!(SD$gt_t_shape == "constant" &
             SD$gt_t_theta != 0), ]

SD <- SD[!(SD$AvgTau_t == 1 &
             SD$tau_t_shape != "constant"), ]

## add stuff that are known by simulation design
# taut
SD$AvgTau_w <- SD$AvgTau_t
SD$tau_w_shape <- SD$tau_t_shape
# gt
#SD$gt_w_shape <- SD$gt_t_shape
#SD$gt_w_theta <- SD$gt_t_theta
SD$ASPNC_w <- SD$ASPNC_t 
# ft
SD$ft_w_shape <- SD$ft_t_shape
SD$ft_w_theta <- SD$ft_t_theta
SD$ft_w_theta2 <- SD$ft_t_theta2
SD$ATE1_w <- SD$ATE1_t
SD$ATE2_w <- SD$ATE2_t


##### fill the SD #####

options(warn = 2)
for (i in 1:nrow(SD)) {
  
  ## gamma, b, m
  gamma <- SD$gamma[i]
  b <- SD$b[i]
  m <- SD$m[i]
  
  ## pt
  pt <- construct_pt(m, SD$pt_shape[i])
  
  ## taut_t
  taut_t <- construct_taut_theta(m , SD$AvgTau_t[i], SD$tau_t_shape[i])
  
  ## gt_t, alpha_t
  gt_t <- construct_ftgt(m, SD$gt_t_shape[i])
  SD$q[i] <- ncol(gt_t)
  alpha_t <- solve_alpha(SD$gt_t_shape[i], m, SD$ASPNC_t[i], gt_t, taut_t, 
                         theta = SD$gt_t_theta[i])
  SD[i, c("alpha_t0", "alpha_t1", "alpha_t2")[1:SD$q[i]]] <- alpha_t
  SD$ASPNC_t_empirical[i] <- compute_ASPNC(gt_t, alpha_t, taut_t)
  stopifnot(all.equal(SD$ASPNC_t[i], SD$ASPNC_t_empirical[i]))
  
  ## ft_t, beta1_t, beta2_t
  ft_t <- construct_ftgt(m, SD$ft_t_shape[i])
  SD$p[i] <- ncol(ft_t)
  # trt 1 effect
  beta1_t <- solve_beta(SD$ft_t_shape[i], m, SD$ATE1_t[i], ft_t, taut_t, theta = SD$ft_t_theta[i])
  SD[i, c("beta1_t0", "beta1_t1", "beta1_t2")[1:SD$p[i]]] <- beta1_t
  SD$ATE1_t_empirical[i] <- compute_ATE(ft_t, beta1_t, taut_t)
  stopifnot(all.equal(SD$ATE1_t[i], SD$ATE1_t_empirical[i]))
  # trt 2 effect
  beta2_tstar <- solve_beta_star(SD$ft_t_shape[i], m, SD$ATE2_t[i], ft_t, taut_t, theta = SD$ft_t_theta2[i], beta1 = beta1_t)
  beta2_t <- beta1_t +beta2_tstar 
  SD[i, c("beta2_t0", "beta2_t1", "beta2_t2")[1:SD$p[i]]] <- beta2_t
  SD$ATE2_t_empirical[i] <- compute_ATE(ft_t, beta2_t, taut_t)
  stopifnot(all.equal(SD$ATE2_t[i], SD$ATE2_t_empirical[i]))
  
  ## taut_w
  taut_w <- construct_taut_theta(m, SD$AvgTau_w[i], SD$tau_w_shape[i])
  
  ## gt_w, alpha_w
  gt_w <- construct_ftgt(m, SD$gt_w_shape[i])
  SD$q_w[i] <- ncol(gt_w)
  alpha_w <- solve_alpha(SD$gt_w_shape[i], m, SD$ASPNC_w[i], gt_w, taut_w,
                         theta = SD$gt_w_theta[i])
  SD[i, c("alpha_w0", "alpha_w1", "alpha_w2")[1:SD$q_w[i]]] <- alpha_w
  SD$ASPNC_w_empirical[i] <- compute_ASPNC(gt_w, alpha_w, taut_w)
  stopifnot(all.equal(SD$ASPNC_w[i], SD$ASPNC_w_empirical[i]))
  
  ## ft_w
  ft_w <- construct_ftgt(m, SD$ft_w_shape[i])
  SD$p_w[i] <- ncol(ft_w)
  # trt 1 eff
  beta1_w <- solve_beta(SD$ft_w_shape[i], m, SD$ATE1_w[i], ft_w,
                        taut_w, theta = SD$ft_w_theta[i])
  SD[i, c("beta1_w0", "beta1_w1", "beta1_w2")[1:SD$p_w[i]]] <- beta1_w
  SD$ATE1_w_empirical[i] <- compute_ATE(ft_w, beta1_w, taut_w)
  stopifnot(all.equal(SD$ATE1_w[i], SD$ATE1_w_empirical[i]))
  # trt 2 eff
  beta2_wstar <- solve_beta_star(SD$ft_w_shape[i], m, SD$ATE2_w[i], ft_w,
                        taut_w, theta = SD$ft_w_theta2[i], beta1 = beta1_w)
  beta2_w <- beta1_w + beta2_wstar
  SD[i, c("beta2_w0", "beta2_w1", "beta2_w2")[1:SD$p_w[i]]] <- beta2_w
  SD$ATE2_w_empirical[i] <- compute_ATE(ft_w, beta2_w, taut_w)
  stopifnot(all.equal(SD$ATE2_w[i], SD$ATE2_w_empirical[i]))
  
  ## sample size n
  tryCatch({
    beta_w <- c(beta1_w, beta2_w)
    L <- construct_L(SD$ft_w_shape[i])
    SD$n[i] <- mrt_mult_trt_ss(ft_w, gt_w, beta_w, alpha_w, pt, mean_var = 1, gamma, b,
                               L, tau = taut_w)
    SD$p_star[i] <- dim(L)[1]
  },
  error = function(cond){
    message("iteration ", i, ": ", cond)
    return(NA)
  })
  
  if (!is.na(SD$n[i])) {
    ## generate a data set to make sure there are no probabilities outside [0,1]
    dta <- dgm(SD$n[i], total_T = m, ft = ft_t, beta_1 = beta1_t, beta_2 = beta2_t, gt = gt_t, alpha = alpha_t, tau = taut_t, pt = pt)
  }
}

summary(SD)
dim(SD)
saveRDS(SD, file = "simulation/3. WA-b violated/sim 4.5a/Simulation_design.RDS")
