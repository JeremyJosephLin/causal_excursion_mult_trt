construct_ftgt <- function(m, shape) {
  # linear_theta is for ft and beta, see goodnotes "dgm_noDE"
  match.arg(shape, c("constant", "linear", "quadratic", 
                     "linear_theta", "quadratic_theta"))
  if (shape == "constant") {
    return(matrix(1, nrow = m, ncol = 1))
  } else if (shape %in% c("linear", "linear_theta")) {
    return(cbind(1, 1:m))
  } else if (shape %in% c("quadratic", "quadratic_theta")) {
    return(cbind(1, 1:m, (1:m)^2))
  }
}

add_analysis_vars <- function(dta, n, ft, gt) {
  stopifnot(nrow(ft) * n == nrow(dta))
  stopifnot(nrow(gt) * n == nrow(dta))
  
  stopifnot(all(ft[, 1] == 1))
  stopifnot(all(gt[, 1] == 1))
  
  if (ncol(ft) >= 2) {
    for (ivar in 1:(ncol(ft) - 1)) {
      dta[, paste0("mod", ivar)] <- rep(ft[, ivar + 1], n)
    }
  }
  if (ncol(gt) >= 2) {
    for (ivar in 1:(ncol(gt) - 1)) {
      dta[, paste0("ctr", ivar)] <- rep(gt[, ivar + 1], n)
    }
  }
  
  return(dta)
}

get_control_varname <- function(gt) {
  if (ncol(gt) == 1) {
    return(NULL)
  } else {
    return(paste0("ctr", 1:(ncol(gt) - 1)))
  }
}

get_moderator_varname <- function(ft) {
  if (ncol(ft) == 1) {
    return(NULL)
  } else {
    return(paste0("mod", 1:(ncol(ft) - 1)))
  }
}


construct_taut_theta <- function(m, AvgTau, shape, theta = 0) {
  match.arg(shape, c("constant", "linear", "sine"))
  stopifnot(AvgTau + theta <= 1 & AvgTau + theta >= 0 &
              AvgTau - theta <= 1 & AvgTau - theta >= 0)
  if (shape == "constant") {
    taut <- rep(AvgTau, m)
  } else if (shape == "linear") {
    taut <- seq(from = AvgTau + theta, to = AvgTau - theta, length.out = m)
  } else if (shape == "sine") {
    taut <- AvgTau + theta * sin(1:m)
  }
  stopifnot(all(taut >= 0 & taut <= 1))
  return(taut)
}


construct_pt <- function(m, shape) {
  match.arg(shape, c("constant equal pk", "constant p1 high", "constant p2 high"))
  if (shape == "constant equal pk") {
    return(matrix(rep(c(0.33, 0.33, 0.33),times = m ), ncol = 3, byrow = TRUE))
  } else if (shape == "constant p1 high") {
    return(matrix(rep(c(0.3, 0.4, 0.3),times = m ), ncol = 3, byrow = TRUE))
  } else if (shape == "constant p2 high") {
    return(matrix(rep(c(0.3, 0.3, 0.4),times = m ), ncol = 3, byrow = TRUE))
  }
}


construct_L <- function(shape){
  match.arg(shape, c("constant", "linear_theta", "quadratic_theta"))
  if (shape == "constant") {
    return(matrix(c(1,-1), nrow = 1))
  } else if (shape == "linear_theta") {
    return(matrix(cbind(diag(1, nrow = 2), - diag(1, nrow = 2)), nrow = 2))
  } else if (shape == "quadratic_theta") {
    return(NA)
  }
}

# function to create linear variance 
construct_var_theta_time <- function(m, mean_var, shape, theta = 0) {
  match.arg(shape, c("constant", "linear"))
  stopifnot(mean_var + theta >= 0 & mean_var - theta >= 0)
  
  if (shape == "constant") {
    exp_var <- rep(mean_var, m)
  } else if (shape == "linear") {
    exp_var <- seq(from = 1 + theta/mean_var, to = 1 - theta/mean_var, length.out = m)
  }
  stopifnot(all(exp_var >= 0))
  return(exp_var)
}

# function to create variance depends on trt
construct_var_theta_trt <- function(m, pt, theta){
  n_trt <- dim(pt)[2]
  var_trt <- h_t <- matrix(rep(NA, times = (m * n_trt) ), ncol = n_trt)
  for (i in 1:m) {
    # store the coefficient of h(At)
    # column are by treatment so at = 0, 1, 2
    h_t[i, 2] = theta
    h_t[i, 3] =  -1 * (pt[i,2] - 1)/(pt[i,3] - 1) * theta
    h_t[i, 1] = 1 - h_t[i,2] - h_t[i,3]
    # store the expected variance to make sure they are all positive
    var_trt[i,1] = h_t[i,1]
    var_trt[i,2] = h_t[i,1] + h_t[i,2]
    var_trt[i,3] = h_t[i,1] + h_t[i,3]
  }
  stopifnot(all(var_trt >= 0))
  return(h_t)
}


##### solve beta #####
solve_beta <- function(shape, m, ATE, ft, taut, ...) {
  match.arg(shape, c("constant", "linear_theta", "quadratic_theta"))
  if (shape == "constant") {
    return(solve_beta_constant(m, ATE, ft, taut))
  } else if (shape == "linear_theta") {
    return(solve_beta_linear_theta(m, ATE, ft, taut, ...))
  } else if (shape == "quadratic_theta") {
    return(solve_beta_quadratic_theta(m, ATE, ft, taut, ...))
  }
}

solve_beta_star <- function(shape, m, ATE, ft, taut, ...) {
  match.arg(shape, c("constant", "linear_theta", "quadratic_theta"))
  if (shape == "constant") {
    return(solve_beta_star_constant_theta(m, ATE, ft, taut,...))
  } else if (shape == "linear_theta") {
    return(solve_beta_star_linear_theta(m, ATE, ft, taut, ...))
  } else if (shape == "quadratic_theta") {
    return(solve_beta_quadratic_theta(m, ATE, ft, taut, ...))
  }
}

compute_ATE <- function(ft, beta, taut) {
  num <- sum((ft %*% beta) * taut)
  denom <- sum(taut)
  return(num/denom)
}

solve_beta_constant <- function(m, ATE, ft, taut) {
  return(ATE)
}


solve_beta_linear_theta <- function(m, ATE, ft, taut, theta) {
  # see goodnotes "dgm_noDE"
  stopifnot(ncol(ft) == 2)
  stopifnot(theta >= -1 & theta <= 1)
  
  if (theta == 1) {
    from_beta1_to_ATE <- function(beta1) {
      beta0 <- -m * beta1
      beta <- c(beta0, beta1)
      return(compute_ATE(ft, beta, taut) - ATE)
    }
    beta1 <- uniroot(from_beta1_to_ATE, interval = c(-0.5, 0.5), extendInt = "yes",
                     tol = 1e-12)$root
    beta0 <- -m * beta1
    return(c(beta0, beta1))
  } else if (theta == 0) {
    beta1 <- 0
    from_beta0_to_ATE <- function(beta0) {
      beta <- c(beta0, 0)
      return(compute_ATE(ft, beta, taut) - ATE)
    }
    beta0 <- uniroot(from_beta0_to_ATE, interval = c(-0.5, 0.5), extendInt = "yes",
                     tol = 1e-12)$root
    return(c(beta0, beta1))
  } else {
    ratio <- (1+theta) / (1-theta)
    from_beta1_to_ATE <- function(beta1) {
      beta0 <- (ratio * m - 1) / (1 - ratio) * beta1
      beta <- c(beta0, beta1)
      return(compute_ATE(ft, beta, taut) - ATE)
    }
    beta1 <- uniroot(from_beta1_to_ATE, interval = c(-0.5, 0.5), extendInt = "yes",
                     tol = 1e-12)$root
    beta0 <- (ratio * m - 1) / (1 - ratio) * beta1
    return(c(beta0, beta1))
  }
}

solve_beta_star_constant_theta <- function(m, ATE, ft, taut, theta, beta1){
  beta1_star <- 0
  from_beta0_to_ATE <- function(beta0_star) {
    beta <- c(beta0_star + beta1[[1]])
    return(compute_ATE(ft, beta, taut) - ATE)
  }
  beta0_star <- uniroot(from_beta0_to_ATE, interval = c(-0.5, 0.5), extendInt = "yes",
                        tol = 1e-12)$root
  return(beta0_star)
}

# solve beta2 with parameterization theta2 % for violation 2.4 see notes 
solve_beta_star_linear_theta <- function(m, ATE, ft, taut, theta, beta1) {
  # see goodnotes "dgm_noDE"
  stopifnot(ncol(ft) == 2)
  stopifnot(length(beta1) == 2)
  stopifnot(theta >= -1 & theta <= 1)
  
  if (theta == 1) {
    from_beta1_to_ATE <- function(beta1_star) {
      beta0_star <- -m * beta1_star
      beta <- c(beta0_star + beta1[[1]], beta1_star + beta1[[2]])
      return(compute_ATE(ft, beta, taut) - ATE)
    }
    beta1_star <- uniroot(from_beta1_to_ATE, interval = c(-0.5, 0.5), extendInt = "yes",
                          tol = 1e-12)$root
    beta0_star <- -m * beta1_star
    return(c(beta0_star, beta1_star))
  } else if (theta == 0) {
    beta1_star <- 0
    from_beta0_to_ATE <- function(beta0_star) {
      beta <- c(beta0_star + beta1[[1]], beta1[[2]])
      return(compute_ATE(ft, beta, taut) - ATE)
    }
    beta0_star <- uniroot(from_beta0_to_ATE, interval = c(-0.5, 0.5), extendInt = "yes",
                          tol = 1e-12)$root
    return(c(beta0_star, 0))
  } else {
    ratio <- (1+theta) / (1-theta)
    from_beta1_to_ATE <- function(beta1_star) {
      beta0_star <- (ratio * m - 1) / (1 - ratio) * beta1_star
      beta <- c(beta0_star + beta1[[1]], beta1_star + beta1[[2]])
      return(compute_ATE(ft, beta, taut) - ATE)
    }
    beta1_star <- uniroot(from_beta1_to_ATE, interval = c(-0.5, 0.5), extendInt = "yes",
                          tol = 1e-12)$root
    beta0_star <- (ratio * m - 1) / (1 - ratio) * beta1_star
    return(c(beta0_star, beta1_star))
  }
}

##### solve alpha #####
solve_alpha <- function(shape, m, ASPNC, gt, taut, ...) {
  match.arg(shape, c("constant", "linear_theta", "quadratic_theta"))
  if (shape == "constant") {
    return(solve_alpha_constant(m, ASPNC, gt, taut))
  } else if (shape == "linear_theta") {
    return(solve_alpha_linear_theta(m, ASPNC, gt, taut, ...))
  } else if (shape == "quadratic_theta") {
    return(solve_alpha_quadratic_theta(m, ASPNC, gt, taut, ...))
  }
}

compute_ASPNC <- function(gt, alpha, taut) {
  num <- sum(gt %*% alpha * taut)
  denom <- sum(taut)
  return(num / denom)
}

solve_alpha_constant <- function(total_T, ASPNC, gt, taut) {
  return(ASPNC)
}

solve_alpha_linear_theta <- function(m, ASPNC, gt, taut, theta) {
  # see goodnotes "dgm_noDE" solve_beta_linear_theta
  # the math for alpha and beta are exactly the same
  stopifnot(ncol(gt) == 2)
  stopifnot(theta >= -1 & theta <= 1)
  
  if (theta == 1) {
    from_alpha1_to_ASPNC <- function(alpha1) {
      alpha0 <- -m * alpha1
      alpha <- c(alpha0, alpha1)
      return(compute_ASPNC(gt, alpha, taut) - ASPNC)
    }
    alpha1 <- uniroot(from_alpha1_to_ASPNC, interval = c(-0.5, 0.5), extendInt = "yes",
                      tol = 1e-12)$root
    alpha0 <- -m * alpha1
    return(c(alpha0, alpha1))
  } else if (theta == 0) {
    alpha1 <- 0
    from_alpha0_to_ASPNC <- function(alpha0) {
      alpha <- c(alpha0, 0)
      return(compute_ASPNC(gt, alpha, taut) - ASPNC)
    }
    alpha0 <- uniroot(from_alpha0_to_ASPNC, interval = c(-0.5, 0.5), extendInt = "yes",
                      tol = 1e-12)$root
    return(c(alpha0, alpha1))
  } else {
    ratio <- (1+theta) / (1-theta)
    from_alpha1_to_ASPNC <- function(alpha1) {
      alpha0 <- (ratio * m - 1) / (1 - ratio) * alpha1
      alpha <- c(alpha0, alpha1)
      return(compute_ASPNC(gt, alpha, taut) - ASPNC)
    }
    alpha1 <- uniroot(from_alpha1_to_ASPNC, interval = c(-0.5, 0.5), extendInt = "yes",
                      tol = 1e-12)$root
    alpha0 <- (ratio * m - 1) / (1 - ratio) * alpha1
    return(c(alpha0, alpha1))
  }
}

solve_alpha_quadratic_theta <- function(m, ASPNC, gt, taut, theta) {
  # see goodnotes "dgm_noDE" solve_beta_quadratic_theta
  # the math for alpha and beta are exactly the same
  stopifnot(ncol(gt) == 3)
  stopifnot(theta >= -1 & theta <= 1)
  
  if (theta == 1) {
    from_alpha2_to_ASPNC <- function(alpha2) {
      alpha0 <- m * alpha2
      alpha1 <- -(m+1) * alpha2
      alpha <- c(alpha0, alpha1, alpha2)
      return(compute_ASPNC(gt, alpha, taut) - ASPNC)
    }
    alpha2 <- uniroot(from_alpha2_to_ASPNC, interval = c(-0.5, 0.5), extendInt = "yes",
                      tol = 1e-12)$root
    alpha0 <- m * alpha2
    alpha1 <- -(m+1) * alpha2
    return(c(alpha0, alpha1, alpha2))
  } else if (theta == 0) {
    alpha1 <- 0
    alpha2 <- 0
    from_alpha0_to_ASPNC <- function(alpha0) {
      alpha <- c(alpha0, 0, 0)
      return(compute_ASPNC(gt, alpha, taut) - ASPNC)
    }
    alpha0 <- uniroot(from_alpha0_to_ASPNC, interval = c(-0.5, 0.5), extendInt = "yes",
                      tol = 1e-12)$root
    return(c(alpha0, alpha1, alpha2))
  } else {
    ratio <- (1+theta) / (1-theta)
    from_alpha2_to_ASPNC <- function(alpha2) {
      alpha0 <- ((m+1)^2/4 - ratio*(m+1)+ratio) / (1-ratio) * alpha2
      alpha1 <- -(m+1) * alpha2
      alpha <- c(alpha0, alpha1, alpha2)
      return(compute_ASPNC(gt, alpha, taut) - ASPNC)
    }
    alpha2 <- uniroot(from_alpha2_to_ASPNC, interval = c(-1, 1), extendInt = "yes",
                      tol = 1e-12)$root
    alpha0 <- ((m+1)^2/4 - ratio*(m+1)+ratio) / (1-ratio) * alpha2
    alpha1 <- -(m+1) * alpha2
    return(c(alpha0, alpha1, alpha2))
  }
}
