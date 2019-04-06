# Helper function for RCSEs

vcovCluster <- function(model, cluster){
  if(nrow(model.matrix(model))!=length(cluster)){
    stop("check your data: cluster variable has different N than model")
  }
  M <- length(unique(cluster))
  N <- length(cluster)           
  K <- model$rank   
  if(M<50){
    warning("Fewer than 50 clusters, variances may be unreliable (could try block bootstrap instead).")
  }
  dfc <- (M/(M - 1)) * ((N - 1)/(N - K))
  uj  <- apply(estfun(model), 2, function(x) tapply(x, cluster, sum));
  rcse.cov <- dfc * sandwich(model, meat = crossprod(uj)/N)
  return(rcse.cov)
}



# Arguments to function ---------------------------------------------------

# N = population size
# n = sample size
# m = number of arms
# prob_assign = vector of units assigned to each arm (control first)
# mu_Y0 = expectation of untreated potential outcome
# ATEs = vector of ATEs for each treatment arm (length prob_each - 1)
# noise_scale = rescale the noise for each experimental group if desired
# coef_X = coefficient on covariate X
# location_scale_X = mean and SD of normally-distributed X covariate
# cov_adjustment = use covariate adjustment in estimation
# block_var_probs = distribution of multinomial blocking variable (probability per level)
# blocked_RA = use blocked random assignment
# block_prob_each = matrix of probabilities if different blocks have different probs of assignment to treatment
# n_clust_pop = number of clusters in population
# n_clust_samp = number of clusters in sample

m_arm_template <- function(N, 
                           n = N, 
                           m = 2,
                           prob_each = rep(1/m, m), 
                           mu_Y0 = 0, 
                           ATEs = 0, 
                           noise_scale = 1, 
                           coef_X = 0,
                           location_scale_X = c(0, 1),
                           cov_adjustment = F,
                           block_var_probs = c(.5, .5),
                           blocked_RA = F,
                           block_prob_each = NULL,
                           n_clust_pop = NULL,
                           clust_noise_coef = 1,
                           n_clust_samp = NULL,
                           robust_SEs = F){
  
  if(n > N){stop("Sample (n) cannot be larger than population (N)!")}
  if(m < 2){stop("Experiment requires that m > 1.")}
  if(length(ATEs) != (length(prob_each) - 1)) {stop("Error: Vector of ATEs is not of length prob_each-1!")}
  if(length(noise_scale) != m  &  length(noise_scale)!= 1)
  {stop("Error: Vector of noise scale must be of length m or 1!")}
  if(sum(noise_scale < 0) > 0){stop("Error: All noise_scale inputs must be positive!")}
  if(location_scale_X[2] <= 0){stop("Error: Standard deviation of X must be positive!")}
  if(length(coef_X)> 1 | !is.numeric(coef_X)){stop("Error: coef_X must be one number!")}
  if(sum(block_var_probs) != 1){stop("Block_var_probs must sum to 1!")}
  if(is.numeric(n_clust_pop) & is.numeric(n_clust_samp)){
    if(n_clust_pop< n_clust_samp){
      stop("Number of clusters in sample cannot exceed number of clusters in population!")}}
  if(!is.null(block_prob_each)){
    if(sum(rowSums(block_prob_each) != 1) != 0){
      stop("Rows of block_prob_each must sum to 1!")}}
  
  if(is.null(n_clust_pop)){  
    population    <- declare_population(
      noise = declare_variable(location_scale = c(0, 1)),
      X = declare_variable(location_scale = location_scale_X),
      block_var = declare_variable(type = "multinomial",
                                   probabilities = block_var_probs, outcome_categories = 1:length(block_var_probs)),
      size = N, super_population = T)
    
    write_PO_formula <- function(ATEs, noise_scale){
      scale <- ifelse(length(noise_scale) == 1, noise_scale, noise_scale[-1])
      tab <- cbind(ATEs, scale)
      
      POs <- paste0("Y ~ ", mu_Y0, " + ", noise_scale[1], " * noise * (Z == 'control') +", 
                    paste(sapply(1:nrow(tab), function(j) {
                      paste0(tab[j, 1], " * (Z == 'treatment",  j, "') + ", 
                             tab[j, 2], " * noise * (Z == 'treatment", j, "')")}), 
                      collapse = " + "), " + X * ", coef_X)
      
      return(POs)
    }
    
    
    potential_outcomes <- declare_potential_outcomes(
      condition_names = c("control", paste0("treatment", 1 : length(ATEs))),
      formula = as.formula(write_PO_formula(ATEs = ATEs, noise_scale = noise_scale)))
    
    if(N == n){
      sampling <- declare_sampling(sampling = F)
    }
    if(N != n){
      sampling <- declare_sampling(n = n)
    }
    
    if(!blocked_RA){
      assignment   <- declare_assignment(probability_each = prob_each, 
                                         potential_outcomes = potential_outcomes)
    }
    
    if(blocked_RA){
      assignment   <- declare_assignment(block_variable_name = "block_var",
                                         block_probabilities = block_prob_each,
                                         potential_outcomes = potential_outcomes)
    }
    estimand_function <- function(data, treatment_num){
      return(mean(data[,paste0("Y_Z_treatment", treatment_num)] - data$Y_Z_control))
    }
    
    if(!cov_adjustment){estimator_formula = "Y ~ Z"}
    if(cov_adjustment){estimator_formula = "Y ~ Z + X"}
    
    linear_estimator <- function(data, treatment_num){
      model <- lm(as.formula(estimator_formula), weights = Z_assignment_probabilities, data = data)
      if(!robust_SEs){
        return(get_regression_coefficient(model, coefficient_name = paste0("Ztreatment", treatment_num)))}
      if(robust_SEs){
        return(get_regression_coefficient_robust(model, coefficient_name = paste0("Ztreatment", treatment_num)))}  
    }
  }
  
  if(!is.null(n_clust_pop)){
    population    <- declare_population(
      individual = list(noise_ind = declare_variable(),
                        X = declare_variable(location_scale = location_scale_X),
                        block_var_ind = declare_variable(type = "multinomial",
                                                         probabilities = block_var_probs, 
                                                         outcome_categories = 1:length(block_var_probs))),
      cluster = list(noise_clust = declare_variable(),
                     block_var_clust = declare_variable(type = "multinomial",
                                                        probabilities = block_var_probs, 
                                                        outcome_categories = 1:length(block_var_probs))),
      size = c(N, n_clust_pop))
    
    
    write_PO_formula <- function(ATEs, noise_scale, clust_noise_coef){
      scale <- ifelse(length(noise_scale) == 1, noise_scale, noise_scale[-1])
      tab <- cbind(ATEs, scale)
      
      POs <- paste0("Y ~ ", mu_Y0, " + ", noise_scale[1], " * noise_ind * (Z == 'control') +", 
                    paste(sapply(1:nrow(tab), function(j) {
                      paste0(tab[j, 1], " * (Z == 'treatment",  j, "') + ", 
                             tab[j, 2], " * noise_ind * (Z == 'treatment", j, "')")}), 
                      collapse = " + "), " + X * ", coef_X, " + noise_clust * ", clust_noise_coef)
      
      return(POs)
    }
    
    potential_outcomes <- declare_potential_outcomes(
      condition_names = c("control", paste0("treatment", 1 : length(ATEs))),
      formula = as.formula(write_PO_formula(ATEs = ATEs, noise_scale = noise_scale, 
                                            clust_noise_coef = clust_noise_coef)))
    
    if(n_clust_pop == n_clust_samp){
      sampling <- declare_sampling(sampling = F)
    }
    
    if(n_clust_pop != n_clust_samp){
      sampling <- declare_sampling(n = n_clust_samp,
                                   cluster_variable_name = "cluster_ID")
    }
    
    
    if(!blocked_RA){
      assignment   <- declare_assignment(cluster_variable_name = "cluster_ID", 
                                         potential_outcomes = potential_outcomes)
    }
    
    if(blocked_RA){
      assignment   <- declare_assignment(block_variable_name = "block_var_clust",
                                         cluster_variable_name = "cluster_ID", 
                                         potential_outcomes = potential_outcomes)
    }
    
    
    estimand_function <- function(data, treatment_num){
      return(mean(data[,paste0("Y_Z_treatment", treatment_num)] - data$Y_Z_control))
    }
    
    if(!cov_adjustment){estimator_formula = "Y ~ Z"}
    if(cov_adjustment){estimator_formula = "Y ~ Z + X"}
    
    get_regression_coefficient_clustered <- function (model, 
                                                      formula = NULL, 
                                                      coefficient_name, 
                                                      statistics = c("est", "se", "p", "ci_lower", "ci_upper", "df"), 
                                                      label = ""){ 
      require(sandwich)
      coef_num <- which(names(coef(model)) %in% coefficient_name)
      df <- df.residual(model)
      est <- coef(model)[coef_num]
      se <- sqrt(diag(vcovCluster(model = model, cluster = data$cluster_ID)))[coef_num]
      p <- 2 * pt(abs(est/se), df = df, lower.tail = FALSE)
      conf_int <- est + se %o% qt(c(0.025, 0.975), summary(model)$df[2])
      output <- matrix(c(est, se, p, conf_int, df), 
                       dimnames = list(c("est", "se", "p", "ci_lower", "ci_upper", "df"), 
                                       paste0(summary(model)$terms[[2]], 
                                              "~", paste(all.vars(summary(model)$terms[[3]]), collapse = "+"), 
                                              "_", label)))
      return(output[which(rownames(output) %in% statistics), , drop = FALSE])
    }
    
    linear_estimator <- function(data, treatment_num){
      model <- lm(as.formula(estimator_formula), weights = Z_assignment_probabilities, data = data)
      return(get_regression_coefficient_clustered(model, coefficient_name = paste0("Ztreatment", treatment_num)))
    }
    
  }
  
  
  make_estimator <- function(treatment_num){
    estimand  <- declare_estimand(estimand_function  = estimand_function, 
                                  treatment_num = treatment_num,
                                  potential_outcomes = potential_outcomes)
    
    estimator <- declare_estimator(estimates = linear_estimator, 
                                   treatment_num = treatment_num, 
                                   estimand = estimand,
                                   labels = paste0("ATE_hat_treatment", treatment_num))
    return(estimator)
  }
  
  
  
  estimator_list <- lapply(X = 1 : length(ATEs), FUN = make_estimator)
  
  
  my_design <- declare_design(
    population         = population, 
    potential_outcomes = potential_outcomes, 
    sampling           = sampling, 
    assignment         = assignment, 
    estimator          = estimator_list)
  
  return(my_design)
}



# 2x2 Factorial Design ----------------------------------------------------


# N = population size
# n = sample size
# m = number of arms
# prob_assign = vector of units assigned to each arm (control first)
# mu_Y0 = expectation of untreated potential outcome
# ATEs = vector of ATEs for each treatment arm (length prob_each - 1)
# noise_scale = rescale the noise for each experimental group if desired
# coef_X = coefficient on covariate X
# location_scale_X = mean and SD of normally-distributed X covariate
# cov_adjustment = use covariate adjustment in estimation
# block_var_probs = distribution of multinomial blocking variable (probability per level)
# blocked_RA = use blocked random assignment
# block_prob_each = matrix of probabilities if different blocks have different probs of assignment to treatment
# n_clust_pop = number of clusters in population
# n_clust_samp = number of clusters in sample

show_condition_names <- function(factors = c(2, 2)){
  k <- prod(factors)  
  perm <- function(v) {
    sapply(1:length(v), function(x) {
      rep(rep(1:v[x], each=prod(v[x:length(v)]) / v[x]),
          length.out=prod(v))
    } ) - 1
  }
  
  ts    <-  perm(factors)
  condition_names  <- sapply(1:nrow(ts), function(j) paste(c("Z",ts[j,]), collapse = ""))
  return(list(condition_names, ts))
}

show_condition_names()[[1]] # use to provide guide for inputting condition means



factorial_template <- function(N, 
                               n, 
                               cond_means = c(0, 0, 0, 0),
                               prob_each = rep(1/4, 4),
                               noise_scale = 1,
                               coef_X = 0,
                               location_scale_X = c(0, 1),
                               cov_adjustment = F,
                               block_var_probs = c(.5, .5),
                               blocked_RA = F,
                               block_prob_each = NULL,
                               n_clust_pop = NULL,
                               clust_noise_coef = 1,
                               n_clust_samp = NULL){
  
  if(n > N){stop("Sample (n) cannot be larger than population (N)!")}
  if(length(cond_means) != 4){stop("Error: Template creates 2x2 design with 4 treament conditions!")}
  if(length(noise_scale) != 4  &  length(noise_scale)!= 1)
  {stop("Error: Vector of noise scale must be of length m or 1!")}
  if(sum(noise_scale < 0) > 0){stop("Error: All noise_scale inputs must be positive!")}
  if(location_scale_X[2] <= 0){stop("Error: Standard deviation of X must be positive!")}
  if(length(coef_X)> 1 | !is.numeric(coef_X)){stop("Error: coef_X must be one number!")}
  if(sum(block_var_probs) != 1){stop("Block_var_probs must sum to 1!")}
  if(is.numeric(n_clust_pop) & is.numeric(n_clust_samp)){
    if(n_clust_pop< n_clust_samp){
      stop("Number of clusters in sample cannot exceed number of clusters in population!")}}
  if(!is.null(block_prob_each)){
    if(sum(rowSums(block_prob_each) != 1) != 0){
      stop("Rows of block_prob_each must sum to 1!")}}
  
  # Primitives --------------------------------------------------------------
  factors <- c(2, 2)
  
  condition_names <- show_condition_names(factors = factors)[[1]]
  
  ts <- show_condition_names(factors = factors)[[2]]
  
  main_conditions <- sapply(1:length(factors), function(f) (
    sapply(1:(factors[f]-1), function(t)
      condition_names[ts[,f]== t])))
  
  # Design ------------------------------------------------------------------
  if(is.null(n_clust_pop)){
    population    <- declare_population(
      noise = declare_variable(location_scale = c(0, 1)),
      X = declare_variable(location_scale = location_scale_X),
      block_var = declare_variable(type = "multinomial",
                                   probabilities = block_var_probs, outcome_categories = 1:length(block_var_probs)),
      size = N, super_population = T)
    
    write_PO_formula <- function(cond_means, noise_scale){
      tab <- cbind(cond_means, noise_scale)
      
      POs <- paste0("Y ~ ", paste(sapply(1:nrow(tab), function(j) {
        paste0(tab[j, 1], " * (Z == '",  condition_names[j], "') + ", 
               tab[j, 2], " * noise * (Z == '", condition_names[j], "')")}), 
        collapse = " + "), " +  X  * ", coef_X)
      return(POs)
    }
    
    potential_outcomes <- declare_potential_outcomes(
      condition_names = condition_names,
      formula = as.formula(write_PO_formula(cond_means = cond_means, noise_scale = noise_scale)))
    
    sampling <- declare_sampling(n = n)
    
    if(N == n){
      sampling <- declare_sampling(sampling = F)
    }
    if(N != n){
      sampling <- declare_sampling(n = n)
    }
    
    if(!blocked_RA){
      assignment   <- declare_assignment(probability_each = prob_each, 
                                         potential_outcomes = potential_outcomes) 
    }
    
    if(blocked_RA){
      assignment   <- declare_assignment(block_variable_name = "block_var",
                                         block_probabilities = block_prob_each,
                                         potential_outcomes = potential_outcomes)
    }
    
  }
  
  if(!is.null(n_clust_pop)){
    population    <- declare_population(
      individual = list(noise_ind = declare_variable(),
                        X = declare_variable(location_scale = location_scale_X),
                        block_var_ind = declare_variable(type = "multinomial",
                                                         probabilities = block_var_probs, 
                                                         outcome_categories = 1:length(block_var_probs))),
      cluster = list(noise_clust = declare_variable(),
                     block_var_clust = declare_variable(type = "multinomial",
                                                        probabilities = block_var_probs, 
                                                        outcome_categories = 1:length(block_var_probs))),
      size = c(N, n_clust_pop))
    
    write_PO_formula <- function(cond_means, noise_scale){
      tab <- cbind(cond_means, noise_scale)
      
      POs <- paste0("Y ~ ", paste(sapply(1:nrow(tab), function(j) {
        paste0(tab[j, 1], " * (Z == '",  condition_names[j], "') + ", 
               tab[j, 2], " * noise_ind * (Z == '", condition_names[j], "')")}), 
        collapse = " + "), " + X * ", coef_X, " + noise_clust * ", clust_noise_coef)
      return(POs)
    }
    
    potential_outcomes <- declare_potential_outcomes(
      condition_names = condition_names,
      formula = as.formula(write_PO_formula(cond_means = cond_means, 
                                            noise_scale = noise_scale)))
    
    if(n_clust_pop == n_clust_samp){
      sampling <- declare_sampling(sampling = F)
    }
    
    if(n_clust_pop != n_clust_samp){
      sampling <- declare_sampling(n = n_clust_samp,
                                   cluster_variable_name = "cluster_ID")
    }
    
    
    if(!blocked_RA){
      assignment   <- declare_assignment(cluster_variable_name = "cluster_ID", 
                                         potential_outcomes = potential_outcomes)
    }
    
    if(blocked_RA){
      assignment   <- declare_assignment(block_variable_name = "block_var_clust",
                                         cluster_variable_name = "cluster_ID", 
                                         potential_outcomes = potential_outcomes)
    }
    
  }  
  
  
  estimand_function <- function(data, estimand_num){ 
    estimand_weights <- matrix(c( -1, 0 , 1, 0,
                                  -1, 1, 0, 0,
                                  0, -1, 0, 1,
                                  0, 0, -1, 1,
                                  -.5, -.5, .5, .5,
                                  -0.5, .5, -.5, .5,
                                  1, -1, -1, 1), nrow = 7, ncol = 4, byrow = T)
    exp <- paste0(estimand_weights[estimand_num,], "*", "data[,'Y_Z_", 
                  condition_names, "']", collapse = "+")
    return(mean(eval(parse(text = exp))))
  }
  
  
  estimator_function <- function(data, estimand_num){
    for(j in 1:nrow(main_conditions)){data[paste0("T", j)] <- 1 * (data$Z %in% main_conditions[,j])}
    
    if(!cov_adjustment){estimator_formula = "Y ~ T1 * T2"}
    if(cov_adjustment){estimator_formula = "Y ~ T1 * T2 + X"}
    
    model <- lm(estimator_formula, data = data, weights = Z_assignment_weights)
    if(!is.null(n_clust_pop)){
      rcse_vcov <- (vcovCluster(model, cluster = data$cluster_ID))}
    
    get_regression_coefficient_clustered <- function (model, 
                                                      formula = NULL, 
                                                      coefficient_name, 
                                                      statistics = c("est", "se", "p", "ci_lower", "ci_upper", "df"), 
                                                      label = ""){ 
      require(sandwich)
      coef_num <- which(names(coef(model)) %in% coefficient_name)
      df <- df.residual(model)
      est <- coef(model)[coef_num]
      se <- sqrt(diag(vcovCluster(model = model, cluster = data$cluster_ID)))[coef_num]
      p <- 2 * pt(abs(est/se), df = df, lower.tail = FALSE)
      conf_int <- est + se %o% qt(c(0.025, 0.975), summary(model)$df[2])
      output <- matrix(c(est, se, p, conf_int, df), 
                       dimnames = list(c("est", "se", "p", "ci_lower", "ci_upper", "df"), 
                                       paste0(summary(model)$terms[[2]], 
                                              "~", paste(all.vars(summary(model)$terms[[3]]), collapse = "+"), 
                                              "_", label)))
      return(output[which(rownames(output) %in% statistics), , drop = FALSE])
    }
    
    if(estimand_num == 1){ifelse(is.null(n_clust_pop),
                                 estimate <- get_regression_coefficient(model, coefficient_name = "T1"),
                                 estimate <- get_regression_coefficient_clustered(model, coefficient_name = "T1"))}
    
    if(estimand_num == 2){ifelse(is.null(n_clust_pop),
                                 estimate <- get_regression_coefficient(model, coefficient_name = "T2"),
                                 estimate <- get_regression_coefficient_clustered(model, coefficient_name = "T2"))}
    
    if(estimand_num == 3){
      est      <- summary(model)$coef["T1",1] + summary(model)$coef["T1:T2", 1]
      ifelse(is.null(n_clust_pop),
             se       <- sqrt(vcov(model)["T1", "T1"] + vcov(model)["T1:T2", "T1:T2"] + 
                                2 * vcov(model)["T1", "T1:T2"]),
             se       <- sqrt(rsce_vcov["T1", "T1"] + rsce_vcov["T1:T2", "T1:T2"] + 
                                2 * rsce_vcov["T1", "T1:T2"]))
      df       <- model$df.residual
      p        <- pt(q = abs(est/se), df = df, lower.tail = F)
      ci_lower <- est - 1.96*se
      ci_upper <- est + 1.96*se
      estimate <- matrix(c(est, se, p, ci_lower, ci_upper, df), 
                         dimnames = list(c("est", "se", "p", "ci_lower", "ci_upper", "df"), "regression"))}
    
    if(estimand_num == 4){
      est      <- summary(model)$coef["T2",1] + summary(model)$coef["T1:T2", 1]
      ifelse(is.null(n_clust_pop),
             se       <- sqrt(vcov(model)["T2", "T2"] + vcov(model)["T1:T2", "T1:T2"] + 
                                2 * vcov(model)["T2", "T1:T2"]),
             se       <- sqrt(rsce_vcov["T2", "T2"] + rsce_vcov["T1:T2", "T1:T2"] + 
                                2 * rsce_vcov["T2", "T1:T2"]))
      df       <- model$df.residual
      p        <- pt(q = abs(est/se), df = df, lower.tail = F)
      ci_lower <- est - 1.96*se
      ci_upper <- est + 1.96*se
      estimate <- matrix(c(est, se, p, ci_lower, ci_upper, df), 
                         dimnames = list(c("est", "se", "p", "ci_lower", "ci_upper", "df"), "regression")) }
    
    if(estimand_num == 5){
      est      <- summary(model)$coef["T1",1] + 0.5 * summary(model)$coef["T1:T2", 1]
      ifelse(is.null(n_clust_pop),
             se       <- sqrt(vcov(model)["T1", "T1"] + 0.25 * vcov(model)["T1:T2", "T1:T2"] + 
                                vcov(model)["T1", "T1:T2"]),
             se       <- sqrt(rsce_vcov["T1", "T1"] + 0.25 * rsce_vcov["T1:T2", "T1:T2"] + 
                                rsce_vcov["T1", "T1:T2"]))
      df       <- model$df.residual
      p        <- pt(q = abs(est/se), df = df, lower.tail = F)
      ci_lower <- est - 1.96*se
      ci_upper <- est + 1.96*se
      estimate <- matrix(c(est, se, p, ci_lower, ci_upper, df), 
                         dimnames = list(c("est", "se", "p", "ci_lower", "ci_upper", "df"), "regression"))}
    if(estimand_num == 6){
      est      <- summary(model)$coef["T2",1] + 0.5 * summary(model)$coef["T1:T2", 1]
      ifelse(is.null(n_clust_pop),
             se       <- sqrt(vcov(model)["T2", "T2"] + 0.25 * vcov(model)["T1:T2", "T1:T2"] + 
                                vcov(model)["T2", "T1:T2"]),
             se       <- sqrt(rsce_vcov["T2", "T2"] + 0.25 * rsce_vcov["T1:T2", "T1:T2"] + 
                                rsce_vcov["T2", "T1:T2"]))
      df       <- model$df.residual
      p        <- pt(q = abs(est/se), df = df, lower.tail = F)
      ci_lower <- est - 1.96*se
      ci_upper <- est + 1.96*se
      estimate <- matrix(c(est, se, p, ci_lower, ci_upper, df), 
                         dimnames = list(c("est", "se", "p", "ci_lower", "ci_upper", "df"), "regression"))}
    
    if(estimand_num == 7){ifelse(is.null(n_clust_pop),
                                 estimate <- get_regression_coefficient(model, coefficient_name = "T1:T2"),
                                 estimate <- get_regression_coefficient_clustered(model, coefficient_name = "T1:T2"))}
    
    return(estimate)
  }
  
  make_estimator <- function(estimand_num){
    estimand  <- declare_estimand(estimand_function = estimand_function, 
                                  estimand_num = estimand_num,
                                  potential_outcomes = potential_outcomes)
    labs = c("ME_T1|T2 = 0", "ME_T2|T1 = 0", "ME_T1|T2 = 1", "ME_T2|T1 = 1",
             "Avg_ME_T1", "Avg_ME_T2", "CME_T1_T2")
    estimator <- declare_estimator(estimates = estimator_function, 
                                   estimand_num = estimand_num, 
                                   estimand = estimand,
                                   labels = labs[estimand_num])
    return(estimator)
  }
  
  
  
  
  estimator_list <- lapply(X = 1 : 7, FUN = make_estimator)
  
  # Declare Design ----------------------------------------------------------
  
  my_design <- declare_design(
    population         = population, 
    potential_outcomes = potential_outcomes, 
    sampling           = sampling, 
    assignment         = assignment, 
    estimator          = estimator_list)
  
  return(my_design)
  
}
