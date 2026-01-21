eval_pois_llhood_homeaway = function(X1, X2, att, def, home, X_mid) {
  match_idx = which(X_mid == 1, arr.ind = TRUE)   # Get all match indices where a game was played
  
  i = match_idx[, 1]  # home teams
  j = match_idx[, 2]  # away teams
  
  # Home team goals
  lambda_H = exp(home + att[i] + def[j])
  loglik_H = dpois(X1[cbind(i, j)], lambda_H, log = TRUE)
  
  # Away team goals
  lambda_A = exp(att[j] + def[i])
  loglik_A = dpois(X2[cbind(i, j)], lambda_A, log = TRUE)
  
  total_loglik = sum(loglik_H) + sum(loglik_A)
  return(total_loglik)
}


eval_pois_llhood_home = function(X1, att, def, home, X_mid) {
  match_idx = which(X_mid == 1, arr.ind = TRUE)   # Get all match indices where a game was played
  
  i = match_idx[, 1]  # home teams
  j = match_idx[, 2]  # away teams
  
  lambda_H = exp(home + att[i] + def[j])
  
  loglik_H = dpois(X1[cbind(i, j)], lambda_H, log = TRUE)
  
  total_loglik = sum(loglik_H, na.rm = TRUE)
  
  return(total_loglik)
}

MH_Pois = function(X1, X2, att_0, def_0, home_0, X_mid = FALSE, iter=100,  
                       sd_prop_att = 1, att_mean_prior = 0, att_sd_prior = 1,
                       sd_prop_def = 1, def_mean_prior = -0, def_sd_prior = 1,
                       sd_prop_home = 1, home_mean_prior = 0, home_sd_prior = 1,
                       verbosity = 0, print_by = 1000, fix_idx = NA){
  
  N = nrow(X1)
  
  if(!is.matrix(X_mid) || any(dim(X_mid) != N)) { #Check that X_mid matrix is correct
    X_mid = matrix(1, N, N)
  }
  diag(X_mid) = 0 #Force no diagonals (self-matches) for safety
  
  #Initialize Posterior Chains
  att_post = matrix(NA_real_, nrow = iter, ncol = N)
  def_post = matrix(NA_real_, nrow = iter, ncol = N)
  home_post = numeric(iter)
  
  att_post[1,] = att_0
  def_post[1,] = def_0
  home_post[1] = home_0
  
  #Star variables indicate currently "active" state of chain
  att_star = att_0
  def_star = def_0
  home_star = home_0
  
  #Acceptance history to compute rates
  acc_history = matrix(0, nrow = iter, ncol = 3)
  colnames(acc_history) = c("att", "def", "home")
  
  for (t in 2:iter){
    
    #### Att update block
    ########################################
    
    att_prop = att_star + rnorm(N, 0, sd_prop_att)
    
    
    if (!is.na(fix_idx)){
      stopifnot(fix_idx %% 1 == 0 & fix_idx>= 0 & fix_idx<= N)
      
      att_prop[fix_idx] = -sum(att_prop[-fix_idx])
    } else {
      att_prop = att_prop - mean(att_prop)
    }
    
    
    log_lik_prop = eval_pois_llhood_homeaway(X1, X2, att_prop, def_star, home_star, X_mid)
    log_prior_prop = sum(dnorm(att_prop, att_mean_prior, att_sd_prior, log = TRUE))
    
    log_lik_prev = eval_pois_llhood_homeaway(X1, X2, att_star, def_star, home_star, X_mid)
    log_prior_prev = sum(dnorm(att_star, att_mean_prior, att_sd_prior, log = TRUE))
    
    alpha = min(1, exp(log_lik_prop - log_lik_prev + log_prior_prop - log_prior_prev))
    u = fast_runif(1)
    if (u <= alpha){
      att_star = att_prop
      acc_history[t, "att"] = 1
    }      
    
    att_post[t,] = att_star
    
    #### Def update block
    ######################################################
    
    def_prop = def_star + rnorm(N, 0 , sd_prop_def)
    if (!is.na(fix_idx)){
      stopifnot(fix_idx %% 1 == 0 & fix_idx>= 0 & fix_idx<= N)
      
      def_prop[fix_idx] = -sum(def_prop[-fix_idx])
    } else {
      def_prop = def_prop - mean(def_prop)
    }
    
    log_lik_prop = eval_pois_llhood_homeaway(X1, X2, att_star, def_prop, home_star, X_mid)
    log_prior_prop = sum(dnorm(def_prop, def_mean_prior, def_sd_prior, log = TRUE))
    
    log_lik_prev = eval_pois_llhood_homeaway(X1, X2, att_star, def_star, home_star, X_mid)
    log_prior_prev = sum(dnorm(def_star, def_mean_prior, def_sd_prior, log = TRUE))
    
    alpha = min(1, exp(log_lik_prop - log_lik_prev + log_prior_prop - log_prior_prev)) 
    u = fast_runif(1)
    if (u <= alpha){
      def_star = def_prop
      acc_history[t, "def"] = 1
    }     
    
    def_post[t,] = def_star
    
    #### Home update block
    ######################################################
    
    home_prop = rnorm(1, home_star, sd_prop_home)
    
    log_lik_prop = eval_pois_llhood_home(X1, att_star, def_star, home_prop, X_mid)
    log_prior_prop = dnorm(home_prop, home_mean_prior, home_sd_prior, log = TRUE)
    
    log_lik_prev = eval_pois_llhood_home(X1, att_star, def_star, home_star, X_mid)
    log_prior_prev = dnorm(home_star, home_mean_prior, home_sd_prior, log = TRUE)
    
    alpha = min(1, exp(log_lik_prop - log_lik_prev + log_prior_prop - log_prior_prev))
    u = fast_runif(1)
    if (u <= alpha){
      home_star = home_prop
      acc_history[t, "home"] = 1
    }
    home_post[t] = home_star
    
    #Verbosity settings (additional comments / visuals during model fit, for debugging)
    if (verbosity >= 1 && t == 2) {pb = txtProgressBar(min = 0, max = iter, style = 3)}
    
    if (verbosity >= 1) {
      setTxtProgressBar(pb, t)
      if (t == iter) {close(pb)}
    }
    
    if (verbosity >= 2 && t %% print_by == 0) {
      start_idx = max(1, t - print_by - 1)
      recent_acc = colMeans(acc_history[start_idx:t, , drop = FALSE])
      
      cat(sprintf(
        "Iter %d | acc_last_1000: att=%.2f, def=%.2f, home=%.2f\n",
        t, recent_acc["att"], recent_acc["def"], recent_acc["home"]
      ))
    }
    
    if (verbosity >= 3 && t %% print_by == 0) {
      start_idx = max(1, t - print_by - 1)
      end_idx   = t
      
      # get last {print_by} iterations means for each parameter
      att_mean  = colMeans(att_post[start_idx:end_idx, ])
      def_mean  = colMeans(def_post[start_idx:end_idx, ])
      home_mean = mean(home_post[start_idx:end_idx])
      
      cat("\n=== Iteration", t, "(mean over last", print_by, "iterations) ===\n")
      cat(sprintf("home = %.3f\n", home_mean))
      cat("ATT:", paste(sprintf("%.3f", att_mean), collapse = " "), "\n")
      cat("DEF:", paste(sprintf("%.3f", def_mean), collapse = " "), "\n\n")
    }
    
  }
  out_list = list(att_post = att_post,
                  def_post = def_post,
                  home_post = home_post,
                  iterations = t,
                  priors = list(
                    att  = list(mean = att_mean_prior,  sd = att_sd_prior),
                    def  = list(mean = def_mean_prior,  sd = def_sd_prior),
                    home = list(mean = home_mean_prior, sd = home_sd_prior)
                  ),
                  proposals = list(
                    att  = sd_prop_att,
                    def  = sd_prop_def,
                    home = sd_prop_home
                  ),
                  fix_idx = fix_idx,
                  X_mid = X_mid,
                  X = list(X1 = X1, X2 = X2),
                  team_names = rownames(X1),
                  distr_type = 'P',
                  constraint = 'STZ',
                  league = league_acro)
  class(out_list) = "MH_posterior"
  return(out_list)
}
