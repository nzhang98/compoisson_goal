eval_pois_llhood_team = function(i, X1, X2, att, def, home, X_mid) {
  N = length(att)
  j_idx = base::setdiff(1:N, i)
  
  # Home matches (i at home vs j away)
  played_home_idx = (X_mid[i, j_idx] == 1) #Only consider played matches
  j_h = j_idx[played_home_idx]
  lambda_H = exp(home + att[i] + def[j_h])
  loglik_H = dpois(X1[i, j_h], lambda_H, log = TRUE)
  
  # Away matches (i away vs j home)
  played_away_idx = (X_mid[j_idx, i] == 1) #Only consider played matches
  j_a = j_idx[played_away_idx]
  lambda_A = exp(att[i] + def[j_a])
  loglik_A = dpois(X2[j_a, i], lambda_A, log = TRUE)
  
  total_loglik = sum(loglik_H, na.rm = TRUE) + sum(loglik_A, na.rm = TRUE)
  
  return(total_loglik)
}

eval_cmp_llhood_team = function(i, X1, X2, att, def, home, nu, X_mid) {
  N = length(att)
  j_idx = base::setdiff(1:N, i)
  
  # Home matches (i at home vs j away)
  played_home_idx = (X_mid[i, j_idx] == 1) #Only consider played matches
  j_h = j_idx[played_home_idx]
  mu_H = exp(home + att[i] + def[j_h])
  loglik_H = eval_cmp_llhood_lookup(X1[i, j_h], mu_H, nu[i])
  
  # Away matches (i away vs j home)
  played_away_idx = (X_mid[j_idx, i] == 1) #Only consider played matches
  j_a = j_idx[played_away_idx]
  mu_A = exp(att[i] + def[j_a])
  loglik_A = eval_cmp_llhood_lookup(X2[j_a, i], mu_A, nu[i])
  
  total_loglik = sum(loglik_H, na.rm = TRUE) + sum(loglik_A, na.rm = TRUE)
  
  return(total_loglik)
}

eval_log_qf_block = function(block_idx, X1, X2, att, def, home, nu, X_mid) {
  N = length(att)
  
  total_log_qf = 0
  
  for (i in block_idx){
    j_idx = base::setdiff(1:N, i)
    
    # Home matches (i at home vs j away)
    played_home_idx = (X_mid[i, j_idx] == 1) #Only consider played matches
    j_h = j_idx[played_home_idx]
    mu_H = exp(home + att[i] + def[j_h])
    log_qf_H = log_qf_cmp(X1[i, j_h], mu_H, nu[i])
    
    # Away matches (i away vs j home)
    played_away_idx = (X_mid[j_idx, i] == 1) #Only consider played matches
    j_a = j_idx[played_away_idx]
    mu_A = exp(att[i] + def[j_a])
    log_qf_A = log_qf_cmp(X2[j_a, i], mu_A, nu[i])
    
    total_log_qf = total_log_qf + sum(log_qf_H, na.rm = TRUE) + sum(log_qf_A, na.rm = TRUE)
  }
  
  return(total_log_qf)
}

eval_log_qf_team = function(i, X1, X2, att, def, home, nu, X_mid) {
  N = length(att)
  j_idx = base::setdiff(1:N, i)
  
  # Home matches (i at home vs j away)
  played_home_idx = (X_mid[i, j_idx] == 1) #Only consider played matches
  j_h = j_idx[played_home_idx]
  mu_H = exp(home + att[i] + def[j_h])
  log_qf_H = log_qf_cmp(X1[i, j_h], mu_H, nu[i])
  
  # Away matches (i away vs j home)
  played_away_idx = (X_mid[j_idx, i] == 1) #Only consider played matches
  j_a = j_idx[played_away_idx]
  mu_A = exp(att[i] + def[j_a])
  log_qf_A = log_qf_cmp(X2[j_a, i], mu_A, nu[i])
  
  total_log_qf = sum(log_qf_H, na.rm = TRUE) + sum(log_qf_A, na.rm = TRUE)
  
  return(total_log_qf)
}

eval_log_qf_home = function(X1, att, def, home, nu, X_mid) {
  match_idx = which(X_mid == 1, arr.ind = TRUE)   # Get all match indices where a game was played
  
  i = match_idx[, 1]  # home teams
  j = match_idx[, 2]  # away teams
  
  mu_H = exp(home + att[i] + def[j])
  
  log_qf_H = dpois(X1[cbind(i, j)], mu_H, nu[i])
  
  total_loglik = sum(log_qf_H, na.rm = TRUE)
  
  return(total_loglik)
}


eval_log_qf_homeaway = function(X1, X2, att, def, home, nu, X_mid) {
  match_idx = which(X_mid == 1, arr.ind = TRUE)   # Get all match indices where a game was played
  
  i = match_idx[, 1]  # home teams
  j = match_idx[, 2]  # away teams
  
  # Home team goals
  mu_H = exp(home + att[i] + def[j])
  log_qf_H = log_qf_cmp(X1[cbind(i, j)], mu_H, nu[i])
  
  # Away team goals
  mu_A = exp(att[j] + def[i])
  log_qf_A = log_qf_cmp(X2[cbind(i, j)], mu_A, nu[j])
  
  total_log_qf = sum(log_qf_H) + sum(log_qf_A)
  return(total_log_qf)
}


eval_cmp_llhood_lookup = function(x, mu, nu) {
  # sanity: nu must be scalar for this function
  if (length(nu) != 1L) stop("This function expects scalar 'nu'.")
  mu = as.numeric(mu)
  x = as.numeric(x)
  if (length(x) != length(mu)) x = rep(x, length.out = length(mu))
  n = length(mu)
  if (n == 0) return(numeric(0))
  
  # preallocate
  Z_vals = numeric(n)
  
  # in-range tests
  mu_in_range = mu >= min(mu_vals) & mu <= max(mu_vals)
  nu_in_range = (nu >= min(nu_vals) & nu <= max(nu_vals))
  
  # compute indices (row per mu, single col for scalar nu)
  row_idx = nearest_index(mu, mu_vals)    # length n
  col_idx = nearest_index(nu, nu_vals)[1] # scalar
  
  # 1) use lookup for observations where both mu and nu are in-range
  in_idx = which(mu_in_range & nu_in_range)
  if (length(in_idx) > 0) {
    # elementwise matrix lookup using cbind(rows, cols)
    Z_vals[in_idx] = Z_lookup[cbind(row_idx[in_idx], rep(col_idx, length(in_idx)))]
  }
  
  # 2) compute Z on-the-fly for out-of-range observations
  out_idx = which(!mu_in_range | !nu_in_range)
  if (length(out_idx) > 0) {
    Z_vals[out_idx] = vapply(out_idx, function(id) compute_Z(mu[id], nu), numeric(1))
  }
  
  # compute log-likelihood
  log_qf = nu * (x * log(mu) - lgamma(x + 1))
  log_qf - log(Z_vals)
}

log_qf_cmp = function(x, mu, nu){
  return(nu * (x * log(mu) - lgamma(x + 1) ) )
}

#### Z lookup functions
#########################################################################
compute_Z = function(mu, nu, j_max = 1000, tol = 1e-20) {
  Z = 0
  for (j in 0:j_max) {
    log_term = (nu * j) * log(mu) - nu * lfactorial(j)
    term = exp(log_term)
    if (is.nan(term) || term < tol) break
    Z = Z + term
  }
  return(Z)
}

nearest_index = function(x, grid) {
  idx = findInterval(x, grid)
  idx[idx == 0] = 1
  too_high = idx == length(grid)
  idx[too_high] = length(grid)
  
  # adjust to whichever endpoint is closer
  left  = grid[idx]
  right = grid[pmin(idx + 1, length(grid))]
  idx + (abs(right - x) < abs(left - x))
}

#### Aux data generation
#########################################################################

generate_full_league = function(N, att_vector, def_vector, nu_vector, X_mid,
                                 home_coef = 0, home_fl = TRUE) {
  
  X_sim = matrix(NA, nrow = N, ncol = N)
  
  for (i in 1:N) {       # row = team i
    for (j in 1:N) {     # column = team j
      if (i == j) next                 
      if (X_mid[i,j] == 0) next        
      
      if (home_fl) {
        # Home team is i, away team is j
        mu = exp(att_vector[i] + def_vector[j] + home_coef)
        nu_att = nu_vector[i]  # attacking team
      } else {
        # Away team is j, home team is i
        mu = exp(att_vector[j] + def_vector[i])
        nu_att = nu_vector[j]  # attacking team
      }
      
      # Draw goals using rejection sampler
      X_sim[i,j] = rejection_sampler_draws(1, mu = mu, nu = nu_att)
    }
  }
  
  return(X_sim)
}

generate_block_goals = function(team_block, N, att_vector, def_vector, nu_vector, X_mid, 
                                home_coef = 0, home_fl = TRUE) {
  # simulated league
  X_sim = matrix(NA, nrow = N, ncol = N)
  
  for (team_i in team_block) {
    j_idx = setdiff(1:N, team_i)
    
    if (home_fl) {
      # Team_i at home against opponents j
      played_idx = (X_mid[team_i, j_idx] == 1)
      j_played = j_idx[played_idx]
      
      mu_vec = exp(att_vector[team_i] + def_vector[j_played] + home_coef)
      for (k in seq_along(j_played)) {
        X_sim[team_i, j_played[k]] = rejection_sampler_draws(
          1, mu = mu_vec[k], nu = nu_vector[team_i]
        )
      }
      
    } else {
      # Team_i away against opponents j
      played_idx = (X_mid[j_idx, team_i] == 1)
      j_played = j_idx[played_idx]
      
      mu_vec = exp(att_vector[team_i] + def_vector[j_played])
      
      for (k in seq_along(j_played)) {
        X_sim[j_played[k], team_i] = rejection_sampler_draws(
          1, mu = mu_vec[k], nu = nu_vector[team_i]
        )
      }
    }
  }
  
  return(X_sim)
}

generate_team_goals = function(team_i, N, att_vector, def_vector, nu_vector, X_mid, 
                               home_coef = 0, home_fl = TRUE) {
  X_sim = matrix(NA, nrow = N, ncol = N)
  
  j_idx = base::setdiff(1:N, team_i) 
  
  if (home_fl) {
    # Home matches: team_i at home
    played_idx = (X_mid[team_i, j_idx] == 1)
    j_played = j_idx[played_idx]
    
    mu_vec = exp(att_vector[team_i] + def_vector[j_played] + home_coef)
    for (k in seq_along(j_played)) {
      X_sim[team_i, j_played[k]] = rejection_sampler_draws(1, mu = mu_vec[k], nu = nu_vector[team_i])
    }
    
  } else {
    # Away matches: team_i away
    played_idx = (X_mid[j_idx, team_i] == 1)
    j_played = j_idx[played_idx]
    
    mu_vec = exp(att_vector[team_i] + def_vector[j_played])
    for (k in seq_along(j_played)) {
      X_sim[j_played[k], team_i] = rejection_sampler_draws(1, mu = mu_vec[k], nu = nu_vector[team_i])
    }
    
  }
  return(X_sim)
}


MH_CMP_SAS = function(X1, X2, att_0, def_0, home_0, Z_0, p_0, nu_0, X_mid = FALSE, iter=100, 
                      sd_prop_att = 1, att_mean_prior = 0, att_sd_prior = 1,
                      sd_prop_def = 1, def_mean_prior = 0, def_sd_prior = 1,
                      sd_prop_home = 1, home_mean_prior = 0, home_sd_prior = 1,
                      sd_prop_nu = 0.5, nu_lmean_prior = 0, nu_sd_prior = 1, rho = 0.5, 
                      p_alpha_prior = 1, p_beta_prior = 1, 
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
  nu_post = matrix(NA_real_, nrow = iter, ncol = N)
  nu_store_mat = matrix(NA_real_, nrow = iter, ncol = N)
  nu_fix_mat = matrix(NA_real_, nrow = iter, ncol = N)
  Z_post = matrix(NA_real_, nrow = iter, ncol = N)
  p_post = matrix(NA_real_, nrow = iter, ncol = N)
  
  att_post[1,] = att_0
  def_post[1,] = def_0
  home_post[1] = home_0
  nu_post[1,] = nu_0
  nu_store_mat[1,] = nu_0
  Z_post[1,] = Z_0
  p_post[1,] = p_0
  
  #Star variables indicate currently "active" state of chain
  att_star = att_0
  def_star = def_0
  home_star = home_0
  nu_star = nu_0
  nu_store = nu_0
  Z_star = Z_0
  p_star = p_0
  
  #Pre-comput (fixed) covariance matrix of the correlated att-nu update
  Sigma = matrix(c(sd_prop_att^2, rho * sd_prop_att * sd_prop_nu,
                   rho * sd_prop_att * sd_prop_nu, sd_prop_nu^2),
                 nrow = 2, 
                 byrow = TRUE
                 )
  sd_phi_prop_fix = sqrt( (1 - rho**2) * sd_prop_nu**2 )
  
  #Acceptance history to compute rates
  acc_history = matrix(0, nrow = iter, ncol = 3)
  colnames(acc_history) = c("att", "def", "home")
  
  # att_acc_history = matrix(0, iter, N)
  nu_store_acc_history = matrix(0L, nrow = iter, ncol = N)
  nu_acc_history = matrix(0L, nrow = iter, ncol = N)
  c = 0
  
  for (t in 2:iter){
    
    #### Z update
    ########################################
    for (i in 1:N){
      
      # Compute w0_i (Pois) and w1_i (CMP) weights
      
      pois_lhood = log(1-p_star[i]) + eval_pois_llhood_team(i, X1, X2, att_star, def_star, home_star, X_mid)
      cmp_lhood = log(p_star[i]) + eval_cmp_llhood_team(i, X1, X2, att_star, def_star, home_star, nu_store, X_mid)
      
      p_z0 = plogis(pois_lhood - cmp_lhood)
      
      u = fast_runif(1)
      if (u <= p_z0) {Z_star[i] = 0} else {Z_star[i] = 1}
      
    }
    # Force Z = 1 for all i for a fully-CMP model
    # Z_star = rep(1, N)
    Z_post[t,] = Z_star
    
    #### p update
    ########################################
    
    #Conjugate (block) update
    p_star = rbeta(N, p_alpha_prior + Z_star, p_beta_prior + 1 - Z_star)
    p_post[t,] = p_star
    
    
    #### Joint ATT and NU updates
    ########################################
    
    for (i in base::setdiff(1:20, fix_idx)){ #iterate over i, except fixed team
      att_prop = att_star
      nu_prop = nu_store
      phi_prop = log(nu_prop)
      
      prop = MASS::mvrnorm(1,
                           c(att_star[i], phi_prop[i]),
                           Sigma,
      )
      att_prop[i] = prop[1]
      phi_prop[i] = prop[2]
      
      att_prop[fix_idx] = -sum(att_prop[-fix_idx])
      delta_att_prop = att_prop[fix_idx] - att_star[fix_idx]
      
      mean_phi_prop_fix = phi_prop[fix_idx] + rho * (sd_prop_nu / sd_prop_att) * delta_att_prop
      phi_prop[fix_idx] = rnorm(1, mean_phi_prop_fix, sd_phi_prop_fix)
      
      nu_prop = exp(phi_prop)
      
      mask_nu_store = exp(Z_star * log(nu_store))
      mask_nu_prop = exp(Z_star * phi_prop)
      
      block_idx = c(i, fix_idx)
      
      aux_X1 = generate_block_goals(block_idx, N, att_prop, def_star, mask_nu_prop, X_mid, home_star, home_fl = TRUE)
      aux_X2 = generate_block_goals(block_idx, N, att_prop, def_star, mask_nu_prop, X_mid, home_fl = FALSE)
      
      log_lik_prop = eval_log_qf_block(block_idx, X1, X2, att_prop, def_star, home_star, mask_nu_prop, X_mid)
      log_prior_prop = sum(dlnorm(nu_prop[block_idx], nu_lmean_prior, nu_sd_prior, log = TRUE) + phi_prop[block_idx] ) +
        sum(dnorm(att_prop[block_idx], att_mean_prior, def_sd_prior, log = TRUE))
      log_lik_prev_aux = eval_log_qf_block(block_idx, aux_X1, aux_X2, att_star, def_star, home_star, mask_nu_store, X_mid)
      
      numerator = log_lik_prop + log_prior_prop + log_lik_prev_aux
      
      log_lik_prev = eval_log_qf_block(block_idx, X1, X2, att_star, def_star, home_star, mask_nu_store, X_mid)
      log_prior_prev = sum(dlnorm(nu_store[block_idx], nu_lmean_prior, nu_sd_prior, log = TRUE) + log(nu_store[block_idx]) ) +
        sum(dnorm(att_star[block_idx], att_mean_prior, def_sd_prior, log = TRUE))
      log_lik_prop_aux = eval_log_qf_block(block_idx, aux_X1, aux_X2, att_prop, def_star, home_star, mask_nu_prop, X_mid)
      
      denominator = log_lik_prev + log_prior_prev + log_lik_prop_aux
      
      alpha = min(1, exp(numerator - denominator))
      
      u = fast_runif(1)
      if (u <= alpha){
        att_star[block_idx] = att_prop[block_idx]
        nu_store[block_idx] = nu_prop[block_idx]
        nu_acc_history[t, i] = 1
      }
      nu_fix_mat[t, i] = nu_store[fix_idx]
    }

    att_post[t,] = att_star
    nu_store_mat[t,] = nu_store
    nu_post[t,] = exp(log(nu_store) * Z_star)

    #### Def update block
    #########################################################################
    
    def_prop = def_star + rnorm(N, 0 , sd_prop_def)
    # 
    if (!is.na(fix_idx)){
      stopifnot(fix_idx %% 1 == 0 & fix_idx>= 0 & fix_idx<= N)
      
      def_prop[fix_idx] = -sum(def_prop[-fix_idx])
    } else {
      def_prop = def_prop - mean(def_prop)
    }
    
    aux_X1 = generate_full_league(N, att_star, def_prop, nu_star, X_mid, home_star, home_fl = TRUE)
    aux_X2 = generate_full_league(N, att_star, def_prop, nu_star, X_mid, home_fl = FALSE)
    
    log_lik_prop = eval_log_qf_homeaway(X1, X2, att_star, def_prop, home_star, nu_star, X_mid)
    log_prior_prop = sum(dnorm(def_prop, def_mean_prior, def_sd_prior, log = TRUE))
    log_lik_prev_aux = eval_log_qf_homeaway(aux_X1, aux_X2, att_star, def_star, home_star, nu_star, X_mid)
    
    numerator = log_lik_prop + log_prior_prop + log_lik_prev_aux
    
    log_lik_prev = eval_log_qf_homeaway(X1, X2, att_star, def_star, home_star, nu_star, X_mid)
    log_prior_prev = sum(dnorm(def_star, def_mean_prior, def_sd_prior, log = TRUE))
    log_lik_prop_aux = eval_log_qf_homeaway(aux_X1, aux_X2, att_star, def_prop, home_star, nu_star, X_mid)
    
    denominator = log_lik_prev + log_prior_prev + log_lik_prop_aux
    
    alpha = min(1, exp(numerator - denominator))
    
    u = fast_runif(1)
    if (u <= alpha){
      def_star = def_prop
      acc_history[t, "def"] = 1
    }
    def_post[t,] = def_star
    
    #### Home update block
    #########################################################################    
    
    home_prop = rnorm(1, home_star, sd = sd_prop_home)
    
    aux_X1 = generate_full_league(N, att_star, def_star, nu_star, X_mid, home_prop, home_fl = TRUE)
    
    log_lik_prop = eval_log_qf_home(X1, att_star, def_star, home_prop, nu_star, X_mid)
    log_prior_prop = dnorm(home_prop, home_mean_prior, home_sd_prior, log = TRUE)
    log_lik_prev_aux = eval_log_qf_home(aux_X1, att_star, def_star, home_star, nu_star, X_mid)
    
    numerator = log_lik_prop + log_prior_prop + log_lik_prev_aux
    
    log_lik_prev = eval_log_qf_home(X1, att_star, def_star, home_star, nu_star, X_mid)
    log_prior_prev = dnorm(home_star, home_mean_prior, home_sd_prior, log = TRUE)
    log_lik_prop_aux = eval_log_qf_home(aux_X1, att_star, def_star, home_prop, nu_star, X_mid)
    
    denominator = log_lik_prev + log_prior_prev + log_lik_prop_aux
    
    alpha = min(1, exp(numerator - denominator))
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
      Z_1_hist = Z_post[start_idx:t,]
      # recent_acc_store = colMeans(nu_store_acc_history[start_idx:t, , drop = FALSE])/(1-Z_1_hist)
      recent_acc_nu_post = colSums(nu_acc_history[start_idx:t, , drop = FALSE]*Z_1_hist)/colSums(Z_1_hist)
      recent_acc_nu_total = colMeans(nu_acc_history[start_idx:t, , drop = FALSE])
      # recent_acc_att = colMeans(att_acc_history[start_idx:t, , drop = FALSE])
      
      print(paste0("\nNu_total acc: ", recent_acc_nu_total, " Nu_post acc: ", recent_acc_nu_post))
      # , " Att acc:", recent_acc_att))
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
      nu_mean   = colMeans(nu_post[start_idx:end_idx, ])
      home_mean = mean(home_post[start_idx:end_idx])
      z1_mean = colMeans(Z_post[start_idx:end_idx, ])
      
      cat("\n=== Iteration", t, "(mean over last", print_by, "iterations) ===\n")
      cat("home = %.3f\n", home_mean))
      cat("ATT:", paste(sprintf("%.3f", att_mean), collapse = " "), "\n")
      cat("DEF:", paste(sprintf("%.3f", def_mean), collapse = " "), "\n")
      cat("NU :", paste(sprintf("%.3f", nu_mean), collapse = " "), "\n")
      cat("Z1 :", paste(sprintf("%.3f", z1_mean), collapse = " "), "\n")
      c = c + 1
      if (c == 5){c = 1}
      plot_MCMC_diagnostics_series(att_post[1:t,], nu_post[1:t,], ((c-1)*5+1):((c)*5))
    }
    
  }
  
  print('done')
  
  out_list = list(att_post = att_post,
                  def_post = def_post,
                  home_post = home_post,
                  nu_post = nu_post,
                  nu_store = nu_store_mat,
                  nu_fix_mat = nu_fix_mat,
                  Z_post = Z_post,
                  p_post = p_post,
                  iterations = t,
                  priors = list(
                    att  = list(mean = att_mean_prior,  sd = att_sd_prior),
                    def  = list(mean = def_mean_prior,  sd = def_sd_prior),
                    home = list(mean = home_mean_prior, sd = home_sd_prior),
                    nu   = list(mean = nu_lmean_prior,  sd = nu_sd_prior),
                    p    = list(alpha= p_alpha_prior, beta = p_beta_prior)
                  ),
                  proposals = list(
                    att  = sd_prop_att,
                    def  = sd_prop_def,
                    home = sd_prop_home,
                    nu   = sd_prop_nu
                  ),
                  fix_idx = fix_idx,
                  X = list(X1 = X1, X2 = X2),
                  X_mid = X_mid,
                  team_names = rownames(X1),
                  distr_type = 'CMP-SAS',
                  constraint = 'STZ',
                  league = league_acro)
  class(out_list) = "MH_posterior"
  return(out_list)
}
