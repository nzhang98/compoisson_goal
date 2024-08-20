eval_log_compois_att = function(x, j, t, att_coef, def, att_intercept, def_intercept, home_coef = 0, home_fl = TRUE, nu_coef = 1){
  if (home_fl){
    mu = exp(att_intercept + att_coef + def_intercept + def[[t-1]][j] + home_coef)
  } else {
    mu = exp(att_intercept + att_coef + def_intercept + def[[t-1]][j])
  }
  return(nu_coef* (x*log(mu) - log(factorial(x)) ))
}

eval_likelihood_att = function(X, i, t, att_coef, def, att_intercept, def_intercept, home_coef = 0, home_fl = TRUE, nu_coef = 1){
  log_sum = 0
  N = dim(X)[1]
  for (j in 1:N){
    if (i == j){next}
    if (home_fl) {if (X_mid[i,j] == 0){next} 
      } else {if (X_mid[j,i] == 0) {next}}
    log_sum = log_sum + eval_log_compois_att(X[i,j], j, t, att_coef, def, att_intercept, def_intercept, home_coef, home_fl, nu_coef[i])
  }
  return(log_sum)
}

eval_log_compois_def = function(x, j, t, def_coef, att, att_intercept, def_intercept, home_coef = 0, home_fl = TRUE, nu_coef = 1){
  if (home_fl){
    mu = exp(att_intercept + att[[t]][j] + def_intercept + def_coef + home_coef)
  } else {
    mu = exp(att_intercept + att[[t]][j] + def_intercept + def_coef)
  }
  return(nu_coef*(x*log(mu) - log(factorial(x))))
}

eval_likelihood_def = function(X, i, t, def_coef, att, att_intercept, def_intercept, home_coef = 0, home_fl = TRUE, nu_coef = 1){
  log_sum = 0
  N = dim(X)[1]
  for (j in 1:N){
    if (i == j){next}
    if (home_fl) {if (X_mid[i,j] == 0){next} 
    } else {if (X_mid[j,i] == 0) {next}}
    log_sum = log_sum + eval_log_compois_def(X[i,j], j, t, def_coef, att, att_intercept, def_intercept, home_coef, home_fl, nu_coef[j])
  }
  return(log_sum)
}

eval_log_compois_home = function(x, i, j, t, att, def, att_intercept, def_intercept, home_coef, nu_coef){
  mu = exp(att_intercept + att[[t]][i] + def_intercept + def[[t]][j] + home_coef)
  return(nu_coef*(x*log(mu) - log(factorial(x))))
}

eval_likelihood_home = function(X, t, att, def, att_intercept, def_intercept, home_coef, nu_coef){
  log_sum = 0
  N = dim(X)[1]
  for (i in 1:N){
    for (j in 1:N){
      if (i == j){next}
      if (X_mid[i,j] == 0){next}
      log_sum = log_sum + eval_log_compois_home(X[i,j], i, j, t, att, def, att_intercept, def_intercept, home_coef, nu_coef[i])
    }
  }
  return(log_sum)
}

eval_log_compois_nu = function(x, i, j, t, att, def, att_intercept, def_intercept, home_coef = 0, nu_coef = 1, home_fl = TRUE){
  if (home_fl){
    mu = exp(att_intercept + att[[t]][i] + def_intercept + def[[t]][j] + home_coef)
  } else {
    mu = exp(att_intercept + att[[t]][i] + def_intercept + def[[t]][j])
  }
  return(nu_coef*(x*log(mu) - log(factorial(x))))
}

eval_likelihood_nu = function(X, i, t, att, def, att_intercept, def_intercept, home_coef = 0, nu = 1, home_fl = TRUE){
  log_sum = 0
  N = dim(X)[1]
  for (j in 1:N){
    if (i == j){next}
    if (home_fl) {if (X_mid[i,j] == 0){next} 
    } else {if (X_mid[j,i] == 0) {next}}
      log_sum = log_sum + eval_log_compois_nu(X[i,j], i, j, t, att, def, att_intercept, def_intercept, home_coef, nu_coef = nu[i], home_fl)
    }
  return(log_sum)
}

generate_auxiliary_compois_league_home = function(N, att_vector, def_vector, att_intercept, def_intercept, home_coef, nu, i = FALSE, update_def = FALSE){
  X1_sim = matrix(NA, nrow = N, ncol = N)
  if (is.numeric(i)){
    if (update_def){
      for (j in 1:N){
        if (i == j){next}
        if (X_mid[j,i] == 0){next}
        X1_sim[j,i] = rejection_sampler_draws(1, 
                                              mu = exp(att_vector[j] + def_vector[i] + att_intercept + def_intercept + home_coef),
                                              nu = nu[j])
      }
    } else {
      for (j in 1:N){
        if (i == j){next}
        if (X_mid[i,j] == 0){next} 
        X1_sim[i,j] = rejection_sampler_draws(1, 
                                              mu = exp(att_vector[i] + def_vector[j] + att_intercept + def_intercept + home_coef),
                                              nu = nu[i])
      }
    }
  } else {
    for (i in 1:N){
      for (j in 1:N){
        if (i == j){next}
        if (X_mid[i,j] == 0){next} 
        X1_sim[i,j] = rejection_sampler_draws(1, 
                                              mu = exp(att_vector[i] + def_vector[j] + att_intercept + def_intercept + home_coef),
                                              nu = nu[i])
      }
    }
  }
  return(X1_sim)
}

generate_auxiliary_compois_league_away = function(N, att_vector, def_vector, att_intercept, def_intercept, nu, i = FALSE, update_def = FALSE){
  X2_sim = matrix(NA, nrow = N, ncol = N)
  if (is.numeric(i)){
    if (update_def){
      for (j in 1:N){
        if (i == j){next}
        if (X_mid[i,j] == 0) {next}
        X2_sim[i,j] = rejection_sampler_draws(1, 
                                              mu = exp(att_vector[j] + def_vector[i] + att_intercept + def_intercept),
                                              nu = nu[j])
      }
    } else {
      for (j in 1:N){
        if (i == j){next}
        if (X_mid[j,i] == 0) {next}
        X2_sim[j,i] = rejection_sampler_draws(1, 
                                              mu = exp(att_vector[i] + def_vector[j] + att_intercept + def_intercept),
                                              nu = nu[i])
      }
    }
  } else {
    for (i in 1:N){
      for (j in 1:N){
        if (i == j){next}
        if (X_mid[j,i] == 0) {next}
        X2_sim[j,i] = rejection_sampler_draws(1, 
                                              mu = exp(att_vector[i] + def_vector[j] + att_intercept + def_intercept),
                                              nu = nu[i])
      }
    }
  }
  return(X2_sim)
}

MH_COMPois_fix = function(X1, X2, att_0, def_0, home_0, att_intercept, fix_att_idx, 
                          def_intercept, fix_def_idx, nu_0, 
                          X_mid = FALSE, iter=100, 
                          sd_prop_att = 1, att_mean_prior = 0, att_sd_prior = 1,
                          sd_prop_def = 1, def_mean_prior = 0, def_sd_prior = 1,
                          sd_prop_home = 1, home_mean_prior = 0, home_sd_prior = 1,
                          sd_prop_nu = 0.5, nu_a_prior = 1, nu_b_prior = 1){
  N = dim(X1)[1]
  att = vector("list", length = iter)
  def = vector("list", length = iter)
  home = numeric(iter)
  nu = vector("list", length = iter)
  
  if(typeof(X_mid) != 'double'){
    X_mid = matrix(1, N, N)
  }
  
  att_0[fix_att_idx] = 0
  def_0[fix_def_idx] = 0
  
  att[[1]] = att_0
  def[[1]] = def_0
  home[1] = home_0
  nu[[1]] = nu_0
  
  att_star = att_0
  def_star = def_0
  home_star = home_0
  nu_star = nu_0
  for (t in 2:iter){
    # if(t %% 100 == 0){print(t)}
    # print(t)
    
    #Update attack
    # print('attack')
    
    for (i in 1:N){
      
      if (i == fix_att_idx){
        att_star[i] = 0
      } else {
        att_prop = att_star
        att_prop[i] = rnorm(1, att_star[i], sd = sd_prop_att)
        
        aux_X1 = generate_auxiliary_compois_league_home(N, att_prop, def_star, att_intercept, def_intercept, home_star, nu_star, i = i)
        aux_X2 = generate_auxiliary_compois_league_away(N, att_prop, def_star, att_intercept, def_intercept, nu_star, i = i)
        
        log_lik_prop = eval_likelihood_att(X1, i, t, att_prop[i], def, att_intercept, def_intercept, home_coef = home_star, nu_coef = nu_star)
        log_lik_prop = eval_likelihood_att(t(X2), i, t, att_prop[i], def, att_intercept, def_intercept, home_fl = FALSE, nu_coef = nu_star) + log_lik_prop
        log_prior_prop = dnorm(att_prop[i], att_mean_prior, att_sd_prior, log = TRUE)
        log_lik_prev_aux = eval_likelihood_att(aux_X1, i, t, att_star[i], def, att_intercept, def_intercept, home_coef = home_star, nu_coef = nu_star)
        log_lik_prev_aux = eval_likelihood_att(t(aux_X2), i, t, att_star[i], def, att_intercept, def_intercept, home_fl = FALSE, nu_coef = nu_star) + log_lik_prev_aux
        numerator = log_lik_prop + log_prior_prop + log_lik_prev_aux
        
        log_lik_prev = eval_likelihood_att(X1, i, t, att_star[i], def, att_intercept, def_intercept, home_coef = home_star, nu_coef = nu_star)
        log_lik_prev = eval_likelihood_att(t(X2), i, t, att_star[i], def, att_intercept, def_intercept, home_fl = FALSE, nu_coef = nu_star) + log_lik_prev
        log_prior_prev = dnorm(att_star[i], att_mean_prior, att_sd_prior, log = TRUE)
        log_lik_prop_aux = eval_likelihood_att(aux_X1, i, t, att_prop[i], def, att_intercept, def_intercept, home_coef = home_star, nu_coef = nu_star)
        log_lik_prop_aux = eval_likelihood_att(t(aux_X2), i, t, att_prop[i], def, att_intercept, def_intercept, home_fl = FALSE, nu_coef = nu_star) + log_lik_prop_aux
        denominator = log_lik_prev + log_prior_prev + log_lik_prop_aux
        
        alpha = min(1, exp(numerator - denominator))
        
        u = runif(1)
        if (u <= alpha){
          att_star[i] = att_prop[i]
        }
      }
    }
    
    att[[t]] = att_star
    # browser()
    #Update defense
    # print('defense')
    
    for (i in 1:N){
      if (i == fix_def_idx){
        def_star[i] = 0
      } else {
        def_prop = def_star
        def_prop[i] = rnorm(1, def_star[i], sd = sd_prop_def)
        
        aux_X1 = generate_auxiliary_compois_league_home(N, att_star, def_prop, att_intercept, def_intercept, home_star, nu_star, i = i, update_def = TRUE)
        aux_X2 = generate_auxiliary_compois_league_away(N, att_star, def_prop, att_intercept, def_intercept, nu_star, i = i, update_def = TRUE)
        
        log_lik_prop = eval_likelihood_def(X2, i, t, def_prop[i], att, att_intercept, def_intercept, home_coef = home_star, nu_coef = nu_star)
        log_lik_prop = eval_likelihood_def(t(X1), i, t, def_prop[i], att, att_intercept, def_intercept, home_fl = FALSE, nu_coef = nu_star) + log_lik_prop
        log_prior_prop = dnorm(def_prop[i], def_mean_prior, def_sd_prior, log = TRUE)
        log_lik_prev_aux = eval_likelihood_def(aux_X2, i, t, def_star[i], att, att_intercept, def_intercept, home_coef = home_star, nu_coef = nu_star)
        log_lik_prev_aux = eval_likelihood_def(t(aux_X1), i, t, def_star[i], att, att_intercept, def_intercept, home_fl = FALSE, nu_coef = nu_star) + log_lik_prev_aux
        numerator = log_lik_prop + log_prior_prop + log_lik_prev_aux
        
        log_lik_prev = eval_likelihood_def(X2, i, t, def_star[i], att, att_intercept, def_intercept, home_coef = home_star, nu_coef = nu_star)
        log_lik_prev = eval_likelihood_def(t(X1), i, t, def_star[i], att, att_intercept, def_intercept, home_fl = FALSE, nu_coef = nu_star) + log_lik_prev
        log_prior_prev = dnorm(def_star[i], def_mean_prior, def_sd_prior, log = TRUE)
        log_lik_prop_aux = eval_likelihood_def(aux_X2, i, t, def_prop[i], att, att_intercept, def_intercept, home_coef = home_star, nu_coef = nu_star)
        log_lik_prop_aux = eval_likelihood_def(t(aux_X1), i, t, def_prop[i], att, att_intercept, def_intercept, home_fl = FALSE, nu_coef = nu_star) + log_lik_prop_aux
        denominator = log_lik_prev + log_prior_prev + log_lik_prop_aux
        
        alpha = min(1, exp(numerator - denominator))
        u = runif(1)
        if (u <= alpha){
          def_star[i] = def_prop[i]
        }
      }
    }
    def[[t]] = def_star
    # browser()
    
    #Update home
    # print('home')
    home_prop = rnorm(1, home_star, sd = sd_prop_home)
    # print(home_prop)
    
    aux_X1 = generate_auxiliary_compois_league_home(N, att_star, def_star, att_intercept, def_intercept, home_prop, nu_star)
    
    log_lik_prop = eval_likelihood_home(X1, t, att, def, att_intercept, def_intercept, home_prop, nu_star)
    log_prior_prop = dnorm(home_prop, home_mean_prior, home_sd_prior, log = TRUE)
    log_lik_prev_aux = eval_likelihood_home(aux_X1, t, att, def, att_intercept, def_intercept, home_star, nu_star)
    numerator = log_lik_prop + log_prior_prop + log_lik_prev_aux
    
    log_lik_prev = eval_likelihood_home(X1, t, att, def, att_intercept, def_intercept, home_star, nu_star)
    log_prior_prev = dnorm(home_star, home_mean_prior, home_sd_prior, log = TRUE)
    log_lik_prop_aux = eval_likelihood_home(aux_X1, t, att, def, att_intercept, def_intercept, home_prop, nu_star)
    denominator = log_lik_prev + log_prior_prev + log_lik_prop_aux
    
    alpha = min(1, exp(numerator - denominator))
    u = runif(1)
    if (u <= alpha){
      home_star = home_prop
    }
    home[t] = home_star
    
    #Update nu home
    
    for (i in 1:N){
      
      nu_prop = nu_star
      nu_prop[i] = rnorm(1, nu_prop[i], sd = sd_prop_nu)
      
      if (nu_prop[i] < 0){
        next
      } else {
        aux_X1 = generate_auxiliary_compois_league_home(N, att_star, def_star, att_intercept, def_intercept, home_star, nu_prop, i = i)
        aux_X2 = generate_auxiliary_compois_league_away(N, att_star, def_star, att_intercept, def_intercept, nu_prop, i = i)
        
        log_lik_prop = eval_likelihood_nu(X1, i, t, att, def, att_intercept, def_intercept, home_coef = home_star, nu = nu_prop)
        log_lik_prop = eval_likelihood_nu(t(X2), i, t, att, def, att_intercept, def_intercept, home_fl = FALSE, nu = nu_prop) + log_lik_prop
        log_prior_prop = dnorm(att_prop[i], att_mean_prior, att_sd_prior, log = TRUE)
        log_lik_prev_aux = eval_likelihood_nu(aux_X1, i, t, att, def, att_intercept, def_intercept, home_coef = home_star, nu = nu_star)
        log_lik_prev_aux = eval_likelihood_nu(t(aux_X2), i, t, att, def, att_intercept, def_intercept, home_fl = FALSE, nu = nu_star) + log_lik_prev_aux
        numerator = log_lik_prop + log_prior_prop + log_lik_prev_aux
        
        log_lik_prev = eval_likelihood_nu(X1, i, t, att, def, att_intercept, def_intercept, home_coef = home_star, nu = nu_star)
        log_lik_prev = eval_likelihood_nu(t(X2), i, t, att, def, att_intercept, def_intercept, home_fl = FALSE, nu = nu_star) + log_lik_prev
        log_prior_prev = dnorm(att_star[i], att_mean_prior, att_sd_prior, log = TRUE)
        log_lik_prop_aux = eval_likelihood_nu(aux_X1, i, t, att, def, att_intercept, def_intercept, home_coef = home_star, nu = nu_prop)
        log_lik_prop_aux = eval_likelihood_nu(t(aux_X2), i, t, att, def, att_intercept, def_intercept, home_fl = FALSE, nu = nu_prop) + log_lik_prop_aux
        denominator = log_lik_prev + log_prior_prev + log_lik_prop_aux
        
        alpha = min(1, exp(numerator - denominator))
        
        u = runif(1)
        if (u <= alpha){
          nu_star[i] = nu_prop[i]
        }
      }
    }
    nu[[t]] = nu_star
    
    if(t %% 1000 == 0){cat("t:", t, "home:", home_star, "\nnu:", nu_star,"\nATT:",att_star, "\nDEF:",def_star ,"\n")}
    
  }
  out_list = list(att_post = att,
                  def_post = def,
                  home_post = home,
                  nu_post = nu,
                  iterations = t,
                  att_fix_idx = att_fix_idx,
                  att_fix = att_fix,
                  def_fix_idx = def_fix_idx,
                  def_fix = def_fix,
                  X_mid = X_mid,
                  team_names = rownames(X1),
                  distr_type = 'CP_D',
                  league = league_acro)
  class(out_list) = "MH_mid"
  return(out_list)
}
