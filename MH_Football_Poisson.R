# #Simple Metropolis Hastings with 1 fixed attack parameter to solve identification problem
# eval_log_pois_att = function(x, j, t, att_coef, def, att_intercept, home_coef = 0, home_fl = TRUE){
#   if (home_fl){
#     return(log(dpois(x, lambda = exp(att_intercept + att_coef + def[[t-1]][j] + home_coef))))
#   } else {
#     return(log(dpois(x, lambda = exp(att_intercept + att_coef + def[[t-1]][j]))))
#   }
# }
# 
# eval_likelihood_att = function(X, i, t, att_coef, def, att_intercept, home_coef = 0, home_fl = TRUE){
#   log_sum = 0
#   N = dim(X)[1]
#   for (j in 1:N){
#     if (i == j){next}
#     log_sum = log_sum + eval_log_pois_att(X[i,j], j, t, att_coef, def, att_intercept, home_coef, home_fl)
#   }
#   return(log_sum)
# }
# 
# eval_log_pois_def = function(x, j, t, def_coef, att, att_intercept, home_coef = 0, home_fl = TRUE){
#   if (home_fl){
#     return(log(dpois(x, lambda = exp(att_intercept + att[[t]][j] + def_coef + home_coef))))
#   } else {
#     return(log(dpois(x, lambda = exp(att_intercept + att[[t]][j] + def_coef))))
#   }
# }
# 
# eval_likelihood_def = function(X, i, t, def_coef, att, att_intercept, home_coef = 0, home_fl = TRUE){
#   log_sum = 0
#   N = dim(X)[1]
#   for (j in 1:N){
#     if (i == j){next}
#     log_sum = log_sum + eval_log_pois_def(X[i,j], j, t, def_coef, att, att_intercept, home_coef, home_fl)
#   }
#   return(log_sum)
# }
# 
# eval_log_pois_home = function(x, i, j, t, att, def, att_intercept, home_coef){
#   return(log(dpois(x, lambda = exp(att_intercept + att[[t]][i] + def[[t]][j] + home_coef))))
# }
# 
# eval_likelihood_home = function(X, t, att, def, att_intercept, home_coef){
#   log_sum = 0
#   N = dim(X)[1]
#   for (i in 1:N){
#     for (j in 1:N){
#       if (i == j){next}
#       log_sum = log_sum + eval_log_pois_home(X[i,j], i, j, t, att, def, att_intercept, home_coef)
#     }
#   }
#   return(log_sum)
# }
# 
# MH_fix = function(X1, X2, att_0, def_0, home_0, att_intercept, fix_att_idx, iter=100, sd_prop = 1){
#   N = dim(X1)[1]
#   att = vector("list", length = iter)
#   def = vector("list", length = iter)
#   home = rep(NULL, iter)
#   
#   att_0[fix_att_idx] = 0
#   
#   att[[1]] = att_0
#   def[[1]] = def_0
#   home[1] = home_0
#   
#   att_star = att_0
#   def_star = def_0
#   home_star = home_0
#   for (t in 2:iter){
#     att_prop = rep(NULL, N)
#     for (i in 1:N){
#       
#       if (i == fix_att_idx){
#         att_star[i] = 0
#       } else {
#         att_prop[i] = rnorm(1, att_star[i], sd = sd_prop)
#         
#         log_lik_prop = eval_likelihood_att(X1, i, t, att_prop[i], def, att_intercept, home_coef = home_star)
#         log_lik_prop = eval_likelihood_att(t(X2), i, t, att_prop[i], def, att_intercept, home_fl = FALSE) + log_lik_prop
#         
#         log_lik_prev = eval_likelihood_att(X1, i, t, att_star[i], def, att_intercept, home_coef = home_star)
#         log_lik_prev = eval_likelihood_att(t(X2), i, t, att_star[i], def, att_intercept, home_fl = FALSE) + log_lik_prev
#         
#         
#         alpha = min(1, exp(log_lik_prop - log_lik_prev))
#         
#         u = runif(1)
#         if (u <= alpha){
#           att_star[i] = att_prop[i]
#         }
#       }
#     }
#     att[[t]] = att_star
#     
#     def_prop = rep(NULL, N)
#     for (i in 1:N){
#       def_prop[i] = rnorm(1, def_star[i], sd = sd_prop)
#       
#       log_lik_prop = eval_likelihood_def(X2, i, t, def_prop[i], att, att_intercept, home_coef = home_star)
#       log_lik_prop = eval_likelihood_def(t(X1), i, t, def_prop[i], att, att_intercept, home_fl = FALSE) + log_lik_prop
#       
#       log_lik_prev = eval_likelihood_def(X2, i, t, def_star[i], att, att_intercept, home_coef = home_star)
#       log_lik_prev = eval_likelihood_def(t(X1), i, t, def_star[i], att, att_intercept, home_fl = FALSE) + log_lik_prev
#       
#       alpha = min(1, exp(log_lik_prop - log_lik_prev))
#       u = runif(1)
#       if (u <= alpha){
#         def_star[i] = def_prop[i]
#       }
#     }
#     def[[t]] = def_star
#     
#     home_prop = rnorm(1, home_star, sd = sd_prop)
#     
#     log_lik_prop = eval_likelihood_home(X1, t, att, def, att_intercept, home_prop)
#     
#     log_lik_prev = eval_likelihood_home(X1, t, att, def, att_intercept, home_star)
#     
#     alpha = min(1, exp(log_lik_prop - log_lik_prev))
#     u = runif(1)
#     if (u <= alpha){
#       home_star = home_prop
#     }
#     home[t] = home_star
#   }
#   out_list = list(att_post = att,
#                   def_post = def,
#                   home_post = home)
#   class(out_list) = "MH"
#   return(out_list)
# }


#Simple Metropolis Hastings with 1 fixed attack and defense parameter to solve identification problem
#Problem to fix: we are not evaluating the prior on the acceptance rate!!
eval_log_pois_att = function(x, j, t, att_coef, def, att_intercept, def_intercept, home_coef = 0, home_fl = TRUE){
  if (home_fl){
    return(log(dpois(x, lambda = exp(att_intercept + att_coef + def_intercept + def[[t-1]][j] + home_coef))))
  } else {
    return(log(dpois(x, lambda = exp(att_intercept + att_coef + def_intercept + def[[t-1]][j]))))
  }
}

eval_likelihood_att = function(X, i, t, att_coef, def, att_intercept, def_intercept, home_coef = 0, home_fl = TRUE){
  log_sum = 0
  N = dim(X)[1]
  for (j in 1:N){
    if (i == j){next}
    if (home_fl) {if (X_mid[i,j] == 0){next} 
    } else {if (X_mid[j,i] == 0) {next}}
    log_sum = log_sum + eval_log_pois_att(X[i,j], j, t, att_coef, def, att_intercept, def_intercept, home_coef, home_fl)
  }
  return(log_sum)
}

eval_log_pois_def = function(x, j, t, def_coef, att, att_intercept, def_intercept, home_coef = 0, home_fl = TRUE){
  if (home_fl){
    return(log(dpois(x, lambda = exp(att_intercept + att[[t]][j] + def_intercept + def_coef + home_coef))))
  } else {
    return(log(dpois(x, lambda = exp(att_intercept + att[[t]][j] + def_intercept + def_coef))))
  }
}

eval_likelihood_def = function(X, i, t, def_coef, att, att_intercept, def_intercept, home_coef = 0, home_fl = TRUE){
  log_sum = 0
  N = dim(X)[1]
  for (j in 1:N){
    if (i == j){next}
    if (home_fl) {if (X_mid[i,j] == 0){next} 
    } else {if (X_mid[j,i] == 0) {next}}
    log_sum = log_sum + eval_log_pois_def(X[i,j], j, t, def_coef, att, att_intercept, def_intercept, home_coef, home_fl)
  }
  return(log_sum)
}

eval_log_pois_home = function(x, i, j, t, att, def, att_intercept, def_intercept, home_coef){
  return(log(dpois(x, lambda = exp(att_intercept + att[[t]][i] + def_intercept + def[[t]][j] + home_coef))))
}

eval_likelihood_home = function(X, t, att, def, att_intercept, def_intercept, home_coef){
  log_sum = 0
  N = dim(X)[1]
  for (i in 1:N){
    for (j in 1:N){
      if (i == j){next}
      if (X_mid[i,j] == 0){next}
      log_sum = log_sum + eval_log_pois_home(X[i,j], i, j, t, att, def, att_intercept, def_intercept, home_coef)
    }
  }
  return(log_sum)
}

MH_fix = function(X1, X2, att_0, def_0, home_0, att_intercept, fix_att_idx, 
                  def_intercept, fix_def_idx, iter=100, X_mid = FALSE, 
                  sd_prop_att = 1, sd_prop_def = 1, sd_prop_home = 1,
                  att_mean_prior = 0, att_sd_prior = 1,
                  def_mean_prior = 0, def_sd_prior = 1,
                  home_mean_prior = 0, home_sd_prior = 1){
  N = dim(X1)[1]
  att = vector("list", length = iter)
  def = vector("list", length = iter)
  home = rep(NULL, iter)
  
  if(typeof(X_mid) != 'double'){
    X_mid = matrix(1, N, N)
  }
  
  att_0[fix_att_idx] = 0
  def_0[fix_def_idx] = 0
  
  att[[1]] = att_0
  def[[1]] = def_0
  home[1] = home_0
  
  att_star = att_0
  def_star = def_0
  home_star = home_0
  for (t in 2:iter){
    if (t %% 1000 == 0){
      cat("t:", t, "home:", home_star)
    }
    # print(t)
    att_prop = rep(NULL, N)
    for (i in 1:N){
      
      if (i == fix_att_idx){
        att_star[i] = 0
      } else {
        att_prop[i] = rnorm(1, att_star[i], sd = sd_prop_att)
        
        log_lik_prop = eval_likelihood_att(X1, i, t, att_prop[i], def, att_intercept, def_intercept, home_coef = home_star)
        log_lik_prop = eval_likelihood_att(t(X2), i, t, att_prop[i], def, att_intercept, def_intercept, home_fl = FALSE) + log_lik_prop
        
        log_prior_prop = dnorm(att_prop[i], att_mean_prior, att_sd_prior, log = TRUE)
        
        log_lik_prev = eval_likelihood_att(X1, i, t, att_star[i], def, att_intercept, def_intercept, home_coef = home_star)
        log_lik_prev = eval_likelihood_att(t(X2), i, t, att_star[i], def, att_intercept, def_intercept, home_fl = FALSE) + log_lik_prev
        
        log_prior_prev = dnorm(att_star[i], att_mean_prior, att_sd_prior, log = TRUE)
        
        
        alpha = min(1, exp(log_lik_prop - log_lik_prev + log_prior_prop - log_prior_prev))
        
        u = runif(1)
        if (u <= alpha){
          att_star[i] = att_prop[i]
        }
      }
    }
    att[[t]] = att_star
    
    def_prop = rep(NULL, N)
    for (i in 1:N){
      if (i == fix_def_idx){
        def_star[i] = 0
      } else {
        def_prop[i] = rnorm(1, def_star[i], sd = sd_prop_def)
        
        log_lik_prop = eval_likelihood_def(X2, i, t, def_prop[i], att, att_intercept, def_intercept, home_coef = home_star)
        log_lik_prop = eval_likelihood_def(t(X1), i, t, def_prop[i], att, att_intercept, def_intercept, home_fl = FALSE) + log_lik_prop
        
        log_prior_prop = dnorm(def_prop[i], def_mean_prior, def_sd_prior, log = TRUE)
        
        log_lik_prev = eval_likelihood_def(X2, i, t, def_star[i], att, att_intercept, def_intercept, home_coef = home_star)
        log_lik_prev = eval_likelihood_def(t(X1), i, t, def_star[i], att, att_intercept, def_intercept, home_fl = FALSE) + log_lik_prev
        
        log_prior_prev = dnorm(def_star[i], def_mean_prior, def_sd_prior, log = TRUE)
        
        alpha = min(1, exp(log_lik_prop - log_lik_prev + log_prior_prop - log_prior_prev))
        u = runif(1)
        if (u <= alpha){
          def_star[i] = def_prop[i]
        }
      }
    }
    def[[t]] = def_star
    
    home_prop = rnorm(1, home_star, sd = sd_prop_home)
    
    log_lik_prop = eval_likelihood_home(X1, t, att, def, att_intercept, def_intercept, home_prop)
    
    log_prior_prop = dnorm(home_prop, home_mean_prior, home_sd_prior, log = TRUE)
    
    log_lik_prev = eval_likelihood_home(X1, t, att, def, att_intercept, def_intercept, home_star)
    
    log_prior_prev = dnorm(home_star, home_mean_prior, home_sd_prior, log = TRUE)
    
    alpha = min(1, exp(log_lik_prop - log_lik_prev + log_prior_prop - log_prior_prev))
    u = runif(1)
    if (u <= alpha){
      home_star = home_prop
    }
    home[t] = home_star
    if(t %% 1000 == 0){cat("t:", t, "home:", home_star,"\nATT:",att_star, "\nDEF:",def_star ,"\n")}
  }
  out_list = list(att_post = att,
                  def_post = def,
                  home_post = home,
                  iterations = t,
                  att_fix_idx = att_fix_idx,
                  att_fix = att_fix,
                  def_fix_idx = def_fix_idx,
                  def_fix = def_fix,
                  X_mid = X_mid,
                  team_names = rownames(X1),
                  distr_type = 'P',
                  league = league_acro)
  class(out_list) = "MH"
  return(out_list)
}

