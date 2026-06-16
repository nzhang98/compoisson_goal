### Generic functions
library(ggplot2)
library(reshape2)
library(purrr)
library(gridExtra)
library(dplyr)
library(tidyr)
library(Rcpp)

Mode = function(x) { #compute the mode of a discrete list x
  ux = unique(x)
  ux[which.max(tabulate(match(x, ux)))]
}

### Read data functions

read_data = function(season, league = 'Premier'){
  # League list: "Premier", "SerieA", "Liga", "Ligue1", "Bundes" 
  Results=read.csv(paste0("Data//Results//Result_",league,"_",season, ".csv"), header = TRUE, sep=",")
  
  rownames(Results)=Results$Home...Away
  Results=Results[,-1]
  
  library(stringi)
  N = nrow(Results)
  
  X1=matrix(NA,N,N)
  X2=matrix(NA,N,N)
  for(i in 1:N){
    for(j in 1:N){
      if(i==j || Results[i,j]=="" || is.na(Results[i,j])){
        X1[i,j]=NA
      }else{
        b=as.character(Results[i,j])
        a=stri_split_fixed(b,"~")
        X1[i,j] = as.numeric(a[[1]][1])
        X2[i,j] = as.numeric(a[[1]][2])
      }
    }
  }
  rownames(X1)=rownames(Results)
  colnames(X1)=colnames(Results)
  rownames(X2)=rownames(Results)
  colnames(X2)=colnames(Results)
  
  return(list(X1,X2))
}

generate_season_string = function(start_year, end_year) {
  if (end_year <= start_year) {
    stop("end_year must be greater than start_year")
  }
  
  years = start_year:(end_year - 1)
  next_years = years + 1
  
  format_two_digits = function(x) sprintf("%02d", x %% 100)
  
  formatted = paste0(format_two_digits(years), format_two_digits(next_years))
  return(formatted)
}

read_historics = function(season, league_acro = 'PL'){
  # league acros: "PL", "SA", "LL", "LC", "BL"
  results = read.csv(paste0("Data//Historics//",league_acro,"_History_",season, ".csv"), header = TRUE, sep=",")
  return(results)
}

get_team_code = function(team_name){
  return(team_codes[which(team_names == team_name)])
}

get_sample_ids = function(n_samples, t, from_t = FALSE, replace_flag = FALSE){
  
  if(typeof(from_t) != 'double'){from_t = t/2}
  
  id = sample(seq(from=from_t, to=t),n_samples, replace = replace_flag)
  return(id)
}

row_to_matrix = function(row, goals = goal_vals) {
  n = length(goals)
  mat = matrix(NA, nrow = n, ncol = n,
                dimnames = list(goals, goals))
  names(row) = scores
  
  for (score in names(row)) {
    parts = str_split(score, ":")[[1]]
    home = as.integer(parts[1])
    away = as.integer(parts[2])
    
    # Convert 0-based to 1-based indexing
    row_idx = home + 1
    col_idx = away + 1
    
    mat[row_idx, col_idx] = as.numeric(row[[score]])
  }
  return(mat)
}

read_odds_matrix = function(df_odds){
  # df_odds = read.csv(file = paste0("Data/Historics/",league_acro,"_Odds_",season,".csv")
  #                    , check.names = FALSE
  # )
  # Select only columns ending with "_odd"
  filtered_odds = df_odds %>%
    select(matches("\\d+:\\d+_odd"))
  
  # filtered_odds = final_df %>%
  # select(matches("X\\d+\\.\\d+_odd"))
  
  
  # Extract score labels
  scores = colnames(filtered_odds) %>%
    str_remove("_odd") %>%
    # str_remove("X") %>%
    unique()
  
  # Get all unique goal values to define matrix size
  goal_vals = sort(unique(as.integer(unlist(str_split(scores, ":", simplify = TRUE)))))
  # goal_vals = sort(unique(as.integer(unlist(str_split(scores, ".", simplify = TRUE)))))
  
  odds_matrices = lapply(1:nrow(filtered_odds), function(i) row_to_matrix(filtered_odds[i, ]))
  return(odds_matrices)
}

merge_hist_odds_df = function(df_odds, df_hist){
  df1 = df_odds
  df1 = df1 %>%
    rename(home_team = home_team, away_team = away_team) %>%
    mutate(date = dmy(date))
  
  df2 = df_hist
  df2 = df2 %>%
    rename(home_team = HomeTeam, away_team = AwayTeam, time = Time, date = Date) %>%
    mutate(date = dmy(date))
  
  # --- Step 1: Create a match identifier for both dataframes ---
  df1 = df1 %>%
    mutate(match_key = paste0(home_team, "_", away_team, "_", date))
  
  df2 = df2 %>%
    mutate(match_key_raw = paste0(home_team, "_", away_team, "_", date))
  
  # --- Step 2: Fuzzy match df1$match_key to df2$match_key_raw ---
  # Define a custom match function
  match_games = function(row_key, all_keys) {
    distances = stringdist::stringdist(row_key, all_keys, method = "jw")
    best_match_index = which.min(distances)
    return(best_match_index)
  }
  
  # Apply matching
  match_indices = map_int(df1$match_key, ~match_games(.x, df2$match_key_raw))
  
  # Create matched dataframe
  df1_matched = df1 %>%
    mutate(matched_index = match_indices)
  
  df2_matched = df2 %>%
    mutate(row_id = row_number())
  
  # --- Step 3: Join based on matched indices ---
  final_df = df2_matched %>%
    left_join(df1_matched, by = c("row_id" = "matched_index"))
  
  return(final_df)
}

load_MH = function(season, distr, n_games = FALSE, league_acro = 'PL'){
  if (is.numeric(n_games)){
    load(file = paste0("Data//MH_Results//Mid//",league_acro,"_",as.character(season),
                       "_n",as.character(n_games),"_",as.character(distr),".RData"),
         envir = .GlobalEnv)
  }
  else {
    load(file = paste0("Data//MH_Results//",league_acro,"_",as.character(season),
                       "_",as.character(distr),".RData"),
         envir = .GlobalEnv)
  }
}

# Rejection Sampler for multiple draws from COM-Poisson as in Benson(2021)

cppFunction('
NumericVector fast_runif(int n) {
  NumericVector result(n);
  for(int i = 0; i < n; ++i) {
    result[i] = R::runif(0.0, 1.0);  // Generate random numbers between 0 and 1
  }
  return result;
}')

geom_envelope_draws = function(n, mu, nu, return_ndraws = FALSE) {
  samples = numeric(n)
  n_draws = numeric(n)
  log_Bfgs = numeric(n)
  
  p = (2*nu)/(2*mu*nu + 1 + nu)
  
  xm = floor(mu/((1-p)**(1/nu)))
  
  log_B = -log(p) + nu*xm*log(mu) - xm*log(1-p) - nu*lfactorial(xm)
  
  for(i in 1:n){
    c = 0
    while (TRUE){
      c = c + 1
      
      u_0 = fast_runif(1)
      x = floor(log(u_0)/log(1-p))
      
      log_alpha = nu*( x*log(mu) - lfactorial(x) ) - log_B - x*log(1-p) - log(p)
      
      u = fast_runif(1)
      
      if (u <= exp(log_alpha)){
        samples[i] = x
        n_draws[i] = c
        log_Bfgs[i] = log_B
        break
      }
    }
  }
  if (return_ndraws) {return(list("samples" = samples, "n_draws" = n_draws, "log_Bfg" = log_Bfgs))} 
  else {return(samples)}
}

pois_envelope_draws = function(n, mu, nu, return_ndraws = FALSE){
  samples = numeric(n)
  n_draws = numeric(n)
  log_Bfgs = numeric(n)
  
  log_B = (nu-1)*( floor(mu)*log(mu) - lfactorial(floor(mu)) )
  
  for(i in 1:n){
    c = 0
    while(TRUE){
      c = c + 1
      
      x = rpois(1, mu)
      
      log_alpha = nu*(x*log(mu) - lfactorial(x)) - log_B - x*log(mu) + lfactorial(x)
      
      u = fast_runif(1)
      
      if (u <= exp(log_alpha)){
        samples[i] = x
        n_draws[i] = c
        log_Bfgs[i] = log_B
        break
      }
    }
  }
  
  if (return_ndraws) {return(list("samples" = samples, "n_draws" = n_draws, "log_Bfg" = log_Bfgs))} 
  else {return(samples)}
}

rejection_sampler_draws = function(n, mu, nu, return_ndraws = FALSE){
  if (nu < 1){
    x = geom_envelope_draws(n, mu, nu, return_ndraws)
  } else { #nu >= 1
    x = pois_envelope_draws(n, mu, nu, return_ndraws)
  }
  return(x)
}

rejection_sampler_draws_vec = function(n, mu, nu, return_ndraws = FALSE){
  x = mapply(function(n_val, mu_val, nu_val, return_ndraws_val) {
    if (nu_val < 1) {
      return(geom_envelope_draws(n, mu, nu_val, return_ndraws))
    } else {
      return(pois_envelope_draws(n, mu, nu_val, return_ndraws))
    }
  }, n, mu, nu, return_ndraws, SIMPLIFY = FALSE)
  return(x)
}
### End sampler for COM-Poisson

### Functions to plot and summarize results

retrieve_acceptance_rates = function(MH_obj, param = 'all', burnin = 'default'){
  N = length(MH_obj$att_post[[1]])
  if (burnin == 'default'){
    t_burn = MH_obj$iterations / 2
  } else {
    t_burn = burnin
  }
  
  iters = MH_obj$iterations - t_burn
  
  # print(iters)
  
  if (param == 'all'){
    output = c() #ADAPT THIS
    att = retrieve_post_longlists(MH_obj, 'att', burnin = burnin)
    def = retrieve_post_longlists(MH_obj, 'def', burnin = burnin)
    if(MH_obj$distr_type == 'CP_D'){print(1)
      nu = retrieve_post_longlists(MH_obj, 'nu', burnin = burnin)}
    
    
    for (i in 1:N){
      temp = att[[i]]
      c = 0
      for (t in 2:iters){
        if (temp[t] != temp[t-1]){c = c+1}
      }
      output = append(output, setNames(c/iters, paste0('att',i)))
    }
    for (i in 1:N){
      temp = def[[i]]
      c = 0
      for (t in 2:iters){
        if (temp[t] != temp[t-1]){c = c+1}
      }
      output = append(output, setNames(c/iters, paste0('def',i)))
    }
    
    temp = MH_obj$home_post[t_burn:MH_obj$iterations]
    # print('home')
    # print(length(temp))
    # print(t_burn)
    # print(MH_obj$iterations)
    c = 0
    for (t in 2:iters){
      if (temp[t] != temp[t-1]){c = c+1}
    }
    output = append(output, c('home' = (c/iters)))
    
    if(MH_obj$distr_type == 'CP_D'){
      for (i in 1:N){
        temp = nu[[i]]
        # print(temp)
        # print(length(temp))
        c = 0
        for (t in 2:iters){
          if (temp[t] != temp[t-1]){c = c+1}
        }
        output = append(output, setNames(c/iters, paste0('nu',i)))
      }
    } else if (MH_obj$distr_type == 'CP'){
      temp = MH_obj$nu_home[t_burn:MH_obj$iterations]
      c = 0
      for (t in 2:iters){
        if (temp[t] != temp[t-1]){c = c+1}
      }
      output = append(output, c('nu_home' = (c/iters)))
      
      temp = MH_obj$nu_away[t_burn:MH_obj$iterations]
      c = 0
      for (t in 2:iters){
        if (temp[t] != temp[t-1]){c = c+1}
      }
      output = append(output, c('nu_away' = (c/iters)))
    }
  }
  
  return(output)
}

retrieve_post_longlists = function(MH_obj, param, burnin = 'default'){
  if (param == 'att'){
    post = MH_obj$att_post
  }
  if (param == 'def'){
    post = MH_obj$def_post
  }
  if (param == 'nu'){
    post = MH_obj$nu
  }
  # if (param == 'all'){
  #   
  # }
  
  if (burnin == 'default'){
    post = post[(MH_obj$iterations/2):MH_obj$iterations]
  } else {
    stopifnot(is.numeric(burnin))
    post = post[burnin:MH_obj$iterations]
  }
  return(lapply(transpose(post), unlist))
}

retrieve_summary = function(MH_obj, burnin = 'default'){
  N = length(MH_obj$att_post[[1]])
  if (burnin == 'default'){
    t_burn = MH_obj$iterations / 2
  } else {
    t_burn = burnin
  }
  
  att_chains = retrieve_post_longlists(MH_obj, 'att', burnin)
  def_chains = retrieve_post_longlists(MH_obj, 'def', burnin)
  home_chain = MH_obj$home_post[t_burn:MH_obj$iterations]
  
  home = c(mean = mean(home_chain), sd = sd(home_chain),  quantile(home_chain, probs = c(0.025, 0.25, 0.5, 0.75, 0.975)))
  
  if (MH_obj$distr_type == 'CP'){
    nu_home_chain = MH_obj$nu_home_post[t_burn:MH_obj$iterations]
    nu_away_chain = MH_obj$nu_away_post[t_burn:MH_obj$iterations]
    nu_home = c(mean = mean(nu_home_chain), sd = sd(nu_home_chain),  quantile(nu_home_chain, probs = c(0.025, 0.25, 0.5, 0.75, 0.975)))
    nu_away = c(mean = mean(nu_away_chain), sd = sd(nu_away_chain),  quantile(nu_away_chain, probs = c(0.025, 0.25, 0.5, 0.75, 0.975)))
    out_df = rbind(home, nu_home, nu_away)
  } else {
    out_df = home
  }
  
  
  for (i in 1:N){
    att_i = att_chains[[i]]
    att_vec = c(mean = mean(att_i), sd = sd(att_i),  quantile(att_i, probs = c(0.025, 0.25, 0.5, 0.75, 0.975)))
    out_df = rbind(out_df, att_vec)
    rownames(out_df)[nrow(out_df)] = paste0('att_',i)
  }
  
  for (i in 1:N){
    def_i = def_chains[[i]]
    def_vec = c(mean = mean(def_i), sd = sd(def_i),  quantile(def_i, probs = c(0.025, 0.25, 0.5, 0.75, 0.975)))
    out_df = rbind(out_df, def_vec)
    rownames(out_df)[nrow(out_df)] = paste0('def_',i)
  }
  
  if (MH_obj$distr_type == 'CP_D'){
    nu_chains = retrieve_post_longlists(MH_obj, 'nu', burnin)
    
    for (i in 1:N){
      nu_i = nu_chains[[i]]
      nu_vec = c(mean = mean(nu_i), sd = sd(nu_i),  quantile(nu_i, probs = c(0.025, 0.25, 0.5, 0.75, 0.975)))
      out_df = rbind(out_df, nu_vec)
      rownames(out_df)[nrow(out_df)] = paste0('nu_',i)
    }
  }
  
  return(out_df)
}

retrieve_pars_summ = function(MH_obj, param_list = NULL, burn_in = 0.5, probs = c(0.05, 0.25, 0.5, 0.75, 0.95)){
  N = length(MH_obj$team_names)
  T_ = MH_obj$iterations
  
  t_start = round(T_*burn_in)
  
  if (is.null(param_list)){
    if (MH_obj$distr_type == 'P'){
      param_list = list('home', 'att', 'def')
    } else{
      param_list = list('home', 'att', 'def', 'nu')
    }
  }
  
  out = list()
  
  for (i in seq_along(param_list)) {
    mat = as.matrix(MH_obj[[paste0(param_list[[i]],'_post')]])
    mat = mat[t_start:T_, ,drop = FALSE]
    
    summaries = apply(mat, 2, function(x) {
      c(mean = mean(x), quantile(x, probs = probs))
    })
    
    summaries = t(summaries)
    
    rownames(summaries) = paste0(param_list[i], "_", seq_len(ncol(mat)))
    
    out[[i]] = as.data.frame(summaries)
  }
  out_df = do.call(rbind, out)
  out_df
}

retrieve_nu_sas_summ = function(Z_post, nu_post, start, end, nu_true = NULL){
  binary_mat = Z_post[start:end,]
  values_mat = nu_post[start:end,]
  
  masked_1 = values_mat * (binary_mat == 1)
  masked_0 = values_mat * (binary_mat == 0)
  count_1 = colSums(binary_mat == 1)
  count_0 = colSums(binary_mat == 0)
  
  colSums(masked_1) / count_1
  colSums(masked_0) / count_0
  
  N = ncol(values_mat)
  medians = numeric(N)
  q10 = numeric(N)
  q90 = numeric(N)
  
  if (is.null(nu_true)){
    nu_true = rep(NA_real_, N)
  }
  
  # Loop through columns
  for (i in 1:N) {
    vals = values_mat[, i]
    mask = binary_mat[, i] == 1
    masked_vals = vals[mask]
    
    if (length(masked_vals) > 0) {
      medians[i] = median(masked_vals)
      q10[i]     = quantile(masked_vals, probs = 0.10, names = FALSE)
      q90[i]    = quantile(masked_vals, probs = 0.90, names = FALSE)
    } else {
      medians[i] = NA
      q10[i] = NA
      q90[i] = NA
    }
  }
  
  
  summary_df = data.frame(
    q0.10 = q10,
    q0.50 = medians,
    q0.90 = q90
  )
  sas_df = t(round(cbind(p_z1 = count_1/(end-start), nu_true, summary_df), 4))
  
  return(sas_df)
}


thin_mcmc = function(MH_object, thin_by){
  t = MH_object$iterations
  
  stopifnot(thin_by %% 1 == 0)
  
  MH_object$att_post = MH_object$att_post[seq(1, t, by = thin_by), ]
  MH_object$def_post = MH_object$def_post[seq(1, t, by = thin_by), ]
  MH_object$home_post = MH_object$home_post[seq(1, t, by = thin_by)]
  MH_object$eta_post = MH_object$eta_post[seq(1, t, by = thin_by)]
  
  
  if (MH_object$distr_type == 'CMP-SAS'){
    MH_object$nu_post = MH_object$nu_post[seq(1, t, by = thin_by),]
    MH_object$Z_post = MH_object$Z_post[seq(1, t, by = thin_by),]
    MH_object$p_post = MH_object$p_post[seq(1, t, by = thin_by),]
  }
  
  if (MH_object$distr_type == 'CMP-Full'){
    MH_object$nu_post = MH_object$nu_post[seq(1, t, by = thin_by),]
    MH_object$Z_post = MH_object$Z_post[seq(1, t, by = thin_by),]
    MH_object$p_post = MH_object$p_post[seq(1, t, by = thin_by),]
  }
  
  MH_object$iterations = t/thin_by
  
  return(MH_object)
}

### Functions to generate matrices for mid-league inferences

generate_X_mid = function(df_hist, n_games){
  # Generate a binary matrix for that holds 1 for the 'n_games' that have been played, 
  # taken in chronological order from df_hist.
  N = length(unique(df_hist[,"HomeTeam"]))
  X_mid = matrix(0L, N, N)
  
  for(t in 1:n_games){
    hometeam = df_hist[t, "HomeTeam"]
    awayteam = df_hist[t, "AwayTeam"]
    i = get_team_code(hometeam)
    j = get_team_code(awayteam)
    X_mid[i,j] = 1
  }
  
  stopifnot(sum(X_mid) == n_games) #sanity check
  return(X_mid)
}

generate_inverse_X_mid = function(X_mid){
  # Generate a binary matrix from 'X_mid', which holds 1 for games that have yet to be played.
  inv_X_mid = +(!X_mid)
  diag(inv_X_mid) = 0
  return(inv_X_mid)}

### Evaluation and simulation functions

evaluate_standings_league = function(X_sim, goal_diff = FALSE){
  X1 = X_sim[[1]]
  X2 = X_sim[[2]]
  
  team_names = rownames(X1)
  N = nrow(X1)
  points = setNames(rep(0, N), team_names) 
  gd = setNames(rep(0, N), team_names)
  for (i in 1:N){
    for (j in 1:N){
      if (i == j || is.na(X1[i,j])){next}
      
      if (X1[i,j] > X2[i,j]) {
        points[i] = points[i] + 3
      }
      if (X1[i,j] == X2[i,j]) {
        points[i] = points[i] + 1
        points[j] = points[j] + 1
      }
      if (X1[i,j] < X2[i,j]) {
        points[j] = points[j] + 3
      }
      
      gd[i] = gd[i] + X1[i,j] - X2[i,j]
      gd[j] = gd[j] + X2[i,j] - X1[i,j]
    }
  }
  if (goal_diff){return(list(points = points, gd = gd))}
  return(points)
}

generate_league_sim_p = function(att, def, home, inv_X_mid = FALSE, X1 = FALSE, X2 = FALSE){
  if(typeof(inv_X_mid) != 'double'){
    X1_sim = matrix(NA, nrow = N, ncol = N)
    X2_sim = matrix(NA, nrow = N, ncol = N)
    for (i in 1:N){
      for (j in 1:N){
        if (i == j){next}
        X1_sim[i,j] = rpois(1, lambda = exp(att[i] + def[j] + home))
        X2_sim[j,i] = rpois(1, lambda = exp(att[i] + def[j]))
      }
    }
  }
  else{
    stopifnot(typeof(X1) == "double" && typeof(X2) == "double")
    
    X1_sim = X1
    X2_sim = X2
    for (i in 1:N){
      for (j in 1:N){
        if (inv_X_mid[i,j] == 0){next}
        X1_sim[i,j] = rpois(1, lambda = exp(att[i] + def[j] + home))
        X2_sim[j,i] = rpois(1, lambda = exp(att[i] + def[j]))
      }
    }
  }
  return(list(X1_sim, X2_sim))
}

generate_league_sim_cp = function(att, def, home, nu_home, nu_away, 
                                  inv_X_mid = FALSE, X1 = FALSE, X2 = FALSE){
  if(typeof(inv_X_mid) != 'double'){
    X1_sim = matrix(NA, nrow = N, ncol = N)
    X2_sim = matrix(NA, nrow = N, ncol = N)
    for (i in 1:N){
      for (j in 1:N){
        if (i == j){next}
        X1_sim[i,j] = rejection_sampler_draws(1, 
                                              mu = exp(att[i] + def[j] + home),
                                              nu = nu_home)
        X2_sim[j,i] = rejection_sampler_draws(1, 
                                              mu = exp(att[i] + def[j]),
                                              nu = nu_away)
      }
    }
  }
  else{
    stopifnot(typeof(X1) == "double" && typeof(X2) == "double")
    
    X1_sim = X1
    X2_sim = X2
    for (i in 1:N){
      for (j in 1:N){
        if (inv_X_mid[i,j] == 0){next}
        X1_sim[i,j] = rejection_sampler_draws(1, 
                                              mu = exp(att[i] + def[j] + home),
                                              nu = nu_home)
        X2_sim[j,i] = rejection_sampler_draws(1, 
                                              mu = exp(att[i] + def[j]),
                                              nu = nu_away)
      }
    }
  }
  return(list(X1_sim, X2_sim))
}

generate_league_sim_cp_d = function(att, def, home, nu, 
                                    inv_X_mid = FALSE, X1 = FALSE, X2 = FALSE){
  if(typeof(inv_X_mid) != 'double'){
    X1_sim = matrix(NA, nrow = N, ncol = N)
    X2_sim = matrix(NA, nrow = N, ncol = N)
    for (i in 1:N){
      for (j in 1:N){
        if (i == j){next}
        X1_sim[i,j] = rejection_sampler_draws(1, 
                                              mu = exp(att[i] + def[j] + home),
                                              nu = nu[i])
        X2_sim[j,i] = rejection_sampler_draws(1, 
                                              mu = exp(att[i] + def[j]),
                                              nu = nu[i])
      }
    }
  }
  else{
    stopifnot(typeof(X1) == "double" && typeof(X2) == "double")
    print('not implemented!')
    
    
    # X1_sim = X1
    # X2_sim = X2
    # for (i in 1:N){
    #   for (j in 1:N){
    #     if (inv_X_mid[i,j] == 0){next}
    #     X1_sim[i,j] = rejection_sampler_draws(1, 
    #                                           mu = exp(att[i] + def[j] + home),
    #                                           nu = nu_home)
    #     X2_sim[j,i] = rejection_sampler_draws(1, 
    #                                           mu = exp(att[i] + def[j]),
    #                                           nu = nu_away)
    #   }
    # }
  }
  return(list(X1_sim, X2_sim))
}


create_heatmatrix = function(home, away){
  df = data.frame(home,away)
  nmax = 8
  heatmatrix = matrix(0, nrow = nmax, ncol =nmax)
  for (n in 1:nrow(df)){
    i = min(df[n,1]+1, 8)
    j = min(df[n,2]+1, 8)
    heatmatrix[i,j] = heatmatrix[i,j] + 1
  }
  m.hmat = melt(heatmatrix)
  m.hmat[,1:2] = m.hmat[,1:2]-1
  return(list(m.hmat, heatmatrix))
}

sim_pair_results_pois = function(i, j, MH_object, id, n_sims = 10){
  att_sample = MH_object$att_post[id]
  def_sample = MH_object$def_post[id]
  home_sample = MH_object$home_post[id]
  att_fix = MH_object$att_fix
  def_fix = MH_object$def_fix
  gh = c()
  ga = c()
  n = length(id)
  for (s in 1:n){
    gh = append(gh, rpois(n_sims, lambda = exp(att_sample[[s]][i] + att_fix + def_fix + def_sample[[s]][j] + home_sample[s])))
    ga = append(ga, rpois(n_sims, lambda = exp(att_sample[[s]][j] + att_fix + def_fix + def_sample[[s]][i])))
  }
  heatmatrix = create_heatmatrix(gh, ga)
  return(heatmatrix)
}

sim_pair_results_compois = function(i, j, MH_object, id, n_sims = 10){
  att_sample = MH_object$att_post[id]
  def_sample = MH_object$def_post[id]
  home_sample = MH_object$home_post[id]
  if (MH_object$distr_type == 'CP'){
    nu_home_sample = MH_object$nu_home_post[id]
    nu_away_sample = MH_object$nu_away_post[id]
  } else {
    nu_sample = MH_object$nu_post[id]
  }
  
  att_fix = MH_object$att_fix
  def_fix = MH_object$def_fix
  gh = c()
  ga = c()
  n = length(id)
  for (s in 1:n){
    if (MH_object$distr_type == 'CP'){
      gh = append(gh, rejection_sampler_draws(n_sims, mu = exp(att_sample[[s]][i] + att_fix + def_fix + def_sample[[s]][j] + home_sample[s]),
                                              nu = nu_home_sample[s]))
      ga = append(ga, rejection_sampler_draws(n_sims, mu = exp(att_sample[[s]][j] + att_fix + def_fix + def_sample[[s]][i]),
                                              nu = nu_away_sample[s]))
    } else {
      gh = append(gh, rejection_sampler_draws(n_sims, mu = exp(att_sample[[s]][i] + att_fix + def_fix + def_sample[[s]][j] + home_sample[s]),
                                              nu = nu_sample[[s]][i]))
      ga = append(ga, rejection_sampler_draws(n_sims, mu = exp(att_sample[[s]][j] + att_fix + def_fix + def_sample[[s]][i]),
                                              nu = nu_sample[[s]][j]))
    }
  }
  heatmatrix = create_heatmatrix(gh, ga)
  return(heatmatrix)
}

compute_implied_probability = function(quotes){ #Distribute margin evenly across bets
  book_margin = sum(1/quotes) - 1
  implied_prob = (1/quotes)/(book_margin+1)
  # print(book_margin)
  return(list(implied_prob,
              book_margin))
}

plot_heatmap = function(i, j, melt_mat, hometeam, awayteam, distr){
  print(ggplot(melt_mat, aes(Var1, Var2)) +
          geom_tile(aes(fill = round(value/sum(value),4)*100)) +
          geom_text(aes(label = round(value/sum(value),4)*100)) +
          scale_fill_gradient(low = "gray95", high = "red4", guide = 'none') +
          scale_x_continuous("Home", labels = as.character(seq(0, 7)), breaks = seq(0, 7)) +
          scale_y_continuous("Away", labels = as.character(seq(0, 7)), breaks = seq(0, 7)) +
          labs(x = "Home", y = "Away", title = paste0(hometeam,"-",awayteam,", Result: ",X1[i,j],"-",X2[i,j],", ",distr," Model")))
}

create_bet = function(game_id,hometeam, awayteam, id_bet, quote, result, book_prob, pred_prob, book_margin){
  bet = switch(as.character(id_bet), "1" = "1", "2" = "X", "3" = "2")
  result = switch(result, "H" = "1", "D"="X", "A"="2")
  out = setNames(list( game_id,
                       hometeam,
                       awayteam,
                       quote,
                       bet,
                       result,
                       bet == result,
                       book_prob,
                       pred_prob,
                       pred_prob - book_prob,
                       book_margin), 
                 c("GameNumber","HomeTeam", "AwayTeam", "AvgQuote", "Bet", 
                   "Result", "BetWin", "ImpliedProb","ModelProb","OurMargin","BookMargin"))
  return(out)
}

bet_evaluator = function(df_row,stake){
  if(df_row["BetWin"]){return(as.numeric(df_row["AvgQuote"])*(stake))}
  else(return(-stake))
}

evaluate_games_prob = function(MH_object_P, MH_object_CP, game_id, id_sample, 
                               bets = c("AvgH", "AvgD", "AvgA"), plot_heatmap_fl = FALSE, bet_margin = 0){
  quotes = df_hist[t,bets]
  implied_probability = compute_implied_probability(quotes)[[1]]
  book_margin = compute_implied_probability(quotes)[[2]]
  
  hometeam = df_hist[t, "HomeTeam"]
  awayteam = df_hist[t, "AwayTeam"]
  result = df_hist[t, "FTR"]
  i = get_team_code(hometeam)
  j = get_team_code(awayteam)
  
  m.hmat = sim_pair_results_pois(i, j, MH_object_P, id_sample)
  melt_mat = m.hmat[[1]]
  if(plot_heatmap_fl){plot_heatmap(i, j, melt_mat, hometeam, awayteam, "Poisson")}
  
  freq_mat = m.hmat[[2]]
  freq_mat = freq_mat/(sum(freq_mat))
  
  home_prob = sum(freq_mat[lower.tri(freq_mat)])
  away_prob = sum(freq_mat[upper.tri(freq_mat)])
  draw_prob = sum(diag(freq_mat))
  predicted_probabilities_P = c(home_prob, draw_prob, away_prob)
  
  df_bets_P = setNames(data.frame(matrix(ncol = 11, nrow = 0)), 
                       c("GameNumber","HomeTeam", "AwayTeam", "AvgQuote", "Bet", 
                         "Result", "BetWin", "ImpliedProb","ModelProb","OurMargin","BookMargin"))
  diff1 = predicted_probabilities_P - implied_probability
  idx_bets = which((predicted_probabilities_P/(1+book_margin)) - implied_probability > bet_margin)
  for(id_bet in idx_bets){
    bet = create_bet(game_id,hometeam, awayteam, id_bet, quotes[[id_bet]], result, 
                     implied_probability[[id_bet]], predicted_probabilities_P[[id_bet]], book_margin)
    df_bets_P = rbind(df_bets_P, bet)
  }
  
  m.hmat = sim_pair_results_compois(i, j, MH_object_CP, id_sample)
  melt_mat = m.hmat[[1]]
  if(plot_heatmap_fl){plot_heatmap(i, j, melt_mat, hometeam, awayteam, "COM-Poisson")}
  
  freq_mat = m.hmat[[2]]
  freq_mat = freq_mat/(sum(freq_mat))
  
  home_prob = sum(freq_mat[lower.tri(freq_mat)])
  away_prob = sum(freq_mat[upper.tri(freq_mat)])
  draw_prob = sum(diag(freq_mat))
  predicted_probabilities_CP = c(home_prob, draw_prob, away_prob)
  
  df_bets_CP = setNames(data.frame(matrix(ncol = 11, nrow = 0)), 
                        c("GameNumber","HomeTeam", "AwayTeam", "AvgQuote", "Bet", 
                          "Result", "BetWin", "ImpliedProb","ModelProb","OurMargin","BookMargin"))
  diff2 = predicted_probabilities_CP - implied_probability
  # idx_bets = which(book_margin+bet_margin < diff2)
  idx_bets = which((predicted_probabilities_CP/(1+book_margin)) - implied_probability > bet_margin)
  for(id_bet in idx_bets){1
    bet = create_bet(game_id,hometeam, awayteam, id_bet, quotes[[id_bet]], result, 
                     implied_probability[[id_bet]], predicted_probabilities_CP[[id_bet]], book_margin)
    df_bets_CP = rbind(df_bets_CP, bet)
  }
  
  out_df = rbind(implied_probability, predicted_probabilities_P, predicted_probabilities_CP,
                 diff1, diff2)
  rownames(out_df) = c("BetAvg", "Poisson", "Com-Poisson", "Diff-P", "Diff-CP")
  colnames(out_df) = c("Home", "Draw", "Away")
  
  out_list = list(probabilities = out_df,
                  bets_P = df_bets_P,
                  bets_CP = df_bets_CP)
  return(out_list)
}


truncate_forecast = function(P, truncate_n = 2) {
  n = nrow(P)
  keep_n = n - truncate_n
  
  # Accumulate truncated rows and columns into last kept row/column
  P[keep_n, 1:keep_n] = P[keep_n, 1:keep_n] + colSums(P[(keep_n+1):n, 1:keep_n])
  P[1:keep_n, keep_n] = P[1:keep_n, keep_n] + rowSums(P[1:keep_n, (keep_n+1):n])
  P[keep_n, keep_n] = P[keep_n, keep_n] + sum(P[(keep_n+1):n, (keep_n+1):n])
  
  # Return truncated matrix, normalized
  # P[1:keep_n, 1:keep_n] / sum(P[1:keep_n, 1:keep_n])
  return(P[1:keep_n, 1:keep_n])
}

safe_log = function(p, eps = 1e-15) {
  p = pmax(p, eps)   
  return(-log2(p))             
}

rotate_ccw = function(m) {
  m[, ncol(m):1] |> t()
}

compute_goaldiff_probs = function(preds){
  n = nrow(preds)
  
  # Initialize vector for goal difference probabilities
  goal_diff_probs = numeric(2*n - 1)  # from +7 to -7 including 0
  names(goal_diff_probs) =  (n-1): - (n-1)
  
  for (k in -(n-1):(n-1)) {
    goal_diff_probs[as.character(k)] = sum(preds[row(preds) - col(preds) == -k])
  }
  
  return(goal_diff_probs)
}

outcome_rps = function(probs, outcome) {
  # probs: vector of length 3 with predicted probabilities
  # outcome: observed outcome (1, 2, or 3)
  
  K = length(probs)
  # Cumulative predicted probabilities
  cum_probs = cumsum(probs)
  
  # Cumulative observed outcome vector
  outcome_vec = rep(0, K)
  outcome_vec[outcome] = 1
  cum_obs = cumsum(outcome_vec)
  
  # Compute RPS
  score = sum((cum_probs[-K] - cum_obs[-K])^2) / (K - 1)
  return(score)
}

gd_rps = function(probs, cum_obs){
  K = length(probs)
  
  cum_probs = cumsum(probs)
  
  score = sum((cum_probs[-K] - cum_obs[-K])^2) / (K - 1)
  return(score)
}

ES_rps2d = function(P, ag_true, hg_true) {
  # P: predicted probability matrix (rows=x, cols=y)
  # x_true, y_true: true outcome coordinates (row, col)
  
  nrow_P = nrow(P)
  ncol_P = ncol(P)
  
  # Observed outcome matrix
  O = matrix(0, nrow=nrow_P, ncol=ncol_P)
  O[ag_true, hg_true] = 1
  
  # Initialize cumulative matrices
  F_cum = O_cum = matrix(0, nrow=nrow_P, ncol=ncol_P)
  
  for(i in 1:nrow_P) {
    for(j in 1:ncol_P) {
      F_cum[i,j] = sum(P[1:i, 1:j])
      O_cum[i,j] = sum(O[1:i, 1:j])
    }
  }
  
  # 2D cumulative RPS score
  score = sum((F_cum - O_cum)^2) / ((nrow_P - 1) * (ncol_P - 1))
  
  return(score)
}

# Manhattan-distance 2D score
ES_rps2d_manhattan = function(P, ag_true, hg_true, exponent=2) {
  # P: predicted probability matrix (rows=x, cols=y)
  # x_true, y_true: true outcome coordinates (row, col)
  # exponent: power applied to distance (default squared)
  
  nrow_P = nrow(P)
  ncol_P = ncol(P)
  
  # Coordinates matrices
  x_coords = matrix(rep(1:nrow_P, times=ncol_P), nrow=nrow_P)
  y_coords = matrix(rep(1:ncol_P, each=nrow_P), nrow=nrow_P)
  
  # Manhattan distance from true outcome
  manhattan_dist = abs(x_coords - ag_true) + abs(y_coords - hg_true)
  
  # Weighted score
  score = sum(P * (manhattan_dist^exponent))
  
  return(score)
}

simulation_league_posterior = function(MH_object, id, distr_type){
  att_sample = MH_object$att_post[id]
  def_sample = MH_object$def_post[id]
  home_sample = MH_object$home_post[id]
  n = length(id)
  # if(MH_object$distr_type == 'poisson'){
  points = data.frame(matrix(ncol=20, nrow=0))
  if(distr_type == 'P'){
    league_sim = generate_league_sim_p(att_sample[[1]] + att_fix, def_sample[[1]] + def_fix, home_sample[1])
    points = evaluate_standings_league(league_sim)
    for (s in 2:n){
      league_sim = generate_league_sim_p(att_sample[[s]] + att_fix, def_sample[[s]] + def_fix, home_sample[s])
      points_sim = evaluate_standings_league(league_sim)
      points = rbind(points, points_sim)
    }
  }
  # if(MH_object$distr_type == 'com_poisson'){
  if(distr_type == 'CP'){
    nu_home_sample = MH_object$nu_home_post[id]
    nu_away_sample = MH_object$nu_away_post[id]
    
    league_sim = generate_league_sim_cp(att_sample[[1]] + att_fix, def_sample[[1]] + def_fix,
                                        home_sample[1], nu_home = nu_home_sample[1], nu_away = nu_away_sample[1])
    points = evaluate_standings_league(league_sim)
    for (s in 2:n){
      league_sim = generate_league_sim_cp(att_sample[[s]] + att_fix, def_sample[[s]] + def_fix, 
                                          home_sample[s], nu_home = nu_home_sample[s], nu_away = nu_away_sample[s])
      points_sim = evaluate_standings_league(league_sim)
      points = rbind(points, points_sim)
    }
  }
  return(points)
}

generate_ci_plot_comparison = function(X, MH_obj_P, MH_obj_CP, id_sample, text_main = "", quantiles = c(0.1, 0.5, 0.9)){
  post_points_P = simulation_league_posterior(MH_obj_P, id_sample, 'P')
  post_points_CP = simulation_league_posterior(MH_obj_CP, id_sample, 'CP')
  
  # print(post_points_P)
  
  true_points = sort(evaluate_standings_league(X), decreasing = TRUE)
  names_standings = names(true_points)
  # print(names_standings)
  lower_ci = numeric(20)
  upper_ci = numeric(20)
  mode = numeric(20)
  
  lower_ci_cp = numeric(20)
  upper_ci_cp = numeric(20)
  mode_cp = numeric(20)
  for(i in 1:20){
    quants = quantile(post_points_P[,names_standings[i]], probs = quantiles)
    lower_ci[i] = quants[1]
    upper_ci[i] = quants[3]
    mode[i] = Mode(post_points_P[,names_standings[i]])
    
    quants = quantile(post_points_CP[,names_standings[i]], probs = quantiles)
    lower_ci_cp[i] = quants[1]
    upper_ci_cp[i] = quants[3]
    mode_cp[i] = Mode(post_points_CP[,names_standings[i]])
  }
  
  df_plot = data.frame(team = names_standings, 
                       true_points = true_points, 
                       lower_ci = lower_ci, 
                       lower_ci_cp = lower_ci_cp,
                       upper_ci = upper_ci,
                       upper_ci_cp = upper_ci_cp,
                       mode = mode,
                       mode_cp = mode_cp)
  df_plot$team = factor(df_plot$team, levels = df_plot$team)
  
  p = ggplot(df_plot) +
    geom_point(aes(y = as.numeric(team) + 0.2, x = mode), color = "blue") +       # First set of points
    geom_errorbarh(aes(y = as.numeric(team) + 0.2, xmin = lower_ci, xmax = upper_ci), 
                   height = 0.2, color = "blue") +                                        # First set of error bars
    geom_point(aes(y = as.numeric(team) - 0.2, x = mode_cp), color = "red") +        # Second set of points
    geom_errorbarh(aes(y = as.numeric(team) - 0.2, xmin = lower_ci_cp, xmax = upper_ci_cp), 
                   height = 0.2, color = "red") +  
    geom_point(aes(y = as.numeric(team) + 0.2, x = true_points)) +
    geom_point(aes(y = as.numeric(team) - 0.2, x = true_points)) +# Second set of error bars
    scale_y_continuous(breaks = 1:20, labels = levels(df_plot$team)) +                      # Correct y-axis labels
    labs(title = paste("CI for final Points", text_main),
         x = "Points",
         y = "Team") +                                                                    # Add labels
    coord_cartesian(xlim = c(0, 105)) +
    theme_minimal()  
  
  return(p)
}

plot_MCMC_diagnostics = function(MH_object_par, cols = 1:5, bins = 30) {
  n = length(cols)
  
  # set up plotting layout (n rows, 2 columns)
  par(mfrow = c(n, 2), mar = c(3, 3, 2, 1))
  
  for (i in cols) {
    # trace plot
    plot(MH_object_par[, i], type = "l", 
         main = paste("Trace:", i), 
         xlab = "Iteration", ylab = "Value")
    
    # histogram
    hist(MH_object_par[, i], breaks = bins, col = "grey", border = "white",
         main = paste("Histogram:", i), xlab = "Value")
  }
}

plot_MCMC_diagnostics_series = function(par1, par2 , cols = 1:5) {
  n = length(cols)
  
  # set up plotting layout (n rows, 2 columns)
  par(mfrow = c(n, 2), mar = c(3, 3, 2, 1))
  
  for (i in cols) {
    # trace plot
    plot(par1[, i], type = "l", 
         main = paste("Att:", i), 
         xlab = "Iteration", ylab = "Value")
    
    plot(par2[, i], type = "l", 
         main = paste("Nu:", i), 
         xlab = "Iteration", ylab = "Value")
    abline(h = 1, col = "red2", lty = "dashed", lwd = 1)
  }
}

### Functions to retrieve results from MH objects

simulate_game_probs = function(MH_object_P, MH_object_CP, game_id, id_sample, df_hist, n_sims = 1,
                               plot_heatmap_fl = FALSE){
  
  hometeam = df_hist[game_id, "HomeTeam"]
  awayteam = df_hist[game_id, "AwayTeam"]
  result = df_hist[game_id, c("FTHG", "FTAG", "FTR", "HomeTeam", "AwayTeam")]
  i = get_team_code(hometeam)
  j = get_team_code(awayteam)
  
  m.hmat = sim_pair_results_pois(i, j, MH_object_P, id_sample, n_sims)
  melt_mat = m.hmat[[1]]
  if(plot_heatmap_fl){plot_heatmap(i, j, melt_mat, hometeam, awayteam, "Poisson")}
  
  freq_mat = m.hmat[[2]]
  freq_mat = freq_mat/(sum(freq_mat))
  
  freq_mat_P = freq_mat
  
  home_prob = sum(freq_mat[lower.tri(freq_mat)])
  away_prob = sum(freq_mat[upper.tri(freq_mat)])
  draw_prob = sum(diag(freq_mat))
  predicted_probabilities_P = c(home_prob, draw_prob, away_prob)
  
  m.hmat = sim_pair_results_compois(i, j, MH_object_CP, id_sample, n_sims)
  melt_mat = m.hmat[[1]]
  if(plot_heatmap_fl){plot_heatmap(i, j, melt_mat, hometeam, awayteam, "COM-Poisson")}
  
  freq_mat = m.hmat[[2]]
  freq_mat = freq_mat/(sum(freq_mat))
  freq_mat_CP = freq_mat
  
  home_prob = sum(freq_mat[lower.tri(freq_mat)])
  away_prob = sum(freq_mat[upper.tri(freq_mat)])
  draw_prob = sum(diag(freq_mat))
  predicted_probabilities_CP = c(home_prob, draw_prob, away_prob)
  
  out_df = rbind(predicted_probabilities_P, predicted_probabilities_CP)
  
  out_list = list(probabilities = out_df,
                  freq_P = freq_mat_P,
                  freq_CP = freq_mat_CP,
                  result = result)
  return(out_list)
}

simulate_game_probs_d = function(MH_object_P, MH_object_CP, MH_object_CPd, game_id, id_sample, df_hist, n_sims = 1,
                                 plot_heatmap_fl = FALSE){
  
  hometeam = df_hist[game_id, "HomeTeam"]
  awayteam = df_hist[game_id, "AwayTeam"]
  result = df_hist[game_id, c("FTHG", "FTAG", "FTR", "HomeTeam", "AwayTeam")]
  i = get_team_code(hometeam)
  j = get_team_code(awayteam)
  
  m.hmat = sim_pair_results_pois(i, j, MH_object_P, id_sample, n_sims)
  melt_mat = m.hmat[[1]]
  if(plot_heatmap_fl){plot_heatmap(i, j, melt_mat, hometeam, awayteam, "Poisson")}
  
  freq_mat = m.hmat[[2]]
  freq_mat = freq_mat/(sum(freq_mat))
  
  freq_mat_P = freq_mat
  
  home_prob = sum(freq_mat[lower.tri(freq_mat)])
  away_prob = sum(freq_mat[upper.tri(freq_mat)])
  draw_prob = sum(diag(freq_mat))
  predicted_probabilities_P = c(home_prob, draw_prob, away_prob)
  
  m.hmat = sim_pair_results_compois(i, j, MH_object_CP, id_sample, n_sims)
  melt_mat = m.hmat[[1]]
  if(plot_heatmap_fl){plot_heatmap(i, j, melt_mat, hometeam, awayteam, "COM-Poisson")}
  
  freq_mat = m.hmat[[2]]
  freq_mat = freq_mat/(sum(freq_mat))
  
  freq_mat_CP = freq_mat
  
  home_prob = sum(freq_mat[lower.tri(freq_mat)])
  away_prob = sum(freq_mat[upper.tri(freq_mat)])
  draw_prob = sum(diag(freq_mat))
  predicted_probabilities_CP = c(home_prob, draw_prob, away_prob)
  
  m.hmat = sim_pair_results_compois(i, j, MH_object_CPd, id_sample, n_sims)
  melt_mat = m.hmat[[1]]
  if(plot_heatmap_fl){plot_heatmap(i, j, melt_mat, hometeam, awayteam, "COM-Poisson-ID")}
  
  freq_mat = m.hmat[[2]]
  freq_mat = freq_mat/(sum(freq_mat))
  freq_mat_CPd = freq_mat
  
  home_prob = sum(freq_mat[lower.tri(freq_mat)])
  away_prob = sum(freq_mat[upper.tri(freq_mat)])
  draw_prob = sum(diag(freq_mat))
  predicted_probabilities_CPd = c(home_prob, draw_prob, away_prob)
  
  out_df = rbind(predicted_probabilities_P, predicted_probabilities_CP,predicted_probabilities_CPd)
  
  out_list = list(probabilities = out_df,
                  freq_P = freq_mat_P,
                  freq_CP = freq_mat_CP,
                  freq_CPd = freq_mat_CPd,
                  result = result)
  return(out_list)
}

compute_scores_match = function(df_prob){
  hg = min(df_prob$result[[1]], 7)
  ag = min(df_prob$result[[2]], 7)
  result_to_index = list('H' = 1, 'D' = 2, 'A' = 3)
  res_idx = result_to_index[[df_prob$result[[3]]]]
  
  df_result_prob = df_prob[[1]]
  
  p_result_prob = df_result_prob[[1, res_idx]]
  p_probs = df_prob[[2]]
  
  p_exact_prob = p_probs[hg+1, ag+1]
  p_probs[hg+1, ag+1] = p_probs[hg+1, ag+1] - 1
  p_exact_brier = sum(p_probs**2)
  
  cp_result_prob = df_result_prob[[2, res_idx]]
  cp_probs = df_prob[[3]]
  
  cp_exact_prob = cp_probs[hg+1, ag+1]
  cp_probs[hg+1, ag+1] = cp_probs[hg+1, ag+1] - 1
  cp_exact_brier = sum(cp_probs**2)
  
  df_result_prob[,res_idx] = df_result_prob[,res_idx] - 1
  p_result_brier = rowSums(df_result_prob**2)[[1]]
  cp_result_brier = rowSums(df_result_prob**2)[[2]]
  
  out_list = list(result_prob_P = p_result_prob,
                  result_prob_CP = cp_result_prob,
                  result_brier_P = p_result_brier,
                  result_brier_CP = cp_result_brier,
                  exact_prob_P = p_exact_prob,
                  exact_prob_CP = cp_exact_prob,
                  exact_brier_P = p_exact_brier,
                  exact_brier_CP = cp_exact_brier)
  
  return(out_list)
}

#### WAIC

compute_Z = function(mu, nu, j_max = 100, tol = 1e-20) {
  Z = 0
  for (j in 0:j_max) {
    term = (mu^(nu * j)) / (factorial(j)^nu)
    Z = Z + term
    if (term < tol) {
      # print(j)
      break
    }
  }
  return(Z)
}

# Create grid of mu and nu
mu_vals = seq(0.01, 10, by = 0.01)
nu_vals = seq(0.01, 10, by = 0.01)

# Initialize lookup table
Z_lookup = matrix(NA, nrow = length(mu_vals), ncol = length(nu_vals),
                   dimnames = list(paste0("mu=", mu_vals), paste0("nu=", nu_vals)))

# Fill table
for (i in seq_along(mu_vals)) {
  for (j in seq_along(nu_vals)) {
    Z_lookup[i, j] = compute_Z(mu_vals[i], nu_vals[j])
  }
}

eval_cmp_llhood_cmp_vec = function(x, mu, nu){
  
  row = nearest_index(mu, mu_vals)
  col = nearest_index(nu, nu_vals)
  
  Z_vals = Z_lookup[cbind(row, col)]
  
  qf = nu * (x * log(mu) - lgamma(x + 1))
  return(qf - log(Z_vals))
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

evaluate_WAIC = function(X, MH_obj, n_samples = 2500, rseed = 1){
  T_ = MH_obj$iterations
  
  set.seed(rseed)
  
  sample_ids = get_sample_ids(n_samples, T_, (T_/5)+1)
  
  att = MH_obj$att_post[sample_ids,]
  def = MH_obj$def_post[sample_ids,]
  home = MH_obj$home_post[sample_ids]
  
  
  if (MH_obj$distr_type == 'P'){
    WAIC = eval_WAIC_pois(X, att, def, home)
  }
  
  else {
    nu = MH_obj$nu_post[sample_ids,]  
    WAIC = eval_WAIC_cmp(X, att, def, home, nu)
  }
  
  return(WAIC)
}


eval_WAIC_cmp = function(X, att, def, home, nu){
  X1 = X[[1]]
  X2 = X[[2]]
  N = nrow(X1)
  
  # WAIC = matrix(NA_real_, N, N)
  lppd = matrix(NA_real_, N, N)
  p_waic = matrix(NA_real_, N, N)
  
  for (i in 1:N){
    for (j in 1:N){
      if (i == j){next}
      
      mu_home = exp(att[,i] + def[,j] + home)
      nu_home = nu[,i]
      
      mu_away = exp(att[,j] + def[,i])
      nu_away = nu[,j]
      
      log_llhood_home = eval_cmp_llhood_cmp_vec(X1[i,j], mu_home, nu_home)
      # lppd_home = log(mean(exp(log_llhood_home)))
      p_waic_home = var(log_llhood_home)
      
      log_llhood_away = eval_cmp_llhood_cmp_vec(X2[i,j], mu_away, nu_away)
      # lppd_away = log(mean(exp(log_llhood_away)))
      p_waic_away = var(log_llhood_away)
      
      lppd_home = log_sum_exp(log_llhood_home)
      lppd_away = log_sum_exp(log_llhood_away)
      
      
      # WAIC[i,j] = (lppd_home - p_waic_home) + (lppd_away - p_waic_away)
      lppd[i,j] = lppd_home + lppd_away
      p_waic[i,j] = p_waic_home + p_waic_away
    }
  }
  # return(WAIC)
  # return(-2*sum(WAIC, na.rm = TRUE))
  return(list(lppd = sum(lppd, na.rm = TRUE),
              p_waic = sum(p_waic, na.rm = TRUE)))
}



eval_WAIC_pois = function(X, att, def, home){
  X1 = X[[1]]
  X2 = X[[2]]
  N = nrow(X1)
  
  # WAIC = matrix(NA_real_, N, N)
  lppd = matrix(NA_real_, N, N)
  p_waic = matrix(NA_real_, N, N)
  
  for (i in 1:N){
    for (j in 1:N){
      if (i == j){next}
      
      lam_home = exp(att[,i] + def[,j] + home)
      
      lam_away = exp(att[,j] + def[,i])
      
      log_llhood_home = dpois(X1[i,j], lam_home, log = TRUE)
      # lppd_home = log(mean(exp(log_llhood_home)))
      p_waic_home = var(log_llhood_home)
      
      log_llhood_away = dpois(X2[i,j], lam_away, log = TRUE)
      # lppd_away = log(mean(exp(log_llhood_away)))
      p_waic_away = var(log_llhood_away)
      
      lppd_home = log_sum_exp(log_llhood_home)
      lppd_away = log_sum_exp(log_llhood_away)
      
      
      # WAIC[i,j] = (lppd_home - p_waic_home) + (lppd_away - p_waic_away)
      lppd[i,j] = lppd_home + lppd_away
      p_waic[i,j] = p_waic_home + p_waic_away
    }
  }
  
  # return(WAIC)
  
  # return(-2*sum(WAIC, na.rm = TRUE))
  return(list(lppd = sum(lppd, na.rm = TRUE),
              p_waic = sum(p_waic, na.rm = TRUE)))
} 

log_sum_exp = function(x) {
  m = max(x)
  m + log(mean(exp(x - m)))
}

### Prediction Evaluation

truncate_forecast = function(P, truncate_n = 2) {
  n = nrow(P)
  keep_n = n - truncate_n
  
  # Accumulate truncated rows and columns into last kept row/column
  P[keep_n, 1:keep_n] = P[keep_n, 1:keep_n] + colSums(P[(keep_n+1):n, 1:keep_n])
  P[1:keep_n, keep_n] = P[1:keep_n, keep_n] + rowSums(P[1:keep_n, (keep_n+1):n])
  P[keep_n, keep_n] = P[keep_n, keep_n] + sum(P[(keep_n+1):n, (keep_n+1):n])
  
  # Return truncated matrix, normalized
  # P[1:keep_n, 1:keep_n] / sum(P[1:keep_n, 1:keep_n])
  return(P[1:keep_n, 1:keep_n])
}

prob_matrix_to_outcome = function(prob_matrix){
  home_win = sum(prob_matrix[upper.tri(prob_matrix)])
  away_win = sum(prob_matrix[lower.tri(prob_matrix)])
  draw = sum(diag(prob_matrix))
  
  return(list(home = home_win/100,
              draw = draw/100,
              away = away_win/100))
}

safe_log = function(p, eps = 1e-15) {
  p = pmax(p, eps)   
  return(-log2(p))             
}

rotate_ccw = function(m) {
  m[, ncol(m):1] |> t()
}

compute_goaldiff_probs = function(preds){
  n = nrow(preds)
  
  # Initialize vector for goal difference probabilities
  goal_diff_probs = numeric(2*n - 1)  # from +7 to -7 including 0
  names(goal_diff_probs) = - (n-1):(n-1)
  
  for (k in -(n-1):(n-1)) {
    goal_diff_probs[as.character(k)] = sum(preds[row(preds) - col(preds) == -k])
  }
  
  return(goal_diff_probs)
}

compute_overunder_probs = function(preds, threshold = 2.5){
  n = nrow(preds)
  
  overunder_vec = numeric(2)
  names(overunder_vec) = c(paste0('under',as.character(threshold)),
                           paste0('over', as.character(threshold)))
  
  goals = 0:(n - 1)
  
  # Create matrix of total goals for each (i, j)
  total_goals = outer(goals, goals, "+")
  
  # Logical masks for under / over threshold
  under_mask = total_goals < threshold
  over_mask  = total_goals > threshold
  
  # Since entries are probabilities *100, divide by 100
  overunder_vec[1] = sum(preds[under_mask]) / 100
  overunder_vec[2]  = sum(preds[over_mask]) / 100
  
  return(overunder_vec)
}


outcome_rps = function(probs, outcome) {
  # probs: vector of length 3 with predicted probabilities
  # outcome: observed outcome (1, 2, or 3)
  
  K = length(probs)
  # Cumulative predicted probabilities
  cum_probs = cumsum(probs)
  
  # Cumulative observed outcome vector
  outcome_vec = rep(0, K)
  outcome_vec[outcome] = 1
  cum_obs = cumsum(outcome_vec)
  
  # Compute RPS
  score = sum((cum_probs[-K] - cum_obs[-K])^2) / (K - 1)
  return(score)
}


gd_rps = function(probs, cum_obs){
  K = length(probs)
  
  cum_probs = cumsum(probs)
  
  score = sum((cum_probs[-K] - cum_obs[-K])^2) / (K - 1)
  return(score)
}

overunder_rps = function(probs, cum_obs){
  K = length(probs)
  
  cum_probs = cumsum(probs)
  
  return(sum((cum_probs[-K] - cum_obs[-K])^2) / (K - 1))
}

compute_outcome_eval = function(league_acro, season_start, season_end, 
                                regime = 'insample', start_id = 1, filter = 0, threshold = 0.5){
  ftr_map = c(H = 1, D = 2, A = 3)
  
  seasons_strvec = generate_season_string(season_start, season_end)
  
  table_summ = data.frame(matrix(NA_real_, nrow = 6, length(seasons_strvec)))
  colnames(table_summ) = seasons_strvec
  
  for (season in seasons_strvec){
    if (league_acro == 'LC' && season == '1920'){next}
    df_hist = read_historics(season, league_acro)
    
    pois_pred = readRDS(file = paste0("Data/Predictions/",league_acro, season,"_Pois_", regime, ".rds"))$list_preds
    
    sas_pred = readRDS(file = paste0("Data/Predictions/",league_acro, season,"_SAS_", regime, ".rds"))$list_preds
    
    cmp_pred = readRDS(file = paste0("Data/Predictions/",league_acro, season,"_CMP_", regime, ".rds"))$list_preds
    
    tot_games = nrow(df_hist)
    pois_rps = rep(NA_real_, tot_games)  
    pois_ign = rep(NA_real_, tot_games)  
    sas_rps = rep(NA_real_, tot_games)  
    sas_ign = rep(NA_real_, tot_games) 
    cmp_rps = rep(NA_real_, tot_games)
    cmp_ign = rep(NA_real_, tot_games)
    
    if (filter == 0){
      match_ids = seq(1:tot_games)
    } 
    if (filter == 1){
      MH_SAS = readRDS(file = paste0("Data//MCMC_Outputs//SAS_FullLeague//",league_acro,"_",season,"_SAS.rds"))
      filter_idx = which(retrieve_nu_sas_summ(MH_SAS$Z_post, MH_SAS$nu_post, MH_SAS$iterations/10, MH_SAS$iterations)[1,] > threshold)
      team_filter = MH_SAS$team_names[filter_idx]
      print(team_filter)
      
      match_ids = as.numeric(rownames(subset(df_hist, HomeTeam %in% team_filter | AwayTeam %in% team_filter)))
    }
    if (filter == 2){
      MH_SAS = readRDS(file = paste0("Data//MCMC_Outputs//SAS_FullLeague//",league_acro,"_",season,"_SAS.rds"))
      filter_idx = which(retrieve_nu_sas_summ(MH_SAS$Z_post, MH_SAS$nu_post, MH_SAS$iterations/10, MH_SAS$iterations)[1,] > threshold)
      team_filter = MH_SAS$team_names[filter_idx]
      print(team_filter)
      
      match_ids = as.numeric(rownames(subset(df_hist, !(HomeTeam %in% team_filter | AwayTeam %in% team_filter))))
    }
    
    # for (game_id in start_id:tot_games){
    for (game_id in match_ids){
      if (regime == 'oos'){if (game_id <= tot_games/2){next}}
      
      game_pred_pois = pois_pred[[game_id]]
      game_pred_sas  = sas_pred[[game_id]]
      game_pred_cmp  = cmp_pred[[game_id]]
      
      home_win_p = sum(game_pred_pois[upper.tri(game_pred_pois)])
      away_win_p = sum(game_pred_pois[lower.tri(game_pred_pois)])
      draw_p = sum(diag(game_pred_pois))
      
      home_win_sas = sum(game_pred_sas[upper.tri(game_pred_sas)])
      away_win_sas = sum(game_pred_sas[lower.tri(game_pred_sas)])
      draw_sas = sum(diag(game_pred_sas))
      
      home_win_cmp = sum(game_pred_cmp[upper.tri(game_pred_cmp)])
      away_win_cmp = sum(game_pred_cmp[lower.tri(game_pred_cmp)])
      draw_cmp = sum(diag(game_pred_cmp))
      
      outcome_prob_p = c(home_win_p, draw_p, away_win_p)/100
      outcome_prob_sas = c(home_win_sas, draw_sas, away_win_sas)/100
      outcome_prob_cmp = c(home_win_cmp, draw_cmp, away_win_cmp)/100
      
      outcome = ftr_map[[ df_hist[game_id,'FTR'] ]]
      
      pois_rps[[game_id]] = outcome_rps(outcome_prob_p, outcome)
      sas_rps[[game_id]] = outcome_rps(outcome_prob_sas, outcome)
      cmp_rps[[game_id]] = outcome_rps(outcome_prob_cmp, outcome)
      
      pois_ign[[game_id]] = safe_log(outcome_prob_p[[outcome]])
      sas_ign[[game_id]] = safe_log(outcome_prob_sas[[outcome]])
      cmp_ign[[game_id]] = safe_log(outcome_prob_cmp[[outcome]])
      
      
    }
    
    table_summ[1, season] = mean(pois_rps, na.rm = TRUE)
    table_summ[2, season] = mean(sas_rps, na.rm = TRUE)
    table_summ[3, season] = mean(cmp_rps, na.rm = TRUE)
    table_summ[4, season] = mean(pois_ign, na.rm = TRUE)
    table_summ[5, season] = mean(sas_ign, na.rm = TRUE)
    table_summ[6, season] = mean(cmp_ign, na.rm = TRUE)
  }
  rownames(table_summ) = c('RPS-Pois', 'RPS-SAS', 'RPS-CMP', 'IGN-Pois', 'IGN-SAS', 'IGN-CMP')
  table_summ = cbind(table_summ, setNames(data.frame(rowMeans(table_summ, na.rm = TRUE)), "Means"))
  return(table_summ)
}

compute_goaldiff_eval = function(league_acro, season_start, season_end, 
                                 regime = 'insample', start_id = 1, filter = 0, threshold = 0.5){
  seasons_strvec = generate_season_string(season_start, season_end)
  table_summ = data.frame(matrix(NA_real_, nrow = 4, length(seasons_strvec)))
  colnames(table_summ) = seasons_strvec
  
  
  for (season in seasons_strvec){
    if (league_acro == 'LC' && season == '1920'){next}
    df_hist = read_historics(season, league_acro)
    pois_pred = readRDS(file = paste0("Data/Predictions/",league_acro, season,"_Pois_", regime, ".rds"))$list_preds
    
    sas_pred = readRDS(file = paste0("Data/Predictions/",league_acro, season,"_SAS_", regime, ".rds"))$list_preds
    
    cmp_pred = readRDS(file = paste0("Data/Predictions/",league_acro, season,"_CMP_", regime, ".rds"))$list_preds
    
    tot_games = nrow(df_hist)
    pois_rps = rep(NA_real_, tot_games)  
    pois_ign = rep(NA_real_, tot_games)  
    sas_rps = rep(NA_real_, tot_games)  
    sas_ign = rep(NA_real_, tot_games) 
    cmp_rps = rep(NA_real_, tot_games)
    cmp_ign = rep(NA_real_, tot_games)
    
    if (filter == 0){
      match_ids = seq(1:tot_games)
    } 
    if (filter == 1){
      MH_SAS = readRDS(file = paste0("Data//MCMC_Outputs//SAS_FullLeague//",league_acro,"_",season,"_SAS.rds"))
      filter_idx = which(retrieve_nu_sas_summ(MH_SAS$Z_post, MH_SAS$nu_post, MH_SAS$iterations/10, MH_SAS$iterations)[1,] > threshold)
      team_filter = MH_SAS$team_names[filter_idx]
      print(team_filter)
      
      match_ids = as.numeric(rownames(subset(df_hist, HomeTeam %in% team_filter | AwayTeam %in% team_filter)))
    }
    if (filter == 2){
      MH_SAS = readRDS(file = paste0("Data//MCMC_Outputs//SAS_FullLeague//",league_acro,"_",season,"_SAS.rds"))
      filter_idx = which(retrieve_nu_sas_summ(MH_SAS$Z_post, MH_SAS$nu_post, MH_SAS$iterations/10, MH_SAS$iterations)[1,] > threshold)
      team_filter = MH_SAS$team_names[filter_idx]
      print(team_filter)
      
      match_ids = as.numeric(rownames(subset(df_hist, !(HomeTeam %in% team_filter | AwayTeam %in% team_filter))))
    }
    
    gd_outcomes = 7:-7
    # for (game_id in start_id:tot_games){
    for (game_id in match_ids){
      if (regime == 'oos'){if (game_id <= tot_games/2){next}}
      hg = min(as.numeric(df_hist[game_id, 'FTHG']), 7)
      ag = min(as.numeric(df_hist[game_id, 'FTAG']), 7)
      goal_diff_obs = hg - ag
      
      cum_gd_obs = as.integer(gd_outcomes <= goal_diff_obs)
      
      # pois_goaldiff = compute_goaldiff_probs(truncate_forecast(pois_pred[[game_id]], 2)/100)
      # sas_goaldiff = compute_goaldiff_probs(truncate_forecast(sas_pred[[game_id]], 2)/100)
      # cmp_goaldiff = compute_goaldiff_probs(truncate_forecast(cmp_pred[[game_id]], 2)/100)
      
      pois_goaldiff = compute_goaldiff_probs(pois_pred[[game_id]]/100)
      sas_goaldiff = compute_goaldiff_probs(sas_pred[[game_id]]/100)
      cmp_goaldiff = compute_goaldiff_probs(cmp_pred[[game_id]]/100)
      
      pois_rps[[game_id]] = gd_rps(pois_goaldiff, cum_gd_obs)
      sas_rps[[game_id]] = gd_rps(sas_goaldiff, cum_gd_obs)
      cmp_rps[[game_id]] = gd_rps(cmp_goaldiff, cum_gd_obs)
      
      pois_ign[[game_id]] = safe_log(pois_goaldiff[as.character(goal_diff_obs)])
      sas_ign[[game_id]] = safe_log(sas_goaldiff[as.character(goal_diff_obs)])
      cmp_ign[[game_id]] = safe_log(cmp_goaldiff[as.character(goal_diff_obs)])
    }
    
    table_summ[1, season] = mean(pois_rps, na.rm = TRUE)
    table_summ[2, season] = mean(sas_rps, na.rm = TRUE)
    table_summ[3, season] = mean(cmp_rps, na.rm = TRUE)
    table_summ[4, season] = mean(pois_ign, na.rm = TRUE)
    table_summ[5, season] = mean(sas_ign, na.rm = TRUE)
    table_summ[6, season] = mean(cmp_ign, na.rm = TRUE)
  }
  
  
  rownames(table_summ) = c('RPS-Pois', 'RPS-SAS', 'RPS-CMP', 'IGN-Pois', 'IGN-SAS', 'IGN-CMP')
  table_summ = cbind(table_summ, setNames(data.frame(rowMeans(table_summ, na.rm = TRUE)), "Means"))
}

compute_overunder_eval = function(league_acro, season_start, season_end, overunder = 2.5,
                                  regime = 'insample', start_id = 1, filter = 0, threshold = 0.5){
  seasons_strvec = generate_season_string(season_start, season_end)
  table_summ = data.frame(matrix(NA_real_, nrow = 6, length(seasons_strvec)))
  colnames(table_summ) = seasons_strvec
  
  for (season in seasons_strvec){
    if (league_acro == 'LC' && season == '1920'){next}
    df_hist = read_historics(season, league_acro)
    
    pois_pred = readRDS(file = paste0("Data/Predictions/",league_acro, season,"_Pois_", regime, ".rds"))$list_preds
    
    sas_pred = readRDS(file = paste0("Data/Predictions/",league_acro, season,"_SAS_", regime, ".rds"))$list_preds
    
    cmp_pred = readRDS(file = paste0("Data/Predictions/",league_acro, season,"_CMP_", regime, ".rds"))$list_preds
    
    
    tot_games = nrow(df_hist)
    
    pois_rps = rep(NA_real_, tot_games)  
    pois_ign = rep(NA_real_, tot_games)  
    sas_rps = rep(NA_real_, tot_games)  
    sas_ign = rep(NA_real_, tot_games)  
    cmp_rps = rep(NA_real_, tot_games)  
    cmp_ign = rep(NA_real_, tot_games)  
    
    if (filter == 0){
      match_ids = seq(1:tot_games)
    } 
    if (filter == 1){
      MH_SAS = readRDS(file = paste0("Data//MCMC_Outputs//SAS_FullLeague//",league_acro,"_",season,"_SAS.rds"))
      filter_idx = which(retrieve_nu_sas_summ(MH_SAS$Z_post, MH_SAS$nu_post, MH_SAS$iterations/10, MH_SAS$iterations)[1,] > threshold)
      team_filter = MH_SAS$team_names[filter_idx]
      print(team_filter)
      
      match_ids = as.numeric(rownames(subset(df_hist, HomeTeam %in% team_filter | AwayTeam %in% team_filter)))
    }
    if (filter == 2){
      MH_SAS = readRDS(file = paste0("Data//MCMC_Outputs//SAS_FullLeague//",league_acro,"_",season,"_SAS.rds"))
      filter_idx = which(retrieve_nu_sas_summ(MH_SAS$Z_post, MH_SAS$nu_post, MH_SAS$iterations/10, MH_SAS$iterations)[1,] > threshold)
      team_filter = MH_SAS$team_names[filter_idx]
      print(team_filter)
      
      match_ids = as.numeric(rownames(subset(df_hist, !(HomeTeam %in% team_filter | AwayTeam %in% team_filter))))
    }
    
    # for (game_id in start_id:tot_games){
    for (game_id in match_ids){
      if (regime == 'oos'){if (game_id <= tot_games/2){next}}
      hg = min(as.numeric(df_hist[game_id, 'FTHG']), 7)
      ag = min(as.numeric(df_hist[game_id, 'FTAG']), 7)
      tot_goals = hg + ag
      
      obs_vec = c(as.numeric(tot_goals < overunder), as.numeric(tot_goals > overunder))
      cum_vec = c(as.numeric(tot_goals < overunder), 1)
      
      pois_overunder = compute_overunder_probs(pois_pred[[game_id]], overunder)
      sas_overunder = compute_overunder_probs(sas_pred[[game_id]], overunder)
      cmp_overunder = compute_overunder_probs(cmp_pred[[game_id]], overunder)
      
      
      pois_rps[[game_id]] = overunder_rps(pois_overunder, cum_vec)
      sas_rps[[game_id]] = overunder_rps(sas_overunder, cum_vec)
      cmp_rps[[game_id]] = overunder_rps(cmp_overunder, cum_vec)
      
      pois_ign[[game_id]] = safe_log(sum(pois_overunder*obs_vec))
      sas_ign[[game_id]] = safe_log(sum(sas_overunder*obs_vec))
      cmp_ign[[game_id]] = safe_log(sum(cmp_overunder*obs_vec))
    }
    
    table_summ[1, season] = mean(pois_rps, na.rm = TRUE)
    table_summ[2, season] = mean(sas_rps, na.rm = TRUE)
    table_summ[3, season] = mean(cmp_rps, na.rm = TRUE)
    table_summ[4, season] = mean(pois_ign, na.rm = TRUE)
    table_summ[5, season] = mean(sas_ign, na.rm = TRUE)
    table_summ[6, season] = mean(cmp_ign, na.rm = TRUE)
  }
  
  
  rownames(table_summ) = c('RPS-Pois', 'RPS-SAS', 'RPS-CMP', 'IGN-Pois', 'IGN-SAS', 'IGN-CMP')
  table_summ = cbind(table_summ, setNames(data.frame(rowMeans(table_summ, na.rm = TRUE)), "Means"))
}
