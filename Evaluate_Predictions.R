source('utils.R')
library(dplyr)
library(kableExtra)

league_acros = c('PL', 'SA', 'LL', 'LC', 'BL', 'WSL')
mcmc_out_dir = "Data/MCMC_Outputs/"
# L = 1
for(L in 1){
  league_acro = league_acros[L]
  
  start_year = 2020
  end_year = 2025
  seasons_strvec = generate_season_string(start_year, end_year)
  
  season_seq = 1
  for (season_seq in 1:5){
    season = seasons_strvec[season_seq]
    
    file.rename(paste0("Data/Predictions/",league_acro, season,"_SAS_insample.RData"), 
                paste0("Data/Predictions/",league_acro, season,"_SAS_insample.rds"))
    # 
    # load(file = paste0("Data/Predictions/",league_acro, season,"_Pois_insample.RData"))
    # saveRDS(predictions, file = paste0("Data/Predictions/",league_acro, season,"_Pois_insample.RData"))
  }
  
  
  
  # if (league_acro == 'LC' && season == '1920'){next}
  # 
  # # saveRDS(predictions, file = paste0("Data/Predictions/",league_acro, season,"_Pois_oos.RData"))
  # # saveRDS(predictions, file = paste0("Data/Predictions/",league_acro, season,"_SAS_oos.RData"))
  # # saveRDS(predictions, file = paste0("Data/Predictions/",league_acro, season,"_CMP_oos.RData"))
  # # 
  # # saveRDS(predictions, file = paste0("Data/Predictions/",league_acro, season,"_Pois_insample.RData"))
  # # saveRDS(predictions, file = paste0("Data/Predictions/",league_acro, season,"_SAS_insample.RData"))
  # # saveRDS(predictions, file = paste0("Data/Predictions/",league_acro, season,"_CMP_insample.RData"))
  
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
    merge_df = read_historics(season, league_acro)
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
      hg = min(as.numeric(merge_df[game_id, 'FTHG']), 5)
      ag = min(as.numeric(merge_df[game_id, 'FTAG']), 5)
      goal_diff_obs = hg - ag
      
      cum_gd_obs = as.integer(gd_outcomes <= goal_diff_obs)
      
      pois_goaldiff = compute_goaldiff_probs(truncate_forecast(pois_pred[[game_id]], 2)/100)
      sas_goaldiff = compute_goaldiff_probs(truncate_forecast(sas_pred[[game_id]], 2)/100)
      cmp_goaldiff = compute_goaldiff_probs(truncate_forecast(cmp_pred[[game_id]], 2)/100)
      
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
    merge_df = read_historics(season, league_acro)
    
    pois_pred = readRDS(file = paste0("Data/Predictions/",league_acro, season,"_Pois_", regime, ".rds"))$list_preds
    
    sas_pred = readRDS(file = paste0("Data/Predictions/",league_acro, season,"_SAS_", regime, ".rds"))$list_preds
    
    cmp_pred = readRDS(file = paste0("Data/Predictions/",league_acro, season,"_CMP_", regime, ".rds"))$list_preds

    
    tot_games = nrow(merge_df)
    
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
      hg = min(as.numeric(merge_df[game_id, 'FTHG']), 7)
      ag = min(as.numeric(merge_df[game_id, 'FTAG']), 7)
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


league_acros = c('PL', 'SA', 'LL', 'LC', 'BL', 'WSL')
L = 1
league_acro = league_acros[L]
game_start = 191
seas_start = 2020
seas_end = 2025
p_id = ''
threshold = 0.49
filter_code = 0

outcome_oos = compute_outcome_eval(league_acro = league_acro, seas_start, seas_end, 'oos',
                                   filter = filter_code, threshold = threshold)
round(outcome_oos,3)

overunder_oos = compute_overunder_eval(league_acro = league_acro, seas_start, seas_end, overunder = 2.5, 'oos', 
                                       filter = filter_code, threshold = threshold)
round(overunder_oos,3)

gd_oos = compute_goaldiff_eval(league_acro = league_acro, seas_start, seas_end, 'oos', 
                               filter = filter_code, threshold = threshold)
round(gd_oos,3)

out_table = round(cbind(outcome_oos[,1:5], overunder_oos[,1:5], gd_oos[,1:5]),4)
out_table = round(rbind(outcome_oos[4:6,1:5], 
                        overunder_oos[4:6,1:5], 
                        gd_oos[4:6,1:5]),3)

out_table = round(rbind(outcome_oos[1:3,1:5], 
                        overunder_oos[1:3,1:5], 
                        gd_oos[1:3,1:5]),3)





################################
############# WAIC #############
################################

compute_Z <- function(mu, nu, j_max = 100, tol = 1e-20) {
  Z <- 0
  for (j in 0:j_max) {
    term <- (mu^(nu * j)) / (factorial(j)^nu)
    Z <- Z + term
    if (term < tol) {
      # print(j)
      break
    }
  }
  return(Z)
}

# Create grid of mu and nu
mu_vals <- seq(0.01, 10, by = 0.01)
nu_vals <- seq(0.01, 10, by = 0.01)

# Initialize lookup table
Z_lookup <- matrix(NA, nrow = length(mu_vals), ncol = length(nu_vals),
                   dimnames = list(paste0("mu=", mu_vals), paste0("nu=", nu_vals)))

# Fill table
for (i in seq_along(mu_vals)) {
  for (j in seq_along(nu_vals)) {
    Z_lookup[i, j] <- compute_Z(mu_vals[i], nu_vals[j])
  }
}

eval_cmp_llhood_cmp_vec = function(x, mu, nu){
  
  row <- nearest_index(mu, mu_vals)
  col <- nearest_index(nu, nu_vals)
  
  Z_vals <- Z_lookup[cbind(row, col)]
  
  qf <- nu * (x * log(mu) - lgamma(x + 1))
  return(qf - log(Z_vals))
}

nearest_index <- function(x, grid) {
  idx <- findInterval(x, grid)
  idx[idx == 0] <- 1
  too_high <- idx == length(grid)
  idx[too_high] <- length(grid)
  
  # adjust to whichever endpoint is closer
  left  <- grid[idx]
  right <- grid[pmin(idx + 1, length(grid))]
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
      
      lppd_home <- log_sum_exp(log_llhood_home)
      lppd_away <- log_sum_exp(log_llhood_away)
      
      
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
      
      lppd_home <- log_sum_exp(log_llhood_home)
      lppd_away <- log_sum_exp(log_llhood_away)
      
      
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

log_sum_exp <- function(x) {
  m <- max(x)
  m + log(mean(exp(x - m)))
}



n_samples = 5000

league_acros = c('PL', 'SA', 'LL', 'LC', 'BL', 'WSL')

leagues = c('Premier', 'SerieA', 'Liga', 'Ligue', 'Bundes', 'WSL')

L = 1
league = leagues[L]
league_acro = league_acros[L]

seasons_strvec = generate_season_string(2015, 2020)

n_seeds = 1

IC_table = data.frame(matrix(NA_real_, nrow = 3, length(seasons_strvec)))
lppd_table = data.frame(matrix(NA_real_, nrow = 3, length(seasons_strvec)))
pwaic_table = data.frame(matrix(NA_real_, nrow = 3, length(seasons_strvec)))
colnames(IC_table) = seasons_strvec
colnames(lppd_table) = seasons_strvec
colnames(pwaic_table) = seasons_strvec
for (season in seasons_strvec){
  if (L == 4 && season == '1920'){next}
  print(season)
  X = read_data(season, league)
  X1 = X[[1]]
  X2 = X[[2]]
  
  for (n_run in 1){
    MH_P = readRDS(file = paste0("Data//MCMC_Outputs//Pois_FullLeague//",league_acro,"_",season,"_Pois.rds"))
    MH_SAS = readRDS(file = paste0("Data//MCMC_Outputs//SAS_FullLeague//",league_acro,"_",season,"_SAS.rds"))
    MH_CMP = readRDS(file = paste0("Data//MCMC_Outputs//CMP_FullLeague//",league_acro,"_",season,"_CMP.rds"))
    
    stopifnot(all(MH_P$team_names == rownames(X1)))
    stopifnot(all(MH_SAS$team_names == rownames(X1)))
    stopifnot(MH_P$att_fix_idx == MH_SAS$att_fix_idx)

    P_WAIC = evaluate_WAIC(X, MH_P, n_samples)
    SAS_WAIC = evaluate_WAIC(X, MH_SAS, n_samples)
    CMP_WAIC = evaluate_WAIC(X, MH_CMP, n_samples)
    
    lppd_table[1, season] = P_WAIC$lppd
    lppd_table[2, season] = SAS_WAIC$lppd
    lppd_table[3, season] = CMP_WAIC$lppd
    
    pwaic_table[1, season] = P_WAIC$p_waic
    pwaic_table[2, season] = SAS_WAIC$p_waic
    pwaic_table[3, season] = CMP_WAIC$p_waic
  }
}

IC_table = -2*(lppd_table - pwaic_table)
IC_table$Sums = rowSums(IC_table, na.rm = TRUE)
lppd_table$Sums = rowSums(lppd_table, na.rm = TRUE)
pwaic_table$Sums = rowSums(pwaic_table, na.rm = TRUE)

# load(file = paste0("Data//MH_Results//CMP_FullLeague//",league_acro,"_",season,"_Run",n_run,"_CMP",".RData"))
IC_table
lppd_table
pwaic_table
out = cbind(lppd_table, pwaic_table, IC_table)
kable(round(IC_table,1), 'latex')
rowMeans(lppd_table - pwaic_table)






































