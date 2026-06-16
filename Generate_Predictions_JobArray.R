source('utils.R')

create_heatmatrix_fast <- function(home, away) {
  nmax <- 8
  i <- pmin(home, nmax - 1) + 1
  j <- pmin(away, nmax - 1) + 1
  
  heatmatrix <- matrix(0, nrow = nmax, ncol = nmax)
  tab <- table(factor(i, levels = 1:nmax), factor(j, levels = 1:nmax))
  heatmatrix[as.matrix(expand.grid(1:nmax, 1:nmax)[, 2:1])] <- as.vector(tab)
  
  m.hmat <- reshape2::melt(heatmatrix)
  m.hmat[, 1:2] <- m.hmat[, 1:2] - 1
  return(list(m.hmat, heatmatrix))
}


sim_halfleague_results_pois = function(season, league_acro, m = 10, n_samples = 5000, games_tot = 37){
  df_hist = read_historics(season, league_acro)
  if (league_acro == 'LC' && season == '1920'){next}
  team_names = sort(unique(df_hist[,"HomeTeam"]))
  N = length(team_names)
  team_codes = setNames(1:N, team_names)
  tot_matches = nrow(df_hist)
  
  list_preds = vector("list", length = tot_matches)
  if(tot_matches == 380){
    games_seq = seq(190, 370, by = 10)
    n_game_week = 10}
  if(tot_matches == 306){
    games_seq = seq(153, 297, by = 9)
    n_game_week = 9}
  # games_seq = 10*seq(19,games_tot, by=1)
  for (n_games in games_seq){
    load(file = paste0("Data//MH_Results//Pois_MidLeague//",league_acro,"_",season,"_n",n_games,"_Pois",".RData"))
    
    MH_object = MH_P
    
    t = MH_object$iterations
    
    ids = get_sample_ids(n_samples, t, t/5)
    
    att_samples = MH_object$att_post[ids,]
    def_samples = MH_object$def_post[ids,]
    eta_samples = MH_object$eta_post[ids]
    home_samples = MH_object$home_post[ids]
    
    for (game_id in (n_games+1):(n_games+n_game_week)){
      if (game_id > tot_matches){break}
      # game_id = 189
      i = team_codes[df_hist[game_id, 'HomeTeam']]
      j = team_codes[df_hist[game_id, 'AwayTeam']]
      
      print(paste0(i, ' vs ', j, '; game-id: ', game_id))
      
      mu_home = exp(att_samples[,i] + def_samples[,j] + home_samples + eta_samples)
      mu_away = exp(att_samples[,j] + def_samples[,i] + eta_samples)
      
      home_goals = matrix(NA_real_, nrow = n_samples, ncol = m)
      away_goals = matrix(NA_real_, nrow = n_samples, ncol = m)
      for (s in 1:n_samples){
        home_goals[s,] = rpois(m, mu_home[s])
        away_goals[s,] = rpois(m, mu_away[s])
      }
      
      list_preds[[game_id]] = create_heatmatrix_fast(c(home_goals), c(away_goals))[[2]]/(m*n_samples/100)  
    }
  }
  return(list_preds)
}


sim_halfleague_results_sas = function(season, league_acro, m = 10, n_samples = 5000, games_tot = 37, ext = ''){
  df_hist = read_historics(season, league_acro)
  if (league_acro == 'LC' && season == '1920'){next}
  team_names = sort(unique(df_hist[,"HomeTeam"]))
  N = length(team_names)
  team_codes = setNames(1:N, team_names)
  tot_matches = nrow(df_hist)
  
  list_preds = vector("list", length = tot_matches)
  if(tot_matches == 380){
    games_seq = seq(190, 370, by = 10)
    n_game_week = 10}
  if(tot_matches == 306){
    games_seq = seq(153, 297, by = 9)
    n_game_week = 9}
  for (n_games in games_seq){
    MH_object = readRDS(paste0(mcmc_out_dir, "SAS_MidLeague/", league_acro, "_", season,"_n",n_games,"_SAS.rds"))
    
    t = MH_object$iterations
    
    ids = get_sample_ids(n_samples, t, t/5)
    
    att_samples = MH_object$att_post[ids,]
    def_samples = MH_object$def_post[ids,]
    nu_samples = MH_object$nu_post[ids,]
    home_samples = MH_object$home_post[ids]
    
    for (game_id in (n_games+1):(n_games+n_game_week)){
      if (game_id > tot_matches){break}
      # game_id = 189
      i = team_codes[df_hist[game_id, 'HomeTeam']]
      j = team_codes[df_hist[game_id, 'AwayTeam']]
      
      print(paste0(i, ' vs ', j, '; game-id: ', game_id))
      
      mu_home = exp(att_samples[,i] + def_samples[,j] + home_samples)
      nu_home = nu_samples[,i]
      mu_away = exp(att_samples[,j] + def_samples[,i])
      nu_away = nu_samples[,j]
      
      home_goals = matrix(NA_real_, nrow = n_samples, ncol = m)
      away_goals = matrix(NA_real_, nrow = n_samples, ncol = m)
      for (s in 1:n_samples){
        home_goals[s,] = rejection_sampler_draws(m, mu_home[s], nu_home[s])
        away_goals[s,] = rejection_sampler_draws(m, mu_away[s], nu_away[s])
      }
      
      list_preds[[game_id]] = create_heatmatrix_fast(c(home_goals), c(away_goals))[[2]]/(m*n_samples/100)  
    }
  }
  return(list_preds)
}

sim_halfleague_results_cmp = function(season, league_acro, m = 10, n_samples = 5000, games_tot = 37){
  df_hist = read_historics(season, league_acro)
  if (league_acro == 'LC' && season == '1920'){next}
  team_names = sort(unique(df_hist[,"HomeTeam"]))
  N = length(team_names)
  team_codes = setNames(1:N, team_names)
  tot_matches = nrow(df_hist)
  
  list_preds = vector("list", length = tot_matches)
  if(tot_matches == 380){
    games_seq = seq(190, 370, by = 10)
    n_game_week = 10}
  if(tot_matches == 306){
    games_seq = seq(153, 297, by = 9)
    n_game_week = 9}
  for (n_games in games_seq){
    load(file = paste0("Data//MH_Results//CMP_MidLeague//",league_acro,"_",season,"_n",n_games,"_CMP",".RData"))
    
    # MH_object = thin_mcmc(MH_SAS, 10)
    
    MH_object = MH_CMP
    
    t = MH_object$iterations
    
    ids = get_sample_ids(n_samples, t, t/5)
    
    att_samples = MH_object$att_post[ids,]
    def_samples = MH_object$def_post[ids,]
    nu_samples = MH_object$nu_post[ids,]
    eta_samples = MH_object$eta_post[ids]
    home_samples = MH_object$home_post[ids]
    
    for (game_id in (n_games+1):(n_games+n_game_week)){
      if (game_id > tot_matches){break}
      # game_id = 189
      i = team_codes[df_hist[game_id, 'HomeTeam']]
      j = team_codes[df_hist[game_id, 'AwayTeam']]
      
      print(paste0(i, ' vs ', j, '; game-id: ', game_id))
      
      mu_home = exp(att_samples[,i] + def_samples[,j] + home_samples)
      nu_home = nu_samples[,i]
      mu_away = exp(att_samples[,j] + def_samples[,i])
      nu_away = nu_samples[,j]
      
      home_goals = matrix(NA_real_, nrow = n_samples, ncol = m)
      away_goals = matrix(NA_real_, nrow = n_samples, ncol = m)
      for (s in 1:n_samples){
        home_goals[s,] = rejection_sampler_draws(m, mu_home[s], nu_home[s])
        away_goals[s,] = rejection_sampler_draws(m, mu_away[s], nu_away[s])
      }
      
      list_preds[[game_id]] = create_heatmatrix_fast(c(home_goals), c(away_goals))[[2]]/(m*n_samples/100)  
    }
  }
  return(list_preds)
}


sim_insample_results_pois = function(season, league_acro, m = 10, n_samples = 5000){
  df_hist = read_historics(season, league_acro)
  team_names = sort(unique(df_hist[,"HomeTeam"]))
  N = length(team_names)
  team_codes = setNames(1:N, team_names)
  tot_matches = nrow(df_hist)
  
  list_preds = vector("list", length = tot_matches)
  # games_seq = 10*seq(19,games_tot, by=1)
  for (game_id in 1:tot_matches){
    load(file = paste0("Data//MH_Results//Pois_FullLeague//",league_acro,"_",season,"_Run1_Pois_CC1",".RData"))
    
    MH_object = MH_P
    
    t = MH_object$iterations
    
    ids = get_sample_ids(n_samples, t, t/5)
    
    att_samples = MH_object$att_post[ids,]
    def_samples = MH_object$def_post[ids,]
    eta_samples = MH_object$eta_post[ids]
    home_samples = MH_object$home_post[ids]
    
    if (game_id > tot_matches){
      print("game_id exceeds number of matches in history")
      break
    }
    
    i = team_codes[df_hist[game_id, 'HomeTeam']]
    j = team_codes[df_hist[game_id, 'AwayTeam']]
    
    print(paste0(i, ' vs ', j, '; game-id: ', game_id))
    
    mu_home = exp(att_samples[,i] + def_samples[,j] + home_samples + eta_samples)
    mu_away = exp(att_samples[,j] + def_samples[,i] + eta_samples)
    
    home_goals = matrix(NA_real_, nrow = n_samples, ncol = m)
    away_goals = matrix(NA_real_, nrow = n_samples, ncol = m)
    for (s in 1:n_samples){
      home_goals[s,] = rpois(m, mu_home[s])
      away_goals[s,] = rpois(m, mu_away[s])
    }
    
    list_preds[[game_id]] = create_heatmatrix_fast(c(home_goals), c(away_goals))[[2]]/(m*n_samples/100)  
    
  }
  return(list_preds)
}

sim_insample_results_sas = function(season, league_acro, m = 10, n_samples = 5000){
  df_hist = read_historics(season, league_acro)
  team_names = sort(unique(df_hist[,"HomeTeam"]))
  N = length(team_names)
  team_codes = setNames(1:N, team_names)
  tot_matches = nrow(df_hist)
  
  list_preds = vector("list", length = tot_matches)
  for (game_id in 1:tot_matches){
    load(file = paste0("Data//MH_Results//SAS_FullLeague//",league_acro,"_",season,"_Run1_SAS_N_Final",".RData"))
    
    MH_object = MH_SAS
    
    t = MH_object$iterations
    
    ids = get_sample_ids(n_samples, t, t/5)
    
    att_samples = MH_object$att_post[ids,]
    def_samples = MH_object$def_post[ids,]
    nu_samples = MH_object$nu_post[ids,]
    eta_samples = MH_object$eta_post[ids]
    home_samples = MH_object$home_post[ids]
    
    
    if (game_id > tot_matches){
      print("game_id exceeds number of matches in history")
      break
    }
    
    i = team_codes[df_hist[game_id, 'HomeTeam']]
    j = team_codes[df_hist[game_id, 'AwayTeam']]
    
    print(paste0(i, ' vs ', j, '; game-id: ', game_id))
    
    mu_home = exp(att_samples[,i] + def_samples[,j] + home_samples + eta_samples)
    nu_home = nu_samples[,i]
    mu_away = exp(att_samples[,j] + def_samples[,i] + eta_samples)
    nu_away = nu_samples[,j]
    
    home_goals = matrix(NA_real_, nrow = n_samples, ncol = m)
    away_goals = matrix(NA_real_, nrow = n_samples, ncol = m)
    for (s in 1:n_samples){
      home_goals[s,] = rejection_sampler_draws(m, mu_home[s], nu_home[s])
      away_goals[s,] = rejection_sampler_draws(m, mu_away[s], nu_away[s])
    }
    
    list_preds[[game_id]] = create_heatmatrix_fast(c(home_goals), c(away_goals))[[2]]/(m*n_samples/100)  
    
  }
  return(list_preds)
}

sim_insample_results_cmp = function(season, league_acro, m = 10, n_samples = 5000){
  df_hist = read_historics(season, league_acro)
  team_names = sort(unique(df_hist[,"HomeTeam"]))
  N = length(team_names)
  team_codes = setNames(1:N, team_names)
  tot_matches = nrow(df_hist)
  
  list_preds = vector("list", length = tot_matches)
  for (game_id in 1:tot_matches){
    load(file = paste0("Data//MH_Results//CMP_FullLeague//",league_acro,"_",season,"_Run1_CMP",".RData"))
    
    MH_object = MH_CMP
    
    t = MH_object$iterations
    
    ids = get_sample_ids(n_samples, t, t/5)
    
    att_samples = MH_object$att_post[ids,]
    def_samples = MH_object$def_post[ids,]
    nu_samples = MH_object$nu_post[ids,]
    eta_samples = MH_object$eta_post[ids]
    home_samples = MH_object$home_post[ids]
    
    
    if (game_id > tot_matches){
      print("game_id exceeds number of matches in history")
      break
    }
    
    i = team_codes[df_hist[game_id, 'HomeTeam']]
    j = team_codes[df_hist[game_id, 'AwayTeam']]
    
    print(paste0(i, ' vs ', j, '; game-id: ', game_id))
    
    mu_home = exp(att_samples[,i] + def_samples[,j] + home_samples + eta_samples)
    nu_home = nu_samples[,i]
    mu_away = exp(att_samples[,j] + def_samples[,i] + eta_samples)
    nu_away = nu_samples[,j]
    
    home_goals = matrix(NA_real_, nrow = n_samples, ncol = m)
    away_goals = matrix(NA_real_, nrow = n_samples, ncol = m)
    for (s in 1:n_samples){
      home_goals[s,] = rejection_sampler_draws(m, mu_home[s], nu_home[s])
      away_goals[s,] = rejection_sampler_draws(m, mu_away[s], nu_away[s])
    }
    
    list_preds[[game_id]] = create_heatmatrix_fast(c(home_goals), c(away_goals))[[2]]/(m*n_samples/100)  
    
  }
  return(list_preds)
}

league_acros = c('PL', 'SA', 'LL', 'LC', 'BL', 'WSL')
mcmc_out_dir = "Data/MCMC_Outputs/"
# L = 1
for(L in 1){
  league_acro = league_acros[L]
  
  start_year = 2023
  end_year = 2024
  seasons_strvec = generate_season_string(start_year, end_year)
  
  # args <- commandArgs(trailingOnly = TRUE)
  # season_seq <- as.numeric(args[1])
  season_seq = 1
  season = seasons_strvec[season_seq]
  
  m = 50
  n_samples = 5000
  
  games_tot = 37

  if (league_acro == 'LC' && season == '1920'){next}
  
  list_preds = sim_halfleague_results_sas(season, league_acro, m = m, n_samples = n_samples, games_tot, ext = '')
  #list_preds = sim_halfleague_results_pois(season, league_acro, m = m, n_samples = n_samples, games_tot)
  #list_preds = sim_halfleague_results_cmp(season, league_acro, m = m, n_samples = n_samples, games_tot)

  
  #list_preds = sim_insample_results_pois(season, league_acro, m = m, n_samples = n_samples)
  #list_preds = sim_insample_results_sas(season, league_acro, m = m, n_samples = n_samples)
  #list_preds = sim_insample_results_cmp(season, league_acro, m = m, n_samples = n_samples)

  
  predictions = list(
    #dist = 'Pois',
    dist = 'SAS',
    #ist = 'CMP',
    replications = m,
    post_samples = n_samples,
    out_of_sample = TRUE,
    #out_of_sample = FALSE,
    season = season,
    league = league_acro,
    list_preds = list_preds)
  #save(predictions, file = paste0("Data/Predictions/",league_acro, season,"_Pois_oos.RData"))
  saveRDS(predictions, file = paste0("Data/Predictions/",league_acro, season,"_SAS_oos.RData"))
  #save(predictions, file = paste0("Data/Predictions/",league_acro, season,"_SAS_oos.RData"))
  #save(predictions, file = paste0("Data/Predictions/",league_acro, season,"_Pois_insample.RData"))
  #save(predictions, file = paste0("Data/Predictions/",league_acro, season,"_SAS_insample.RData"))
  #save(predictions, file = paste0("Data/Predictions/",league_acro, season,"_CMP_oos.RData"))

}

