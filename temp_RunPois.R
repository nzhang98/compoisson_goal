source('utils.R')
source('MH_Poisson.R')

start_year = 2015
end_year = 2025
seasons_strvec = generate_season_string(start_year, end_year)

n_seeds = 1

league_acros = c('PL', 'SA', 'LL', 'LC', 'BL', 'WSL')

leagues = c('Premier', 'SerieA', 'Liga', 'Ligue', 'Bundes', 'WSL')

for (L in 1){
  league = leagues[L]
  league_acro = league_acros[L]
  print(league)
  print(league_acro)
  
  args <- commandArgs(trailingOnly = TRUE)
  season_seq <- as.numeric(args[1])
  
  season = seasons_strvec[season_seq]
  
  if (L == 4 && season == '1920'){next}
  
  cat("Running model for season =", season, "\n")
  
  X = read_data(season, league)
  X1 = X[[1]]
  X2 = X[[2]]
  
  N = dim(X1)[1]
  
  df_hist = read_historics(season, league_acro)
  team_names = sort(unique(df_hist[,"HomeTeam"]))
  rownames(X1) = team_names
  rownames(X2) = team_names
  team_codes = setNames(1:N, team_names)
  # print(team_names)
  # print(team_codes)
  
  { 
    verbosity = 3
    
    sd_prop_att = 0.1
    att_mean_prior = 0
    att_sd_prior = 10
    
    sd_prop_def = 0.04
    def_mean_prior = -0
    def_sd_prior = 10
    
    sd_prop_home = 0.08
    home_mean_prior = 0
    home_sd_prior = 10
    
    # sd_prop_eta = 0.4
    # eta_mean_prior = 0
    # eta_sd_prior = 1
    
    # rho = 0.85
    
    p_alpha_prior = 1
    p_beta_prior = 1
    
    att_0 = rep(0, N)
    def_0 = rep(-0, N)
    # eta_0 = rep(0, N)
    home_0 = 0
    # Z_0 = rep(1, N)
    # p_0 = rep(0.5,N)
    
    print_freq = 10000
    
    X_mid = matrix(1L, N, N)
    diag(X_mid) = 0L
    
    iter = 500000
    
    fixed_i = N
  } 
  
  for (n_run in 1){
    set.seed(n_run)
    
    MH_P = MH_Pois(X1, X2, att_0, def_0, home_0,
                        X_mid, iter,
                        sd_prop_att, att_mean_prior, att_sd_prior, 
                        sd_prop_def, def_mean_prior, def_sd_prior,
                        sd_prop_home, home_mean_prior, home_sd_prior,
                        verbosity = verbosity, print_by = print_freq, fix_idx = fixed_i, 
                        league_acro = league_acro, season = season)
    
    MH_P = thin_mcmc(MH_P, 5)
    mcmc_out_dir = "Data/MCMC_Outputs/"
    saveRDS(MH_SAS, paste0(mcmc_out_dir, "Pois_FullLeague/", league_acro, "_", season, "_Pois.rds"))
  }
}


source('utils.R')
source('MH_Poisson.R')

start_year = 2020
end_year = 2025
seasons_strvec = generate_season_string(start_year, end_year)

set.seed(1)

league_acros = c('PL', 'SA', 'LL', 'LC', 'BL')

leagues = c('Premier', 'SerieA', 'Liga', 'Ligue', 'Bundes')

for (L in 1){
  league = leagues[L]
  league_acro = league_acros[L]
  print(league)
  print(league_acro)
  
  for (season in seasons_strvec){
    
    if (L == 4 && season == '1920'){next}
    
    cat("Running model for season =", season, "\n")
    print(season)
    
    X = read_data(season, league)
    X1 = X[[1]]
    X2 = X[[2]]
    
    N = dim(X1)[1]
    
    df_hist = read_historics(season, league_acro)
    team_names = sort(unique(df_hist[,"HomeTeam"]))
    rownames(X1) = team_names
    rownames(X2) = team_names
    team_codes = setNames(1:N, team_names)
    
    { 
      verbosity = 3
      
      sd_prop_att = 0.1
      att_mean_prior = 0
      att_sd_prior = 10
      
      sd_prop_def = 0.04
      def_mean_prior = -0
      def_sd_prior = 10
      
      sd_prop_home = 0.08
      home_mean_prior = 0
      home_sd_prior = 10
      
      att_0 = rep(0, N)
      def_0 = rep(-0, N)
      home_0 = 0
      
      print_freq = 10000
      
      X_mid = matrix(1L, N, N)
      diag(X_mid) = 0L
      
      iter = 250000
      
      fixed_i = N
    } 
    
    games_seq = 10*seq(19,37, by=1)
    args <- commandArgs(trailingOnly = TRUE)
    games_seq_id <- as.numeric(args[1])
    
    n_games = games_seq[games_seq_id]
    
    tot_games = nrow(df_hist)
    #if(tot_games == 380){games_seq = seq(190, 370, by = 10)}
    #if(tot_games == 306){games_seq = seq(153, 288, by = 9)}
    
    print(n_games)
    X_mid = generate_X_mid(df_hist, n_games)
    
    MH_P = MH_Pois(X1, X2, att_0, def_0, home_0, 
                        X_mid, iter,
                        sd_prop_att, att_mean_prior, att_sd_prior, 
                        sd_prop_def, def_mean_prior, def_sd_prior,
                        sd_prop_home, home_mean_prior, home_sd_prior,
                        verbosity = verbosity, print_by = print_freq, fix_idx = fixed_i, 
                        league_acro = league_acro, season = season)
    
    MH_P = thin_mcmc(MH_P, 5)
    mcmc_out_dir = "Data/MCMC_Outputs/"
    saveRDS(MH_P, paste0(mcmc_out_dir, "Pois_MidLeague/", league_acro, "_", season,"_n",n_games,"_Pois.rds"))  
  }
}


source('utils.R')

source('MH_Poisson.R')

{ 
  N = 20
  verbosity = 3
  
  sd_prop_att = 0.1
  att_mean_prior = 0
  att_sd_prior = 10
  
  sd_prop_def = 0.04
  def_mean_prior = -0
  def_sd_prior = 10
  
  sd_prop_home = 0.08
  home_mean_prior = 0
  home_sd_prior = 10
  
  att_0 = rep(0, N)
  def_0 = rep(-0, N)
  home_0 = 0
  
  print_freq = 10000
  
  X_mid = matrix(1L, N, N)
  diag(X_mid) = 0L
  
  iter = 250000
  
  fixed_i = N
  
  league_acro = 'Sim'
}

args <- commandArgs(trailingOnly = TRUE)
nu_seq <- as.numeric(args[1])

nus_underdispersed = seq(1.2, 5, by = 0.4)
nus_overdispersed = seq(0.1, 0.90, by = 0.1)

nu_scalar = nus_overdispersed[nu_seq]
#nu_scalar = nus_underdispersed [nu_seq]

load(file = paste0("Data//Simulations//SAS_sim_nubase_",gsub("\\.", "_", as.character(nu_scalar)),".RData"))
n_seeds = length(simulations_list)
for (run in 1:n_seeds){
  print(paste0("Nu_scalar: ",nu_scalar, "  Run: ", run))
  set.seed(1)
  
  Xsim_list = simulations_list[[run]]
  X = Xsim_list$X_sim
  X1 = X[[1]]
  X2 = X[[2]]
  
  MH_P = MH_Pois(X1, X2, att_0, def_0, home_0, 
                 X_mid, iter,
                 sd_prop_att, att_mean_prior, att_sd_prior, 
                 sd_prop_def, def_mean_prior, def_sd_prior,
                 sd_prop_home, home_mean_prior, home_sd_prior,
                 verbosity = verbosity, print_by = print_freq, fix_idx = fixed_i, 
                 league_acro = league_acro, season = season)
  
  MH_P = thin_mcmc(MH_P, 5)
  mcmc_out_dir = "Data/MCMC_Outputs/"
  saveRDS(MH_P, paste0(mcmc_out_dir, "Simulations/SAS_SIM_Nu",gsub("\\.", "_", as.character(nu_scalar)),"_Run",run,"_Pois.rds"))  
}


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


sim_halfleague_results_pois = function(season, league_acro, m = 10, n_samples = 5000, out_dir = "Data/MCMC_Outputs/"){
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
    MH_object = readRDS(paste0(out_dir, "Pois_MidLeague/", league_acro, "_", season, "_n", n_games, "_Pois.rds"))
    # load(file = paste0("Data//MH_Results//Pois_MidLeague//",league_acro,"_",season,"_n",n_games,"_Pois",".RData"))
    
    # MH_object = MH_P
    
    t = MH_object$iterations
    
    ids = get_sample_ids(n_samples, t, t/5)
    
    att_samples = MH_object$att_post[ids,]
    def_samples = MH_object$def_post[ids,]
    home_samples = MH_object$home_post[ids]
    
    for (game_id in (n_games+1):(n_games+n_game_week)){
      if (game_id > tot_matches){break}
      # game_id = 189
      i = team_codes[df_hist[game_id, 'HomeTeam']]
      j = team_codes[df_hist[game_id, 'AwayTeam']]
      
      print(paste0(i, ' vs ', j, '; game-id: ', game_id))
      
      mu_home = exp(att_samples[,i] + def_samples[,j] + home_samples)
      mu_away = exp(att_samples[,j] + def_samples[,i])
      
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


sim_halfleague_results_sas = function(season, league_acro, m = 10, n_samples = 5000, ext = '', out_dir = "Data/MCMC_Outputs/"){
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
    MH_object = readRDS(paste0(out_dir, "SAS_MidLeague/", league_acro, "_", season,"_n",n_games,"_SAS.rds"))
    
    t = MH_object$iterations
    
    ids = get_sample_ids(n_samples, t, t/5)
    
    att_samples = MH_object$att_post[ids,]
    def_samples = MH_object$def_post[ids,]
    nu_samples = MH_object$nu_post[ids,]
    home_samples = MH_object$home_post[ids]
    
    for (game_id in (n_games+1):(n_games+n_game_week)){
      if (game_id > tot_matches){break}
      
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

sim_halfleague_results_cmp = function(season, league_acro, m = 10, n_samples = 5000, out_dir = "Data/MCMC_Outputs/"){
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
    MH_object = readRDS(paste0(out_dir, "CMP_MidLeague/", league_acro, "_", season,"_n",n_games,"_CMP.rds"))
    
    
    MH_object = MH_CMP
    
    t = MH_object$iterations
    
    ids = get_sample_ids(n_samples, t, t/5)
    
    att_samples = MH_object$att_post[ids,]
    def_samples = MH_object$def_post[ids,]
    nu_samples = MH_object$nu_post[ids,]
    home_samples = MH_object$home_post[ids]
    
    for (game_id in (n_games+1):(n_games+n_game_week)){
      if (game_id > tot_matches){break}

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


sim_insample_results_pois = function(season, league_acro, m = 10, n_samples = 5000, out_dir = "Data/MCMC_Outputs/"){
  df_hist = read_historics(season, league_acro)
  team_names = sort(unique(df_hist[,"HomeTeam"]))
  N = length(team_names)
  team_codes = setNames(1:N, team_names)
  tot_matches = nrow(df_hist)
  
  list_preds = vector("list", length = tot_matches)

  for (game_id in 1:tot_matches){
    MH_object = readRDS(paste0(out_dir, "Pois_FullLeague/", league_acro, "_", season, "_n", n_games, "_Pois.rds"))
    
    t = MH_object$iterations
    
    ids = get_sample_ids(n_samples, t, t/5)
    
    att_samples = MH_object$att_post[ids,]
    def_samples = MH_object$def_post[ids,]
    home_samples = MH_object$home_post[ids]
    
    if (game_id > tot_matches){
      print("game_id exceeds number of matches in history")
      break
    }
    
    i = team_codes[df_hist[game_id, 'HomeTeam']]
    j = team_codes[df_hist[game_id, 'AwayTeam']]
    
    print(paste0(i, ' vs ', j, '; game-id: ', game_id))
    
    mu_home = exp(att_samples[,i] + def_samples[,j] + home_samples)
    mu_away = exp(att_samples[,j] + def_samples[,i])
    
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

sim_insample_results_sas = function(season, league_acro, m = 10, n_samples = 5000, out_dir = "Data/MCMC_Outputs/"){
  df_hist = read_historics(season, league_acro)
  team_names = sort(unique(df_hist[,"HomeTeam"]))
  N = length(team_names)
  team_codes = setNames(1:N, team_names)
  tot_matches = nrow(df_hist)
  
  list_preds = vector("list", length = tot_matches)
  for (game_id in 1:tot_matches){
    MH_object = readRDS(paste0(out_dir, "SAS_FullLeague/", league_acro, "_", season, "_n", n_games, "_SAS.rds"))
    
    MH_object = MH_SAS
    
    t = MH_object$iterations
    
    ids = get_sample_ids(n_samples, t, t/5)
    
    att_samples = MH_object$att_post[ids,]
    def_samples = MH_object$def_post[ids,]
    nu_samples = MH_object$nu_post[ids,]
    home_samples = MH_object$home_post[ids]
    
    
    if (game_id > tot_matches){
      print("game_id exceeds number of matches in history")
      break
    }
    
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
  return(list_preds)
}

sim_insample_results_cmp = function(season, league_acro, m = 10, n_samples = 5000, out_dir = "Data/MCMC_Outputs/"){
  df_hist = read_historics(season, league_acro)
  team_names = sort(unique(df_hist[,"HomeTeam"]))
  N = length(team_names)
  team_codes = setNames(1:N, team_names)
  tot_matches = nrow(df_hist)
  
  list_preds = vector("list", length = tot_matches)
  for (game_id in 1:tot_matches){
    MH_object = readRDS(paste0(out_dir, "CMP_FullLeague/", league_acro, "_", season, "_n", n_games, "_CMP.rds"))
    
    MH_object = MH_CMP
    
    t = MH_object$iterations
    
    ids = get_sample_ids(n_samples, t, t/5)
    
    att_samples = MH_object$att_post[ids,]
    def_samples = MH_object$def_post[ids,]
    nu_samples = MH_object$nu_post[ids,]
    home_samples = MH_object$home_post[ids]
    
    
    if (game_id > tot_matches){
      print("game_id exceeds number of matches in history")
      break
    }
    
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
  return(list_preds)
}

league_acros = c('PL', 'SA', 'LL', 'LC', 'BL', 'WSL')
mcmc_out_dir = "Data/MCMC_Outputs/"
# L = 1
for(L in 1){
  league_acro = league_acros[L]
  
  start_year = 2020
  end_year = 2025
  seasons_strvec = generate_season_string(start_year, end_year)
  
  # args <- commandArgs(trailingOnly = TRUE)
  # season_seq <- as.numeric(args[1])
  season_seq = 1
  season = seasons_strvec[season_seq]
  
  m = 50
  n_samples = 5000
  
  if (league_acro == 'LC' && season == '1920'){next}
  
  #list_preds = sim_halfleague_results_sas(season, league_acro, m = m, n_samples = n_samples)
  list_preds = sim_halfleague_results_pois(season, league_acro, m = m, n_samples = n_samples)
  #list_preds = sim_halfleague_results_cmp(season, league_acro, m = m, n_samples = n_samples)
  
  
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
  saveRDS(predictions, file = paste0("Data/Predictions/",league_acro, season,"_Pois_oos.RData"))
  saveRDS(predictions, file = paste0("Data/Predictions/",league_acro, season,"_SAS_oos.RData"))
  saveRDS(predictions, file = paste0("Data/Predictions/",league_acro, season,"_CMP_oos.RData"))
  
  saveRDS(predictions, file = paste0("Data/Predictions/",league_acro, season,"_Pois_insample.RData"))
  saveRDS(predictions, file = paste0("Data/Predictions/",league_acro, season,"_SAS_insample.RData"))
  saveRDS(predictions, file = paste0("Data/Predictions/",league_acro, season,"_CMP_insample.RData"))
  
}


