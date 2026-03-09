source('utils.R')
source('MH_CMP.R')

start_year = 2023
end_year = 2024
seasons_strvec = generate_season_string(start_year, end_year)

n_seeds = 1

league_acros = c('PL', 'SA', 'LL', 'LC', 'BL', 'WSL')

leagues = c('Premier', 'SerieA', 'Liga', 'Ligue', 'Bundes', 'WSL')

for (L in 1){
  league = leagues[L]
  league_acro = league_acros[L]
  print(league)
  print(league_acro)
  
  for (season in seasons_strvec){
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
      
      sd_prop_eta = 0.4
      eta_mean_prior = 0
      eta_sd_prior = 1.15
      
      rho = 0.85
      
      p_alpha_prior = 1
      p_beta_prior = 1
      
      att_0 = rep(0, N)
      def_0 = rep(-0, N)
      eta_0 = 0
      home_0 = 0
      nu_0 = rep(1, N)
      Z_0 = rep(1, N)
      p_0 = rep(0.5,N)
      print_freq = 1000
      
      X_mid = matrix(1L, N, N)
      diag(X_mid) = 0L
      
      iter = 50000
      
      fixed_i = N
    } 
    
    
    for (n_run in 1){
      set.seed(n_run)
      MH_SAS = MH_SAS_STZ(X1, X2, att_0, def_0, eta_0, home_0, Z_0, p_0, nu_0, X_mid, iter = iter,
                          sd_prop_att, att_mean_prior, att_sd_prior, sd_prop_def, def_mean_prior, def_sd_prior,
                          sd_prop_home, home_mean_prior, home_sd_prior, sd_prop_eta, eta_mean_prior, eta_sd_prior,
                          sd_prop_eta, eta_mean_prior, eta_sd_prior, rho, p_alpha_prior, p_beta_prior, verbosity, 
                          print_by = print_freq, fix_idx = fixed_i, league_acro = league_acro, season = season)
      # save(MH_SAS, file = paste0("Data//MH_Results//SAS_FullLeague//",league_acro,"_",season,"_Run",n_run,"_SAS_NewZ",".RData"))
    }
  }
}
