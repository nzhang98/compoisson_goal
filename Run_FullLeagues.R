source('utils.R')
source('MH_Poisson.R')
source('MH_CMP_SAS.R')
source('MH_CMP_Full.R')

start_year = 2020
end_year = 2025
seasons_strvec = generate_season_string(start_year, end_year)

league_acros = c('PL', 'SA', 'LL', 'LC', 'BL')

leagues = c('Premier', 'SerieA', 'Liga', 'Ligue', 'Bundes')

# E.g. L = 1, Premier League
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
      
      mcmc_out_dir = "Data/MCMC_Outputs/"
      
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
      eta_sd_prior = 1
      
      rho = 0.85
      
      p_alpha_prior = 1
      p_beta_prior = 1
      
      att_0 = rep(0, N)
      def_0 = rep(-0, N)
      eta_0 = rep(0, N)
      home_0 = 0
      Z_0 = rep(1, N)
      p_0 = rep(0.5,N)
      print_freq = 10000
      
      X_mid = matrix(1L, N, N) #Fail safe
      diag(X_mid) = 0L
      
      iter = 250000
      
      fixed_i = N #Use the last team as the fixed anchor team for the att and def parameters
    } 

    set.seed(1)
    
    ### Poisson Model
    
    MH_P = MH_Pois(X1, X2, att_0, def_0, home_0,
                   X_mid, iter,
                   sd_prop_att, att_mean_prior, att_sd_prior, 
                   sd_prop_def, def_mean_prior, def_sd_prior,
                   sd_prop_home, home_mean_prior, home_sd_prior,
                   verbosity = verbosity, print_by = print_freq, fix_idx = fixed_i, 
                   league_acro = league_acro, season = season)
    
    MH_P = thin_mcmc(MH_P, 5) # Optional: thin MCMC, keeping 1 every 5 samples to save storage space
    
    saveRDS(MH_P, paste0(mcmc_out_dir, "Pois_FullLeague/", league_acro, "_", season, "_Pois.rds"))
    
    ### CMP-SAS Model
    MH_SAS = MH_CMP_SAS(X1, X2, att_0, def_0, home_0, Z_0, p_0, eta_0, 
                        X_mid, iter,
                        sd_prop_att, att_mean_prior, att_sd_prior, 
                        sd_prop_def, def_mean_prior, def_sd_prior,
                        sd_prop_home, home_mean_prior, home_sd_prior,
                        sd_prop_eta, eta_mean_prior, eta_sd_prior, rho,
                        p_alpha_prior, p_beta_prior, 
                        verbosity = verbosity, print_by = print_freq, fix_idx = fixed_i, 
                        league_acro = league_acro, season = season)
    
    MH_SAS = thin_mcmc(MH_SAS, 5) 
    
    saveRDS(MH_SAS, paste0(mcmc_out_dir, "SAS_FullLeague/", league_acro, "_", season, "_SAS.rds"))
    
    ### CMP-Full model
    
    MH_CMP = MH_CMP_Full(X1, X2, att_0, def_0, home_0, Z_0, p_0, eta_0, 
                         X_mid, iter,
                         sd_prop_att, att_mean_prior, att_sd_prior, 
                         sd_prop_def, def_mean_prior, def_sd_prior,
                         sd_prop_home, home_mean_prior, home_sd_prior,
                         sd_prop_eta, eta_mean_prior, eta_sd_prior, rho,
                         p_alpha_prior, p_beta_prior, 
                         verbosity = verbosity, print_by = print_freq, fix_idx = fixed_i, 
                         league_acro = league_acro, season = season)
    
    MH_CMP = thin_mcmc(MH_CMP, 5)
    
    saveRDS(MH_CMP, paste0(mcmc_out_dir, "CMP_FullLeague/", league_acro, "_", season, "_CMP.rds"))
  }
}






