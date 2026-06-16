source('utils.R')
source('MH_Poisson.R')
source('MH_CMP_SAS.R')
source('MH_CMP_Full.R')


{ 
  N = 20
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
  
  print_freq = 1000
  
  X_mid = matrix(1L, N, N)
  diag(X_mid) = 0L
  
  iter = 250000
  
  fixed_i = N
  
  league_acro = 'Sim'
}

nus_underdispersed = seq(1.2, 4, by = 0.4)
nus_overdispersed = seq(0.1, 0.90, by = 0.1)

nus_full_list = c(nus_overdispersed, nus_underdispersed)

for (nu_scalar in nus_full_list){
  
  simulations_list = readRDS(file = paste0("Data//Simulations//SAS_sim_nubase_",gsub("\\.", "_", as.character(nu_scalar)),".rds"))
  
  n_seeds = length(simulations_list) # Number of different seeds generated in the simulations
  
  for (run in 1:n_seeds){
    print(paste0("Nu_scalar: ",nu_scalar, "  Run: ", run))
    set.seed(1)
    
    Xsim_list = simulations_list[[run]]
    X = Xsim_list$X_sim
    X1 = X[[1]]
    X2 = X[[2]]
    
    ### Poisson Model
    
    MH_P = MH_Pois(X1, X2, att_0, def_0, home_0, 
                   X_mid, iter,
                   sd_prop_att, att_mean_prior, att_sd_prior, 
                   sd_prop_def, def_mean_prior, def_sd_prior,
                   sd_prop_home, home_mean_prior, home_sd_prior,
                   verbosity = verbosity, print_by = print_freq, fix_idx = fixed_i, 
                   league_acro = league_acro)
    
    MH_P = thin_mcmc(MH_P, 5)

    saveRDS(MH_P, paste0(mcmc_out_dir, "Simulations/SAS_SIM_Nu",gsub("\\.", "_", as.character(nu_scalar)),"_Run",run,"_Pois.rds"))  
    
    ### CMP-SAS Model
    
    MH_SAS = MH_CMP_SAS(X1, X2, att_0, def_0, home_0, Z_0, p_0, eta_0, 
                        X_mid, iter,
                        sd_prop_att, att_mean_prior, att_sd_prior, 
                        sd_prop_def, def_mean_prior, def_sd_prior,
                        sd_prop_home, home_mean_prior, home_sd_prior,
                        sd_prop_eta, eta_mean_prior, eta_sd_prior, rho,
                        p_alpha_prior, p_beta_prior, 
                        verbosity = verbosity, print_by = print_freq, fix_idx = fixed_i, 
                        league_acro = league_acro)
    
    MH_SAS = thin_mcmc(MH_SAS, 5)
    
    saveRDS(MH_SAS, paste0(mcmc_out_dir, "Simulations/SAS_SIM_Nu",gsub("\\.", "_", as.character(nu_scalar)),"_Run",run,"_SAS.rds")) 
    
    ### CMP-Full Model
    
    MH_CMP = MH_CMP_Full(X1, X2, att_0, def_0, home_0, Z_0, p_0, eta_0, 
                         X_mid, iter,
                         sd_prop_att, att_mean_prior, att_sd_prior, 
                         sd_prop_def, def_mean_prior, def_sd_prior,
                         sd_prop_home, home_mean_prior, home_sd_prior,
                         sd_prop_eta, eta_mean_prior, eta_sd_prior, rho,
                         p_alpha_prior, p_beta_prior, 
                         verbosity = verbosity, print_by = print_freq, fix_idx = fixed_i, 
                         league_acro = league_acro)
    
    MH_CMP = thin_mcmc(MH_CMP, 5)

    saveRDS(MH_CMP, paste0(mcmc_out_dir, "Simulations/SAS_SIM_Nu",gsub("\\.", "_", as.character(nu_scalar)),"_Run",run,"_CMP.rds"))  
  }
}

