source('utils.R')

### Generate Data

att = seq(-1, 1, length.out = 20)
def = seq(1, -1, length.out = 20)
home = 0.2
sim_distr = 'CMP-SAS'

nus_underdispersed = seq(1.2, 4, by = 0.4)
for (nu_scalar in nus_underdispersed){
  nu_base = c(rep(1, 10), rep(nu_scalar, 10))
  simulations_list = list()
  
  n_seeds = 5
  
  for (rseed in 1:n_seeds){
    set.seed(rseed+10)
    nu = sample(nu_base)
    X_sim = generate_league_sim_cp_d(att, def, home, nu)
    
    Xsim_list = list(sim_distr_type = sim_distr,
                     att_true = att,
                     def_true = def,
                     home_true = home,
                     rseed = rseed,
                     X_sim = X_sim,
                     nu_true = nu)
    
    simulations_list[[rseed]] = Xsim_list
  }
  saveRDS(simulations_list, file = paste0("Data//Simulations//SAS_sim_nubase_",gsub("\\.", "_", as.character(nu_scalar)),".rds"))
}

nus_overdispersed = seq(0.20, 0.90, by = 0.1)
for (nu_scalar in nus_overdispersed){
  nu_base = c(rep(1, 10), rep(nu_scalar, 10))
  simulations_list = list()
  
  n_seeds = 5
  
  for (rseed in 1:n_seeds){
    set.seed(rseed+10)
    nu = sample(nu_base)
    X_sim = generate_league_sim_cp_d(att, def, home, nu)
    
    Xsim_list = list(sim_distr_type = sim_distr,
                     att_true = att,
                     def_true = def,
                     home_true = home,
                     rseed = rseed,
                     X_sim = X_sim,
                     nu_true = nu)
    
    simulations_list[[rseed]] = Xsim_list
  }
  saveRDS(simulations_list, file = paste0("Data//Simulations//SAS_sim_nubase_",gsub("\\.", "_", as.character(nu_scalar)),".rds"))
}

#### Varying N

# N_vec = c(6, 8, 10, 12, 14, 16, 18, 20, 22, 24, 26, 28, 30, 40, 50, 60, 80)
N_vec = c(12, 14, 16, 18, 20, 30, 40, 50, 60, 70)

for (N in N_vec){
  att = seq(-1, 1, length.out = N)
  def = seq(1, -1, length.out = N)
  home = 0.2
  sim_distr = 'CMP'
  n_disp = floor(N/3)+1
  nu_base = c(rep(0.6, n_disp), rep(2, n_disp), rep(1, N-2*n_disp))
  
  simulations_list = list()
  
  n_seeds = 10
  
  for (rseed in 1:n_seeds){
    set.seed(rseed+10)
    nu = sample(nu_base)
    X_sim = generate_league_sim_cp_d(att, def, home, nu)
    
    Xsim_list = list(sim_distr_type = sim_distr,
                     att_true = att,
                     def_true = def,
                     home_true = home,
                     rseed = rseed,
                     X_sim = X_sim,
                     nu_true = nu)
    
    simulations_list[[rseed]] = Xsim_list
  }
  saveRDS(simulations_list, file = paste0("Data//Simulations//SAS_sim_N",N,".rds"))
}
