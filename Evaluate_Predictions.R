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
    
    file.rename(paste0("Data/Predictions/",league_acro, season,"_Pois_oos.RData"), 
                paste0("Data/Predictions/",league_acro, season,"_Pois_oos.rds"))
    
    file.rename(paste0("Data/Predictions/",league_acro, season,"_SAS_oos.RData"), 
                paste0("Data/Predictions/",league_acro, season,"_SAS_oos.rds"))
    
    file.rename(paste0("Data/Predictions/",league_acro, season,"_CMP_oos.RData"), 
                paste0("Data/Predictions/",league_acro, season,"_CMP_oos.rds"))
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
                        gd_oos[4:6,1:5]),4)

out_table
out_table = round(rbind(outcome_oos[1:3,1:5], 
                        overunder_oos[1:3,1:5], 
                        gd_oos[1:3,1:5]),3)





################################
############# WAIC #############
################################





n_samples = 5000

league_acros = c('PL', 'SA', 'LL', 'LC', 'BL', 'WSL')

leagues = c('Premier', 'SerieA', 'Liga', 'Ligue', 'Bundes', 'WSL')

L = 1
league = leagues[L]
league_acro = league_acros[L]

seasons_strvec = generate_season_string(2015, 2025)

n_seeds = 1

set.seed(1)
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






































