### Generic functions
library(ggplot2)
library(reshape2)
library(purrr)
library(gridExtra)

Mode <- function(x) { #compute the mode of a discrete list x
  ux <- unique(x)
  ux[which.max(tabulate(match(x, ux)))]
}

### Read data functions

read_data = function(season, league = 'Premier'){
  # League list: "Premier", "SerieA", "Liga", "Ligue1", "Bundes" 
  Results<-read.csv(paste0("Data//Results//Result_",league,"_",season, ".csv"), header = TRUE, sep=",")
  
  rownames(Results)=Results$Home...Away
  Results=Results[,-1]
  
  library(stringi)
  N = nrow(Results)
  
  X1<-matrix(NA,N,N)
  X2<-matrix(NA,N,N)
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

get_available_seasons = function(league = 'PL'){
  if (league == 'PL'){
    return(list(list('9596', 4), list('9697', 14), list('9798', 3),  #3
                list('9899', 15),list('9900', 18), list('0001', 3),  #6
                list('0102', 12),list('0203', 17), list('0304', 20), #9
                list('0405', 16),list('0506', 16), list('0607', 18), #12
                list('0708', 7), list('0809', 18), list('0910', 19), #15
                list('1011', 18),list('1112', 20), list('1213', 11), #18
                list('1314', 7), list('1415', 13), list('1516', 2),  #21
                list('1617', 15),list('1718', 19), list('1819', 10), #24
                list('1920', 14),list('2021', 15), list('2122', 15), #27
                list('2223', 17),list('2324', 17)))
  }
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

geom_envelope_draws = function(n, mu, nu){
  p = (2*nu)/(2*mu*nu + 1 + nu)
  
  xm = floor(mu/((1-p)**(1/nu)))
  
  B = (1/p) * (mu**(nu*xm)) /
    ((1-p)**xm * (factorial(xm)**nu))
  
  samples = numeric(n)
  
  for(i in 1:n){
    while (TRUE){
      #attempts
      u_0 = runif(1)
      x = floor(log(u_0)/log(1-p))
      
      alpha = (mu**x/factorial(x))**nu/(p*B*(1-p)**x)
      if (is.nan(alpha)){
        next
      }
      
      u = runif(1)
      
      if (u <= alpha){
        samples[i] = x
        break
      }
    }
  }
  return(samples)
}

pois_envelope_draws = function(n, mu, nu){
  B = ((mu**floor(mu))/factorial(floor(mu)))**(nu-1)
  
  samples = numeric(n)
  
  for(i in 1:n){
    while(TRUE){
      x = rpois(1, mu)
      
      alpha = (mu**x/factorial(x))**nu/(B * (mu**x)/factorial(x))
      
      if (is.nan(alpha)){
        next
      }
      u = runif(1)
      
      if (u <= alpha){
        samples[i] = x
        break
      }
    }
  }
  return(samples)
  
}

rejection_sampler_draws = function(n, mu, nu){
  if (nu < 1){
    x = geom_envelope_draws(n, mu, nu)
  } else {
    x = pois_envelope_draws(n, mu, nu)
  }
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
  
  print(iters)
  
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
    print('home')
    print(length(temp))
    print(t_burn)
    print(MH_obj$iterations)
    c = 0
    for (t in 2:iters){
      if (temp[t] != temp[t-1]){c = c+1}
    }
    output = append(output, c('home' = (c/iters)))
    
    if(MH_obj$distr_type == 'CP_D'){
      for (i in 1:N){
        temp = nu[[i]]
        print(temp)
        print(length(temp))
        c = 0
        for (t in 2:iters){
          if (temp[t] != temp[t-1]){c = c+1}
        }
        output = append(output, setNames(c/iters, paste0('nu',i)))
      }
    } else {
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
  
  att_chains = retrieve_post_longlists(MH_obj, 'att')
  def_chains = retrieve_post_longlists(MH_obj, 'def')
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
    nu_chains = retrieve_post_longlists(MH_obj, 'nu')
    
    for (i in 1:N){
      nu_i = nu_chains[[i]]
      nu_vec = c(mean = mean(nu_i), sd = sd(nu_i),  quantile(nu_i, probs = c(0.025, 0.25, 0.5, 0.75, 0.975)))
      out_df = rbind(out_df, nu_vec)
      rownames(out_df)[nrow(out_df)] = paste0('nu_',i)
    }
  }
  
  return(out_df)
}

### Functions to generate matrices for mid-league inferences

generate_X_mid = function(df_hist, n_games){
  # Generate a binary matrix for that holds 1 for the 'n_games' that have been played, 
  # taken in chronological order from df_hist.
  N = length(unique(df_hist[,"HomeTeam"]))
  X_mid = matrix(0, N, N)
  
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