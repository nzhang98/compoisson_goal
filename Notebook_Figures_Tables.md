# Figures and Tables - Notebook


## Notebook

This notebook reproduces all figures and tables in the paper

#### Preliminary

``` r
source('utils.R') # Contains all utility functions to support modelling and analysis, including functions to import/export data, manipulate the MH_posterior objects from the MCMC routines, and generate predictions
```

    Warning: package 'purrr' was built under R version 4.3.3


    Attaching package: 'dplyr'

    The following object is masked from 'package:gridExtra':

        combine

    The following objects are masked from 'package:stats':

        filter, lag

    The following objects are masked from 'package:base':

        intersect, setdiff, setequal, union

    Warning: package 'tidyr' was built under R version 4.3.3


    Attaching package: 'tidyr'

    The following object is masked from 'package:reshape2':

        smiths

``` r
library(patchwork)
mcmc_out_dir = "Data/MCMC_Outputs/"
```

##### Figure 1

``` r
X = read_data('2324', 'Premier')
X1 = X[[1]]
X2 = X[[2]]

# Remove diagonals
X1_nodiag = X1[upper.tri(X1) | lower.tri(X1)]
X2_nodiag = X2[upper.tri(X2) | lower.tri(X2)]

# Truncate values > 7
X1_nodiag = X1_nodiag[X1_nodiag <= 7]
X2_nodiag = X2_nodiag[X2_nodiag <= 7]

# Convert to factors for exact bins 0 to 7
X1_factor = factor(X1_nodiag, levels = 0:7)
X2_factor = factor(X2_nodiag, levels = 0:7)

# Calculate means for vertical lines
mean_X1 = mean(as.numeric(as.character(X1_factor)))
mean_X2 = mean(as.numeric(as.character(X2_factor)))
var_X1 = var(as.numeric(as.character(X1_factor)))
var_X2 = var(as.numeric(as.character(X2_factor)))

# Create individual plots

p1 = ggplot(data.frame(value = X1_factor), aes(x = value)) +
  geom_bar(fill = "#87CEEB", color = "black", width = 1) +
  geom_vline(xintercept = mean_X1 + 1, linetype = "solid", color = "#0072B2", linewidth = 1) +
  geom_vline(xintercept = var_X1 + 1, linetype = "dashed", color = "#8E44AD", linewidth = 0.8) +
  scale_x_discrete(drop = FALSE) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.1))) +
  labs(x = "Home Goals Scored", y = "Number of Matches") +
  theme_minimal() +
  theme(
    axis.title.x = element_text(size = 12),
    axis.title.y = element_text(size = 12),
    axis.text = element_text(size = 10),
    plot.margin = margin(t = 5, r = 10, b = 5, l = 5)
  )

p2 = ggplot(data.frame(value = X2_factor), aes(x = value)) +
  geom_bar(fill = "#87CEEB", color = "black", width = 1) +
  geom_vline(xintercept = mean_X2 + 1, linetype = "solid", color = "#0072B2", linewidth = 1) +
  geom_vline(xintercept = var_X2 + 1, linetype = "dashed", color = "#8E44AD", linewidth = 0.8) +
  scale_x_discrete(drop = FALSE) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.1))) +
  labs(x = "Away Goals Scored", y = NULL) +
  theme_minimal() +
  theme(
    axis.title.x = element_text(size = 12),
    axis.text = element_text(size = 10),
    axis.ticks.y = element_blank(),
    axis.text.y = element_blank(),
    plot.margin = margin(t = 5, r = 5, b = 5, l = 10)
  )

# Get max y-axis limit to sync y scales
max_y = max(ggplot_build(p1)$data[[1]]$count, ggplot_build(p2)$data[[1]]$count)

# Fix y-axis limits to be the same
p1 = p1 + coord_cartesian(ylim = c(0, max_y))
p2 = p2 + coord_cartesian(ylim = c(0, max_y))

# Combine with patchwork, add space between plots
combined_plot = p1 + p2 + plot_layout(ncol = 2, widths = c(1, 1)) & 
  theme(plot.margin = margin(10, 10, 10, 10)) 

print(combined_plot)
```

![](Notebook_Figures_Tables_files/figure-commonmark/fig1-1.png)

##### Figure 2

``` r
season = '2324'
league_acro = 'PL'
league = 'Premier'

X = read_data(season, league)
X1 = X[[1]]
X2 = X[[2]]

df_hist = read_historics(season, league_acro)
team_names = sort(unique(df_hist[,"HomeTeam"]))
N = length(team_names)
rownames(X1) = team_names
rownames(X2) = team_names
team_codes = setNames(1:N, team_names)

plots = list()
for (j in 1:N){
  goals = unname(append(X1[j,],X2[,j]))
  goals = goals[!is.na(goals)]
  
  if (mean(goals)/var(goals) > 1){disp_color = '#D55E00'}else{disp_color = '#009E73'}
  
  mean_goals = mean(goals)
  var_goals = var(goals)
  
  p = ggplot(data.frame(goals), aes(x = goals)) +
    coord_cartesian(xlim=c(0,7), ylim=c(0,15)) +
    geom_histogram(binwidth = 1, fill = disp_color, color = 'black') +
    geom_vline(xintercept = mean_goals, linetype = "solid", color = "#0072B2", linewidth = 0.8) +
    geom_vline(xintercept = var_goals, linetype = "dashed", color = "#8E44AD", linewidth = 0.6) +
    theme_minimal() +
    theme(axis.title.x = element_blank(), axis.title.y = element_blank()) +
    ggtitle(paste(team_names[j]))
  
  plots[[j]] = p
}

do.call("grid.arrange", c(plots, nrow = 4, ncol = 5, left = 'Number of Matches', bottom = 'Goals'))
```

![](Notebook_Figures_Tables_files/figure-commonmark/fig2-1.png)

##### Figure 4

``` r
compute_Z = function(mu, nu, j_max = 100) {
  j = 0:j_max
  terms = (mu^(nu * j)) / (factorial(j)^nu)
  sum(terms)
}

get_index = function(value, grid) {
  pmax(1, pmin(length(grid), round((value - grid[1]) / (grid[2] - grid[1])) + 1))
}

generate_cmp_contour_plot = function(x, mu_true, nu_true, 
                                     mu_limit_lower = 0.01, mu_limit_upper = 3, mu_by = 0.01,
                                     nu_limit_lower = 0.01, nu_limit_upper = 3, nu_by = 0.01){
  mu_vals = seq(mu_limit_lower, mu_limit_upper, by = mu_by)
  nu_vals = seq(nu_limit_lower, nu_limit_upper, by = nu_by)
  
  Z_lookup = outer(mu_vals, nu_vals, Vectorize(compute_Z))
  LL = matrix(0, nrow = length(mu_vals), ncol = length(nu_vals))

  x_sum = sum(lgamma(x + 1))
  
  for (i in seq_along(mu_vals)) {
    mu = mu_vals[i]
    
    for (j in seq_along(nu_vals)) {
      nu = nu_vals[j]
      
      i_idx = get_index(mu, mu_vals)
      j_idx = get_index(nu, nu_vals)
      
      Z_val = Z_lookup[i_idx, j_idx]
      
      LL[i, j] =
        nu * sum(x) * log(mu) - length(x) * log(Z_val) - nu * x_sum
    }
  }
  
  df = expand.grid(
    mu = mu_vals,
    nu = nu_vals
  )
  df$LL = as.vector(LL)

  z_max = max(LL, na.rm = TRUE)
  z_min = min(LL, na.rm = TRUE)
  
  levels_fine   = seq(z_max - (z_max - z_min)/20, z_max, length.out = 20)
  levels_coarse = seq(z_min, z_max - (z_max - z_min)/10, length.out = 5)
  levels = c(levels_coarse, levels_fine)
  
  # Colors (same length as levels)
  my_colors = colorRampPalette(
    c("lightyellow", "orange", "red", "darkred")
  )(length(levels))
  
  p = ggplot(df, aes(mu, nu, z = LL)) +
    geom_contour(
      aes(colour = after_stat(level)),
      breaks = levels,
      linewidth = 0.4
    ) +
    xlim(floor(mu_limit_lower), mu_limit_upper) + 
    scale_colour_gradientn(colours = my_colors, breaks = levels) +
    geom_point(aes(x = mu_true, y = nu_true), color = "black", size = 2) +
    labs(
      x = expression(mu),
      y = expression(nu),
      colour = "LL"
    ) +
    theme_minimal() +
    theme(legend.position = "none")
  
  return(p)
}

x = rejection_sampler_draws(20000, 2, 0.5)
p1 = generate_cmp_contour_plot(x, 2, 0.5)

x = rejection_sampler_draws(20000, 1, 2)
p2 = generate_cmp_contour_plot(x, 1, 2)

p1+p2
```

![](Notebook_Figures_Tables_files/figure-commonmark/fig4-1.png)

##### Figure 6

``` r
library(ggridges)
library(forcats)
```

    Warning: package 'forcats' was built under R version 4.3.3

``` r
league_acro = 'PL'
season = '2324'
MH_object = readRDS(paste0(mcmc_out_dir, "SAS_FullLeague/", league_acro, "_", season, "_SAS.rds"))
t = MH_object$iterations
N = length(MH_object$team_names)

nus = MH_object$nu_post[(t/10):t, ]
nus[nus == 1] = NA

quantiles = apply(nus, 2, function(x) quantile(x, probs = c(0.1, 0.9), na.rm = TRUE))

classify = function(low, high) {
  if (low > 1) {
    return("Underdispersed")
  } else if (high < 1) {
    return("Overdispersed")
  } else {
    return("Equidispersed")
  }
}

group_colors = c("Underdispersed" = "#D55E00", "Overdispersed" = "#009E73", "Equidispersed" = "gray")

groups = mapply(classify, quantiles[1, ], quantiles[2, ])
names(groups) =  MH_object$team_names

p_z1 = retrieve_nu_sas_summ(MH_object$Z_post, MH_object$nu_post, (t/10)+1, t)[1,]
names(p_z1) = MH_object$team_names

groups[p_z1 <= 0.49] = "Equidispersed"

df_long = as.data.frame(nus)
colnames(df_long) = MH_object$team_names

df_long = df_long |>
  pivot_longer(cols = everything(), names_to = "team", values_to = "value") |>
  filter(!is.na(value)) |>
  mutate(group = groups[team]) |>
  mutate(
    log_value = log(value),
    group = factor(group, levels = c("Underdispersed", "Equidispersed", "Overdispersed"))
  )

summary_df = df_long |>
  group_by(team, group) |>
  summarise(
    ymin = quantile(log(value), 0.10),
    lower = quantile(log(value), 0.25),
    middle = quantile(log(value), 0.5),
    upper = quantile(log(value), 0.75),
    ymax = quantile(log(value), 0.90),
    .groups = "drop"
  )

summary_means = df_long |>
  group_by(team) |>
  summarise(mean_val = median(log(value), na.rm = TRUE))

summary_df = summary_df |>
  left_join(summary_means, by = "team") |>
  mutate(team = fct_reorder(team, mean_val))

team_levels = levels(summary_df$team)
p_z1 = p_z1[team_levels]

df_z1 = data.frame(
  team = factor(names(p_z1), levels = levels(summary_df$team)),
  p_z1 = as.numeric(p_z1)
)

summary_df2 = df_long |>
  group_by(team, group) |>
  mutate(team = factor(team, levels = levels(summary_df$team))) |>
  arrange(team)

p2 = ggplot(summary_df2, aes(
  x = log_value,
  y = factor(team, levels = unique(team)), # reverse to match histogram
  fill = group
)) +
  geom_density_ridges(
    scale = 1,             # reduce ridge height to fit rows
    rel_min_height = 0.01, # minimum ridge height
    alpha = 0.8,
    color = "black",
    size = 0.3
  ) +
  scale_fill_manual(values = group_colors) +
  scale_y_discrete(expand = c(0,0)) +   # remove extra padding
  # scale_x_continuous(limits = c(-3, 2.5)) + 
  coord_cartesian(clip = "off") +
  theme_minimal() +
  theme(
    axis.text.y = element_text(size = 10),
    legend.position = "none",
    plot.title = element_text(hjust = 0.5),
    panel.grid.minor = element_blank()
  ) +
  labs(
    title = "(a)",
    x = expression(log ~ P(nu[i] ~ "|" ~ Z[i] == 1, bold(Y))),
    y = NULL
  ) +
  geom_vline(xintercept = 0, color = "blue", linetype = "dashed", linewidth = 0.7)
```

    Warning in geom_density_ridges(scale = 1, rel_min_height = 0.01, alpha = 0.8, :
    Ignoring unknown parameters: `size`

``` r
offset = 0
p3 = ggplot(df_z1, aes(
  x = p_z1,
  y = factor(team, levels = (unique(team)))  # match order of rows
)) +
  geom_col(fill = "skyblue3", width = 0.7,
           position = position_nudge(y = -offset)) +   # shift bars down
  geom_text(aes(label = sprintf("%.2f", p_z1)),
            hjust = -0.1, size = 3.3,
            position = position_nudge(y = -offset)) +  # shift text too
  geom_segment(
    aes(x = 0.5, xend = 0.5, y = 0.5, yend = 20.87),  # set y range manually
    color = "red3",
    linetype = "dashed",
    linewidth = 0.5
  ) +
  scale_x_continuous(limits = c(0, 1)) +
  scale_y_discrete(expand = c(0,0)) +   # remove padding
  coord_cartesian(clip = "off", ylim = c(1, 20.87)) +
  labs(
    title = "(b)",
    x = expression(P(Z[i] == 1 ~ "|" ~ bold(Y))),
    y = NULL
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_blank(),
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    panel.grid.minor = element_blank(),
    axis.title.x = element_text(margin = margin(b = 2)),
    plot.title = element_text(hjust = 0.5)
  )

p2+p3
```

    Picking joint bandwidth of 0.049

    Warning in geom_segment(aes(x = 0.5, xend = 0.5, y = 0.5, yend = 20.87), : All aesthetics have length 1, but the data has 20 rows.
    ℹ Please consider using `annotate()` or provide this layer with data containing
      a single row.

![](Notebook_Figures_Tables_files/figure-commonmark/fig6-1.png)

##### Figure 7

##### Att and Mean scatter plot comparison

``` r
# library(matrixStats)

league_acro = 'PL'
season = '2324'

load(file = paste0("C:/Users/nakaz/OneDrive/Desktop/Research/SBM-CB/CMP_STZ//Data//MH_Results//Pois_FullLeague//",league_acro,"_",season,"_Run1_Pois",".RData"))

t = MH_P$iterations
df_plot = data.frame(att_means = colMeans(MH_P$att_post[(t/5):t,]),
                     def_means = colMeans(MH_P$def_post[(t/5):t,]),
                     label = MH_P$team_names)

MH_object = readRDS(paste0(mcmc_out_dir, "SAS_FullLeague/", league_acro, "_", season, "_SAS.rds"))

t = MH_object$iterations
df_plot2 = data.frame(att_means = colMeans(MH_object$att_post[(t/5):t,]),
                      def_means = colMeans(MH_object$def_post[(t/5):t,]),
                      label = MH_object$team_names)

team_names = MH_object$team_names

df_colors = data.frame(
  name = team_names,
  color = c(
    "gray",  
    "gray",
    "gray",
    "gray",
    "#009E73",
    "gray",  
    "gray",  
    "gray",  
    "gray",
    "#009E73",
    "gray",
    "gray",
    "gray",
    "gray", 
    "gray",
    "#D55E00",
    "gray", 
    "gray", 
    "gray", 
    "#009E73"  
  )
)
df_plot = cbind(df_plot, colors = df_colors$color)
df_plot2 = cbind(df_plot2, colors = df_colors$color)

library(ggrepel)
```

    Warning: package 'ggrepel' was built under R version 4.3.3

``` r
p1 =  ggplot(df_plot, aes(x = att_means, y = def_means, label = team_names)) +
  geom_jitter(aes(color = colors), shape = 19, size = 2) + geom_text_repel(size = 4, max.overlaps = 12) +
  scale_color_identity() +
  scale_y_reverse() +
  theme_minimal() + 
  geom_vline(xintercept = 0, linetype = 'dashed', alpha = 0.2) +
  geom_hline(yintercept = 0, linetype = 'dashed', alpha = 0.2) +
  coord_cartesian(xlim = c(-2, 0.8)) +
  scale_y_continuous(limits = c(-1, 0.8)) + 
  ylab('') + xlab('') +
  theme(plot.title = element_text(face = "bold")) +
  ggtitle('Poisson Model')
```

    Scale for y is already present.
    Adding another scale for y, which will replace the existing scale.

``` r
p2 = ggplot(df_plot2, aes(x = att_means, y = def_means, label = team_names)) +
  geom_jitter(aes(color = colors), shape = 19, size = 2) + geom_text_repel(size = 4, max.overlaps = 12) +
  scale_color_identity() +
  scale_y_reverse() +
  theme_minimal() + 
  geom_vline(xintercept = 0, linetype = 'dashed', alpha = 0.2) +
  geom_hline(yintercept = 0, linetype = 'dashed', alpha = 0.2) +
  coord_cartesian(xlim = c(-2, 0.8)) +
  scale_y_continuous(limits = c(-1, 0.8)) +
  ylab('') + xlab('') +
  theme(plot.title = element_text(face = "bold")) +
  ggtitle('CMP-SAS Model')
```

    Scale for y is already present.
    Adding another scale for y, which will replace the existing scale.

``` r
plots = list()
plots[[1]] = p1
plots[[2]] = p2


do.call("grid.arrange", c(plots, nrow = 1, ncol = 2, left = 'Mean DEF', bottom = 'Mean ATT'))
```

    Warning: ggrepel: 10 unlabeled data points (too many overlaps). Consider
    increasing max.overlaps

    Warning: ggrepel: 7 unlabeled data points (too many overlaps). Consider
    increasing max.overlaps

![](Notebook_Figures_Tables_files/figure-commonmark/fig7-1.png)

##### Table 3

``` r
library(kableExtra)
```


    Attaching package: 'kableExtra'

    The following object is masked from 'package:dplyr':

        group_rows

``` r
league_acro = 'PL'
season = '2324'

load(file = paste0("C:/Users/nakaz/OneDrive/Desktop/Research/SBM-CB/CMP_STZ//Data//MH_Results//Pois_FullLeague//",league_acro,"_",season,"_Run1_Pois",".RData"))

t = MH_P$iterations
pois_pars = retrieve_pars_summ(MH_P, burn_in = 0.2)[,1]
names(pois_pars) = rownames(retrieve_pars_summ(MH_P))

MH_SAS = readRDS(paste0(mcmc_out_dir, "SAS_FullLeague/", league_acro, "_", season, "_SAS.rds"))

sas_pars = retrieve_pars_summ(MH_SAS, burn_in = 0.2)[,1]
names(sas_pars) = rownames(retrieve_pars_summ(MH_SAS))

df_table = cbind(pois_pars[2:21],
                 pois_pars[22:41],
                 sas_pars[2:21],
                 sas_pars[22:41],
                 (sas_pars[42:61]),
                 retrieve_nu_sas_summ(MH_SAS$Z_post, MH_SAS$nu_post, (t/5), t)[1,])

df_sd = cbind(apply(MH_P$att_post[(t/5):t,], 2, sd),
              apply(MH_P$def_post[(t/5):t,], 2, sd),
              apply(MH_SAS$att_post[(t/5):t,], 2, sd),
              apply(MH_SAS$def_post[(t/5):t,], 2, sd),
              apply((MH_SAS$nu_post[(t/5):t,]), 2, sd))
df_sd = round(df_sd[order(-df_table[,5]),],3)

rownames(df_table) = MH_SAS$team_names

df_table = cbind(round(df_table[,1:5], 3),
      round(df_table[,6], 2))

df_sort = df_table[order(-df_table[,5]),]

df_latex = df_sort |>
  as.data.frame()

for (j in 1:5) {
  df_latex[[j]] = sprintf(
    "%.3f (%.3f)",
    df_latex[[j]],
    df_sd[, j]
  )
}

rn = rownames(df_latex)
rownames(df_latex) = paste0("\\textit{", rn, "}")

# ---- 3. LaTeX table ----
table3 = kable(
  df_latex,
  format = "latex",
  booktabs = TRUE,                        # cleaner table, no row-by-row hlines
  align = c("l", rep("S", ncol(df_latex) - 1)),# row name column is left-aligned
  escape = FALSE                          # allow italic row names
)

print(table3)
```


    \begin{tabular}{llSSSSS}
    \toprule
      & V1 & V2 & V3 & V4 & V5 & V6\\
    \midrule
    \textit{Nott'm Forest} & -0.221 (0.150) & 0.116 (0.125) & 0.054 (0.202) & 0.121 (0.126) & 1.395 (0.570) & 0.53\\
    \textit{Liverpool} & 0.391 (0.112) & -0.419 (0.167) & 0.559 (0.124) & -0.399 (0.164) & 1.261 (0.419) & 0.45\\
    \textit{Bournemouth} & -0.113 (0.145) & 0.124 (0.126) & 0.116 (0.170) & 0.121 (0.129) & 1.241 (0.448) & 0.45\\
    \textit{Brentford} & -0.064 (0.139) & 0.092 (0.129) & 0.110 (0.177) & 0.084 (0.125) & 1.059 (0.261) & 0.31\\
    \textit{Burnley} & -0.433 (0.165) & 0.279 (0.118) & -0.237 (0.398) & 0.273 (0.116) & 1.025 (0.292) & 0.34\\
    \addlinespace
    \textit{Man City} & 0.500 (0.104) & -0.643 (0.183) & 0.624 (0.125) & -0.607 (0.175) & 1.005 (0.161) & 0.24\\
    \textit{Luton} & -0.136 (0.146) & 0.388 (0.112) & 0.025 (0.181) & 0.375 (0.110) & 0.999 (0.205) & 0.29\\
    \textit{Newcastle} & 0.396 (0.111) & 0.065 (0.132) & 0.522 (0.142) & 0.068 (0.127) & 0.986 (0.158) & 0.25\\
    \textit{Aston Villa} & 0.269 (0.118) & 0.039 (0.134) & 0.393 (0.156) & 0.037 (0.128) & 0.985 (0.165) & 0.26\\
    \textit{Crystal Palace} & -0.053 (0.140) & -0.034 (0.135) & 0.089 (0.196) & -0.033 (0.132) & 0.973 (0.184) & 0.28\\
    \addlinespace
    \textit{Man United} & -0.049 (0.141) & -0.036 (0.135) & 0.059 (0.198) & -0.037 (0.135) & 0.968 (0.190) & 0.28\\
    \textit{West Ham} & 0.020 (0.135) & 0.244 (0.121) & 0.132 (0.228) & 0.232 (0.117) & 0.963 (0.194) & 0.28\\
    \textit{Arsenal} & 0.440 (0.109) & -0.854 (0.204) & 0.549 (0.143) & -0.793 (0.202) & 0.962 (0.156) & 0.26\\
    \textit{Everton} & -0.477 (0.165) & -0.202 (0.147) & -0.343 (0.612) & -0.196 (0.140) & 0.958 (0.278) & 0.35\\
    \textit{Tottenham} & 0.242 (0.120) & 0.039 (0.131) & 0.333 (0.199) & 0.031 (0.135) & 0.952 (0.187) & 0.29\\
    \addlinespace
    \textit{Sheffield United} & -0.603 (0.178) & 0.588 (0.102) & -0.606 (0.566) & 0.568 (0.100) & 0.919 (0.301) & 0.37\\
    \textit{Chelsea} & 0.289 (0.118) & 0.075 (0.132) & 0.252 (0.417) & 0.070 (0.128) & 0.819 (0.263) & 0.44\\
    \textit{Wolves} & -0.202 (0.147) & 0.088 (0.130) & -0.475 (0.798) & 0.074 (0.127) & 0.766 (0.301) & 0.50\\
    \textit{Brighton} & -0.097 (0.139) & 0.035 (0.129) & -0.325 (0.472) & 0.028 (0.132) & 0.729 (0.281) & 0.54\\
    \textit{Fulham} & -0.098 (0.141) & 0.014 (0.134) & -1.831 (1.615) & -0.018 (0.130) & 0.323 (0.216) & 0.96\\
    \bottomrule
    \end{tabular}
