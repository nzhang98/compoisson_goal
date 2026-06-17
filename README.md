# compoisson_goal
Repository for the paper on the Conway-Maxwell Poisson Goal model with Spike-and-Slab priors on the dispersion parameters for association football scores.

- 'utils.R': contains all utility functions needed to import data, working with MH objects, creating predictions, evaluating the models, etc.;
- 'MH_Poisson.R': Contains the MH algorithm to infer the parameters for the Poisson Goal model;
- 'MH_CMP_SAS.R': Contains the MWGS to infer the parameters for the COMPoisson Goal model with Spike-and-Slab;
- 'MH_CMP_Full.R': Contains the MWGS for the fully CMP model (i.e. indicators Z = 1 for all i);
- 'Notebook_Figures_Tables.qmd/md': Notebook with example usages and reproducing all the figures and tables in the paper;
- 'Generate_Simulations.R': Script to generate the (seeded) simulation data used in Chapter 3;
- 'Run_FullLeagues.R', 'Run_MidLeagues.R', 'Run_Simulations.R': Script to run the MCMC for the full leagues, partial leagues, and simulation leagues respectively;
- 'Generate_Predictions.R': Script to generate predictions from the posterior predictive distribution for each match.
