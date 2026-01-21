# compoisson_goal
Repository for the paper on the Conway-Maxwell Poisson Goal model with Spike-and-Slab priors on the dispersion parameters for association football scores.

- 'utils.R': contains all miscellaneous functions needed to import data, working with MH objects, creating predictions;
- 'MH_Poisson.R': Contains the MH algorithm to infer the parameters for the Poisson Goal model;
- 'MH_CMP.R': Contains the MWGS to infer the parameters for the COMPoisson Goal model with Spike-and-Slab (NB: A fully-CMP model can be retrieved by forcing indicators Z = 1)
