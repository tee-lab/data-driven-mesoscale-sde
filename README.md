# Data driven mesoscale SDE

## Simulating the flocking model

**Generate data of positions and velocities to calculate group polarisation for given group size and parameters.**

Use Matlab code `/spp_model/sim_data.m` to generate all the required data to calculate group polarisation. Positions, velocities of agents are stored in the .mat file named `n_pw.mat`.

Variables are commented within the code. To reproduce the results, it is advised to run for at least T = 3500 and the number of realisations (no_it) = 15. Interaction rates, number of interacting neighbours (k_alg) and all the parameters can be changed in sim_data.m.

To simulate stochastic pairwise interaction set `k_alg = 1` in `sim_data.m`

