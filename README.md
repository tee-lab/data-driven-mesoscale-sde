# Data driven mesoscale SDE

This repository contains the code for the manuscript:

Data-driven discovery of stochastic dynamical equations of collective motion. 
Nabeel, A., Jadhav, V., Sire, C., Raj, M. D., Theraulaz, G., Escobedo, R., Iyer, S. K., & Guttal, V. (2023).
Physical Biology, 20(5), 056003

DOI: [10.1088/1478-3975/ace22d](https://dx.doi.org/10.1088/1478-3975/ace22d)

The code is divided into two parts:
  - The directory `/spp_model` contains MATLAB code to simulate the agent-based models.
  - The directory `/sde_discovery` contains code to discover mesoscale SDEs from the simulated trajectories.

## Simulating the flocking model

**Generate data of positions and velocities to calculate group polarisation for given group size and parameters.**

Use Matlab code `/spp_model/sim_data.m` to generate all the required data to calculate group polarisation. Positions, velocities of agents are stored in .mat file named `n_pw.mat`.

Variables are commented within the code. To reproduce the results, it is advised to run for at least T = 3500 and the number of realisations (no_it) = 15. Interaction rates, number of interacting neighbours (k_alg) and all the parameters can be changed in `sim_data.m`.

To simulate stochastic pairwise interaction set `k_alg = 1` in `sim_data.m`. Similarly set `k_alg = 2` and `k_alg = n` for ternary and Vicsek like averaging model respectively.

After running the simulation, use the data in `n_pw.mat` to calculate the order parameter, i.e., group polarisation. To do so, run the Matlab code `/figures/grp_pol.m`. This file stores all the required data in the `n_pw.csv` file in the format `[mx, my, m]`.

### Simulation video

To see the collective motion of agents, run the code `/spp_model/simulations.m`. Make sure to load `n_pw.mat`. Variables are defined within the code and can be changed accordingly.

## Data driven SDE discovery

The notebook `/sde_discovery/sde-discovery.ipynb` contains code to analyze the order-parameter time series and discover mesoscale stochastic differential equations. The code uses the CSV file generated from the MATLAB code above. An example CSV file, corresponding to the _ternary local interaction model_ is provided.
