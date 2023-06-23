close all
clear
clc

%%

tic

n = 30; % No.of agents
S0 = 0.2; % speed
dt = 0.05; % Integration time.
T = 200; % Simulation time
n_iter = round(T/dt);
int_rad = 1; % interaction radius
box_length = 3.5; % box length

k_alg = 1; % no of agents to interact with. change k_alg = 2 for ternary and k_alg = n for Vicsek like interaction

r_spon = 0.1; % Spontaneous interaction rate
sigma_theta = pi;

r_align = 0.8; % alignment interaction rate

no_it = 2; % No.of realisations

parfor i = 1:no_it

    [theta_t, pos_t, sum_int] = n_particles(n, r_spon, r_align, sigma_theta, dt, n_iter, ...
        k_alg, S0, box_length, int_rad)

    theta(:,:,i) = theta_t; % orientation 
    pos(:,:,:,i) = pos_t; % position

end

% Store all data in .mat file as structure
n_n = struct('pos_t', pos, 'theta_t', theta, 'S0', S0, 'dt', dt, 'n_iter', n_iter, ...
    'box_length', box_length, 'r_spon', r_spon, 'r_align', r_align, 'n', n, 'int_rad', int_rad, ...
    'k_alg', k_alg, 'no_it', no_it, 'sigma_t', sigma_theta);
save('n_pw.mat', '-struct', 'n_n', '-v7.3')

disp('Simulation complete')

toc