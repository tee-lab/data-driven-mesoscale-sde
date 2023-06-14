close all
clear
clc

%%

tic

grp_size = 250;
S0 = 0.2;
dt = 0.05;
T = 200;
n_iter = round(T/dt);
int_rad = 1;
box_length = 10;

r_spon = 0.1;
sigma_theta = pi;

r_align = 0.8;
k_alg = 1;

no_it = 2;

for gp = 1:length(grp_size)

    n = grp_size(gp);
    
    theta = zeros(n,n_iter,no_it);
    pos = zeros(n,2,n_iter,no_it);
    parfor i = 1:no_it

        [theta_t, pos_t, sum_int] = n_particles(n, r_spon, r_align, sigma_theta, dt, n_iter, ...
            k_alg, S0, box_length, int_rad)

        theta(:,:,i) = theta_t;
        pos(:,:,:,i) = pos_t;

    end

    n_n = struct('pos_t', pos, 'theta_t', theta, 'S0', S0, 'dt', dt, 'n_iter', n_iter, ...
        'box_length', box_length, 'r_spon', r_spon, 'r_align', r_align, 'n', n, 'int_rad', int_rad, ...
        'k_alg', k_alg, 'no_it', no_it, 'sigma_t', sigma_theta);
    file_name = sprintf('n%d_try.mat', n);
    save(file_name, '-struct', 'n_n', '-v7.3')

end

toc

disp('Simulation complete')