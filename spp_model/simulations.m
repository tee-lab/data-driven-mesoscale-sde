close all
clear
clc

%%
load('n_pw.mat') % Load data

iter = 2; % iteration 

pos_t = pos_t(:,:,:,iter);
theta_t = theta_t(:,:,iter);

for t = 1:10:n_iter

    vel_x = cos(theta_t(:,t));
    vel_y = sin(theta_t(:,t));
    pos_x = pos_t(:,1,t);
    pos_y = pos_t(:,2,t);

    quiver(pos_x, pos_y, vel_x, vel_y, 0.2, 'LineWidth', 2.5, 'ShowArrowHead','off',...
        'Color', '#B0E0E6')

    hold all

    plot(pos_x, pos_y, '.', 'Color', '#00A693', 'MarkerSize', 25)

    hold off

    axis('equal')
    axis([0 box_length 0 box_length])

    drawnow('limitrate')

end
