close all
clear
clc

%%

tic

load('n_pw.mat')

st_time = 10;

mx = zeros(no_it, (n_iter - st_time + 1)); % store x component of group polarisation (mx)
my = zeros(no_it, (n_iter - st_time + 1)); % store y component of group polarisation (my)
m = zeros (no_it, (n_iter - st_time + 1)); % group polarisation m = sqrt(mx^2 + my^2)

for iter = 1:no_it

    mx_iter = mean(cos(theta_t(:,st_time:n_iter,iter)),1); % calculate mx
    my_iter = mean(sin(theta_t(:,st_time:n_iter,iter)),1); % calculate my 
    m_iter = sqrt(mx_iter.^2 + my_iter.^2); % calculate m

    mx(iter,:) = mx_iter;
    my(iter,:) = my_iter;
    m(iter,:) = m_iter;

end

% storing mx, my and m
mx = mx.';
mx = mx(:);
my = my.';
my = my(:);
m = m.';
m = m(:);

ord_para = [mx my m];

% storing order paramter data as [mx, my, m] in .csv format
writematrix(ord_para, 'n_pw.csv', 'Delimiter', 'tab')

toc