close all
clear
clc

%%

tic

grp_size = 10:5:20;

for gp = 1:length(grp_size)

    file_name = sprintf('n%d_try.mat', grp_size(gp));
    load(file_name)

    st_time = 10;

    mx = zeros(no_it, (n_iter - st_time + 1));
    my = zeros(no_it, (n_iter - st_time + 1));
    m = zeros (no_it, (n_iter - st_time + 1));

    for iter = 1:no_it

        mx_iter = mean(cos(theta_t(:,st_time:n_iter,iter)),1);
        my_iter = mean(sin(theta_t(:,st_time:n_iter,iter)),1);
        m_iter = sqrt(mx_iter.^2 + my_iter.^2);

        mx(iter,:) = mx_iter;
        my(iter,:) = my_iter;
        m(iter,:) = m_iter;

    end

    mx = mx.';
    mx = mx(:);
    my = my.';
    my = my(:);
    m = m.';
    m = m(:);

    ord_para = [mx my m];
    
    op_file_name = sprintf('n%d_try.csv', grp_size(gp));
    writematrix(ord_para, op_file_name, 'Delimiter', 'tab')

end

toc