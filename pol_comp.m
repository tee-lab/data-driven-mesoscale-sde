close all
clear
clc

% %% Plotting Cohesion Coefficient vs K 
% 
grp_size = 5:5:60;
% 
% figure(1)
% 
% mean_pol = zeros(length(grp_size),1);
% err_pol = zeros(length(grp_size),1);
% 
% for j = 1:length(grp_size)
%     
%     file = sprintf('n%d_pw.csv', grp_size(j));
%     pol_dat = readmatrix(file);
%     pol_dat = pol_dat(:,3);
%     mean_pol(j,1) = mean(pol_dat);
%     std_pol = std(pol_dat);
%     err_pol(j,1) = std_pol/size(pol_dat,1);
%     
% end
% 
% errorbar(grp_size, mean_pol, err_pol, 'Marker', 'o', 'MarkerSize', 10, ...
%         'Color', '#7E2F8E', 'LineStyle', '--', 'Color', '#7E2F8E',...
%         'LineWidth', 2, 'CapSize', 4);
% 
% hold all
% 
% for j = 1:length(grp_size)
%     
%     file = sprintf('n%d_ter.csv', grp_size(j));
%     pol_dat = readmatrix(file);
%     pol_dat = pol_dat(:,3);
%     mean_pol(j,1) = mean(pol_dat);
%     std_pol = std(pol_dat);
%     err_pol(j,1) = std_pol/size(pol_dat,1);
%     
% end
% 
% errorbar(grp_size, mean_pol, err_pol, 'Marker', 'p', 'MarkerSize', 10, ...
%         'Color', '#A2142F', 'LineStyle', ':', 'Color', '#A2142F',...
%         'LineWidth', 2, 'CapSize', 4);
% 
% hold all
% 
% for j = 1:length(grp_size)
%     
%     file = sprintf('n%d_n.csv', grp_size(j));
%     pol_dat = readmatrix(file);
%     pol_dat = pol_dat(:,3);
%     mean_pol(j,1) = mean(pol_dat);
%     std_pol = std(pol_dat);
%     err_pol(j,1) = std_pol/size(pol_dat,1);
%     
% end
% 
% errorbar(grp_size, mean_pol, err_pol, 'Marker', '^', 'MarkerSize', 10, ...
%         'Color', '#0072BD', 'LineStyle', '-.', 'Color', '#0072BD',...
%         'LineWidth', 2, 'CapSize', 4);
% 
% 
% xlabel('N', 'FontName', 'Helvetica', 'FontSize', 14, 'FontWeight', 'bold')
% ylabel('|m|', 'FontName', 'Helvetica', 'FontSize', 14, 'FontWeight', 'bold')
% axis([grp_size(1) grp_size(end) 0.4 1])
% ax = gca;
% ax.XTick = 0:50:max(grp_size);
% ax.YTick = 0.4:0.1:1;
% ax.Box = 'off';
% ax.TickDir = 'out';
% ax.TickLength = [0.01 0.01];
% ax.XMinorTick = 'off';
% ax.YMinorTick = 'off';
% ax.LineWidth = 1.0;
% 
% legend({'Pairwiae', 'Ternary', 'Average'}, 'Location','best')

%% 1/sqrt(N) vs diffusion

% % Pairwise
% drift = -[0.019, 0.019, 0.018, 0.021, 0.016, 0.014, 0.015, 0.016, 0.013, 0.016, 0.016, 0.018];
diffusion_I = [0.089, 0.064, 0.054, 0.046, 0.047, 0.035, 0.033, 0.032, 0.039, 0.034, 0.030, 0.031];
diffusion_m = -[0.079, 0.059, 0.050, 0.044, 0.044, 0.033, 0.031, 0.030, 0.037, 0.032, 0.029, 0.029];

% Ternary
% drift_c = [0.032, 0.045, 0.078, 0.092, 0.129, 0.127, 0.132, 0.182, 0.178, 0.176, 0.230, 0.219];
% drift_m = -[0.058, 0.078, 0.125, 0.143, 0.196, 0.191, 0.195, 0.266, 0.255, 0.253, 0.328, 0.311];
% diffusion_I = [0.026, 0.021, 0.021, 0.029, 0.030, 0.031, 0.033, 0.039, 0.029, 0.028, 0.038, 0.034];
% diffusion_m = -[0.019, 0.016, 0.017, 0.023, 0.024, 0.026, 0.027, 0.032, 0.024, 0.023, 0.031, 0.028];

n_ratio = zeros(length(grp_size), 1);
diffusion_ratio_I = zeros(length(grp_size),1);
diffusion_ratio_m = zeros(length(grp_size), 1);

for i = 1:(length(grp_size))

    diffusion_ratio_I(i) = diffusion_I(i)/diffusion_I(1);
    diffusion_ratio_m(i) = diffusion_m(i)/diffusion_m(1);

end

% diffusion_ratio = 1./diffusion_ratio;

% figure(2)
% 
% scatter(grp_size, drift, 100, 'MarkerEdgeColor',[0 .5 .5],...
%               'MarkerFaceColor',[0 .7 .7],...
%               'LineWidth',1.5)
% 
% xlabel('N', 'FontName', 'Helvetica', 'FontSize', 14, 'FontWeight', 'bold')
% ylabel('f(m)', 'FontName', 'Helvetica', 'FontSize', 14, 'FontWeight', 'bold')

% figure(2)
% 
% plot(grp_size, drift_c, '-o', 'Color', '#7E2F8E',...
%         'LineWidth', 2)
% 
% hold all
% 
% plot(grp_size, drift_m, '-x', 'Color', '#A2142F',...
%         'LineWidth', 2)
% 
% hold off
% 
% xlabel('N', 'FontName', 'Helvetica', 'FontSize', 14, 'FontWeight', 'bold')
% ylabel('Drift', 'FontName', 'Helvetica', 'FontSize', 14, 'FontWeight', 'bold')
% 
% ax = gca;
% ax.Box = 'off';
% 
% legend({'Cons Term', '|m|^2 Term'}, 'Location','best')

figure(3)

plot(grp_size, (1./grp_size)*5, '-o', 'Color', '#7E2F8E',...
        'LineWidth', 2)

hold all

plot(grp_size, diffusion_ratio_I, '--*', 'Color', '#0072BD',...
    'LineWidth', 2)

hold all

plot(grp_size, diffusion_ratio_m, '--x', 'Color', '#A2142F',...
    'LineWidth', 2)
hold off

xlabel('N', 'FontName', 'Helvetica', 'FontSize', 14, 'FontWeight', 'bold')
ylabel('S^2(N)', 'FontName', 'Helvetica', 'FontSize', 14, 'FontWeight', 'bold')
title('Pairwise interaction', 'FontName', 'Helvetica')
axis([grp_size(1) grp_size(end) min((1./grp_size)*5) 1])
ax = gca;
ax.XTick = 0:5:max(grp_size);
ax.YTick = 0.0:0.2:1;
ax.Box = 'off';
ax.TickDir = 'out';
ax.TickLength = [0.01 0.01];
ax.XMinorTick = 'off';
ax.YMinorTick = 'off';
ax.LineWidth = 1.0;

legend({'Expected', 'Cons Term', '|m|^2 Term'}, 'Location','east')