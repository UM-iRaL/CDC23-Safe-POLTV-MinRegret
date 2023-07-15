%% plot varing horizon result
clc; close all; clear;
disturbance.profiles = ["Gaussian" "Uniform" "Gamma" "Exponential " "Bernoulli" "Weibull" "Poisson" "Worst-case"];

% Random open-loop stable system
load("./Data/data_system_1_T_1_30_rho_085.mat");
[horizon_max, ~, noise_num] = size(cost_mean);
T = 2:horizon_max+1;
cost_mean_low = cost_mean - cost_std;
cost_mean_up = cost_mean + cost_std;

row_num = 2;
col_num = 4;
figure(1)
for i = 1:noise_num
    subplot(row_num,col_num,i)
    plot(T,cost_mean(1:horizon_max,1,i),'rs-');hold on;
    plot(T,cost_mean(1:horizon_max,2,i),'go-');hold on;
    plot(T,cost_mean(1:horizon_max,3,i),'*-','Color',[0 0.7 1]);hold on;
    plot(T,cost_mean(1:horizon_max,4,i),'k^-');hold on;

    plot(T,cost_mean_low(1:horizon_max,1,i),'r-');
    plot(T,cost_mean_up(1:horizon_max,1,i),'r-');
    patch([T(:); flipud(T(:))],[cost_mean_low(1:horizon_max,1,i); flipud(cost_mean_up(1:horizon_max,1,i))], 'r', 'FaceAlpha',0.05, 'EdgeColor','none');

    plot(T,cost_mean_low(1:horizon_max,2,i),'g-');
    plot(T,cost_mean_up(1:horizon_max,2,i),'g-');
    patch([T(:); flipud(T(:))],[cost_mean_low(1:horizon_max,2,i); flipud(cost_mean_up(1:horizon_max,2,i))], 'g', 'FaceAlpha',0.05, 'EdgeColor','none');
    
    plot(T,cost_mean_low(1:horizon_max,3,i),'-','Color',[0 0.7 1]);
    plot(T,cost_mean_up(1:horizon_max,3,i),'-','Color',[0 0.7 1]);
    patch([T(:); flipud(T(:))],[cost_mean_low(1:horizon_max,3,i); flipud(cost_mean_up(1:horizon_max,3,i))],[0 0.7 1], 'FaceAlpha',0.05, 'EdgeColor','none');

    plot(T,cost_mean_low(1:horizon_max,4,i),'k-');
    plot(T,cost_mean_up(1:horizon_max,4,i),'k-');
    patch([T(:); flipud(T(:))],[cost_mean_low(1:horizon_max,4,i); flipud(cost_mean_up(1:horizon_max,4,i))], 'k', 'FaceAlpha',0.05, 'EdgeColor','none');

    set(gca,'FontSize', 8)
    set(gca,'TickLabelInterpreter','latex')
    set(gca,'LooseInset', max(get(gca,'TightInset'), 0.3))
    xlabel('$T$','interpreter','latex')
    ylabel('$Cost$','interpreter','latex') 
    title(disturbance.profiles(i))
end

legend('$\mathcal{H}_2$ Control','$\mathcal{H}_{\infty}$ Control','Ours', 'Clairvoyant', 'interpreter','latex', 'FontSize', 5);



% Linearized UAV translational model
load("./Data/data_system_4_T_1_25.mat");
horizon_max = 24; T = 2:horizon_max+1;
cost_mean(:,:,8) = cost_mean(:,:,8)/10;
cost_mean_low = cost_mean - cost_std;
cost_mean_up = cost_mean + cost_std;

figure(2)
for i = 1:noise_num
    subplot(row_num,col_num,i)
    plot(T,cost_mean(1:horizon_max,1,i),'rs-');hold on;
    plot(T,cost_mean(1:horizon_max,2,i),'go-');hold on;
    plot(T,cost_mean(1:horizon_max,3,i),'*-','Color',[0 0.7 1]);hold on;
    plot(T,cost_mean(1:horizon_max,4,i),'k^-');hold on;

    plot(T,cost_mean_low(1:horizon_max,1,i),'r-');
    plot(T,cost_mean_up(1:horizon_max,1,i),'r-');
    patch([T(:); flipud(T(:))],[cost_mean_low(1:horizon_max,1,i); flipud(cost_mean_up(1:horizon_max,1,i))], 'r', 'FaceAlpha',0.05, 'EdgeColor','none');

    plot(T,cost_mean_low(1:horizon_max,2,i),'g-');
    plot(T,cost_mean_up(1:horizon_max,2,i),'g-');
    patch([T(:); flipud(T(:))],[cost_mean_low(1:horizon_max,2,i); flipud(cost_mean_up(1:horizon_max,2,i))], 'g', 'FaceAlpha',0.05, 'EdgeColor','none');
    
    plot(T,cost_mean_low(1:horizon_max,3,i),'-','Color',[0 0.7 1]);
    plot(T,cost_mean_up(1:horizon_max,3,i),'-','Color',[0 0.7 1]);
    patch([T(:); flipud(T(:))],[cost_mean_low(1:horizon_max,3,i); flipud(cost_mean_up(1:horizon_max,3,i))],[0 0.7 1], 'FaceAlpha',0.05, 'EdgeColor','none');

    plot(T,cost_mean_low(1:horizon_max,4,i),'k-');
    plot(T,cost_mean_up(1:horizon_max,4,i),'k-');
    patch([T(:); flipud(T(:))],[cost_mean_low(1:horizon_max,4,i); flipud(cost_mean_up(1:horizon_max,4,i))], 'k', 'FaceAlpha',0.05, 'EdgeColor','none');

    set(gca,'FontSize', 8)
    set(gca,'TickLabelInterpreter','latex')
    set(gca,'LooseInset', max(get(gca,'TightInset'), 0.3))
    xlabel('$T$','interpreter','latex')
    ylabel('$Cost$','interpreter','latex') 
    xlim([0 25])
    title(disturbance.profiles(i))
end

legend('$\mathcal{H}_2$ Control','$\mathcal{H}_{\infty}$ Control','Ours', 'Clairvoyant', 'interpreter','latex', 'FontSize', 5);



% Random open-loop unstable system
load("./Data/data_system_2_T_1_30_rho_105.mat");
[horizon_max, ~, noise_num] = size(cost_mean);
T = 2:horizon_max+1;
cost_mean_low = cost_mean - cost_std;
cost_mean_up = cost_mean + cost_std;

figure(3)
for i = 1:noise_num
    subplot(row_num,col_num,i)
    plot(T,cost_mean(1:horizon_max,1,i),'rs-');hold on;
    plot(T,cost_mean(1:horizon_max,2,i),'go-');hold on;
    plot(T,cost_mean(1:horizon_max,3,i),'*-','Color',[0 0.7 1]);hold on;
    plot(T,cost_mean(1:horizon_max,4,i),'k^-');hold on;

    plot(T,cost_mean_low(1:horizon_max,1,i),'r-');
    plot(T,cost_mean_up(1:horizon_max,1,i),'r-');
    patch([T(:); flipud(T(:))],[cost_mean_low(1:horizon_max,1,i); flipud(cost_mean_up(1:horizon_max,1,i))], 'r', 'FaceAlpha',0.05, 'EdgeColor','none');

    plot(T,cost_mean_low(1:horizon_max,2,i),'g-');
    plot(T,cost_mean_up(1:horizon_max,2,i),'g-');
    patch([T(:); flipud(T(:))],[cost_mean_low(1:horizon_max,2,i); flipud(cost_mean_up(1:horizon_max,2,i))], 'g', 'FaceAlpha',0.05, 'EdgeColor','none');
    
    plot(T,cost_mean_low(1:horizon_max,3,i),'-','Color',[0 0.7 1]);
    plot(T,cost_mean_up(1:horizon_max,3,i),'-','Color',[0 0.7 1]);
    patch([T(:); flipud(T(:))],[cost_mean_low(1:horizon_max,3,i); flipud(cost_mean_up(1:horizon_max,3,i))],[0 0.7 1], 'FaceAlpha',0.05, 'EdgeColor','none');
    
    plot(T,cost_mean_low(1:horizon_max,4,i),'k-');
    plot(T,cost_mean_up(1:horizon_max,4,i),'k-');
    patch([T(:); flipud(T(:))],[cost_mean_low(1:horizon_max,4,i); flipud(cost_mean_up(1:horizon_max,4,i))], 'k', 'FaceAlpha',0.05, 'EdgeColor','none');

    set(gca,'FontSize',8)
    set(gca,'TickLabelInterpreter','latex')
    set(gca,'LooseInset', max(get(gca,'TightInset'), 0.3))
    xlabel('$T$','interpreter','latex')
    ylabel('$Cost$','interpreter','latex') 
    title(disturbance.profiles(i))
end

legend('$\mathcal{H}_2$ Control','$\mathcal{H}_{\infty}$ Control','Ours', 'Clairvoyant', 'interpreter','latex', 'FontSize', 5);
