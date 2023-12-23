%**************************************************************************
% The following code reproduces Figure 4
%**************************************************************************
clear all; close all; clc
RATE_range = [37 43];

load('') % Load firing rate data using gridsearch.m when gK = 2 (gridsearch.m)
t_RATE_inds1 = find(RATE>RATE_range(1) & RATE < RATE_range(2));
[y_RATE_2D_1,x_RATE_2D_1] = ind2sub(size(RATE),t_RATE_inds1);
RATE_sf2D_1 = fit(gsubNa(x_RATE_2D_1)',gL(y_RATE_2D_1)','poly2');

load('') % Load firing rate data using gridsearch.m when gK = 0 (gridsearch.m)
t_RATE_inds2 = find(RATE>RATE_range(1) & RATE < RATE_range(2));
[y_RATE_2D_2,x_RATE_2D_2] = ind2sub(size(RATE),t_RATE_inds2);
RATE_sf2D_2 = fit(gsubNa(x_RATE_2D_2)',gL(y_RATE_2D_2)','poly2');



figure('name','scatter')
RATE_line = plot(RATE_sf2D_1,'--r');
RATE_line.LineWidth = 2;
xlim([0 4]); ylim([1 4])
hold on
RATE_line = plot(RATE_sf2D_2,'-r');
RATE_line.LineWidth = 2;

legend off
xlabel('g_{Na}'); ylabel('g_{leak}')
pbaspect([1 1 1])

%% Generate initial values that lie on the iso-firing rate line
gsubNa_gdist = normrnd(3.3,0.1,[1 300]);
gL_gdist = RATE_sf2D_1(gsubNa_gdist);

%% HOMEOSTATICALLY FOUND SOLUTIONS
load('') % solutions found by gNa compensation(feedback_control.m)
gi_x = gi_traj;
hold on
scatter(pre_g(:,1),pre_g(:,3),'Marker','o','MarkerFaceColor',[0.7 0.7 0.7],'MarkerEdgeColor','k')
hold on
scatter(post_g(:,1),post_g(:,3),'Marker','o','MarkerFaceColor','m','MarkerEdgeColor','k')
hold on

for i=1:n_mdls
    last_it = find (M_traj(:,4,i)~=0,1,'last');
    final_EE_1(i) =  M_traj(last_it,4,i);
end


load('') % solutions found by gleak compensation(feedback_control.m)
gi_y = gi_traj;
scatter(post_g(:,1),post_g(:,3),'MarkerFaceColor','c','MarkerEdgeColor','k')

pbaspect([1 1 1]); axis([0 4 1 3])
set(gca,'TickDir','out','FontSize',15); box off



%% Plot energy efficiency distribution
figure('name','EE distribution')

% find last_it to store final EE values
for i=1:n_mdls
    last_it = find (M_traj(:,4,i)~=0,1,'last');
    final_EE_2(i) =  M_traj(last_it,4,i);
end
h1 = histfit(final_EE_1,20); % fit distribution in 20 bins
set(h1(1),'facecolor','k'); set(h1(2),'color','b')

hold on
h2 = histfit(final_EE_2,20); % fit distribution in 20 bins
set(h2(1),'facecolor','k'); set(h2(2),'color','g')

EE = ; % Stores energy efficiency values for neuron models pre- and post- 
% gNa/gleak compensation (calc_EE)
h3 = histfit(EE,20); % fit distribution in 20 bins
set(h3(1),'facecolor','k'); set(h3(2),'color','m');
ax = gca;
ax.XTick = 0.15:0.05:0.3; xlim([0.15 0.3])
set(ax,'TickDir','out','FontSize',15); box off
xlabel('Energy Efficiency')



%% PLOT TRAJECTORIES
%***NOTE: % set to gi_x for delta gNa, g_y for delta gL****

pre_knockout = 10; % number of iterations pre-knockout for plotting
gi_traj = gi_y; 
figure('name','gi_trajectories')
max_it = size(gi_traj,1);


gi_traj(gi_traj==0) = NaN;
mean_traj = mean(gi_traj,3);

for i_mdl = 1:n_mdls
    x = 0:max_it + pre_knockout-1;
    y = gi_traj(1:max_it,gsubNa_ind,i_mdl);
    y = [ones(pre_knockout,1)*gsubNa_gdist(i_mdl); y];
    
    plot(x,y,'r')
    hold on
    y = gi_traj(1:max_it,gL_ind,i_mdl);
    y = [ones(pre_knockout,1)*gL_gdist(i_mdl); y];
    plot(x,y,'c')
    hold on
    
end

xlim([0 30+pre_knockout]); ylim([0 4])
ylabel('g_{i}');xlabel('# of iterations');
set(gcf,'position',[506   605   420   191])
set(gca,'TickDir','out','FontSize',15); box off