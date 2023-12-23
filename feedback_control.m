%**************************************************************************
% The following code iterates feedback control without  noise (P_control.m)
%**************************************************************************

clear all; close all; clc
tic
graph = 1;

fprintf('>>>>>>>>>  feedback_control.m   <<<<<<<<<\n\n')
n_mdls = 300;
N_inputs = [1 1 0 1 0];
N_outputs = [1 1]; % 0 and 1's determine which outputs are regulated 1 - FR; 2 - EE
EE_lowB_only = 0; % for fig.8

if EE_lowB_only == 1
    fprintf('************   EE - only lower bound   ************\n\n')
end

N_max = 5; % maximum number of tunable inputs (n)
new_signal = 0;
new_ref = 0;
fprintf('PARAMETERS : \n')
fprintf([' - N = ', num2str(sum(N_inputs)),';  M = ', num2str(sum(N_outputs)),'\n'])

%% Conductance indices:
gsubNa_ind = 1;
gsubK_ind = 2;
gL_ind = 3;
gM_ind = 4;
gAHP_ind = 5;

%% REGULATION RATES: (1) gsubNa (2) gsubK (3)gL (4)gM (5)AHP
% NOTE: Must store NON_ZERO values

% e.g. RR's used in Figure 6
rate_t_reg=[-400 100 9999 100 100];
EE_t_reg=[0.15 0.105 9999 0.105 0.105];



%% UPDATE FILENAME
note = ['sample'];

%% LOAD INITIAL CONDITION
% without noise
load('GD_gsubNa=3p5_gsubK=p7_gM=1p75_gAHP=p5_gL=2_stdev=0.1_300mdls.mat')
mu_gi = [3.5 0.7 2 1.75 0.5];

% for Figure4:
% load('GD_gsubNa=3p3_gL=fitted.mat')

if N_inputs(gsubNa_ind) == 1 % varied
    gsubNa_i = gsubNa_gdist(1:n_mdls);
else % fixed
    gsubNa_i = mean(gsubNa_gdist(1:n_mdls));
end

if N_inputs(gsubK_ind) == 1
    gsubK_i = gsubK_gdist(1:n_mdls);
else
    gsubK_i = mean(gsubK_gdist(1:n_mdls));
end
if N_inputs(gL_ind) == 1
    gL_i = gL_gdist(1:n_mdls);
else
    gL_i = 2;
end

if N_inputs(gM_ind) == 1
    gM_i = gM_gdist(1:n_mdls);
else
    gM_i = 1.75;
end

if N_inputs(gAHP_ind) == 1
    gAHP_i = gAHP_gdist(1:n_mdls);
else
    gAHP_i = 0.5;
end


%% TRIAL PARAMETERS:
time = 1500; % (1 stime length in ms)
dt = 0.05;  % time step for forward euler method
loop  = time/dt;   % no. of iterations of euler
t = (1:loop)*dt;
dc = 0; % (uA/cm2) - DC amplitude (NOT USED)

%% Noisy Istim (constant across trials)
if new_signal == 1 % generate new Istim
    isignal_avg = 40;
    tau_isignal = 5; % (ms)
    D_isignal = 10; %1; % sigma(signal) ~ 1 uA/cm2
    i_signal = zeros(1,loop); % Allocate output vector, set initial condition
    for step=1:loop-1
        % Signal (Ornstein-Uhlenbeck process)
        di_signal_dt = -1/tau_isignal*(i_signal(step)-isignal_avg)+ D_isignal/sqrt(dt)*sqrt(2/tau_isignal)*randn;
        i_signal(step+1) = i_signal(step)+dt*di_signal_dt;
    end
    std(i_signal)
    fprintf(' - NEW ')
else % load saved Istim
    load('fast_signal_mu_40_std=10p5_length=1p5sec.mat')
    fprintf(' - SAVED ')
end

fprintf('Fast signal :\n')
fprintf(['         mu_signal = ', num2str(mean(i_signal)),' (uA/cm2)\n'])
fprintf(['         STDEV = ',num2str(std(i_signal)),' (uA/cm2)\n'])

%% SET TARGET OUTPUT VALUES:
load('20201103_target_ranges.mat')

% EE with only a lower bound (fig.8)
if EE_lowB_only == 1
    EE_Trange = 0.22; % fig.8a
%     EE_Trange =  0.27; % fig.8b
%     EE_min = 0.30; % fig.8c
end


%% PRINT TARGET RANGES:
fprintf('\n\n\n - Target outputs :\n')
fprintf(['     Firing rate = ',num2str(RATE_Trange(1)),' ~ ',num2str(RATE_Trange(2)),'\n'])
if numel(EE_Trange) == 2 % with upper and lower bound
    fprintf(['     EE = ',num2str(EE_Trange(1)),' ~ ',num2str(EE_Trange(2)),' (fixed) \n'])
else % only a lower bound
    fprintf(['     EE >= ',num2str(EE_Trange),'\n'])
end
fprintf('*************************************************\n\n\n')
%% INITIALIZE EMPTY VECTORS:
pre_g = zeros(n_mdls, N_max); % pre-gi's
post_g = zeros(n_mdls,N_max); % post-gi's
success = zeros(n_mdls,1);
final_M = zeros(n_mdls,2); % store final output values
% gi_traj = zeros(300,N_max,n_mdls);

% Store pre-gi's (row = number of mdls; cols = gi)
pre_g(:,gsubNa_ind) = gsubNa_i;
pre_g(:,gsubK_ind) = gsubK_i;
pre_g(:,gL_ind) = gL_i;
pre_g(:,gM_ind) = gM_i;
pre_g(:,gAHP_ind) = gAHP_i;


%% HOMEOSTATIC RULE
for i = 1:n_mdls
    fprintf('\n')
    [post_g(i,:), success(i), gi_traj(:,:,i),M_traj(:,:,i),final_M(i,:),n_it(i)] = P_control(N_inputs,N_outputs,RATE_Trange,EE_Trange,rate_t_reg, EE_t_reg, i_signal, pre_g(i,gsubNa_ind),pre_g(i,gsubK_ind),pre_g(i,gL_ind),pre_g(i,gM_ind),pre_g(i,gAHP_ind));
    % print progress
    if mod(i,10)== 0 % for every 10 models,
        P_progress = i/n_mdls*100;
        fprintf(['\n\n\n>>>>>>>>   ',num2str(P_progress),'   <<<<<<<<<\n\n\n'])
    end
end

n_fails = sum(success==0);
fprintf(['Number of failures = ', num2str(n_fails),'\n'])

elapsed_time = toc;



%% SAVE
% FileName=[datestr(now, 'yyyymmdd'),'_N=',num2str(sum(N_inputs)),'_M=',num2str(sum(N_outputs)),'_',num2str(n_mdls),'mdls_',note,'.mat'];
% save(FileName)

%% VISUALIZE << NEEDS REVISION
if graph == 1
    fontsize = 15;
    figure('name','pre&post_homeo_2D')
    
    % Select gi indicies for plotting:
    ax_ind = [gsubNa_ind gsubK_ind gM_ind];
    
    fail_inds = find(success~=1);
    if (numel(fail_inds)>0) % if there are failures
        scatter3(pre_g(fail_inds,ax_ind(1)),pre_g(fail_inds,ax_ind(2)),pre_g(fail_inds,ax_ind(3)),'MarkerEdgeColor','k','MarkerFaceColor','w')
        hold on
        scatter3(post_g(fail_inds,ax_ind(1)),post_g(fail_inds,ax_ind(2)),post_g(fail_inds,ax_ind(3)),'MarkerEdgeColor','k','MarkerFaceColor','y')
        for i=1:numel(fail_inds)
            x = gi_traj(:,ax_ind(1),fail_inds(i));
            y = gi_traj(:,ax_ind(2),fail_inds(i)); 
            z = gi_traj(:,ax_ind(3),fail_inds(i)); 
            
            plot3(x,y,z,'r-')
            hold on
        end
    end % end of if (numel(fail_inds)>0)
    
    % exclude failures
    pre_g(fail_inds,:) = []; post_g(fail_inds,:) = [];
    gi_traj(:,:,fail_inds) = [];
    
    hold on
    for i=1:n_mdls -numel(fail_inds)       
        x = gi_traj(:,ax_ind(1),i);
        y = gi_traj(:,ax_ind(2),i); 
        z = gi_traj(:,ax_ind(3),i); 
        
        plot3(x,y,z,'r-')
        hold on
    end
    
    hold on
    
    scatter3(pre_g(:,ax_ind(1)),pre_g(:,ax_ind(2)),pre_g(:,ax_ind(3)),'MarkerEdgeColor','k','MarkerFaceColor','w')
    hold on
    scatter3(post_g(:,ax_ind(1)),post_g(:,ax_ind(2)),post_g(:,ax_ind(3)),'MarkerEdgeColor','k','MarkerFaceColor','r')
    pbaspect([1 1 1])
    xlabel('g_{Na}');ylabel('g_{K}');zlabel('g_{M}')
    
    
    grid on
    axis([0 4 0 4 0 4]); view([-149 12])
    set(gca,'FontSize',fontsize,'TickDir','out')
    
    %% To overlay surfaces
    %%%%% RATE
    load('') % gridsearch data for firing rate and energy efficiency (gridsearch.m)
    
    t_RATE_inds = find(RATE>RATE_Trange(1) & RATE < RATE_Trange(2));
    [x_RATE_3D,y_RATE_3D,z_RATE_3D] = ind2sub(size(RATE),t_RATE_inds); % target indices
    RATE_sf3D = fit([gsubNa(x_RATE_3D)',gsubK(y_RATE_3D)'],gM(z_RATE_3D)','poly21'); %original
    % Below is x fitted based on y & z (to find the intersection of FR and EE
    % isosurfaces)
    RATE_sf3D_graph = fit([gsubK(y_RATE_3D)',gM(z_RATE_3D)'],gsubNa(x_RATE_3D)','poly23');
    
    
    %%%%% ENERGY EFFICIENCY
    t_EE_inds = find(EE>EE_Trange(1) & EE<EE_Trange(2)); % find indicies for models within range
    [x_EE_3D,y_EE_3D,z_EE_3D] = ind2sub(size(EE),t_EE_inds);
    EE_sf3D = fit([gsubNa(x_EE_3D)',gsubK(y_EE_3D)'],gM(z_EE_3D)','poly21'); % fitting based on x and y
    % Below is x fitted based on y & z (to find the intersection of FR and EE
    % isosurfaces)
    EE_sf3D_graph = fit([gsubK(y_EE_3D)',gM(z_EE_3D)'],gsubNa(x_EE_3D)','poly23'); % fitting based on y and z
    % scatter3(gsubNa(x_EE_3D)',gsubK(y_EE_3D)',gM(z_EE_3D)')
    
    
    %**** Plot surfaces
    %%%%% RATE
    [x1,y1] = meshgrid(linspace(0,4,30),linspace(0,4,30)); %original
    % RATE_sf = surf(x1,y1,RATE_sf3D(x1,y1),'FaceColor','r','FaceAlpha',0.5);
    RATE_sf = surf(RATE_sf3D_graph(x1,y1),x1,y1,'FaceColor','r','FaceAlpha',0.5);
    RATE_sf.EdgeColor = [170 0 0]./255; % red
    
    hold on
    
    %%%%% ENERGY EFFICIENCY
    [y1,z1] = meshgrid(linspace(0,4,30),linspace(0,4,30)); %original
    EE_sf = surf(EE_sf3D_graph(y1,z1),y1,z1,'FaceColor','g','FaceAlpha',0.3);
    EE_sf.EdgeColor = [0 0.7 0]; EE_sf.EdgeAlpha = 0.6;
    
    %%%%% INTERSECTION
    [yL_RATE_EE, zL_RATE_EE,xL_RATE_EE,]= find_inter(RATE_sf3D_graph, EE_sf3D_graph);
    hold on
    line(xL_RATE_EE,yL_RATE_EE,zL_RATE_EE,'Color','y','LineWidth',1);
    
    
end