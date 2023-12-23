%**************************************************************************
% The following code reproduces Figure 6, but can be modified to reproduce 
% 3D plots of iso-output surfaces, conductance trajectories and 
% correlation matrices in Figures 7-10 
%**************************************************************************
clear all; close all; clc

%% TARGET RANGES:
load('output_target_values.mat')

%% LOAD DATA:
load('') % gridsearch data for FR and EE (gridsearch.m)

%% Visualize SOSSs
figure
%%%%% RATE
t_RATE_inds = find(RATE>RATE_Trange(1) & RATE < RATE_Trange(2));
[x_RATE_3D,y_RATE_3D,z_RATE_3D] = ind2sub(size(RATE),t_RATE_inds); % target indices
RATE_sf3D = fit([gsubNa(x_RATE_3D)',gsubK(y_RATE_3D)'],gM(z_RATE_3D)','poly21'); %original
RATE_sf3D_graph = fit([gsubK(y_RATE_3D)',gM(z_RATE_3D)'],gsubNa(x_RATE_3D)','poly23');

%****visualize surface fit 
[x1,y1] = meshgrid(linspace(0,4,30),linspace(0,4,30)); %original
% RATE_sf = surf(x1,y1,RATE_sf3D(x1,y1),'FaceColor','r','FaceAlpha',0.5);
RATE_sf = surf(RATE_sf3D_graph(x1,y1),x1,y1,'FaceColor','r','FaceAlpha',0.5);
RATE_sf.EdgeColor = [170 0 0]./255; % red
hold on


%%%%% ENERGY EFFICIENCY 
t_EE_inds = find(EE>EE_Trange(1) & EE<EE_Trange(2)); % find indicies for models within range
[x_EE_3D,y_EE_3D,z_EE_3D] = ind2sub(size(EE),t_EE_inds);
EE_sf3D = fit([gsubNa(x_EE_3D)',gsubK(y_EE_3D)'],gM(z_EE_3D)','poly21'); % fitting based on x and y
EE_sf3D_graph = fit([gsubK(y_EE_3D)',gM(z_EE_3D)'],gsubNa(x_EE_3D)','poly23'); % fitting based on y and z 
% scatter3(gsubNa(x_EE_3D)',gsubK(y_EE_3D)',gM(z_EE_3D)')

%****visualize surface fit
% [x1,y1] = meshgrid(linspace(0,4,30),linspace(0,4,30)); %original
% EE_sf = surf(x1,y1,EE_sf3D(x1,y1),'FaceColor','g','FaceAlpha',0.5);
[y1,z1] = meshgrid(linspace(0,4,30),linspace(0,4,30)); %original
EE_sf = surf(EE_sf3D_graph(y1,z1),y1,z1,'FaceColor','g','FaceAlpha',0.3);
EE_sf.EdgeColor = [0 0.7 0]; EE_sf.EdgeAlpha = 0.6;

% intersection
[yL_RATE_EE, zL_RATE_EE,xL_RATE_EE,]= find_inter(RATE_sf3D_graph, EE_sf3D_graph); 
hold on
line(xL_RATE_EE,yL_RATE_EE,zL_RATE_EE,'Color','y','LineWidth',1); 

%% HOMEOSTATICALLY FOUND SOLUTIONS
%**** plot solutions for RR1 (e.g. cyan RR's) 
load('') % data generated using feedback_control.m for cyan RR's

% plot starting poisition
scatter3(pre_g(:,1),pre_g(:,2),pre_g(:,4),'MarkerFaceColor','w','MarkerEdgeColor','k')
hold on
scatter3(pre_g(:,1),pre_g(:,2),pre_g(:,4)*0,'marker','.','MarkerEdgeColor',[0.5 0.5 0.5] ) % plot shadows

% plot final position
hold on
scatter3(post_g(:,1),post_g(:,2),post_g(:,4),'MarkerFaceColor','c','MarkerEdgeColor','k')

% plot trajectories
for n=1:n_mdls % find last ind
    Na_ind = max(find(gi_traj(:,1,n)~=0));
    K_ind = max(find(gi_traj(:,2,n)~=0));
    M_ind = max(find(gi_traj(:,4,n)~=0));
    AHP_ind = max(find(gi_traj(:,5,n)~=0));
    last_it = max([Na_ind K_ind M_ind AHP_ind]);
    gi_traj(last_it+1:end,:,n) = NaN;
end

for i=1:n_mdls
    x = gi_traj(:,1,i);
    y = gi_traj(:,2,i);
    z = gi_traj(:,4,i); 
    
    plot3(x,y,z, 'Color','c')
    hold on
end
hold on
% plot shadows
scatter3(post_g(:,1),post_g(:,2),post_g(:,4)*0,'marker','o','MarkerEdgeColor',[0.5 0.5 0.5],'MarkerFaceColor',[0.5 0.5 0.5] )



%**** plot solutions for RR2 (e.g. magenta RR's) 
load('') % data generated using feedback_control.m for magenta RR's

hold on
scatter3(post_g(:,1),post_g(:,2),post_g(:,4),'MarkerFaceColor','m','MarkerEdgeColor','k')

% plot trajectories
for n=1:n_mdls % find last ind
    Na_ind = max(find(gi_traj(:,1,n)~=0));
    K_ind = max(find(gi_traj(:,2,n)~=0));
    M_ind = max(find(gi_traj(:,4,n)~=0));
    AHP_ind = max(find(gi_traj(:,5,n)~=0));
    last_it = max([Na_ind K_ind M_ind AHP_ind]);
    gi_traj(last_it+1:end,:,n) = NaN;
end

for i=1:n_mdls
    x = gi_traj(:,1,i);
    y = gi_traj(:,2,i);
    z = gi_traj(:,4,i);
    
    plot3(x,y,z,'m')
    hold on
end

% plot shadows
hold on
scatter3(post_g(:,1),post_g(:,2),post_g(:,4)*0,'marker','o','MarkerEdgeColor',[0.5 0.5 0.5],'MarkerFaceColor',[0.5 0.5 0.5] )


%****plot settings
set(gca,'TickDir','out','FontSize',15)
axis([0 4 0 4 0 4]); xlabel('g_{Na}');ylabel('g_{K}');zlabel('g_{M}')
pbaspect([1 1 1])
view([-149 12])




%% Plot correlation matrix
clear all

load('') % data generated using feedback_control.m for RR's

% Conductance indicies
gsubNa_ind=1;
gsubK_ind=2;
gL_ind=3;
gM_ind=4;
gAHP_ind=5;

center_mean = 0; % 1 - center points around their means 

N = 4; % number of channels  

% data = pre_g; fprintf('************* INITIAL CONDITION: \n\n') 
data = squeeze(gi_traj(end,:,:))'; % the last it 

data(:,gAHP_ind) = []; % exclude gAHP

if center_mean == 1 % z-scores: standardized = subtracted by its mean, then divided by its stdev
    data(:,gsubNa_ind) =  (data(:,gsubNa_ind)-mean(data(:,gsubNa_ind)))/std(data(:,gsubNa_ind));
    if N == 3
        data(:,gsubK_ind) =  zeros(300,1);
    else
        data(:,gsubK_ind) =  (data(:,gsubK_ind)-mean(data(:,gsubK_ind)))/std(data(:,gsubK_ind));
    end
    data(:,gL_ind) =  (data(:,gL_ind)- mean(data(:,gL_ind)))/std(data(:,gL_ind));
    data(:,gM_ind) =  (data(:,gM_ind)- mean(data(:,gM_ind)))/std(data(:,gM_ind));    
    axis_rng = [-5 5];
else
    axis_rng = [0 4 ]; %[-1 1 ]
end

% Plot settings
scat_color = 'k';
[S,AX,BigAx,H,HAx] = plotmatrix(data,'.k');

[R P] = corrcoef(data)

ylabel(AX(1,1),'gNa')
ylabel(AX(2,1),'gK')
ylabel(AX(3,1),'gL')
ylabel(AX(4,1),'gM')

xlabel(AX(4,1),'gNa')
xlabel(AX(4,2),'gK')
xlabel(AX(4,3),'gL')
xlabel(AX(4,4),'gM')
set(S,'MarkerSize',8)


xlim(AX(1,1),[-1 1 ]); ylim(AX(1,1),[-1 1 ])
xlim(AX(2,1),axis_rng); ylim(AX(2,1),axis_rng)
xlim(AX(3,1),axis_rng); ylim(AX(3,1),axis_rng)
xlim(AX(4,1),axis_rng); ylim(AX(4,1),axis_rng)
xlim(AX(3,2),axis_rng); ylim(AX(3,2),axis_rng)
xlim(AX(4,2),axis_rng); ylim(AX(4,2),axis_rng)
xlim(AX(4,3),axis_rng); ylim(AX(4,3),axis_rng)

% histogram
xlim(HAx(1),[-1 1]);xlim(HAx(2),[-1 1]);xlim(HAx(3),[-1 1]);xlim(HAx(4),[-1 1])
ylim(HAx(1),[0 80]);ylim(HAx(2),[0 80]);ylim(HAx(3),[0 80]);ylim(HAx(4),[0 80])


% remove ticks
for ii = 1:N
    for jj = 1:N
        set(AX(ii,jj),'xtick',[])
        set(AX(ii,jj),'ytick',[])
        pbaspect(AX(ii,jj),[1 1 1])
    end
end

delete(AX(1,1:4))
delete(AX(2,2:4))
delete(AX(3,3:4))
delete(AX(4,4:4))
delete(HAx)

% Add linear fit
line_color = 'r';

hold(AX(2,1),'on') % gNa - gK
x = data(:,gsubNa_ind);
y = data(:,gsubK_ind);

[c,S] = polyfit(x,y,1);
y_est = polyval(c,x);
plot(AX(2,1),x,y_est,line_color,'LineWidth',1)
set(gca,'xtick',[]);set(gca,'ytick',[]) % remove ticks

hold(AX(3,1),'on') % gNa - gL
x = data(:,gsubNa_ind);
y = data(:,gL_ind);

[c,S] = polyfit(x,y,1)
y_est = polyval(c,x);
plot(AX(3,1),x,y_est,line_color,'LineWidth',1)


hold(AX(4,1),'on') % gNa - gM
x = data(:,gsubNa_ind);
y = data(:,gM_ind);

[c,S] = polyfit(x,y,1); % best fit for data in y
y_est = polyval(c,x);
plot(AX(4,1),x,y_est,line_color,'LineWidth',1)
set(gca,'xtick',[]);set(gca,'ytick',[]) % remove ticks


hold(AX(3,2),'on') % gK - gL
x = data(:,gsubK_ind);
y = data(:,gL_ind);


[c,S] = polyfit(x,y,1) 
y_est = polyval(c,x);
plot(AX(3,2),x,y_est,line_color,'LineWidth',1)
set(gca,'xtick',[]);set(gca,'ytick',[]) % remove ticks


hold(AX(4,2),'on') % gK - gM
x = data(:,gsubK_ind);
y = data(:,gM_ind);

c = polyfit(x,y,1); %original
y_est = polyval(c,x);
plot(AX(4,2),x,y_est,line_color,'LineWidth',1)
set(gca,'xtick',[]);set(gca,'ytick',[]) % remove ticks


hold(AX(4,3),'on') % gL - gM
x = data(:,gL_ind);
y = data(:,gM_ind);

c = polyfit(x,y,1); %original
y_est = polyval(c,x);
plot(AX(4,3),x,y_est,line_color,'LineWidth',1)
set(gca,'xtick',[]);set(gca,'ytick',[]) % remove ticks

set(gcf,'position',[680   530   482   448])


%% Plot histograms
binsz = (axis_rng(2) - axis_rng(1))/20;
binedges = axis_rng(1):binsz:axis_rng(2); 
figure(2)
h = histogram(data(:,gsubNa_ind));
set(h,'BinEdges', binedges); set(h,'FaceColor',scat_color );set(h,'FaceAlpha',1 )
ylim([0 300]); pbaspect([1  1   1])
xlabel('gsubNa')
set(gca,'xtick',[]);set(gca,'ytick',[]) % remove ticks

figure
h = histogram(data(:,gsubK_ind));
set(h,'BinEdges', binedges); set(h,'FaceColor',scat_color );set(h,'FaceAlpha',1 )
ylim([0 300]); pbaspect([1  1   1])
xlabel('gsubK')
set(gca,'xtick',[]);set(gca,'ytick',[]) % remove ticks

figure
h = histogram(data(:,gM_ind));
set(h,'BinEdges', binedges); set(h,'FaceColor',scat_color );set(h,'FaceAlpha',1 )
ylim([0 300]); pbaspect([1  1   1])
xlabel('gM')
set(gca,'xtick',[]);set(gca,'ytick',[]) % remove ticks

