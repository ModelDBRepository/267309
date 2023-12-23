%**************************************************************************
% The following code reproduces Figure 5D, which is essentially equal to
% a combination of Figure 5A-C. 
%**************************************************************************
clear all; close all; clc
load('output_target_values.mat')

%% Visualize SOSSs for Nin = 2
load('') % gridsearch data for firing rate, energy efficiency and input resistance for Nin = 2 (gNa&gK)

%%%%% ENERGY EFFICIENCY
t_EE_inds = find(EE>EE_Trange(1) & EE<EE_Trange(2)); 
[y_EE_2D,x_EE_2D] = ind2sub(size(EE),t_EE_inds); 
EE_sf2D = fit(gsubNa(x_EE_2D)',gsubK(y_EE_2D)','poly2');

%%%%%   RIN
t_RIN_inds = find(RIN>RIN_Trange(1) & RIN < RIN_Trange(2)); 
[y_RIN_2D,x_RIN_2D] = ind2sub(size(RIN),t_RIN_inds);
RIN_sf2D = fit(gsubNa(x_RIN_2D)',gsubK(y_RIN_2D)','poly2'); 

%%%%%   RATE
t_RATE_inds = find(RATE>RATE_Trange(1) & RATE < RATE_Trange(2)); 
[y_RATE_2D,x_RATE_2D] = ind2sub(size(RATE),t_RATE_inds);
RATE_sf2D = fit(gsubNa(x_RATE_2D)',gsubK(y_RATE_2D)','poly2'); 


figure('name','2D_fitted')
EE_line = plot(EE_sf2D,'g'); EE_line.LineWidth = 2;
hold on
RIN_line = plot(RIN_sf2D,'b'); RIN_line.LineWidth = 2; 
hold on
RATE_sf = plot(RATE_sf2D,'r');
RATE_sf.LineWidth = 2; 
legend off
axis([2 4 0 3])
xlabel('g_{Na}'); ylabel('g_{K}')
pbaspect([1 1 1])

%% Visualize SOSSs for Nin = 3
load('') % gridsearch data for firing rate, energy efficiency and input resistance for Nin = 3 (gNa&gK&gleak)

%%%%% ENERGY EFFICIENCY (green)
t_EE_inds = find(EE>EE_Trange(1) & EE<EE_Trange(2)); % find indicies for models within range
[x_EE_3D,y_EE_3D,z_EE_3D] = ind2sub(size(EE),t_EE_inds);
EE_sf3D = fit([gsubNa(x_EE_3D)',gsubK(y_EE_3D)'],gL(z_EE_3D)','poly21'); % fitting based on x and y
EE_sf3D_graph = fit([gsubK(y_EE_3D)',gL(z_EE_3D)'],gsubNa(x_EE_3D)','poly23'); % fitting based on y and z 

figure('name','surface_fits')
[y1,z1] = meshgrid(linspace(0,4,30),linspace(1,4,30)); %original
EE_sf = surf(EE_sf3D_graph(y1,z1),y1,z1,'FaceColor','g','FaceAlpha',0.3); %0.3
EE_sf.EdgeColor = 'none'; %[0 0.7 0]; 
EE_sf.EdgeAlpha = 0.6;
EE_sf.FaceLighting = 'gouraud';
EE_sf.BackFaceLighting = 'lit';
hold on


%%%%% RATE (red)
t_RATE_inds = find(RATE>RATE_Trange(1) & RATE < RATE_Trange(2));
[x_RATE_3D,y_RATE_3D,z_RATE_3D] = ind2sub(size(RATE),t_RATE_inds); % target indices
RATE_sf3D = fit([gsubNa(x_RATE_3D)',gsubK(y_RATE_3D)'],gL(z_RATE_3D)','poly23'); %original
RATE_sf3D_graph = fit([gsubK(y_RATE_3D)',gL(z_RATE_3D)'],gsubNa(x_RATE_3D)','poly23');

[x1,y1] = meshgrid(linspace(0,4,30),linspace(0,4,30)); %original
RATE_sf = surf(x1,y1,RATE_sf3D(x1,y1),'FaceColor','r','FaceAlpha',0.5); % 0.5
RATE_sf.EdgeColor = 'none'; %[170 0 0]./255; % red
RATE_sf.FaceLighting = 'gouraud';
hold on

%%%%% INPUT RESISTANCE (blue)
t_RIN_inds = find(RIN>RIN_Trange(1) & RIN < RIN_Trange(2));
[x_RIN_3D,y_RIN_3D,z_RIN_3D] = ind2sub(size(RIN),t_RIN_inds); % target indices
RIN_sf3D = fit([gsubNa(x_RIN_3D)',gsubK(y_RIN_3D)'],gL(z_RIN_3D)','poly21'); %original
RIN_sf3D_graph = fit([gsubK(y_RIN_3D)',gL(z_RIN_3D)'],gsubNa(x_RIN_3D)','poly23');

[x1,y1] = meshgrid(linspace(0,4,30),linspace(0,4,30)); %original
RIN_sf = surf(x1,y1,RIN_sf3D(x1,y1),'FaceColor','b','FaceAlpha',0.5); % 0.5 
RIN_sf.EdgeColor = 'none'; %[0 0 170]./255; % blue
RIN_sf.FaceLighting = 'gouraud';
% RIN_sf.BackFaceLighting = 'unlit';
hold on

%% INTERSECTION
[xL_RIN_RATE, yL_RIN_RATE, zL_RIN_RATE]= find_inter(RATE_sf3D, RIN_sf3D);
[yL_Rinput_EE, zL_Rinput_EE,xL_Rinput_EE,]= find_inter(RATE_sf3D_graph, EE_sf3D_graph); 
[yL_EE_Rheo, zL_EE_Rheo, xL_EE_Rheo]= find_inter(EE_sf3D_graph, RIN_sf3D_graph); 

line(xL_RIN_RATE,yL_RIN_RATE,zL_RIN_RATE,'Color','y'); 
line(xL_Rinput_EE,yL_Rinput_EE,zL_Rinput_EE,'Color','y','LineWidth',1); 
line(xL_EE_Rheo,yL_EE_Rheo,zL_EE_Rheo,'Color','y','LineWidth',1); 

axis([0 4 0 4 1 4]); xlabel('g_{Na}');ylabel('g_{K}');zlabel('g_{L}')
pbaspect([1 1 1]); view([-66 22])
