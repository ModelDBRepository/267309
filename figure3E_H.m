%**************************************************************************
% The following code reproduces Figure 3E-H, but can be modified to 
% reproduce contour plots and isooutput contour vs. other output plots 
% in Figure 3A-D
%**************************************************************************
clear all; close all; clc
%% Load target range: 
load('output_target_values.mat')

%% Fig.2a - Energy consumption rate contours + iso-FR line
load('') % gridsearch data for energy consumption (gridsearch.m)
figure('name','Figure 3E')
[C,h] = contourf(gsubNa,gsubK,totATP,50);
set(h,'LineColor','none')
set(gca,'TickDir','out','FontSize',15); box off
pbaspect([1 1 1])
xlabel('g_{Na}');ylabel('g_{K}')
cb = colorbar;ylabel(cb,'Energy Consumption Rate (10^{4} ATP/cm2 s)')
colormap(jet)

load('')% gridsearch data for firing rate (gridsearch.m)
%**** iso-firing rate contour
t_RATE_inds = find(RATE>RATE_Trange(1) & RATE<RATE_Trange(2)); % find indicies for models within range
[y_RATE_2D,x_RATE_2D] = ind2sub(size(RATE),t_RATE_inds);
RATE_2D = fit(gsubNa(x_RATE_2D)',gsubK(y_RATE_2D)','poly2'); % linear fit 
hold on
isotRATE_line= plot(RATE_2D);
isotRATE_line.LineWidth = 3; 
isotRATE_line.Color = [0.7 0.7 0.7];
legend('Energy consumption rate','Firing Rate = 40 spk/s')
axis([0 4 0 4]); xlabel('g_{Na}');ylabel('g_{K}')


%% Fig.2b - Energy consumption rate vs. FR 
t_EE_inds = find(EE>EE_Trange(1) & EE<EE_Trange(2)); 
[y_EE_2D,x_EE_2D] = ind2sub(size(EE),t_EE_inds); 
EE_sf2D = fit(gsubNa(x_EE_2D)',gsubK(y_EE_2D)','poly2');

figure('name','Figure 3F')
iso_values = zeros(numel(x_RATE_2D),1);
for i = 1: numel(x_RATE_2D)
    iso_values(i) = totATP(y_RATE_2D(i),x_RATE_2D(i));
end
plot(gsubNa(x_RATE_2D),iso_values)


% find min and max values for the same x - to draw thickness
for i = 1:numel(x_RATE_2D)
    same_x = find(x_RATE_2D(i)==x_RATE_2D);
    max_x(i) = x_RATE_2D(i);
    max_y(i) = max(iso_values(same_x));
end
hold on; line(gsubNa(max_x),max_y)

for i = 1:numel(x_RATE_2D)
    same_x = find(x_RATE_2D(i)==x_RATE_2D);
    min_x(i) = x_RATE_2D(i);
    min_y(i) = min(iso_values(same_x));
end
hold on; line(gsubNa(min_x),min_y)

xlim([gsubNa(x_RATE_2D(1)) gsubNa(x_RATE_2D(end))]);ylim([1.5e14 3e14])
xlabel('g_{Na} along firing rate = 40 spk/s'); ylabel('Energy Consumption Rate')
set(gca,'TickDir','out','FontSize',15); box off 
set(gcf, 'position',[760   589   274   231])
box off

%% Fig.2c - EE contours 
figure('name','Figure 3G')
[C,h] = contourf(gsubNa,gsubK,EE,50);
set(h,'LineColor','none')
cb = colorbar;ylabel(cb,'Energy Efficiency')
colormap(flipud(jet))

set(gca,'TickDir','out','FontSize',15); box off
xlabel('g_{Na}');ylabel('g_{K}')
pbaspect([1 1 1])
%% Fig.2d - EE vs. FR vs. Rheobase
figure('name','Figure 3H')
subplot(2,1,1)
iso_values = zeros(numel(x_EE_2D),1);
for i = 1: numel(x_EE_2D)
    iso_values(i) = RATE(y_EE_2D(i),x_EE_2D(i));
end
plot(gsubNa(x_EE_2D),iso_values)


% find min and max values for the same x - to draw thickness
for i = 1:numel(x_EE_2D)
    same_x = find(x_EE_2D(i)==x_EE_2D);
    max_x(i) = x_EE_2D(i);
    max_y(i) = max(iso_values(same_x));
end
hold on; line(gsubNa(max_x),max_y)

for i = 1:numel(x_EE_2D)
    same_x = find(x_EE_2D(i)==x_EE_2D);
    min_x(i) = x_EE_2D(i);
    min_y(i) = min(iso_values(same_x));
end
hold on; line(gsubNa(min_x),min_y)
xlim([gsubNa(x_EE_2D(1)) gsubNa(x_EE_2D(end))]);ylim([0 50])
ylabel('Firing Rate (spk/s)')
set(gca,'TickDir','out','FontSize',15); box off


subplot(2,1,2)
iso_values = zeros(numel(x_EE_2D),1);
for i = 1: numel(x_EE_2D)
    iso_values(i) = RHEO(y_EE_2D(i),x_EE_2D(i));
end
plot(gsubNa(x_EE_2D),iso_values)
ylim([0 50]); xlim([gsubNa(x_EE_2D(1)) gsubNa(x_EE_2D(end))])

% find min and max values for the same x - to draw thickness
for i = 1:numel(x_EE_2D)
    same_x = find(x_EE_2D(i)==x_EE_2D);
    max_x(i) = x_EE_2D(i);
    max_y(i) = max(iso_values(same_x));
end
hold on; line(gsubNa(max_x),max_y)

for i = 1:numel(x_EE_2D)
    same_x = find(x_EE_2D(i)==x_EE_2D);
    min_x(i) = x_EE_2D(i);
    min_y(i) = min(iso_values(same_x));
end
hold on; line(gsubNa(min_x),min_y)
ylim([0 50]); xlim([gsubNa(x_EE_2D(1)) gsubNa(x_EE_2D(end))])

xlabel('gNa along energy Efficiency = 0.235'); ylabel('Rheobase')
set(gcf,'position',[ 477   368   362   450])
set(gca,'TickDir','out','FontSize',15); box off