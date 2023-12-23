% Examples of 2D and 3D grid search

clear all; close all; clc;

tic
%% LOAD TARGET RANGES:
load('output_target_values.mat')

%% LOAD Istim for firing rate:
load('fast_signal_mu_40_std=10p5_length=1p5sec.mat')

%% TRIAL PARAMETERS:
time = 1500; % (1 stime length in ms)
pre_stim = 500; 
dt = 0.05;  % time step for forward euler method

loop  = time/dt;   % no. of iterations of euler
t = (1:loop)*dt;
dc = 0; % (uA/cm2) - DC amplitude (NOT USED)

%% INITIALIZE EMPTY MATRICES
n_pts = 100; % # of points on an axis 

gsubNa = linspace(0,4,n_pts);
gsubK = linspace(0,4,n_pts);
gM = linspace(0,4,n_pts);
gL = linspace(1,4,n_pts);
% gAHP = linspace(0,4,n_pts);

gNa_fixed = 2; 
gM_fixed = 1.75;
gL_fixed = 2;
gAHP_fixed = 0.5; 

%**** for fig.3b
gK_fixed = 0;

[gna,gk] = ndgrid(gsubNa,gsubK); % gNaxgK grid
% [gna,gl] = ndgrid(gsubNa,gL); % gNaxgL grid
% [gna,gk,gm] = ndgrid(gsubNa,gsubK,gM); % gNaxgKxgM grid
% [gna,gk,gl] = ndgrid(gsubNa,gsubK,gL); % gNaxgKxgL grid


RATE = zeros(size(gna));
EE = zeros(size(gna));
RHEO = zeros(size(gna));
RIN = zeros(size(gna));
VREST = zeros(size(gna));

%% 2D grid search (e.g. gNa x gK)
%**** NOTE: x&y are flipped for visualization
for a = 1:length(gsubNa)
    for b = 1:length(gsubK)
        
        % Calculate firing rate 
        [spike] = ML_HH_adapt(time, dt, pre_stim, dc, i_signal, gsubNa(a),gsubK(b),gL_fixed,gM_fixed,gAHP_fixed);
        % replace ML_HH_adapt with ML_HH_adapt_ver2 to calculate FR AND
        % energy consumption rate simultaneously
        
        RATE(b,a) = sum(spike);
        
        % Calculate energy efficiency
        [min_Na,tot_Na] = calc_EE(gsubNa(a),gsubK(b),gL_fixed,gM_fixed,gAHP_fixed);
        EE(b,a) = min_Na/tot_Na;
        
        % Find rheobase
        [rheo, fmin, spike] = Rheo_ML_HH_JY(1100,dt,100,gsubNa(a),gsubK(b),gL_fixed,gM_fixed,gAHP_fixed);
        RHEO(b,a) = rheo;
        
        % Find Vrest and Rinput 
        [Vrest,Rinput]  = Vrest_Rinput_ML_HH(gsubNa(a),gsubK(b),gL_fixed,gM_fixed,gAHP_fixed);
        RIN(b,a) = Rinput;
        VREST(b,a) = Vrest;
    end
    
    % Print progress 
    fprintf(['Progress = ',num2str(a/length(gsubNa)*100),'%%\n'])
end

%% 3D grid search (e.g. gNa x gK x gM)
% for a = 1:length(gsubNa)
%     for b = 1:length(gsubK)
%         for c = 1:length(gL)
%             [spike] = ML_HH_adapt(time, dt, pre_stim, dc, i_signal, gsubNa(a),gsubK(b),gL_fixed,gM(c),gAHP_fixed);
%             RATE(a,b,c) = sum(spike);
%             
%             [min_Na,toat_Na] = calc_EE(gsubNa(a),gsubK(b),gL_fixed,gM(c),gAHP_fixed);
%             EE(a,b,c) = min_Na/tot_Na;
%             
%             [Vrest,Rinput]  = Vrest_Rinput_ML_HH(gsubNa(a),gsubK(b),gL_fixed,gM(c),gAHP_fixed);
%             RIN(a,b,c) = Rinput;
%             VREST(a,b,c) = Vrest;
%             
%             [rheo, fmin, spike] = Rheo_ML_HH_JY(1100,dt,100,gsubNa(a),gsubK(b),gL_fixed,gM(c),gAHP_fixed);
%             RHEO(a,b,c) = rheo;
%         end        
%     end
%     fprintf(['Progress = ',num2str(a/length(gsubNa)*100),'%%\n\n'])
% end


%% SAVE RESULTS
FileName=[datestr(now, 'yyyymmdd'),'_testing.mat'];
% save(FileName)

toc