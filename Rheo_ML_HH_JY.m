% LAST MODIFIED: Nov04-2020
% Called inside gridsearch.m and P_control.m
% Returns rheobase for repetitive spiking (change sum(spikes)> 1 to 0 to
% find rheobase for transient spiking)

function [rheo,fmin,spike] = Rheo_ML_HH_JY(time, dt, pre_stim, gsubNa,gsubK,gleak,gM,gAHP)
%   This function simulates the Morris Lecar model for user specified input
%   current.
% time = 1300; %time length in ms
% dt = 0.05;  % time step for forward euler method
loop  = time/dt;   % no. of iterations of euler
i_stim = zeros(1,loop);
on = pre_stim/dt; %100/dt; %switch on a input after 1000 ms

%initialize constants
gNa = 20;  %same as g fast in the PLoS paper
eNa = 50;
gK = 20;
eK = -100;
gL = gleak;%2;
eL = -70;
C = 2;
phi = 0.15;

% gAHP and gM parameters (from Prescott & Sejnowski, 2008)
tau_z = 100; % [msec] - same for gAHP and gM
beta_z_M = -35; % [mV]
beta_z_AHP = 0; % [mV]
gamma_z = 4; % [mV] - same for gAHP and gM

% Initializing variable vectors
t = (1:loop)*dt;
V = zeros(1,loop);
w = zeros(1,loop);
m = zeros(1,loop);
z_M = zeros(1,loop); % M-type
z_AHP = zeros(1,loop); % AHP
spike = zeros(1,loop);
ref = zeros(1,loop);

% Set initial values for the variables
V(1)= -70;
w(1)= 0.000025;
m(1) = 0; %HH style currents
spike(1) = 0;
ref(1) = 0;



% NEED TO CHANGE RHEO_LIST to RH_RANGE
rh_range = [0 350];

%% Variables
distim = ( rh_range(2) - rh_range(1) );

prev_type = 0; % assume no spikes
curr_type = 0;

sign = 1; % -1 or +1
% initialized as 1. Need to increase istim at initial condition (istim = min rheo)

prev_istim = rh_range(2);
curr_istim = rh_range(1); % initialized as min rheo - 1 bc rheobase range includes min rheo, not min rheo + 1
new_istim = curr_istim;
rheo = 99999;
%%


while ( distim >= 1 && new_istim < rh_range(2)) % until res == 1
    
    
    amp = new_istim;
    curr_istim = new_istim;
    
    
    i_stim(on:end) = amp; % set stim amplitude
%     i_stim(1:end) = amp; % set stim amplitude

    % Euler method
    for step=1:loop-1
        dV_dt = (i_stim(step) - gNa*m_infinity(V(step))*(V(step)-eNa) - gK*w(step)*(V(step)-eK) - gL*(V(step)-eL)...
            -gsubNa*m(step)*(V(step) - eNa) -gsubK*m(step)*(V(step)-eK)...
            -gM*z_M(step)*(V(step)-eK) -gAHP*z_AHP(step)*(V(step)-eK)...
            )/C;
        V(step+1) = V(step) + dt*dV_dt; %forward Euler equation
        
        %Rectifying potassium
        dw_dt = phi*alpha_w(V(step))*(w_infinity(V(step))- w(step));
        w(step+1) = w(step) + dt*dw_dt; %forward Euler equation
        
        %HH style conductances
        dm_dt = (1-m(step))*alpha_HH(V(step)) - m(step)*beta_HH(V(step));
        m(step+1) = m(step) + dt*dm_dt; %forward Euler equation
        
        % M-type
        dz_M_dt = (1/(1+exp((beta_z_M-V(step))/gamma_z))-z_M(step))/tau_z;
        z_M(step+1) = z_M(step) + dt*dz_M_dt; % forward Euler equation
        
        % AHP
        dz_AHP_dt = (1/(1+exp((beta_z_AHP-V(step))/gamma_z))-z_AHP(step))/tau_z;
        z_AHP(step+1) = z_AHP(step) + dt*dz_AHP_dt; % forward Euler equation
        
        spike(step) = (V(step) > 0).*(~ref(step));
        ref(step+1) = (V(step) > 0);
    end
    
    
    if sum(spike)>0 && curr_istim == rh_range(1) %spontaneous spiking
        rheo = 0;
        break;
    end
    
    % throw away the first 50 msec 
%     spike = spike(pre_stim/dt:end); %stabilized spikes
    
    % Store activity pattern (spks vs. no spks)
%     curr_type = sum(spike)>1; % repetitive spiking
        curr_type = sum(spike)>0; % transient spiking (1 spike)
    
    % Store switch in activity
    d_type = prev_type - curr_type;
    
    % 1) Specify changes in spking pattern(d_type)
    % 2) if resolution = 1 uA/cm2, then set rheobase based on 1),
    % and terminate loop
    
    
    if d_type ~= 0 % spks to no spikes (or vice versa)
        if distim == 1 % found rheo at res of 1 uA/cm2
            switch d_type
                case -1 % spks -> no spks
                    rheo = curr_istim;
%                                         fprintf(['curr_stim = ',num2str(sum(spike))])
                case 1 % no spks -> spks
                    rheo = prev_istim;
%                                         fprintf(['prev_stim = ',num2str(sum(prev_spike))])
                                        spike = prev_spike;
            end
            
            break; % exit loop
        end % closes distim == 1
        sign = -sign;
    end % end of if d_type ~= 0
    % NOTE: if d_type = 0, no need to change sign.
    
    % Update istim using (updated) sign
    distim = ceil(distim/2); % decrement distim ( never equal to 0)
    new_istim = curr_istim + sign*distim;
    
    prev_istim = curr_istim; % store curr_istim
    prev_type = curr_type; % store curr_spks
    prev_spike = spike;  % store previous spike train
    
    
    
end

[ISI,fmin] = calc_fmin(dt,spike);
end

% steady state and decay functions for the gating variables
function [minf] = m_infinity(V)
beta_m = -1.2;
gamma_m = 14; % Ratte et al. 2014
% gamma_m = 18; % (Prescott & Sejnowski, 2008)
minf = 0.5*(1+tanh((V-beta_m)/gamma_m));
end

function [winf] = w_infinity(V)
beta_w = -10; % Ratte et al. 2014
% beta_w = 0; % (Prescott & Sejnowski, 2008)
gamma_w = 10;
winf = 0.5*(1+tanh((V-beta_w)/gamma_w));
end

function [aw] = alpha_w(V)
beta_w = -10;
gamma_w = 10;
aw = cosh((V-beta_w)/(2*gamma_w));
end

function [a] = alpha_HH(V)
V_a = -24;
s_a = -17;
k_a = 1;
a = (V-V_a)/s_a;
a = k_a*a/(exp(a) - 1);
end

function [b] = beta_HH(V)
V_b = -24;
s_b = -17;
k_b = 1;
b = (V-V_b)/s_b;
b = k_b*exp(b);
end

