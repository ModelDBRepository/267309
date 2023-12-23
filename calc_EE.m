% LAST MODIFIED: Nov04-2020
% Called inside gridsearch.m and P_control.m
function [min_Na,tot_Na] = calc_EE(gsubNa,gsubK,gleak,gM, gAHP)
%   Calculates energy efficiency = theoretical min : total Na+ load
%   theoretical min = Na+ load required for charging a pure capacitor 
%   to Vdiff = Vpeak - Vrest 

graph = 0;

time = 500;
dt = 0.05;  % time step for forward euler method
loop  = time/dt;   % no. of iterations of euler
i_stim = zeros(1,loop);
stim_onset = 200; %%100;
stim_length = 0.05;
t_on = stim_onset/dt; %switch on a input after 100 ms
reset_V = -40;

i_stim(t_on:t_on + stim_length/dt) = 0; % no current injection; voltage reset is used instead
  
%% PARAMETERS
gNa = 20;  %same as g fast in the PLoS paper
eNa = 50;
gK = 20;
eK = -100;
gL = gleak; %2;
eL = -70;
C = 2;
phi = 0.15;

% gAHP and gM parameters (from Prescott & Sejnowski, 2008)
tau_z = 100; % [msec] - same for gAHP and gM
beta_z_M = -35; % [mV]
beta_z_AHP = 0; % [mV]
gamma_z = 4; % [mV] - same for gAHP and gM
%% INITIALIZE EMPTY VECTORS
t = (1:loop)*dt;
V = zeros(1,loop);
w = zeros(1,loop);
m = zeros(1,loop);
z_M = zeros(1,loop); % M-type  
z_AHP = zeros(1,loop); % AHP

spike = zeros(1,loop);
ref = zeros(1,loop);
DATA = zeros(6,loop);

Q_tot = zeros(1,loop);
Q_Na = zeros(1,loop);
Q_NaS = zeros(1,loop);

%% SET INITIAL VALUES 
V(1)= -70;
w(1)= 0.000025;
m(1) = 0; %HH style currents
spike(1) = 0;
ref(1) = 0;
Q_Na(1) = 0; 
Q_NaS(1) = 0;
Q_tot(1) = 0; 

% Euler method
for step=1:loop-1
    if(step == t_on)
        V(step) = reset_V; % reset to 0
    else
%         dV_dt = (i_stim(step) - gNa*m_infinity(V(step))*(V(step)-eNa) - gK*w(step)*(V(step)-eK) - gL*(V(step)-eL)  -gsubNa*m(step)*(V(step) - eNa) -gsubK*m(step)*(V(step)-eK))/C;
        dV_dt = (i_stim(step) - gNa*m_infinity(V(step))*(V(step)-eNa) - gK*w(step)*(V(step)-eK) - gL*(V(step)-eL)...
            -gsubNa*m(step)*(V(step) - eNa) -gsubK*m(step)*(V(step)-eK)...
            -gM*z_M(step)*(V(step)-eK) -gAHP*z_AHP(step)*(V(step)-eK)...
            )/C;
        
        V(step+1) = V(step) + dt*dV_dt; %forward Euler equation
    end
    
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
    
    INa = gNa*m_infinity(V(step))*(V(step)-eNa);
    INas = gsubNa*m(step)*(V(step) - eNa);
    IK = gK*w(step)*(V(step)-eK);
    IKs = gsubK*m(step)*(V(step)-eK);
    
    IKm = gM*z_M(step)*(V(step)-eK);
    IKahp = gAHP*z_AHP(step)*(V(step)-eK);
    
    spike(step) = (V(step) > 0).*(~ref(step));
    ref(step+1) = (V(step) > 0);
    DATA(:,step) = [INa;INas;IK;IKs;IKm;IKahp];
    Q_tot(step+1) = Q_tot(step)+1e-9*abs(dt*sum(INa+INas)); % row - 1(Na) + 2(SubNa) at each time point 

end

if graph == 1
    figure('name','voltage_trace')
    t_range = (stim_onset-10)/dt: (stim_onset+10)/dt;
    
    subplot(2,1,1)
    plot(t(t_range),V(t_range))

    subplot(2,1,2) 
    Q_tot = Q_tot/1e-6; % unit converision (C -> uC)
    plot(t(t_range),Q_tot(t_range),'r') 
    hold on
    plot(t(t_range),ones(numel(t_range))*0.8975,'k')
    xlabel('Time (msec)'); ylabel('Q_{Na} (uC/cm^{2})')
end


if sum(spike)~=1 % spontaneous spiking
    min_Na = NaN;
    tot_Na = NaN;
else % not spontaneous spiking
    % Calculate theoretical minimum Na+ load:
    V_rest = V(t_on-1);
    V_peak = max(V);
    AP_amp = V_peak - V_rest;
    min_Na = 1e-9 * C * AP_amp; %[C/cm2]
    
    % Calculate total Na+ load
    DATA = DATA(:,t_on:end);
    V = V(t_on:end);
    pl = 10/dt; % 10msec cutoff
    V = V(1:pl);
    t = t(t_on:end);
    t = t(1:pl);
    DATA = DATA(:,1:pl);
    
    %For graphing:
    Q_tot = Q_tot(t_on:end);
    Q_tot = Q_tot(1:pl);

    charge = dt*sum(DATA,2); %charge moved across each channel during one AP
    chargeNa = abs(charge(1) + charge(2)); %sodium flux in nC/cm^2
    tot_Na = 1e-9*(chargeNa); %convert to C/cm^2
end

% print energy efficiency
EE = min_Na/tot_Na; % print energy efficiency

end

function [minf] = m_infinity(V)
beta_m = -1.2;
gamma_m = 14;
minf = 0.5*(1+tanh((V-beta_m)/gamma_m));
end

function [winf] = w_infinity(V)
beta_w = -10;
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
