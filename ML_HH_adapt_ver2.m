% LAST MODIFIED: 2019-11-25
% Called inside gridsearch.m and P_control.m
% Modified from ML_HH_adapt.m to calculate FR AND energy consumption rate
% in fig2 

function [spike, ATPn, ATPk] = ML_HH_adapt_ver2(time, dt, pre_stim, i_signal, gsubNa,gsubK,gleak,gM, gAHP)

graph = 0;

loop  = time/dt;   % no. of iterations of euler
i_stim = zeros(1,loop);
t_on = pre_stim/dt; %switch on a input after 100 ms

%% PARAMETERS
gNa = 20;  %same as g fast in the PLoS paper
gK = 20;
gL = gleak;

eNa = 50;
eK = -100;
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
w = zeros(1,loop); % Rectifying potassium
m = zeros(1,loop); % Sodium
z_M = zeros(1,loop); % M-type
z_AHP = zeros(1,loop); % AHP

spike = zeros(1,loop);
ref = zeros(1,loop);
DATA = zeros(6,loop);

%% Set initial values for the variables
V(1)= -70;
w(1)= 0.000025;
m(1) = 0; %HH style currents
spike(1) = 0;
ref(1) = 0;

% Euler method
for step=1:loop-1
    dV_dt = (i_signal(step) - gNa*m_infinity(V(step))*(V(step)-eNa) - gK*w(step)*(V(step)-eK) - gL*(V(step)-eL)...
        -gsubNa*m(step)*(V(step) - eNa) -gsubK*m(step)*(V(step)-eK)...
        -gM*z_M(step)*(V(step)-eK) -gAHP*z_AHP(step)*(V(step)-eK)...
        )/C;
    V(step+1) = V(step) + dt*dV_dt; %forward Euler equation
    %     end
    
    % Rectifying potassium
    dw_dt = phi*alpha_w(V(step))*(w_infinity(V(step))- w(step));
    w(step+1) = w(step) + dt*dw_dt; % forward Euler equation
    
    % HH style conductances
    dm_dt = (1-m(step))*alpha_HH(V(step)) - m(step)*beta_HH(V(step));
    m(step+1) = m(step) + dt*dm_dt; % forward Euler equation
    
    % M-type
    dz_M_dt = (1/(1+exp((beta_z_M-V(step))/gamma_z))-z_M(step))/tau_z;
    z_M(step+1) = z_M(step) + dt*dz_M_dt; % forward Euler equation
    
    % AHP
    dz_AHP_dt = (1/(1+exp((beta_z_AHP-V(step))/gamma_z))-z_AHP(step))/tau_z;
    z_AHP(step+1) = z_AHP(step) + dt*dz_AHP_dt; % forward Euler equation
    
    % Calculate currents
    INa = gNa*m_infinity(V(step))*(V(step)-eNa);
    INas = gsubNa*m(step)*(V(step) - eNa);
    IK = gK*w(step)*(V(step)-eK);
    IKs = gsubK*m(step)*(V(step)-eK);
    IKm = gM*z_M(step)*(V(step)-eK);
    IKahp = gAHP*z_AHP(step)*(V(step)-eK);
    
    spike(step) = (V(step) > 0).*(~ref(step));
    ref(step+1) = (V(step) > 0);
    DATA(:,step) = [INa;INas;IK;IKs; IKm; IKahp];
    
 
end
spike = spike(t_on:end);
DATA = DATA(:,t_on:end);

charge = dt*sum(DATA,2); %charge moved across each channel during one AP
chargeNa = abs(charge(1) + charge(2)); %sodium flux in nC/cm^2
chargeK = abs(charge(3) + charge(4) + charge(5) + charge(6)); %potassium flux in nC/cm^2
chargeNa = 1e-9*chargeNa; %convert to C/cm^2
chargeK = 1e-9*chargeK; %convert to C/cm^2
q = 1.60217662*1e-19;
nNa = chargeNa/q; %calculate number of Na+
nK = chargeK/q; %calculate number of K+ 
ATPn = nNa/3; % 3:2 ratio based on Na/K pump
ATPk = nK/2;

if graph == 1
    close all;
    figure(1)
    subplot(2,1,1)
    plot(t,V)
    ylabel('Voltage (mV)')
    ylim([-90 40])
    
    subplot(2,1,2)
    plot(t, i_signal,'b')
    hold on
    plot(t,i_noise,'r')
    ylabel('Current (uA/cm^{2})')
    
    xlabel('Time (msec)');
    set(gcf,'position',[ 3         558        1911         420])
end

end % closes function()

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
