% LAST MODIFIED: Nov04-2020
% Called inside gridsearch.m

function [Vrest,Rinput] =  Vrest_Rinput_ML_HH(gsubNa,gsubK,gleak,gM,gAHP)

graph = 0; 

time = 300; %time length in ms
dt = 0.05;  % time step for forward euler method
loop  = time/dt;   % no. of iterations of euler
i_stim = zeros(1,loop);
on = 100/dt; %switch on a input after 1000 ms
off = 200/dt;
amp = - 1;
i_stim(on:off) = amp; %determined to be minimum amount to evoke AP


%initialize constants
gNa = 20;  %same as g fast in the PLoS paper
eNa = 50;
gK = 20;
eK = -100;
gL = gleak;  % 2; % for 2 random parameters
eL = -70;
C = 2;
phi = 0.15;

% Initializing variable vectors
t = (1:loop)*dt;
V = zeros(1,loop);
w = zeros(1,loop);
m = zeros(1,loop);
z_M = zeros(1,loop); % M-type
z_AHP = zeros(1,loop); % AHP

% gAHP and gM parameters (from Prescott & Sejnowski, 2008)
tau_z = 100; % [msec] - same for gAHP and gM
beta_z_M = -35; % [mV]
beta_z_AHP = 0; % [mV]
gamma_z = 4; % [mV] - same for gAHP and gM

spike = zeros(1,loop);
ref = zeros(1,loop);
%   DATA = zeros(4,loop);

% Set initial values for the variables
V(1)= -70;
w(1)= 0.000025;
m(1) = 0; %HH style currents


spike(1) = 0;
ref(1) = 0;
Vrest = NaN;
Rinput = NaN;

% Euler method
for step=1:loop-1
    dV_dt = (i_stim(step)  - gNa*m_infinity(V(step))*(V(step)-eNa) - gK*w(step)*(V(step)-eK) - gL*(V(step)-eL)...
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
    %       DATA(:,step) = [INa;INas;IK;IKs];
end


if sum(spike) == 0
    Vrest = V(on);
    Rinput = ( V(150/dt) - Vrest )/amp;
end

if graph == 1
    figure(1)
    subplot(2,1,1)
    plot(t,V)
    ylabel('Voltage (mV)')
    
    subplot(2,1,2)
    plot(t,i_stim)
    ylim([-1.5 0.5])
    ylabel('I_{DC} (uA/cm^{2})'); xlabel('Time (msec)'); 
end


end

% steady state and decay functions for the gating variables
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
