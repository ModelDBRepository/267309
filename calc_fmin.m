% Last modified: Jan06-2020 
% Called inside gridsearch.m

function [ISI,fmin,CoV] = calc_fmin(dt,spike)

spike_times = find(spike);
spike_times = spike_times*dt; % [msec]

ISI = diff(spike_times); % [msec]
avgISI = mean(ISI); % [msec]
avgISI = avgISI*1e-3; % [sec]
fmin = 1/avgISI;

CoV = std(ISI)/mean(ISI);

%% UNCOMMENT TO PLOT t vs. ISI
% figure('name','ISI plot')
% scatter(spike_times(1:end-1),ISI)
% ylabel('ISI (msec)')
% xlabel('t (msec)')
end