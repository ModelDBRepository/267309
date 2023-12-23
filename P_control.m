%**************************************************************************
% The following code runs proportional feedback control without noise
%**************************************************************************

function [post_gi, success, gi_traj, M_traj, final_M, n_it] = P_control(N_inputs, N_outputs, rate_Trange, EE_min, rate_t_reg, EE_t_reg,i_signal,i_gsubNa, i_gsubK,i_gL, i_gM, i_gAHP)
graph = 0;
no_EEmax = numel(EE_min)==1; % 1 - only minimum requirement (range); 2 - fixed with tolerance
max_it = 300;  % max number of iterations

%% TRIAL PARAMETERS:
time = 1500; % 1.5 sec
pre_stim = 500; 
dt = 0.05;  % time step for forward euler method
dc = 0; % (uA/cm2) - DC amplitude (NOT USED)
rate_T = mean(rate_Trange);
EE_T = mean(EE_min);

g_lim = [0 4]; % max and min values for conductances

%% Output indices
FR_ind = 1;
EE_ind = 2;

%% Conductance indices
gsubNa_ind = 1;
gsubK_ind = 2;
gL_ind = 3;
gM_ind = 4;
gAHP_ind = 5;



%% Initialize empty arrays
% Conductances (gi)
gsubNa = NaN(max_it,1);
gsubK = NaN(max_it,1);
gleak = NaN(max_it,1);
gM = NaN(max_it,1);
gAHP = NaN(max_it,1);

% Outputs
if N_outputs(EE_ind) ==0 % not regulated
    EE = zeros(max_it,1);
    EE_error = zeros(max_it,1);
else
    EE = NaN(max_it,1);
    EE_error = NaN(max_it,1);
end


if N_outputs(FR_ind) ==0 % not regulated
    rate = zeros(max_it,1);
    rate_error = zeros(max_it,1);
else
    rate = NaN(max_it,1);
    rate_error = NaN(max_it,1);
end



% Output errors



%% Store initial gi
gsubNa(1) = i_gsubNa;
gsubK(1) = i_gsubK;
gleak(1) = i_gL;
gM(1) = i_gM;
gAHP(1) = i_gAHP;

%% Generate Test Trial
[test] = ML_HH_adapt(time, dt, pre_stim, i_signal, i_gsubNa,i_gsubK,i_gL,i_gM,i_gAHP);

% Calculate rate
curr_rate = sum(test);

% Calculate Energy Efficiency
[min_Na,tot_Na] = calc_EE( i_gsubNa,i_gsubK,i_gL,i_gM,i_gAHP);
curr_EE = min_Na/tot_Na;


rate(1) = curr_rate;
EE(1) = curr_EE;

fprintf(' - Initial values:\n')
fprintf(['          gsubNa = ',num2str(i_gsubNa),' ,gsubK = ',num2str(i_gsubK),' ,gleak = ',num2str(i_gL),' ,gM = ',num2str(i_gM),' ,gAHP = ',num2str(i_gAHP),'\n\n'])
fprintf(['          Rate = ', num2str(rate(1)),  '\n'])
fprintf(['          Energy efficiency = ', num2str(EE(1)),  '\n\n'])


n_avg = 5; %10; % number of iterations during which outputs need to be on target
on_target = 0; % counts # of it at which output was within target range
prev_on_target = 0;

n_it = 1; % number of iterations b/f reaching target rheobase

while on_target < n_avg && n_it <= max_it-1 % compensation needed
    % Terminate if EE is NaN
    if isnan(curr_EE) %|| curr_rate == 0
        fprintf('Error in calc_EE. Need to adjust Vreset\n')
        break;
    else
        % Compute error (P)
        if N_outputs(FR_ind) == 1 % only if rate is being regulated
            rate_error(n_it) = curr_rate - rate_T;
        end
        
        if N_outputs(EE_ind) == 1 % Energy Efficiency
            EE_error(n_it) = curr_EE - EE_T;
        end
        
        
        if no_EEmax == 1 % only minimum requirement (>=)
            if EE_error(n_it) >= 0 % if more efficient than minimum
                EE_error(n_it) = 0;  % no error
            end
        end
        
        
        %% Update gi's based on the errors
        if N_inputs(gsubNa_ind) == 0 % Subthreshold Na
            gsubNa(n_it+1) = i_gsubNa;
        else
            dgsubNa_it(n_it) = rate_error(n_it)/rate_t_reg(gsubNa_ind)*N_outputs(FR_ind) + EE_error(n_it)/EE_t_reg(gsubNa_ind)*N_outputs(EE_ind);
            gsubNa(n_it+1) = gsubNa(n_it) + dgsubNa_it(n_it);
        end
        
        if N_inputs(gsubK_ind) == 0 % Subthreshold K
            gsubK(n_it+1) = i_gsubK;
        else
            dgsubK_it(n_it) = rate_error(n_it)/rate_t_reg(gsubK_ind)*N_outputs(FR_ind) + EE_error(n_it)/EE_t_reg(gsubK_ind)*N_outputs(EE_ind);
            gsubK(n_it+1) = gsubK(n_it) + dgsubK_it(n_it);
        end
        
        if N_inputs(gL_ind) == 0
            gleak(n_it+1) = i_gL;
        else
            dgleak_it(n_it) = rate_error(n_it)/rate_t_reg(gL_ind)*N_outputs(FR_ind) + EE_error(n_it)/EE_t_reg(gL_ind)*N_outputs(EE_ind);
            gleak(n_it+1) = gleak(n_it) + dgleak_it(n_it);
        end
        
        if N_inputs(gM_ind) == 0
            gM(n_it+1) = i_gM;
        else
            dgM_it(n_it) = rate_error(n_it)/rate_t_reg(gM_ind)*N_outputs(FR_ind) + EE_error(n_it)/EE_t_reg(gM_ind)*N_outputs(EE_ind);
            gM(n_it+1) = gM(n_it) + dgM_it(n_it);
        end
        
        if N_inputs(gAHP_ind) == 0
            gAHP(n_it+1) = i_gAHP;
        else
            dgAHP_it(n_it) = rate(n_it)/rate_t_reg(gAHP_ind)*N_outputs(FR_ind) + EE_error(n_it)/EE_t_reg(gAHP_ind)*N_outputs(EE_ind);
            gAHP(n_it+1) = gAHP(n_it) + dgAHP_it(n_it);
        end
        
        %% Conductances capped at 0 (no neg. conductances)
        if gsubNa(n_it+1) < 0
            gsubNa(n_it+1) = 0;
        else if gsubNa(n_it+1) > g_lim(2)
                gsubNa(n_it+1) = g_lim(2);
            end
        end
        if gsubK(n_it+1) < 0
            gsubK(n_it+1) = 0;
        else if gsubK(n_it+1) > g_lim(2)
                gsubK(n_it+1) = g_lim(2);
            end
        end
        if gleak(n_it+1) < 1 %**************
            gleak(n_it+1) = 1;
        else if gleak(n_it+1) > g_lim(2)
                gleak(n_it+1) = g_lim(2);
            end
        end
        
        if gM(n_it+1) < 0 %**************
            gM(n_it+1) = 0;
        else if gM(n_it+1) > g_lim(2)
                gM(n_it+1) = g_lim(2);
            end
        end
        
        if gAHP(n_it+1) < 0 %**************
            gAHP(n_it+1) = 0;
        else if gAHP(n_it+1) > g_lim(2)
                gAHP(n_it+1) = g_lim(2);
            end
        end
        
        %% Update output values rate, R, P, and EE
        % Generate Test Trial
        [test] = ML_HH_adapt(time, dt, pre_stim, i_signal, gsubNa(n_it+1),gsubK(n_it+1),gleak(n_it+1),gM(n_it+1),gAHP(n_it+1));
        
        % Calculate firing rate
        curr_rate = sum(test);
        rate(n_it+1) = curr_rate;
        
        % Calculate energy efficiency
        [min_Na,tot_Na] = calc_EE(gsubNa(n_it+1),gsubK(n_it+1), gleak(n_it+1), gM(n_it+1),gAHP(n_it+1));
        curr_EE = min_Na/tot_Na;
        EE(n_it+1) = round_JY(curr_EE,4);
        
        %% Check if current output values meet target
        M_on_target = zeros(1,numel(N_outputs));
        
        if N_outputs(FR_ind) ~= 0 % rate
            M_on_target(FR_ind) = (curr_rate >= rate_Trange(1) && curr_rate <= rate_Trange(2));
        end
        
        if N_outputs(EE_ind) ~= 0 % EE
            if numel(EE_min) == 1 % When EE > = EEmin
                M_on_target(EE_ind) = EE(n_it+1) >= EE_min;
            else % EE is fixed (i.e. with tolerance)
                M_on_target(EE_ind) = (curr_EE >= EE_min(1) && curr_EE <= EE_min(2) );
            end
        end
        
        if isequal(N_outputs, M_on_target) % if all outputs are on target
            on_target = on_target + 1;
        else
            on_target = 0;
        end
        n_it  = n_it + 1;  % increment number of iterations
    end % end of test spontaneous spiking
end % end of while loop


if isnan(curr_EE)
    post_gi = [NaN NaN NaN NaN NaN];
    final_M = [NaN NaN];
else
    % AVERAGED final gi values
    avg_gsubNa = mean(gsubNa(n_it-n_avg:n_it));
    avg_gsubK = mean(gsubK(n_it-n_avg:n_it));
    avg_gleak = mean(gleak(n_it-n_avg:n_it));
    avg_gM = mean(gM(n_it-n_avg:n_it));
    avg_gAHP = mean(gAHP(n_it-n_avg:n_it));
    % NOTE = gsubNa(1) stores initial value. gsubNa(2) stores the value after
    % 1st interation
    
    % AVERAGED final output values
    avg_rate = mean(rate(n_it-n_avg:n_it));
    avg_EE = mean(EE(n_it-n_avg:n_it));
    
    % Store avg values
    post_gi = [avg_gsubNa avg_gsubK avg_gleak avg_gM avg_gAHP];
    final_M = [avg_rate avg_EE];
    
end

fprintf(' - Final values: \n')
fprintf(['     gsubNa = ',num2str(post_gi(gsubNa_ind)),' ,gsubK = ',num2str(post_gi(gsubK_ind)),' ,gleak = ',num2str(post_gi(gL_ind)),', gM = ',num2str(post_gi(gM_ind)),', gAHP = ',num2str(post_gi(gAHP_ind)),'\n'])

fprintf(['     Rate: ' , num2str(curr_rate),'\n'])
fprintf(['     EE: ' , num2str(curr_EE),'\n\n'])


if on_target >= n_avg
    fprintf(['\n\nSuccess (number of iterations = ',num2str(n_it), ')\n'])
    success = 1;
else
    fprintf(['\n\nFailure(number of iterations = ',num2str(n_it), ')\n'])
    success = 0;
end

on_target

fprintf('\n\n*************************************************\n')


%% Store trajectories of gsubNa, gsubK and gleak
gi_traj = [gsubNa gsubK gleak gM gAHP];
M_traj = [rate EE];

%% Graph
if (graph == 1)
    %% settings
    grouping = 1; % plot success or failures together
    success_only = 2; % 2- plot success & failures together
    fontsize = 15;
    %% plots
    if (success_only == 1 && success==1) || (success_only == 0 && success==0 || (success_only == 2))
        if grouping == 1
            figure(1)
        else
            figure %figure(1)
        end
        hold on
        subplot(3,1,1) % conductances
        p1_Na = semilogx(1:max_it,gi_traj(:,1),'Color',[0.9290, 0.6940, 0.1250]); % gNa
        hold on
        p1_K = semilogx(1:max_it,gi_traj(:,2),'Color',[0.4940, 0.1840, 0.5560]); % gK
        hold on
        p1_M = semilogx(1:max_it,gi_traj(:,4),'b'); % gM
        
        % axis settings
        set(gca,'TickDir','out'); box off
        ylabel('g_{i}')
        xlim([1 max_it])
        set(gca,'xscale','log','FontSize',fontsize); % switch to logscale
        
        
        subplot(3,1,2) % Rate
        hold on
        p2 = semilogx(1:max_it,M_traj(:,FR_ind),'r'); % firing rate
        %     xlim([1 n_it])
        line([1 n_it],[rate_Trange(1) rate_Trange(1)],'Color','k','LineWidth',1,'LineStyle','--')
        line([1 n_it],[rate_Trange(2) rate_Trange(2)],'Color','k','LineWidth',1,'LineStyle','--')
        line([1 n_it],[mean(rate_Trange) mean(rate_Trange)],'Color','r','LineWidth',1,'LineStyle','--')
        
        % axis settings
        set(gca,'TickDir','out'); box off
        ylabel('Rate (spk/s)')
        xlim([1 max_it])
        set(gca,'xscale','log','FontSize',fontsize); % switch to logscale
        
        subplot(3,1,3) % energy efficiency
        hold on
        p3 = semilogx(1:max_it,M_traj(:,EE_ind),'g'); % EE
        
        if numel(EE_min)==1
            line([1 n_it],[EE_min EE_min],'Color','k','LineWidth',1,'LineStyle','--')
        else
            line([1 n_it],[EE_min(1) EE_min(1)],'Color','k','LineWidth',1,'LineStyle','--')
            line([1 n_it],[EE_min(2) EE_min(2)],'Color','k','LineWidth',1,'LineStyle','--')
        end
        line([1 n_it],[mean(EE_min) mean(EE_min)],'Color','g','LineWidth',1,'LineStyle','--')
        hold off
        
        % axis settings
        xlabel('# of iterations');ylabel('EE (%)');
        set(gca,'TickDir','out'); box off
        xlim([1 200])
        set(gca,'xscale','log','FontSize',fontsize); % switch to logscale
        set(gcf,'position',[10    86   560   892])
        
        
    end
end % end of if graph==1

end % end of function