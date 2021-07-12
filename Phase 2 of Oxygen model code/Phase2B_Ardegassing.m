%% 2019/08/05 Phase 2B of argon degassing -- Thermal evolution
% This file read "out_phase2a.dat" and select successful solutions
% according to the mantle potential temperature misfit of your choice.
% Successful solutions are stored in file "out_phase2b.dat".
%
% Meng Guo, Yale University
% Summer, 2019

clear all;

rng('shuffle');% to avoid have same random results everytime

% Open output file
fileID = fopen('out_phase2b.dat','w');

%% Read file "out_phase1a.dat"
data = load('out_phase2a.dat');
% the variables in 'out_phase2a.dat' are in following order:
% 1.Krw_factor, 2.kappa_r, 3.kappa_g, 4.Rs, 5.Rp, ...
% 6.ts, 7.RMSE_F, 8.RMSE_S, 9.H_BSE_tp, 10.H_cc_tp,...
% 11.Q_total_tp, 12.Qc_tp, 13.d_Qc, 14.misfit_Tp

%% Define the misfit of mantle potential temperature
RMSE_Tp = data(:,14);

%% Select successful solutions for plotting results
[ndata,nparam] = size(data);% the size of the data file
j = 0;% counter for the number of successful solutions

crit_misfit_Tp = 1;% set the selection criterion for misfit
success_Tp =[];% empty matrix for saving successful solutions

for i = 1:ndata
    if (RMSE_Tp(i)<crit_misfit_Tp)
        j = j+1;
        for k = 1:nparam
            success_Tp(j,k) = data(i,k);
        end % for k = 1:parameter
        
        % Save all the variables and RMSEs
        fprintf(fileID,['%6g %6g %6g %6g %6g %6g %6g %6g %6g %6g %6g %6g %6g %6g\n'], ...
            success_Tp(j,1), success_Tp(j,2), success_Tp(j,3), success_Tp(j,4),...
            success_Tp(j,5), success_Tp(j,6), success_Tp(j,7), success_Tp(j,8),...
            success_Tp(j,9), success_Tp(j,10),success_Tp(j,11), success_Tp(j,12),...
            success_Tp(j,13), success_Tp(j,14));
    end % if (RMSE_Tp(i)<crit_misfit_Tp)
    
end % for i = 1:ndata

fclose(fileID);

%% Set the time period
tmax = 4.567;% the age of solar system, in unit Ga
dt = 0.001;% length of each timestep
t = 0:dt:tmax;
nt = length(t);% number of timesteps
t = t';

%% Herzburg et all.(2010) data for mantle potential temperature
% load the Ti from Herzberg et al., 2010
data_Tp = xlsread('Herz data.xlsx');
% set the anchor points for misfit calculation
[Tp_anchorHerz1,Tp_anchorHerz2,Tp_anchorHerz3,Tp_anchorHerz4,...
    t_anchorHerz1,t_anchorHerz2,t_anchorHerz3,t_anchorHerz4,...
    t_Herz,Tp_Herz] = load_Tp_fun(data_Tp);

%% Histogram of parameter distributions and correlations
figure(1);
subplot(3,4,1);
histogram(success_Tp(:,1)); xlabel('init Krw factor'); ylabel('count');
subplot(3,4,2);
histogram(success_Tp(:,2)); xlabel('Kappa_r');
subplot(3,4,3);
histogram(success_Tp(:,3)); xlabel('Kappa_g');
subplot(3,4,4);
histogram(success_Tp(:,4)); xlabel('Rs');
subplot(3,4,5);
histogram(success_Tp(:,5)); xlabel('Rp');ylabel('Count')
subplot(3,4,6);
histogram(success_Tp(:,6)); xlabel('Ts');
subplot(3,4,7);
histogram(success_Tp(:,9)); xlabel('Hbse(t_p)');
subplot(3,4,8);
histogram(success_Tp(:,10)); xlabel('Hcc(t_p)');
subplot(3,4,9);
histogram(success_Tp(:,11)); xlabel('Q(t_p)');ylabel('Count')
subplot(3,4,10);
histogram(success_Tp(:,12)); xlabel('Q_c(t_p)');
subplot(3,4,11);
histogram(success_Tp(:,13)); xlabel('dQ_c');

% calculate the correlation matrix
correlation_Tp = corrcoef(success_Tp);

figure(2);
subplot(2,3,1);
plot(success_Tp(:,1),success_Tp(:,5),'b.');
xlabel('init Krw fac'); ylabel('Rp');
subplot(2,3,2);
plot(success_Tp(:,2),success_Tp(:,4),'b.');
xlabel('kappa_r'); ylabel('Rs');
subplot(2,3,3);
plot(success_Tp(:,2),success_Tp(:,5),'b.');
xlabel('kappa_r'); ylabel('Rp');
subplot(2,3,4);
plot(success_Tp(:,9),success_Tp(:,11),'b.');
xlabel('Hbse(t_p)'); ylabel('Q(t_p)');
subplot(2,3,5);
plot(success_Tp(:,11),success_Tp(:,12),'b.');
xlabel('Q(t_p)'); ylabel('Q_c(t_p)');

%% Calculate the corresponding thermal evolution
% Constants for performing crustal growth model
Mcp = 2.09e22;% mass of continental crust at present-day, in unit kg

% constants for thermal evolution model
p_K = 2.79e-5;% heat production for K, in unit W/kg
p_U235 = 5.69e-4;% unit heat production for 235U, in unit W/kg
p_U238 = 9.37e-5;% unit heat production for 238U, in unit W/kg
p_Th = 2.69e-5;% unit heat production for Th, in unit W/kg
Th_U = 4;% Th/U ratio
K_U = 1.27e4;% K/U ratio
K40_K = 1.28e-4;% 40K/K ratio
U238_U = 0.9927;% 238U/U ratio
Mmp = 4.015e24;% mass of the mantle at present-day, in unit kg
Na = 6.022E+23;% Avogadro constant
mK_40 = 40;% K's atomic weight
K_factor = K_U*K40_K;% 40K relative to U
heat_factor = U238_U*p_U238 + (1-U238_U)*p_U235 + Th_U*p_Th + K_factor*p_K;
Tkg_atoms = 1e12*1e3*Na/mK_40;% convert unit from Tkg to number of atoms
% Q_Haden = 36; % Q from 0 to 0.5Ga, in unit TW
Ti_tp = 1350; % present-day mantle potential temperature, in unit degree C
V_tp = 5; % present-day plate velocity, in unit cm/yr (Parsons,1981)
rhom = 3300; % the average density of mantle, in unit kg/m3
dTdP = 1.54e-8; % dT/dp, in unit K/Pa (Korenaga et al., 20020)
type = 2;% set the type to be constant Q scaling law

% create empty matrixes to store results
[size_Tp,nparameter ]= size(success_Tp);
Ti = nan(nt,size_Tp);
Q = nan(nt,size_Tp);
H = nan(nt,size_Tp);
V = nan(nt,size_Tp);
Z = nan(nt,size_Tp);

for l=1:size_Tp
    if mod(l,10)==0
        disp(['l=' num2str(l) ' of ' num2str(size_Tp)]);% keep track of the calculation
    end

    % independent variables
    Krw_factor_model = success_Tp(l,1);
    kappa_r_model = success_Tp(l,2);
    kappa_g_model = success_Tp(l,3);
    Rs_model = success_Tp(l,4);
    Rp_model = success_Tp(l,5);
    ts_model = success_Tp(l,6);
    H_BSE_tp_model = success_Tp(l,9);
    H_cc_tp_model = success_Tp(l,10);
    Q_total_tp_model = success_Tp(l,11);
    Qc_tp_model = success_Tp(l,12);
    d_Qc_model = success_Tp(l,13);
    
    % calculate the dependent variables
    Krw_s_model = Rs_model * Krw_factor_model;% initial Krw_factor
    Q_tp_model = Q_total_tp_model - H_cc_tp_model;
    K40_BSE_tp_model = (K_factor * H_BSE_tp_model) / heat_factor * Tkg_atoms;
    K40_CC_observe_model = (K_factor * H_cc_tp_model) / heat_factor * Tkg_atoms;
    
    % calculate the corresponding crustal growth pattern
    [Mc_model,Mdd_model,Mud_model,Krw_model_first] = CC_growth_fun1(t,ts_model,tmax,...
        Mcp,kappa_g_model,Rp_model,Rs_model,kappa_r_model,Krw_s_model);
    
    [Qc_model,Qc_backward] = Qc_backward_fun(Qc_tp_model,d_Qc_model,nt,t);
    
    % calculate the corresponding thermal evolution
    % Since we are calculating H backward in time, Mc need to be backward in time as well
    Mc_backward = flipud(Mc_model);
    [Ti_backward, Q_backward, H_backward,V_backward, Z_backward] = ...
        Thermal_history_fun_test(t,type,Q_tp_model,Qc_backward,Ti_tp,V_tp,...
        rhom,dTdP,Mc_backward,Mcp,H_BSE_tp_model,H_cc_tp_model);
    
    % change the results to be forward in time for later calculation
    Ti_model = flipud(Ti_backward);
    Q_model  = flipud(Q_backward);
    H_model  = flipud(H_backward);
    V_model  = flipud(V_backward);
    Z_model  = flipud(Z_backward);
    
    % save the results
    Ti(:,l) = Ti_model;
    Q(:,l) = Q_model;
    H(:,l) = H_model;
    V(:,l) = V_model;
    Z(:,l) = Z_model;
    
end % for l=1:size_Tp

%% Plot successful solutions of thermal evolution
% calculate the middle 50% and 90% of all successful results
[Ti_5,Ti_25,Ti_50,Ti_75,Ti_95] = calculate_percentile_fun(Ti,nt,t,size_Tp);
[H_5,H_25,H_50,H_75,H_95] = calculate_percentile_fun(H,nt,t,size_Tp);
[Q_5,Q_25,Q_50,Q_75,Q_95] = calculate_percentile_fun(Q,nt,t,size_Tp);
[V_5,V_25,V_50,V_75,V_95] = calculate_percentile_fun(V,nt,t,size_Tp);
[Z_5,Z_25,Z_50,Z_75,Z_95] = calculate_percentile_fun(Z,nt,t,size_Tp);

sz =90;% scatter size
figure(3);
subplot(2,2,1); hold off;
plot(t,Ti_25,'k--',t,Ti_75,'k--',t,Ti_5,'k:',t,Ti_95,'k:',t,Ti_50,'r-'); hold on;
scatter(t_Herz,Tp_Herz,sz,'MarkerEdgeColor',[0 .5 .5],'MarkerFaceColor','c',...
    'LineWidth',1.5);
xlabel('Time (Gyr)'); ylabel('Mantle Potential Temperature (degree C)');
grid on;

subplot(2,2,2); hold off;
plot(t,H_25,'k--',t,H_75,'k--',t,H_5,'k:',t,H_95,'k:',t,H_50,'r-'); hold on;
plot(t,Q_25,'k--',t,Q_75,'k--',t,Q_5,'k:',t,Q_95,'k:',t,Q_50,'b-'); grid on;
xlabel('Time (Gyr)');ylabel('H and Q (TW)');

subplot(2,2,3); hold off;
plot(t,V_25,'k--',t,V_75,'k--',t,V_5,'k:',t,V_95,'k:',t,V_50,'r-');
xlabel('Time (Gyr)');ylabel('Plate velocity'); grid on;

subplot(2,2,4); hold off;
plot(t,Z_25,'k--',t,Z_75,'k--',t,Z_5,'k:',t,Z_95,'k:',t,Z_50,'r-'); hold on;
xlabel('Time (Gyr)');ylabel('Melting depth'); grid on;
