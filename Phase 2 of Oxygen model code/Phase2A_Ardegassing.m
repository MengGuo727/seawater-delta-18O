%% 2019/08/05 Phase 2A of argon degassing -- Thermal evolution
% This file uses the successful solutions "out_phase1b.dat" from crustal
% growth and uses Monte Carlo sampling to test thermal history. The
% results are stored in file "out_phase2a.dat".
%
% Meng Guo, Yale University
% Summer, 2019

clear all;

rng('shuffle');% to avoid have same random results everytime

% Open output file
fileID = fopen('out_phase2a.dat','w');

%% Read selected results from stage 1 (out_phase1b.dat)
data = load('out_phase1b.dat');
% the variables in 'out_phase1a.dat' are in following order:
% 1.Krw_factor, 2.kappa_r, 3.kappa_g, 4.Rs, 5.Rp, 6.ts, 7.RMSE_F, 8.RMSE_S

%% Define the variables and misfits in file "out_phase1b.dat"
Krw_factor = data(:,1);
kappa_r = data(:,2);
kappa_g = data(:,3);
Rs = data(:,4);
Rp = data(:,5);
ts = data(:,6);
RMSE_F = data(:,7);
RMSE_S = data(:,8);

% size of the data
[ndata,nparam ]= size(Krw_factor);

%% Set the time period
tmax = 4.567;% the age of solar system, in unit Ga
dt = 0.001;% length of each timestep
t = 0:dt:tmax;
nt = length(t);% number of timesteps
t = t';

%% Define the a priori ranges for independent variables
% heat production in the BSE, in unit TW
H_BSE_tp_min = 13;
H_BSE_tp_max = 19;
% heat production in the continental crust, in unit TW
H_cc_tp_min = 5;
H_cc_tp_max = 10;
% total heat flux, in unit TW
Q_total_tp_min = 43;
Q_total_tp_max = 49;
% core heat flux, in unit TW (Lay et al.,2008)
Qc_tp_max = 15;
Qc_tp_min = 5;
% difference between the initial and the present-day Qc,
% in unit TW, O'Rourke et al. (2017)
d_Qc_max = 5;
d_Qc_min = 2;

%% Constants used to run thermal history and crustal growth functions
% constants for performing crustal growth model
Mcp = 2.09e22;% mass of continental crust at present-day, in unit kg

% constants for performing thermal evolution model
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

%% Herzburg et all., 2010 data for mantle potential temperature
% load the Ti from Herzberg et al., 2010
data_Tp = xlsread('Herz data.xlsx');
% decide anchor points for misfit calculation
[Tp_anchorHerz1,Tp_anchorHerz2,Tp_anchorHerz3,Tp_anchorHerz4,...
    t_anchorHerz1,t_anchorHerz2,t_anchorHerz3,t_anchorHerz4,...
    t_Herz,Tp_Herz] = load_Tp_fun(data_Tp);

%% Decide the number of iteration
itermax = 3e4;

%% Using Monte Carlo sampling to test different variable combinations
for l=1:itermax
    if mod(l,10) == 0
        disp(['l=' num2str(l)]);% keep track of the calculation
    end % if mod(l,10) == 0
    
    % Randomly sample the independent variables within the a priori distributions
    H_BSE_tp_model = H_BSE_tp_min + rand(1)*(H_BSE_tp_max - H_BSE_tp_min);
    H_cc_tp_model = H_cc_tp_min + rand(1)*(H_cc_tp_max - H_cc_tp_min);
    Q_total_tp_model = Q_total_tp_min + rand(1)*(Q_total_tp_max - Q_total_tp_min);
    Qc_tp_model = Qc_tp_min + rand(1)*(Qc_tp_max - Qc_tp_min);
    d_Qc_model = d_Qc_min + rand(1)*(d_Qc_max - d_Qc_min);
    
    % select one successful CC growth model from the results of stage 1
    r_index = randi(ndata);
    Krw_factor_model = Krw_factor(r_index);
    kappa_r_model = kappa_r(r_index);
    kappa_g_model = kappa_g(r_index);
    Rs_model = Rs(r_index);
    Rp_model = Rp(r_index);
    ts_model = ts(r_index);
    RMSE_F_model = RMSE_F(r_index);
    RMSE_S_model = RMSE_S(r_index);
    
    % Calculate the dependent variables
    Q_tp_model = Q_total_tp_model - H_cc_tp_model;
    K40_BSE_tp_model = (K_factor * H_BSE_tp_model) / heat_factor * Tkg_atoms;
    K40_CC_observe_model = (K_factor * H_cc_tp_model) / heat_factor * Tkg_atoms;
    Krw_s_model = Rs_model * Krw_factor_model;
    
    % Calculate the corresponding crustal growth pattern
    [Mc_model,Mdd_model,Mud_model,Krw_model_first] = CC_growth_fun1(t,ts_model,tmax,...
        Mcp,kappa_g_model,Rp_model,Rs_model,kappa_r_model,Krw_s_model);
    
    [Qc_model,Qc_backward] = Qc_backward_fun(Qc_tp_model,d_Qc_model,nt,t);
    
    % Calculate Thermal evolution
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
    
    % Calculate misfit for mantle potential temperature using 4 anchor points
    misfit_Tp_model = sqrt((((Ti_model(3768)-Tp_anchorHerz1)/70)^2 ...
        +((Ti_model(2698)-Tp_anchorHerz2)/90)^2 ...
        +((Ti_model(1818)-Tp_anchorHerz3)/100)^2 ...
        +((Ti_model(1178)-Tp_anchorHerz4)/130)^2)/4);
    
    % Save all the variables and misfit
    fprintf(fileID,['%6g %6g %6g %6g %6g '...
        '%6g %6g %6g %6g %6g ' ...
        '%6g %6g %6g %6g\n'], ...
        Krw_factor_model, kappa_r_model, kappa_g_model, Rs_model, Rp_model, ...
        ts_model, RMSE_F_model, RMSE_S_model, H_BSE_tp_model, H_cc_tp_model,...
        Q_total_tp_model, Qc_tp_model, d_Qc_model, misfit_Tp_model);
    
end % for l=1:itermax;

fclose(fileID);


