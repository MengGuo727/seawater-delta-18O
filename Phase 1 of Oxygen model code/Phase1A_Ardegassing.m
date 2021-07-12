%% 2019/08/05 Phase 1A of argon degassing -- Crustal growth model
% This file uses Monte Carlo sampling to simulate different patterns of
% crustal growth and calculate their formation age and surface age distributions.
% The simulated parameters and misfits for two age distribution are stored in
% the file "out_phase1a.dat"
% 
% Meng Guo, Yale University
% Summer 2019


clear all;

% Open output file
fileID = fopen('out_phase1a.dat','w');

rng('shuffle');% to avoid have same "rand" results everytime

%% Set the time period
tmax = 4.567;% the age of solar system, in unit Ga
dt = 0.001;% length of each timestep
t = 0:dt:tmax;
nt = length(t);% number of timesteps
t = t';

%% Define the a priori ranges for independent variables
%  decay constant for crustal recycling rate, in unit Gyr-1
kappa_r_min = -3;
kappa_r_max = 3;
% decay constant for crustal growth rate, in unit Gyr-1
kappa_g_min = -1;
kappa_g_max = 30;
% the initial crustal recycling rate (t=ts), in unit kg/Gyr
Rs_min = 0;
Rs_max = 10e22;
% the present-day crustal recycling rate (t=tp), in unit kg/Gyr
Rp_min = 0;
Rp_max = 2e22;
% the onset time for crustal growth and recycle, in unit Ga
ts_min = tmax-4.51;
ts_max = tmax-4;
% the crustal reworking rate factor
Krw_factor_min = 0.1;
Krw_factor_max = 0.8;

%% Constants used in the model
Mcp = 2.09e22;% mass of continental crust at present-day, in unit kg

%% Decide the number of iteration
itermax = 2e5;

%% Load in the observed formation & surface age distributions
% Formation age distribution data from Korenaga (2018a)
data_formationage = load('korenaga18a_Tunmix_orig.dat');
% Surface age distribution data from Roberts & Spencer (2015)
data_zircon_surf = load('korenaga18a_T_U_Pb.dat');
% use load_FandS_fun function to set data in the same dimension as the time series
[F_Jun_same,S_Jun_same] = load_FandS_fun(t,nt,data_formationage,data_zircon_surf);

%% Using Monte Carlo sampling to calculate Mc (crustal growth pattern) and 
% formation age and surface age distributions
for l=1:itermax
    if mod(l,10) == 0
        disp(['l=' num2str(l)]);% keep track of the calculation
    end % if mod(l,10) == 0
    
    % Randomly sample the independent variables within the a priori distributions
    Krw_factor_model = Krw_factor_min + rand(1)*(Krw_factor_max - Krw_factor_min);
    kappa_r_model = kappa_r_min + rand(1)*(kappa_r_max - kappa_r_min);
    kappa_g_model = kappa_g_min + rand(1)*(kappa_g_max - kappa_g_min);
    Rs_model = Rs_min + rand(1)*(Rs_max - Rs_min);
    Rp_model = Rp_min + rand(1)*(Rp_max - Rp_min);
    ts_model = ts_min + rand(1)*(ts_max - ts_min);
    
    % Calculate the dependent variable (crustal reworking rate)
    Krw_s_model = Rs_model * Krw_factor_model;% initial Krw_factor
    
    % Calculate the corresponding crustal growth pattern
    [Mc_model,Mdd_model,Mud_model,Krw_model_first] = CC_growth_fun1(t,ts_model,tmax,...
        Mcp,kappa_g_model,Rp_model,Rs_model,kappa_r_model,Krw_s_model);
    
    % Calculate formation age and surface age distributions and recalculate 
    % crustal reworking rate using Formation_surface_age_fun function
    [F_model,S_model,m_tp,m,Krw_model] = ...
        Formation_surface_age_fun(t,Mud_model,Mdd_model,Mc_model,Krw_model_first);
    
    % Calculate RMSE for formation age and surface age distributions
    RMSE_F_model = sqrt((sum((F_model-F_Jun_same).^2))/(sum((F_Jun_same).^2)));
    RMSE_S_model = sqrt((sum((S_model-S_Jun_same).^2))/(sum((S_Jun_same).^2)));
    
    % Save all the variables and RMSEs
    fprintf(fileID,['%6g %6g %6g ' ...
        '%6g %6g %6g ' ...
        '%6g  %6g\n'], ...
        Krw_factor_model, kappa_r_model, kappa_g_model, ...
        Rs_model, Rp_model, ts_model,...
        RMSE_F_model, RMSE_S_model);
end % for l=1:itermax

fclose(fileID);
