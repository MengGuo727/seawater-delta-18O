%% 2020.03.17 Delta Oxygen Model A -- Monte Carlo simulation
% This file simulates the delta 18O evolution trend in seawater and
% collects the parameters of each Monte Carlo simulations.
close all;
clear all;
clc;

rng('shuffle');% to avoid have same random results everytime

% Open output file
fileID = fopen('out_oxygenmodel1_noFandPHI.dat','w');

%% load in the data for delta 18O 
% data = xlsread('Database_1.xlsx'); % Fe data from Noah science paper
% age = data(:,1); age = 3.5-age;
% deltaO= data(:,2);

data = xlsread('delta 18O of ocean in phanerozoic.xlsx'); % oxygen of seawater from past studies
age = data(:,1); age = 3.24-age;
deltaO= data(:,2);


%% Using crustal evolution of the Earth to constrain the low temperature processes during 18O exchange
% Read selected results from crustal growth models and thermal evolution models
% that satisfy formation and surface age distributions and mantle potential
% temperature
data = load('out_phase2b.dat');
% the variables in 'out_phase2a.dat' are in following order:
% 1.Krw_factor, 2.kappa_r, 3.kappa_g, 4.Rs, 5.Rp, ...
% 6.ts, 7.RMSE_F, 8.RMSE_S, 9.H_BSE_tp, 10.H_cc_tp,...
% 11.Q_total_tp, 12.Qc_tp, 13.d_Qc, 14.misfit_Tp

% Define the variables and misfits in file "out_phase2b.dat"
Krw_factor = data(:,1);
kappa_r = data(:,2);
kappa_g = data(:,3);
Rs = data(:,4);
Rp = data(:,5);
ts = data(:,6);
H_BSE_tp = data(:,9);
H_cc_tp = data(:,10);
Q_total_tp = data(:,11);
Qc_tp = data(:,12);
d_Qc = data(:,13);

% Set up constants to run the model
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

% set the constants to operate argon degassing model
Fpm_tp = 2.17e14;% the total plume flux at present-day, in unit kg/yr
Kmo_tp = 6.7e23;% present-day mantle processing rate to generate oceanic crust, in unit kg/Ga (Korenaga,2006)
Kpm_tp = 2.17e14*1e9;% the total plume flux at present-day, in unit kg/Gyr

%%  Set up constant model parameters, values are from Wallmann (2001) Talbe. 5; and the variable value ranges
% Set the time serie for calculating oxygen model
dt = 0.001;
tspan = 0.001:dt:3.24; 
tspan = tspan';

% The recylcing efficiency
% r_HT = nan(size(tspan)); % The recylcing efficiency for high T altered crust
% r_LT = nan(size(tspan)); % The recylcing efficiency for low T altered cru
% r_W = nan(size(tspan)); % The recylcing efficiency for sediments

r_HT_min = 0.40;
r_HT_max = 0.99;

% The isotope fractionation factors
alpha_w_min  = 1.017; % isotopic fractionation factor for H2O uptake in clays, 1.02
alpha_w_max  = 1.023; % isotopic fractionation factor for H2O uptake in clays, 1.02
alpha_LT_min = 1.012; % isotopic fractionation factor for H2O uptake in low T altered OC, 1.015
alpha_LT_max = 1.018; % isotopic fractionation factor for H2O uptake in low T altered OC, 1.015
alpha_rev_min = 1.022; % isotopic fractionation factor for H2O uptake during reverse weathering, 1.025
alpha_rev_max = 1.028; % isotopic fractionation factor for H2O uptake during reverse weathering
alpha_HT = 1;   % isotopic fractionation factor for H2O uptake in high T altered OC
alpha_m = 1; % isotopic fractionation factor for H2O uptake in mantle

% The kinetic constants for isotope exchange processes, in unit mole 18O/Gyr
k_w = 1.2e20;  % kinetic constant for isotopic exchange during weathering, 0.23e21
k_LT = 1.0e20;  % kinetic constant for isotopic exchange during low T alteration OC, 0.2e21
k_HT = 4.9e20;  % kinetic constant for isotopic exchange during high T alteration OC, was 0.3e21
k_rev = 1.2e20;  % kinetic constant for isotopic exchange during weathering, 0.23e21

% The 18O model fractions
phi_oc = 2.0126e-3;   % 18O mole fraction in fresh crust
phi_m = 2.015e-3;   % 18O mole fraction in mantle water
phi_si = nan(size(tspan));
for i = 1:2740
    phi_si(i) = 2.017e-3;   % 18O mole fraction in weathering silicate rocks, 2.017e-3
end
for i = 2741:3240
    phi_si(i) = 2.017e-3;   % 18O mole fraction in weathering silicate rocks, 2.032e-3;
end

% The masses of oxygen in different reservoirs, in unit mole
M_w_0 = 7e22;   % Mass of oxygen in sedimentary rocks
M_LT_0 = 1e22;   % Mass of oxygen in low T altered OC
M_HT_0 = 1.6e23;   % Mass of oxygen in high T altered OC

% Current water fixation, in unit mol/Gyr
F_w_0 = 7e21;   % Current H2O fixation in weathering products,7e21
F_LT_0 = 6e21;   % Current H2O fixation in low T altered OC, 6e21
F_HT_0_min = 20e21; % Current H2O fixation in high T altered OC, 20-70e21
F_HT_0_max = 70e21; % Current H2O fixation in high T altered OC, 20-70e21
F_m_0 = 3e21; % current H20 released from mantle, 3e21

% Current prodection of oxygen, in unit mole/Gyr
F_mw_0 = 7.8e22;   % Current production of sedimentary rock oxygen
F_mLT_0 = 1.2e23;   % Current production of low T altered OC oxygen
F_mHT_0 = 1.6e24;   % Current production of high T altered OC oxygen

% Current 18O loss, in unit mole/Gyr
F_pw18_0 = 2e17; % Current 18O loss during turnover of porewater

% The 18O/16O ratio in standard mean ocean water
R_smow = 0.0020052; % Hoefs, 1997

% The fraction between reverse weathering and the total LT processes
f_reverse = nan(size(tspan));
f_reverse_min = 0.01;
f_reverse_max = 0.6;

%  The corresponding constants calculated using the phi defined above
% 18O/16O ratios (constant)
R_oc = phi_oc / (1 - phi_oc);
R_si = phi_si ./ (1 - phi_si);
R_m = phi_m / (1 - phi_m);

delta_oc = (R_oc / R_smow - 1) * 1000; % Fresh crust isotope value
delta_si = (R_si / R_smow - 1) * 1000; % Fresh silicate isotope value
delta_m = (R_m / R_smow - 1) * 1000; % Fresh mantle water isotope value

% Set up initial conditions before using s_model
% free water mass and isotope (read from Fig.3 of Wallmann, 2001)
M_f_ini_min= 2.51e24/18; % in unit mole
M_f_ini_max= 3.035e24/18; % in unit mole
delta_f18_ini_min = 2;
delta_f18_ini_max = 4;

% 18O in sedimentary rocks (read from Fig.4 of Wallmann, 2001)
delta_w18_ini_min = 12;%15
delta_w18_ini_max = 18;%15

% 18O in low T altered OC (read from Fig.4 of Wallmann, 2001)
delta_LT18_ini_min = 7;%8
delta_LT18_ini_max = 11;%8

% 18O in high T altered OC (read from Fig.4 of Wallmann, 2001)
delta_HT18_ini_min = 3;
delta_HT18_ini_max = 7;

%% observed mass of ocean at present-day
M_f_observe = 1.46e24;

%% Decide the number of iteration
itermax = 50000;
                                                 
l = 1;
for l=1:itermax   
    if mod(l,10) == 0
        disp(['l=' num2str(l)]);% keep track of the calculation
    end % if mod(l,10) == 0

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
    % Set the time period for running crustal growth and thermal evolution functions
    tmax = 4.567;% the age of solar system, in unit Ga
    dt = 0.001;% length of each timestep
    t = 0:dt:tmax;
    nt = length(t);% number of timesteps
    t = t';
    
    % size of the data
    [ndata,nparam ]= size(Krw_factor);
    
    % Select one of the successful CC growth models and thermal evolutions
    r_index = randi(ndata);% randomly get an index
    Krw_factor_model = Krw_factor(r_index);
    kappa_r_model = kappa_r(r_index);
    kappa_g_model = kappa_g(r_index);
    Rs_model = Rs(r_index);
    Rp_model = Rp(r_index);
    ts_model = ts(r_index);
    H_BSE_tp_model = H_BSE_tp(r_index);
    H_cc_tp_model = H_cc_tp(r_index);
    Q_total_tp_model = Q_total_tp(r_index);
    Qc_tp_model = Qc_tp(r_index);
    d_Qc_model = d_Qc(r_index);
    
    % Calculate the dependent variables
    Q_tp_model = Q_total_tp_model - H_cc_tp_model;
    K40_BSE_tp_model = (K_factor * H_BSE_tp_model) / heat_factor * Tkg_atoms;
    K40_CC_observe_model = (K_factor * H_cc_tp_model) / heat_factor * Tkg_atoms;
    Krw_s_model = Rs_model * Krw_factor_model;

    % Calculate the corresponding crustal growth model
    [Mc_model,Mdd_model,Mud_model,Krw_model_first] = CC_growth_fun1(t,ts_model,tmax,...
        Mcp,kappa_g_model,Rp_model,Rs_model,kappa_r_model,Krw_s_model);
    
    % Calculate the corresponding formation age and surface age distributions,
    % and the crustal reworking rate
    [F_model,S_model,m_tp,m,Krw_model] =...
        Formation_surface_age_fun(t,Mud_model,Mdd_model,Mc_model,Krw_model_first);
    
    [Qc_model,Qc_backward] = Qc_backward_fun(Qc_tp_model,d_Qc_model,nt,t);
    
    % Calculate the corresponding thermal evolution
    % Since we are calculating H backward in time, Mc need to be backward in time as well
    Mc_backward = flipud(Mc_model);
    [Ti_backward, Q_backward, H_backward,V_backward, Z_backward] = ...
        Thermal_history_fun_test(t,type,Q_tp_model,Qc_backward,Ti_tp,V_tp,...
        rhom,dTdP,Mc_backward,Mcp,H_BSE_tp_model,H_cc_tp_model);
    
    % Change the results to be forward in time for later calculation
    Ti_model = flipud(Ti_backward);
    Q_model  = flipud(Q_backward);
    H_model  = flipud(H_backward);
    V_model  = flipud(V_backward);
    Z_model  = flipud(Z_backward);
    
    % Calculate the generation rates of oceanic crust and hotspot islands,
    % which control the rate of high temperature
    % processes and partially contribute to the low temperature processes
    % Use decreasing Q to calcualte spreading rate
    % set the constants to operate argon degassing model
    Fpm_tp = 2.17e14;% the total plume flux at present-day, in unit kg/yr
    Kmo_tp = 6.7e23;% present-day mantle processing rate to generate oceanic crust, in unit kg/Ga (Korenaga,2006)
    Kpm_tp = 2.17e14*1e9;% the total plume flux at present-day, in unit kg/Gyr
    
    % Calculate mantle processing rate at mantle plume
    Kmp_model = nan(nt,1);
    for i = 1:nt
        Kmp_model(i) = Kpm_tp * (Qc_model(i)/Qc_tp_model);
    end
    
    % Calculate generation rate of oceanic crust (Kmo), continental crust (Kmc=Mud),
    % generate rate of CC from mantle in addition to secondary metling from OC
    % (Kmc_c), crustal recycling rate (Krc=Mdd)
    
    [Kmo_model,Kmc_model,Krc_model, Kmc_c_model] = ...
        Kmo_fun(t,Z_model,V_model,Kmo_tp,Mc_model,Mdd_model,Mud_model,Mcp,Mmp,...
        K40_BSE_tp_model,K40_CC_observe_model);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
    % Use the mass transfer rates from thermal and crustal evolution model
    K_highT= V_model(1329:4568,1);
    K_lowT_Krc = Krc_model(1329:4568,1);
    K_lowT_Krw = Krw_model(1329:4568,1);
    K_lowT_Kmo = Kmo_model(1329:4568,1);
    K_lowT_Kmc = Kmc_model(1329:4568,1);% this is Kmc + Kmo
    K_lowT_Kmp = Kmp_model(1329:4568,1);
    K_lowT = K_lowT_Kmc + 100*K_lowT_Krc + 100*K_lowT_Krw + K_lowT_Kmp;
    
    f_V = K_highT./K_highT(3240);
    f_Krc = K_lowT_Krc./K_lowT_Krc(3240);
    f_Kmc = K_lowT_Kmc./K_lowT_Kmc(3240);
    f_Kmp = Kmp_model./Kmp_model(3240);
    
    f_LT = K_lowT./K_lowT(3240);
    n_range = [1,2,3];
    n = n_range(randperm(length(n_range),1));
    f_HT = (f_V).^n;
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%     
    % Determine the variable values for oxygen model
  
    % The recylcing efficiency for high T altered crust
    r_factor = r_HT_min + rand(1)*(r_HT_max - r_HT_min);
    r_HT(1:3240) = r_factor;
    r_LT = r_HT;% The recylcing efficiency for low T altered crust
    r_w = r_HT;% The recylcing efficiency for sediments

    % The isotope fractionation factors
    alpha_w  = alpha_w_min + rand(1)*(alpha_w_max - alpha_w_min); % isotopic fractionation factor for H2O uptake in clays, 1.02
    alpha_LT = alpha_LT_min + rand(1)*(alpha_LT_max - alpha_LT_min); % isotopic fractionation factor for H2O uptake in clays, 1.02
    alpha_rev = alpha_rev_min + rand(1)*(alpha_rev_max - alpha_rev_min); % isotopic fractionation factor for H2O uptake in clays, 1.02
    
    % Current water fixation, in unit mol/Gyr 
    F_HT_0 = F_HT_0_min + rand(1)*(F_HT_0_max - F_HT_0_min);
    
    % reverse weathering factore
    f_reverse_factore = f_reverse_min + rand(1)*(f_reverse_max - f_reverse_min);
    for i = 1:2740
        f_reverse(i) = f_reverse_factore;%0.1
    end
    for i = 2741:3240
        f_reverse(i) = f_reverse_factore;% 0 
    end
    
    % Set up initial conditions before running s_model function 
    % free water mass and isotope (read from Fig.3 of Wallmann, 2001)
    M_f_ini= M_f_ini_min + rand(1)*(M_f_ini_max - M_f_ini_min); % in unit mole
    delta_f18_ini = delta_f18_ini_min + rand(1)*(delta_f18_ini_max - delta_f18_ini_min); 
    R_f18_ini = R_smow * (delta_f18_ini / 1000 + 1);
    M_f18_ini = M_f_ini * (R_f18_ini / (1 + R_f18_ini)); % in unit mole
    
    % 18O in sedimentary rocks (read from Fig.4 of Wallmann, 2001)
    delta_w18_ini = delta_w18_ini_min + rand(1)*(delta_w18_ini_max - delta_w18_ini_min);
    R_w18_ini = R_smow * (delta_w18_ini / 1000 + 1);
    M_w18_ini = M_w_0 * (R_w18_ini / (1 + R_w18_ini)); % in unit mole
    phi_w_ini = R_w18_ini / (1 + R_w18_ini);
    M_w_ini = M_w18_ini / phi_w_ini;
    
    % 18O in low T altered OC (read from Fig.4 of Wallmann, 2001)
    delta_LT18_ini = delta_LT18_ini_min + rand(1)*(delta_LT18_ini_max - delta_LT18_ini_min);
    R_LT18_ini = R_smow * (delta_LT18_ini / 1000 + 1);
    M_LT18_ini = M_LT_0 * (R_LT18_ini / (1 + R_LT18_ini)); % in unit mole
    phi_LT_ini = R_LT18_ini / (1 + R_LT18_ini);
    M_LT_ini = M_LT18_ini / phi_LT_ini ;
    
    % 18O in high T altered OC (read from Fig.4 of Wallmann, 2001)
    delta_HT18_ini = delta_HT18_ini_min + rand(1)*(delta_HT18_ini_max - delta_HT18_ini_min);
    R_HT18_ini = R_smow * (delta_HT18_ini / 1000 + 1);
    M_HT18_ini = M_HT_0 * (R_HT18_ini / (1 + R_HT18_ini)); % in unit mole
    phi_HT_ini = R_HT18_ini / (1 + R_HT18_ini);
    M_HT_ini = M_HT18_ini / phi_HT_ini ;
    
    y_ini  = [M_f_ini, M_f18_ini, M_w_ini, M_w18_ini, ...
        M_LT_ini, M_LT18_ini, M_HT_ini, M_HT18_ini];
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Set the time serie for calculating oxygen model
    dt = 0.001;
    tspan = 0.001:dt:3.24; 
    tspan = tspan';
    
    % Calculate oxygen circulation history
    [M_f, M_f18, M_w, M_w18, M_LT, M_LT18, M_HT, M_HT18,...
        F_w, F_LT, F_HT, F_rew, F_reLT, F_reHT, ...
        F_w18, F_LT18, F_HT18, F_w_ex, ...
        F_LT_ex, F_HT_ex, F_pw18, phi_w, phi_LT, phi_HT, F_rew18, ...
        F_reLT18, F_reHT18, F_m18, F_mw18, F_mLT18, F_mHT18, F_sw18, ...
        F_sLT18, F_sHT18,phi_f, R_f, omega_w, omega_LT , omega_HT,...
        F_reverse, F_reverse18, F_reverse_ex, F_m] = ...
        s_model_withCC2(y_ini, tspan, f_LT, f_HT, f_V, f_Krc, f_Kmc, ...
        f_Kmp,F_w_0, F_LT_0, F_HT_0, r_w,...
        r_LT, r_HT, alpha_w, alpha_LT, alpha_HT, alpha_rev, R_si, k_w, ...
        R_oc, k_LT, k_HT,k_rev,F_pw18_0, phi_m, phi_si, F_mw_0, phi_oc, ...
        F_mLT_0, F_mHT_0, dt,F_m_0, f_reverse);
 
    % First, calculate misfit of present-day's 36Ar in the atmosphere
    misfit_Mf = ((M_f(3240)*18 - M_f_observe)/1.460e+23)^2;

    if (misfit_Mf > 1)
        % adjust recycling rate to better match with M_f_observe
        r_factor = r_factor * (M_f_observe/M_f(3240)/18);
        r_HT(1:3240) = r_factor;
        r_LT = r_HT;% The recylcing efficiency for low T altered crust
        r_w = r_HT;% The recylcing efficiency for sediments
        
        % recalculate the argon degassing history using new initial 36Ar
        [M_f, M_f18, M_w, M_w18, M_LT, M_LT18, M_HT, M_HT18,...
            F_w, F_LT, F_HT, F_rew, F_reLT, F_reHT, ...
            F_w18, F_LT18, F_HT18, F_w_ex, ...
            F_LT_ex, F_HT_ex, F_pw18, phi_w, phi_LT, phi_HT, F_rew18, ...
            F_reLT18, F_reHT18, F_m18, F_mw18, F_mLT18, F_mHT18, F_sw18, ...
            F_sLT18, F_sHT18,phi_f, R_f, omega_w, omega_LT , omega_HT,...
            F_reverse, F_reverse18, F_reverse_ex, F_m] = ...
            s_model_withCC2(y_ini, tspan, f_LT, f_HT, f_V, f_Krc, f_Kmc, ...
            f_Kmp,F_w_0, F_LT_0, F_HT_0, r_w,...
            r_LT, r_HT, alpha_w, alpha_LT, alpha_HT, alpha_rev, R_si, k_w, ...
            R_oc, k_LT, k_HT,k_rev,F_pw18_0, phi_m, phi_si, F_mw_0, phi_oc, ...
            F_mLT_0, F_mHT_0, dt,F_m_0, f_reverse);
        
        misfit_Mf = ((M_f(3240)*18 - M_f_observe)/1.460e+23 )^2;
    end
    
    % Calcualte the delta 18O of each reservoir
    delta_f18_all = phi_to_delta(M_f18 , M_f, R_smow,1);
    delta_w18_all = phi_to_delta(M_w18 , M_w, R_smow,1);
    delta_LT18_all = phi_to_delta(M_LT18 , M_LT, R_smow,1);
    delta_HT18_all = phi_to_delta(M_HT18 , M_HT, R_smow,1);
    
    % Calculate misfits
    misfit_deltaO_lateArchean = ((delta_f18_all(1850)-(-8))/1)^2;
    misfit_deltaO_phanerozoic1 = ((delta_f18_all(2740)-(-3))/1)^2;
    misfit_deltaO_phanerozoic2 = ((delta_f18_all(3240)-(0))/1)^2;
    misfit_deltaO_phanerozoic3 = ((delta_f18_all(2500)-(-6))/1)^2;
    
    % Save all the variables and misfits
    fprintf(fileID,['%6g %6g %6g %6g %6g ' ...
        '%6g %6g %6g %6g %6g '  ...
        '%6g %6g '  ...
        '%6g %6g %6g %6g %6g '  ...
        '%6g %6g %6g %6g '  ...
        '%6g %6g %6g %6g %6g %6g %6g\n'], ...
        Krw_factor_model, kappa_r_model, kappa_g_model, Rs_model, Rp_model, ...
        ts_model,  H_BSE_tp_model, H_cc_tp_model,Q_total_tp_model, Qc_tp_model,...
        d_Qc_model,...
        f_reverse_factore,r_factor, alpha_w, alpha_LT, alpha_rev,...
        F_HT_0, M_f_ini, delta_w18_ini, delta_LT18_ini,...
        delta_HT18_ini, delta_f18_ini, misfit_Mf, misfit_deltaO_lateArchean, misfit_deltaO_phanerozoic1, misfit_deltaO_phanerozoic2, misfit_deltaO_phanerozoic3,n);
    
    
end % for l=1:itermax;

fclose(fileID);

