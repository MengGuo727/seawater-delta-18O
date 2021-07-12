%% 2020.03.17 Delta Oxygen Model C -- plot simulation results
% This file reads the output from Oxygen_model_MonteCarloB and plots the 
% results. 
close all;
clear all;
clc;

rng('shuffle');% to avoid have same random results everytime

%% Read file "out_phase2.dat"
success_oxygen = load('out_oxygenmodel3_scenario1.dat');
% the variables in 'out_phase3a.dat' are in following order:
% 1.Krw_factor, 2.kappa_r, 3.kappa_g, 4.Rs, 5.Rp, ...
% 6.ts, 7.H_BSE_tp, 8.H_cc_tp,...
% 9.Q_total_tp, 10.Qc_tp, 11.d_Qc, ...
% 12.f_reverse_factore, 13.r_factor, 14.alpha_w, 15.alpha_LT, 16.alpha_rev,...
% 17.F_HT_0, 18.M_f_ini, 19.delta_w18_ini, 20.delta_LT18_ini, 21.delta_HT18_ini,...
% 22.delta_f18_ini, 23,misfit_Mf, 24.misfit_deltaO_lateArchean,...
% 25.misfit_deltaO_phanerozoic1, 26.misfit_deltaO_phanerozoic2, 27.misfit_deltaO_phanerozoic3
% 28.n

%% Set up constants to run the model
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

% set constants to run oxygen model
alpha_HT = 1;   % isotopic fractionation factor for H2O uptake in high T altered OC
alpha_m = 1; % isotopic fractionation factor for H2O uptake in mantle
% The kinetic constants for isotope exchange processes, in unit mole 18O/Gyr
k_w = 1.2e20;  % kinetic constant for isotopic exchange during weathering, 0.23e21
k_LT = 1.0e20;  % kinetic constant for isotopic exchange during low T alteration OC, 0.2e21
k_HT = 4.9e20;  % kinetic constant for isotopic exchange during high T alteration OC, was 0.3e21
k_rev = 1.2e20;  % kinetic constant for isotopic exchange during weathering, 0.23e21
% The 18O model fractions
for i = 1:2740
    phi_si(i) = 2.017e-3;   % 18O mole fraction in weathering silicate rocks, 2.017e-3
end
for i = 2741:3240
    phi_si(i) = 2.017e-3;   % 18O mole fraction in weathering silicate rocks, 2.032e-3;
end

phi_oc = 2.0126e-3;   % 18O mole fraction in fresh crust
phi_m = 2.015e-3;   % 18O mole fraction in mantle water
% The masses of oxygen in different reservoirs, in unit mole
M_w_0 = 7e22;   % Mass of oxygen in sedimentary rocks
M_LT_0 = 1e22;   % Mass of oxygen in low T altered OC
M_HT_0 = 1.6e23;   % Mass of oxygen in high T altered OC
% Current water fixation, in unit mol/Gyr
F_w_0 = 7e21;   % Current H2O fixation in weathering products,7e21
F_LT_0 = 6e21;   % Current H2O fixation in low T altered OC, 6e21
F_m_0 = 3e21; % current H20 released from mantle, 3e21
% Current prodection of oxygen, in unit mole/Gyr
F_mw_0 = 7.8e22;   % Current production of sedimentary rock oxygen
F_mLT_0 = 1.2e23;   % Current production of low T altered OC oxygen
F_mHT_0 = 1.6e24;   % Current production of high T altered OC oxygen
% Current 18O loss, in unit mole/Gyr
F_pw18_0 = 2e17; % Current 18O loss during turnover of porewater
% The 18O/16O ratio in standard mean ocean water
R_smow = 0.0020052; % Hoefs, 1997
%  The corresponding constants calculated using the phi defined above
% 18O/16O ratios (constant)
R_oc = phi_oc / (1 - phi_oc);
R_si = phi_si ./ (1 - phi_si);
R_m = phi_m / (1 - phi_m);
delta_oc = (R_oc / R_smow - 1) * 1000; % Fresh crust isotope value
delta_si = (R_si / R_smow - 1) * 1000; % Fresh silicate isotope value
delta_m = (R_m / R_smow - 1) * 1000; % Fresh mantle water isotope value
% observed mass of ocean at present-day
M_f_observe = 1.46e24;

%% set up empty matrixes to store results
nt1 = 3240;
nt2 = 4568;
[size_oxygen, nparameter ]= size(success_oxygen);
M_f_success= nan(nt1,size_oxygen);
M_f_g_success= nan(nt1,size_oxygen);
delta_f18_all_success = nan(nt1,size_oxygen);
delta_w18_all_success = nan(nt1,size_oxygen);
delta_LT18_all_success = nan(nt1,size_oxygen);
delta_HT18_all_success = nan(nt1,size_oxygen);
f_LT_success = nan(nt1,size_oxygen);
f_HT_success = nan(nt1,size_oxygen);
f_Krc_success = nan(nt1,size_oxygen);

Kmc_c_success = nan(nt2,size_oxygen);
Mc_success = nan(nt2,size_oxygen);
Krc_success = nan(nt2,size_oxygen);
Kmc_success = nan(nt2,size_oxygen);
Krw_success = nan(nt2,size_oxygen);
Ti_success = nan(nt2,size_oxygen);
Q_success = nan(nt2,size_oxygen);
H_success = nan(nt2,size_oxygen);
V_success = nan(nt2,size_oxygen);
Z_success = nan(nt2,size_oxygen);
F_success = nan(nt2,size_oxygen);
S_success = nan(nt2,size_oxygen);

%% Calculate oxgen transfer history
for l=1:size_oxygen
    if mod(l,10)==0
        disp(['l=' num2str(l) ' of ' num2str(size_oxygen)]);% keep track of the calculation
    end
    % the independent variables
    Krw_factor_model = success_oxygen(l,1);
    kappa_r_model = success_oxygen(l,2);
    kappa_g_model = success_oxygen(l,3);
    Rs_model = success_oxygen(l,4);
    Rp_model = success_oxygen(l,5);
    ts_model = success_oxygen(l,6);
    H_BSE_tp_model = success_oxygen(l,7);
    H_cc_tp_model = success_oxygen(l,8);
    Q_total_tp_model = success_oxygen(l,9);
    Qc_tp_model = success_oxygen(l,10);
    d_Qc_model = success_oxygen(l,11);

    f_reverse_factore= success_oxygen(l,12);
    r_factor = success_oxygen(l,13);
    alpha_w = success_oxygen(l,14);
    alpha_LT = success_oxygen(l,15);
    alpha_rev = success_oxygen(l,16);
    F_HT_0 = success_oxygen(l,17);
    M_f_ini = success_oxygen(l,18);
    delta_w18_ini = success_oxygen(l,19);
    delta_LT18_ini = success_oxygen(l,20);
    delta_HT18_ini = success_oxygen(l,21);
    delta_f18_ini = success_oxygen(l,22);
    misfit_Mf = success_oxygen(l,23);
    misfit_deltaO_lateArchean = success_oxygen(l,24);
    misfit_deltaO_phanerozoic = success_oxygen(l,25);
    misfit_deltaO_phanerozoic2  = success_oxygen(l,26);
    n = success_oxygen(l,28);

    % calculate the dependent variables
    Krw_s_model = Rs_model * Krw_factor_model;% initial Krw_factor
    Q_tp_model = Q_total_tp_model - H_cc_tp_model;
    K40_BSE_tp_model = (K_factor * H_BSE_tp_model) / heat_factor * Tkg_atoms;
    K40_CC_observe_model = (K_factor * H_cc_tp_model) / heat_factor * Tkg_atoms;
    
    % Set the time period for running crustal growth and thermal evolution functions
    tmax = 4.567;% the age of solar system, in unit Ga
    dt = 0.001;% length of each timestep
    t = 0:dt:tmax;
    nt = length(t);% number of timesteps
    t = t';
    
    % Calculate the corresponding crustal growth model
    [Mc_model,Mdd_model,Mud_model,Krw_model_first] = CC_growth_fun1(t,ts_model,tmax,...
        Mcp,kappa_g_model,Rp_model,Rs_model,kappa_r_model,Krw_s_model);
    
    % Calculate the corresponding formation age and surface age distribution...
    % and the crustal reworking rate
    [F_model,S_model,m_tp,m,Krw_model] = ...
        Formation_surface_age_fun(t,Mud_model,Mdd_model,Mc_model,Krw_model_first);
    
    [Qc_model,Qc_backward] = Qc_backward_fun(Qc_tp_model,d_Qc_model,nt,t);
    
    % Calculate the corresponding thermal evolution
    % Since we are calcualting H backward in time, Mc need to be backward in time as well
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
    
    % Use the mass transfer rates from thermal and crustal evolution model
    K_highT= V_model(1329:4568,1); % ask V or Q
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
    f_HT = (f_V).^n;
    
    % Calculate the corresponding formation age and surface age distribution...
    % and the crustal reworking rate
    [F_model,S_model,m_tp,m,Krw_model] = ...
        Formation_surface_age_fun(t,Mud_model,Mdd_model,Mc_model,Krw_model_first);
    
    
    % recycling rates
    r_HT(1:3240) = r_factor;
    r_LT = r_HT;% The recylcing efficiency for low T altered crust
    r_w = r_HT;% The recylcing efficiency for sediments
    
    % reverse weathering factore
    for i = 1:2740
        f_reverse(i) = f_reverse_factore;
    end
    for i = 2741:3240
        f_reverse(i) = f_reverse_factore;% 0
    end
    
    % Set up initial conditions before running s_model function
    % free water mass and isotope (read from Fig.3 of Wallmann, 2001)
    R_f18_ini = R_smow * (delta_f18_ini / 1000 + 1);
    M_f18_ini = M_f_ini * (R_f18_ini / (1 + R_f18_ini)); % in unit mole
    
    % 18O in sedimentary rocks (read from Fig.4 of Wallmann, 2001)
    R_w18_ini = R_smow * (delta_w18_ini / 1000 + 1);
    M_w18_ini = M_w_0 * (R_w18_ini / (1 + R_w18_ini)); % in unit mole
    phi_w_ini = R_w18_ini / (1 + R_w18_ini);
    M_w_ini = M_w18_ini / phi_w_ini;
    
    % 18O in low T altered OC (read from Fig.4 of Wallmann, 2001)
    R_LT18_ini = R_smow * (delta_LT18_ini / 1000 + 1);
    M_LT18_ini = M_LT_0 * (R_LT18_ini / (1 + R_LT18_ini)); % in unit mole
    phi_LT_ini = R_LT18_ini / (1 + R_LT18_ini);
    M_LT_ini = M_LT18_ini / phi_LT_ini ;
    
    % 18O in high T altered OC (read from Fig.4 of Wallmann, 2001)
    R_HT18_ini = R_smow * (delta_HT18_ini / 1000 + 1);
    M_HT18_ini = M_HT_0 * (R_HT18_ini / (1 + R_HT18_ini)); % in unit mole
    phi_HT_ini = R_HT18_ini / (1 + R_HT18_ini);
    M_HT_ini = M_HT18_ini / phi_HT_ini ;
    
    y_ini  = [M_f_ini, M_f18_ini, M_w_ini, M_w18_ini, ...
        M_LT_ini, M_LT18_ini, M_HT_ini, M_HT18_ini];
    
    % Set the time serie for calculating oxygen model
    dt = 0.001;
    tspan = 0.001:dt:3.24; % the time serie to reproduce Wallmann 2001
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

    % Calcualte the delta 18O of each reservoir
    delta_f18_all = phi_to_delta(M_f18 , M_f, R_smow,1);
    delta_w18_all = phi_to_delta(M_w18 , M_w, R_smow,1);
    delta_LT18_all = phi_to_delta(M_LT18 , M_LT, R_smow,1);
    delta_HT18_all = phi_to_delta(M_HT18 , M_HT, R_smow,1);
        
    % save the results
    M_f_success(:,l) = M_f;
    M_f_g_success(:,l) = M_f * 18;
    delta_f18_all_success(:,l) = delta_f18_all;
    delta_w18_all_success(:,l) = delta_w18_all;
    delta_LT18_all_success(:,l) = delta_LT18_all;
    delta_HT18_all_success(:,l) = delta_HT18_all;
    f_LT_success(:,l) = f_LT;
    f_HT_success(:,l) = f_HT;
    f_Krc_success(:,l) = f_Krc;
    
    Kmc_c_success(:,l) = Kmc_c_model;
    Mc_success(:,l) = Mc_model;
    Krc_success(:,l) = Krc_model;
    Kmc_success(:,l) = Kmc_model;
    Krw_success(:,l) = Krw_model;
    Ti_success(:,l) = Ti_model;
    Q_success(:,l) = Q_model;
    H_success(:,l) = H_model;
    V_success(:,l) = V_model;
    Z_success(:,l) = Z_model;
    F_success(:,l) = F_model;
    S_success(:,l) = S_model;
    
end % for l=1:itermax;


%% plot net crustal growth(Mc), crustal generation rate(Mud), ...
%  and crustal recycling rate (Mdd)
figure(1);
% Set the time period for running crustal growth and thermal evolution functions
tmax = 4.567;% the age of solar system, in unit Ga
dt = 0.001;% length of each timestep
t = 0:dt:tmax;
nt = length(t);% number of timesteps
t = t';
[Mc_5,Mc_25,Mc_50,Mc_75,Mc_95] = calculate_percentile_fun(Mc_success,nt,t,size_oxygen);
[Kmc_c_5,Kmc_c_25,Kmc_c_50,Kmc_c_75,Kmc_c_95] = calculate_percentile_fun(Kmc_c_success,nt,t,size_oxygen);
[Krc_5,Krc_25,Krc_50,Krc_75,Krc_95] = calculate_percentile_fun(Krc_success,nt,t,size_oxygen);
[Krw_5,Krw_25,Krw_50,Krw_75,Krw_95] = calculate_percentile_fun(Krc_success,nt,t,size_oxygen);

subplot(2,2,1); hold off;
plot(t,Mc_25,'k',t,Mc_75,'k',t,Mc_5,'k',t,Mc_95,'k','Linewidth',1);hold on
xlabel('Time (Gyr)'); ylabel('Net crustal growth'); 
xlim([0,4.567]);
fill([t' fliplr(t')],[Mc_25 fliplr(Mc_75)],[0.4660 0.540 0.1880]);
fill([t' fliplr(t')],[Mc_5 fliplr(Mc_25)],'g');
fill([t' fliplr(t')],[Mc_75 fliplr(Mc_95)],'g');

subplot(2,2,2); hold off;
plot(t,Kmc_c_25,'k',t,Kmc_c_75,'k',t,Kmc_c_5,'k',t,Kmc_c_95,'k','Linewidth',1);hold on
xlabel('Time (Gyr)'); ylabel('Crustal generation rate');xlim([0,4.567]);
fill([t' fliplr(t')],[Kmc_c_25 fliplr(Kmc_c_75)],[0.4660 0.540 0.1880]);
fill([t' fliplr(t')],[Kmc_c_5 fliplr(Kmc_c_25)],'g');
fill([t' fliplr(t')],[Kmc_c_75 fliplr(Kmc_c_95)],'g');

subplot(2,2,3); hold off;
plot(t,Krc_25,'k',t,Krc_75,'k',t,Krc_5,'k',t,Krc_95,'k','Linewidth',1);hold on
xlabel('Time (Gyr)'); ylabel('Crustal recycling rate'); xlim([0,4.567]);
fill([t' fliplr(t')],[Krc_25 fliplr(Krc_75)],[0.4660 0.540 0.1880]);
fill([t' fliplr(t')],[Krc_5 fliplr(Krc_25)],'g');
fill([t' fliplr(t')],[Krc_75 fliplr(Krc_95)],'g');

subplot(2,2,4); hold off;
plot(t,Krw_25,'k',t,Krw_75,'k',t,Krw_5,'k',t,Krw_95,'k','Linewidth',1);hold on
xlabel('Time (Gyr)'); ylabel('Crustal reworking rate'); xlim([0,4.567]);
fill([t' fliplr(t')],[Krw_25 fliplr(Krw_75)],[0.4660 0.540 0.1880]);
fill([t' fliplr(t')],[Krw_5 fliplr(Krw_25)],'g');
fill([t' fliplr(t')],[Krw_75 fliplr(Krw_95)],'g');

myfig = figure(1);
myfig.Renderer = 'Painters';

%% plot thermal evolution history
data_Tp = xlsread('Herz data.xlsx');
% set the anchor points for misfit calculation
[Tp_anchorHerz1,Tp_anchorHerz2,Tp_anchorHerz3,Tp_anchorHerz4,...
    t_anchorHerz1,t_anchorHerz2,t_anchorHerz3,t_anchorHerz4,...
    t_Herz,Tp_Herz] = load_Tp_fun(data_Tp);

sz =90;% scatter size

[Ti_5,Ti_25,Ti_50,Ti_75,Ti_95] = calculate_percentile_fun(Ti_success,nt,t,size_oxygen);
[Q_5,Q_25,Q_50,Q_75,Q_95] = calculate_percentile_fun(Q_success,nt,t,size_oxygen);
[H_5,H_25,H_50,H_75,H_95] = calculate_percentile_fun(H_success,nt,t,size_oxygen);
[V_5,V_25,V_50,V_75,V_95] = calculate_percentile_fun(V_success,nt,t,size_oxygen);
[Z_5,Z_25,Z_50,Z_75,Z_95] = calculate_percentile_fun(Z_success,nt,t,size_oxygen);

figure(2);
subplot(2,2,1); hold off;
% plot the mantle potential temperature
plot(t,Ti_5,'g-','LineWidth',0.5);hold on;
plot(t,Ti_95,'g-','LineWidth',0.5);
plot(t,Ti_25,'Color',[0.4660 0.6740 0.1880],'LineWidth',0.5);
plot(t,Ti_75,'Color',[0.4660 0.6740 0.1880],'LineWidth',0.5);
fill([t' fliplr(t')],[Ti_25 fliplr(Ti_75)],[0.4660 0.540 0.1880]);
fill([t' fliplr(t')],[Ti_5 fliplr(Ti_25)],'g');
fill([t' fliplr(t')],[Ti_75 fliplr(Ti_95)],'g');
scatter(t_Herz,Tp_Herz,sz,'MarkerEdgeColor',[0.9290 0.6940 0.1250],...
    'MarkerFaceColor',[0.9100 0.4100 0.1700],...
    'LineWidth',1.5);
xlabel('Time (Gyr)'); ylabel('Mantle potential temperature (\circ C)');
xlim([0,4.567]);
xticklabels({'0','','','','','1','','','','','2',...
    '','','','','3','','','','','4','','',''}); xticks(0:0.2:4.567);
yticklabels({'1200','','1300','','1400','','1500','','1600','','1700',...
    '','1800'}); yticks(1200:50:1800);
box on;

subplot(2,2,2); hold off;
plot(t,H_25,'k--',t,H_75,'k--',t,H_5,'k:',t,H_95,'k:',t,H_50,'r-'); hold on;
plot(t,Q_25,'k--',t,Q_75,'k--',t,Q_5,'k:',t,Q_95,'k:',t,Q_50,'b-'); 
fill([t' fliplr(t')],[H_25 fliplr(H_75)],[0.3 0.5 0.3]);
fill([t' fliplr(t')],[H_5 fliplr(H_25)],'y');
fill([t' fliplr(t')],[H_75 fliplr(H_95)],'y');
fill([t' fliplr(t')],[Q_25 fliplr(Q_75)],[0.4660 0.540 0.1880]);
fill([t' fliplr(t')],[Q_5 fliplr(Q_25)],'g');
fill([t' fliplr(t')],[Q_75 fliplr(Q_95)],'g');
xlabel('Time (Gyr)');ylabel('H and Q (TW)');xlim([0,4.567]);

subplot(2,2,3); hold off;
plot(t,V_25,'k',t,V_75,'k',t,V_5,'k',t,V_95,'k');hold on;
fill([t' fliplr(t')],[V_25 fliplr(V_75)],[0.4660 0.540 0.1880]);
fill([t' fliplr(t')],[V_5 fliplr(V_25)],'g');
fill([t' fliplr(t')],[V_75 fliplr(V_95)],'g');
xlabel('Time (Gyr)');ylabel('Plate velocity'); xlim([0,4.567]);

subplot(2,2,4); hold off;
plot(t,Z_25,'k',t,Z_75,'k',t,Z_5,'k',t,Z_95,'k'); hold on;
fill([t' fliplr(t')],[Z_25 fliplr(Z_75)],[0.4660 0.540 0.1880]);
fill([t' fliplr(t')],[Z_5 fliplr(Z_25)],'g');
fill([t' fliplr(t')],[Z_75 fliplr(Z_95)],'g');
xlabel('Time (Gyr)');ylabel('Melting depth'); xlim([0,4.567]);

myfig = figure(2);
myfig.Renderer = 'Painters';

%% plot the selected f_LT,f_HT, and recycling rates
nt = 3240;
[f_LT_5,f_LT_25,f_LT_50,f_LT_75,f_LT_95] = calculate_percentile_fun(f_LT_success,nt,t,size_oxygen);
[f_HT_5,f_HT_25,f_HT_50,f_HT_75,f_HT_95] = calculate_percentile_fun(f_HT_success,nt,t,size_oxygen);
% [f_LT_Krc_5,f_LT_Krc_25,f_LT_Krc_50,f_LT_Krc_75,f_LT_Krc_95] = calculate_percentile_fun(f_LT_Krc_success,nt,t,size_oxygen);

figure(3);
subplot(2,2,1);
plot(tspan,f_LT_25,'k',tspan,f_LT_75,'k',tspan,f_LT_5,'k',tspan,f_LT_95,'k');hold on;
fill([tspan' fliplr(tspan')],[f_LT_25 fliplr(f_LT_75)],[0.4660 0.540 0.1880]);
fill([tspan' fliplr(tspan')],[f_LT_5 fliplr(f_LT_25)],'g');
fill([tspan' fliplr(tspan')],[f_LT_75 fliplr(f_LT_95)],'g'); 
plot(tspan,f_HT_25,'k',tspan,f_HT_75,'k',tspan,f_HT_5,'k',tspan,f_HT_95,'k');
fill([tspan' fliplr(tspan')],[f_HT_25 fliplr(f_HT_75)],[0.3 0.5 0.3]);
fill([tspan' fliplr(tspan')],[f_HT_5 fliplr(f_HT_25)],'y');
fill([tspan' fliplr(tspan')],[f_HT_75 fliplr(f_HT_95)],'y'); 

xlabel('Time (Gyr)','Fontsize',16); xlim([0,3.24]);
ylabel('f_{LT} and f_{HT}','Fontsize',16); 
legend('f_{LT} (Low T alteration rate)','f_{HT} (High T alteration rate)','Location','Northeast','Fontsize',14);
set(gca,'FontSize',16); box on;

myfig = figure(3);
myfig.Renderer = 'Painters';

%% load in the data for delta 18O 
data = xlsread('delta 18O of ocean in phanerozoic.xlsx'); % oxygen of seawater from past studies
age = data(:,1); age = 3.24-age;
deltaO= data(:,2);

%% Plot the necessary figures
nt = length(tspan);
% Reproduce Fig.3 of Wallmann, 2001
% First, convert mass from mole to g
[M_f_g_5,M_f_g_25,M_f_g_50,M_f_g_75,M_f_g_95] = calculate_percentile_fun(M_f_g_success,nt,t,size_oxygen);

figure(4);
subplot(2,1,1);
plot(tspan,M_f_g_25,'k',tspan,M_f_g_75,'k',tspan,M_f_g_5,'k',tspan,M_f_g_95,'k');hold on;
fill([tspan' fliplr(tspan')],[M_f_g_25 fliplr(M_f_g_75)],[0.4660 0.540 0.1880]);
fill([tspan' fliplr(tspan')],[M_f_g_5 fliplr(M_f_g_25)],'g');
fill([tspan' fliplr(tspan')],[M_f_g_75 fliplr(M_f_g_95)],'g'); 
scatter(tspan(nt), M_f_observe,90,'MarkerEdgeColor',[0.9290 0.6940 0.1250],'MarkerFaceColor',[0.9100 0.4100 0.1700],'LineWidth',1.5);
hold off;
xlabel('Time (Gyr)','Fontsize',16);xlim([0,3.24]);
ylabel('Mass of free water (g)','Fontsize',16);
set(gca,'FontSize',16); box on;

myfig = figure(4);
myfig.Renderer = 'Painters';

%% 
[delta_f18_all_5,delta_f18_all_25,delta_f18_all_50,delta_f18_all_75,delta_f18_all_95] = ...
    calculate_percentile_fun(delta_f18_all_success,nt,t,size_oxygen);
[delta_HT18_all_5,delta_HT18_all_25,delta_HT18_all_50,delta_HT18_all_75,delta_HT18_all_95] = ...
    calculate_percentile_fun(delta_HT18_all_success,nt,t,size_oxygen);
[delta_LT18_all_5,delta_LT18_all_25,delta_LT18_all_50,delta_LT18_all_75,delta_LT18_all_95] = ...
    calculate_percentile_fun(delta_LT18_all_success,nt,t,size_oxygen);
[delta_w18_all_5,delta_w18_all_25,delta_w18_all_50,delta_w18_all_75,delta_w18_all_95] = ...
    calculate_percentile_fun(delta_w18_all_success,nt,t,size_oxygen);

%% Dealta 18O data from previous studies
% Tartese et al., 2017 data
age_1 = age(32:35);
deltaO_1 = deltaO(32:35);
% Came et al., 2007 data
age_2 = age(36:37);
deltaO_2 = deltaO(36:37);
% Dennis et al., 2013 data
age_3 = age(38:39);
deltaO_3 = deltaO(38:39);
% Bergmann et al., 2018 data
age_4 = age(40:43);
deltaO_4 = deltaO(40:43);
% Galili et al., 2019 data
age_5 = age(1:31);
deltaO_5 = deltaO(1:31);
%%
figure(5);
subplot(2,1,1);
% plot(tspan,delta_f18_all_50,'k',tspan,delta_f18_all_5,'k-.',tspan,delta_f18_all_95,'k-.','LineWidth',2);hold on;
plot(tspan,delta_f18_all_25,'k',tspan,delta_f18_all_75,'k',tspan,delta_f18_all_5,'k',tspan,delta_f18_all_95,'k');hold on;
fill([tspan' fliplr(tspan')],[delta_f18_all_25 fliplr(delta_f18_all_75)],[0.4660 0.540 0.1880]);
fill([tspan' fliplr(tspan')],[delta_f18_all_5 fliplr(delta_f18_all_25)],'g');
fill([tspan' fliplr(tspan')],[delta_f18_all_75 fliplr(delta_f18_all_95)],'g'); 
scatter(age_1,deltaO_1,90,'d','MarkerEdgeColor','k','MarkerFaceColor','r','LineWidth',1.5);% Tartese et al., 2017 data
scatter(age_2,deltaO_2,90,'<','MarkerEdgeColor','k','MarkerFaceColor','b','LineWidth',1.5);% Came et al., 2007 data
scatter(age_3,deltaO_3,90,'>','MarkerEdgeColor','k','MarkerFaceColor','m','LineWidth',4);% Dennis et al., 2013 data
scatter(age_4,deltaO_4,90,'v','MarkerEdgeColor','k','MarkerFaceColor','c','LineWidth',1.5);% Bergmann et al., 2018 data
scatter(age_5,deltaO_5,90,'MarkerEdgeColor',[0.9290 0.6940 0.1250],'MarkerFaceColor',[0.9100 0.4100 0.1700],'LineWidth',1.5);% Galili et al., 2019 data
scatter(0,3,90,'MarkerEdgeColor',[0.9290 0.6940 0.1250],'MarkerFaceColor',[0.9100 0.4100 0.1700],'LineWidth',1.5);% Johnson and Wing, 2020 data
errorbar(0,3,1);
xlabel('Time (Gyr)','Fontsize',16);xlim([0,3.24]);
ylabel('\delta^{18}O of free water','Fontsize',16);
set(gca,'FontSize',16); box on;

subplot(2,1,2);
plot(tspan,delta_HT18_all_25,'k',tspan,delta_HT18_all_75,'k',tspan,delta_HT18_all_5,'k',tspan,delta_HT18_all_95,'k');hold on;
fill([tspan' fliplr(tspan')],[delta_HT18_all_25 fliplr(delta_HT18_all_75)],[0.7 0.7 0.7]);
fill([tspan' fliplr(tspan')],[delta_HT18_all_5 fliplr(delta_HT18_all_25)],'y');
fill([tspan' fliplr(tspan')],[delta_HT18_all_75 fliplr(delta_HT18_all_95)],'y'); 

plot(tspan,delta_LT18_all_25,'k',tspan,delta_LT18_all_75,'k',tspan,delta_LT18_all_5,'k',tspan,delta_LT18_all_95,'k');hold on;
fill([tspan' fliplr(tspan')],[delta_LT18_all_25 fliplr(delta_LT18_all_75)],[0.4660 0.540 0.1880]);
fill([tspan' fliplr(tspan')],[delta_LT18_all_5 fliplr(delta_LT18_all_25)],'g');
fill([tspan' fliplr(tspan')],[delta_LT18_all_75 fliplr(delta_LT18_all_95)],'g'); 

plot(tspan,delta_w18_all_25,'k',tspan,delta_w18_all_75,'k',tspan,delta_w18_all_5,'k',tspan,delta_w18_all_95,'k');hold on;
fill([tspan' fliplr(tspan')],[delta_w18_all_25 fliplr(delta_w18_all_75)],[0.6 0.6 0.6]);
fill([tspan' fliplr(tspan')],[delta_w18_all_5 fliplr(delta_w18_all_25)],'c');
fill([tspan' fliplr(tspan')],[delta_w18_all_75 fliplr(delta_w18_all_95)],'c'); 

legend('Deep crust','Upper crust','Clay','Location','Northwest');
hold off;
xlabel('Time (Gyr)','Fontsize',16);xlim([0,3.24]);
ylabel('\delta^{18}O of oceanic crust','Fontsize',16);
set(gca,'FontSize',16); box on;

myfig = figure(5);
myfig.Renderer = 'Painters';

%% plot formation age and surface age distributions
% Set the time period for running crustal growth and thermal evolution functions
tmax = 4.567;% the age of solar system, in unit Ga
dt = 0.001;% length of each timestep
t = 0:dt:tmax;
nt = length(t);% number of timesteps
t = t';

% formation age and surface age distributions
[F_5,F_25,F_50,F_75,F_95] = calculate_percentile_fun(F_success,nt,t,size_oxygen);
[S_5,S_25,S_50,S_75,S_95] = calculate_percentile_fun(S_success,nt,t,size_oxygen);

% Formation age distribution data from Korenaga (2018a)
data_formationage = load('korenaga18a_Tunmix_orig.dat');
t_Jun = data_formationage(:,1);
t_Jun = 4.568-(t_Jun/1000);
F_Jun = data_formationage(:,2);
% Surface age distribution data from Roberts & Spencer (2015)
data_zircon_surf = load('korenaga18a_T_U_Pb.dat');
S_Jun = data_zircon_surf(:,2);

figure(6);
subplot(2,1,1); hold off;
plot(t,F_25,'Color',[0.4660 0.6740 0.1880],'LineWidth',0.5);hold on;
plot(t,F_75,'Color',[0.4660 0.6740 0.1880],'LineWidth',0.5);
plot(t,F_5,'g-','LineWidth',0.5);
plot(t,F_95,'g-','LineWidth',0.5);
% plot(t,F_50,'r-','LineWidth',2);
fill([t' fliplr(t')],[F_25 fliplr(F_75)],[0.4660 0.540 0.1880]);
fill([t' fliplr(t')],[F_5 fliplr(F_25)],'g');
fill([t' fliplr(t')],[F_75 fliplr(F_95)],'g');
plot(t_Jun, F_Jun,'Color',[0.9100 0.4100 0.1700],'LineWidth',2.5);
xlabel('Time (Gyr)','FontSize',16); ylabel('Formation age distribution');
xlim([0,4.567]); ylim([0,1]);
xticklabels({'0','','','','','1','','','','','2',...
    '','','','','3','','','','','4','','',''}); xticks(0:0.2:4.567);
yticklabels({'0','','0.2','','0.4','','0.6','','0.8','','1.0'}); yticks(0:0.1:1);
set(gca,'FontSize',14); box on;

subplot(2,1,2); hold off;
plot(t,S_25,'Color',[0.4660 0.6740 0.1880],'LineWidth',0.5);hold on;
plot(t,S_75,'Color',[0.4660 0.6740 0.1880],'LineWidth',0.5);
plot(t,S_5,'g-','LineWidth',0.5);
plot(t,S_95,'g-','LineWidth',0.5);
fill([t' fliplr(t')],[S_25 fliplr(S_75)],[0.4660 0.540 0.1880]);
fill([t' fliplr(t')],[S_5 fliplr(S_25)],'g');
fill([t' fliplr(t')],[S_75 fliplr(S_95)],'g');
plot(t_Jun,S_Jun,'Color',[0.9100 0.4100 0.1700],'LineWidth',2.5);
% plot(t,S_50,'r-','LineWidth',2);
xlabel('Time (Gyr)','FontSize',16); ylabel('Surface age distribution');
xlim([0,4.567]); ylim([0,1]);
xticklabels({'0','','','','','1','','','','','2',...
    '','','','','3','','','','','4','','',''}); xticks(0:0.2:4.567);
yticklabels({'0','','0.2','','0.4','','0.6','','0.8','','1.0'}); yticks(0:0.1:1);
box on;
set(gca,'FontSize',14); box on;

myfig = figure(6);
myfig.Renderer = 'Painters';
