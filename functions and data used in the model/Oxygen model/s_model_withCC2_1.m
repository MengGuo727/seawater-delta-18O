function [M_f, M_f18, M_w, M_w18, M_LT, M_LT18, M_HT, M_HT18,...
    F_w, F_LT, F_HT, F_rew, F_reLT, F_reHT, ...
    F_w18, F_LT18, F_HT18, F_w_ex, ...
    F_LT_ex, F_HT_ex, F_pw18, phi_w, phi_LT, phi_HT, F_rew18, ...
    F_reLT18, F_reHT18, F_m18, F_mw18, F_mLT18, F_mHT18, F_sw18, ...
     F_sLT18, F_sHT18,phi_f, R_f, omega_w, omega_LT , omega_HT,...
     F_reverse, F_reverse18, F_reverse_ex, F_m] = ...
    s_model_withCC2_1(y_ini, tspan, f_LT,f_HT,f_V, f_Krc, f_Kmc, f_Kmp, F_w_0, F_LT_0, F_HT_0, r_w,...
    r_LT, r_HT, alpha_w, alpha_LT,alpha_HT, alpha_rev,R_si, k_w, ...
    R_oc, k_LT, k_HT,k_rev,F_pw18_0, phi_m, phi_si, F_mw_0, phi_oc, ...
    F_mLT_0, F_mHT_0, dt,F_m_0, f_reverse)

%% The number of time steps
nt = length(tspan);

%% set-up empty matrix to store each variables
% mass for each reservoir
M_f= nan(size(tspan));   % free water mass
M_w = nan(size(tspan));  % 18O in sedimentary rocks
M_LT = nan(size(tspan));   % 18O in upper crust
M_HT = nan(size(tspan)); % 18O in deep crust
M_f18 = nan(size(tspan));   % free H218O mass
M_w18 = nan(size(tspan));  % 18O in sedimentary rocks
M_LT18 = nan(size(tspan));   % 18O in upper crust
M_HT18 = nan(size(tspan)); % 18O in deep crust
% mass flux of H2O
F_w =  nan(size(tspan));   % H2O fixation in weathering products
F_LT = nan(size(tspan));   % H2O uptake in low T altered OC 
F_HT = nan(size(tspan));  % H2O uptake in high T altered OC (at ridge or hot spots)
F_rew = nan(size(tspan));   % H2O release from weathering products
F_reLT = nan(size(tspan));   % H2O release from low T altered OC  at subduction zones
F_reHT = nan(size(tspan));   % H2O release from high T altered OC  at subduction zones
F_m = nan(size(tspan));   % release H20 from mantle

F_mw = nan(size(tspan));   % Formation rate of continental sediments
F_sw = nan(size(tspan));   % subduction rate of sediments
F_mLT = nan(size(tspan));   % Formation rate of low temperature altered OC
F_sLT = nan(size(tspan));   % subduction rate of low temperature altered OC
F_mHT = nan(size(tspan));   % Formation rate of high temperature altered OC
F_sHT = nan(size(tspan));   % subduction rate of high temperature altered OC

% mass flux of 18O
F_w18 = nan(size(tspan));   % H218O uptake in continental weathering products
F_LT18 = nan(size(tspan));   % H218O uptake in low T altered OC
F_HT18 = nan(size(tspan));   % H218O uptake in high T altered OC (at ridge or hot spots)
F_w_ex = nan(size(tspan));   % 18O exchange during weathering
F_LT_ex = nan(size(tspan));   % 18O exchange between low T altered OC and seawater
F_HT_ex = nan(size(tspan));   % 18O exchange between high T altered OC and seawater
F_pw18 = nan(size(tspan));   % 18O loss during turnover of porewater
F_rew18 = nan(size(tspan));  % H218O release from weathering products
F_reLT18 = nan(size(tspan));   % H218O release from low T altered OC at subduction zones
F_reHT18 = nan(size(tspan));   % H218O release from high T altered OC at subduction zones
F_m18 = nan(size(tspan));   % H218O release through mantle degassing
F_mw18 = nan(size(tspan));   % formation of sedimentary rock 18O
F_mLT18 = nan(size(tspan));   % formation of low T altered OC 18O at spreading zones
F_mHT18 = nan(size(tspan));   % formation of high T altered OC 18O at spreading zones
F_sw18 = nan(size(tspan));   % subduction of sedimentary rock 18O
F_sLT18 = nan(size(tspan));   % subduction of low T altered OC 18O
F_sHT18 = nan(size(tspan));   % subduction of high T altered OC 18O
F_reverse = nan(size(tspan));   % reverse weathering
F_reverse18 = nan(size(tspan));   % reverse weathering 18O
F_reverse_ex = nan(size(tspan));   % reverse weathering exchanged 18O
phi_f = nan(size(tspan)); % 18O mole fraction in free water
R_f = nan(size(tspan)); % 18O/16O ratios in free water
omega_w = nan(size(tspan));% isotope saturation index for 18O exchange in silicate weathering
omega_LT = nan(size(tspan));% isotope saturation index for 18O exchange between low T altered OC and seawater
omega_HT = nan(size(tspan));% isotope saturation index for 18O exchange between high T altered OC and seawater
omega_reverse = nan(size(tspan));% isotope saturation index for 18O exchange during reverse weathering
phi_w = nan(size(tspan));% 18O model fraction in weathering product
phi_LT = nan(size(tspan));% 18O model fraction in low T altered OC
phi_HT = nan(size(tspan));% 18O model fraction in high T altered OC

% Set the initial values, cause our for loop starts from 2nd time step
M_f(1)= y_ini(1);   % free water mass
M_f18(1) = y_ini(2);   % free H218O mass
M_w(1) = y_ini(3);   % free water mass in sedimentary rocks
M_w18(1) = y_ini(4);  % 18O in sedimentary rocks
M_LT(1) = y_ini(5);   % free water mass in low T altered OC
M_LT18(1) = y_ini(6);   % 18O in low T altered OC
M_HT(1) = y_ini(7); % free water mass in high T altered OC
M_HT18(1) = y_ini(8); % 18O in high T altered OC


for i = 2:nt
   
    if (M_f(i-1) == 0)
        phi_f(i) = 0;
    else
        phi_f(i) = M_f18(i-1)  / M_f(i-1) ; % 18O mole fraction in free water
    end
    R_f(i)  = phi_f(i)  / (1 - phi_f(i)); % 18O/16O ratios in free water

    omega_w(i)  = (R_si(i) / R_f(i) ) / alpha_w;% isotope saturation index for 18O exchange in silicate weathering
    omega_LT(i)  = (R_oc / R_f(i) ) / alpha_LT;% isotope saturation index for 18O exchange between low T altered OC and seawater
    omega_HT(i)  = (R_oc / R_f(i) ) / alpha_HT;% isotope saturation index for 18O exchange between high T altered OC and seawater
    omega_reverse(i)  = (R_oc / R_f(i) ) / alpha_rev;% isotope saturation index for 18O exchange between DC and seawater
    
    % Calcualte the water fluxes and 18O turnover processes, all the fluxes
    % F are in unit mol/Gyr
    % Water uptake or fixation
    %     F_m(i) = (f_Kmc(i) + f_Kmp(i))/(f_Kmc(3500)+f_Kmp(3500)) * F_m_0;  % H2O release from mantle magmatism
    F_m(i) = f_HT(i)* F_m_0;
    F_LT(i) = f_LT(i) * f_V(i) * F_LT_0 ; % H2O uptake in low T altered OC
    F_HT(i) = f_HT(i) * F_HT_0 ;  % H2O uptake in high T altered OC
    F_w(i) =  f_LT(i) * f_Krc(i) * F_w_0;  % H2O fixation in weathering products
    F_reverse(i) = f_reverse(i) * (F_w(i) + F_LT(i)); % H2O fixation in sediments during reverse weathering
    
    % Water release
    F_rew(i) = r_w(i) * (F_w(i) + F_reverse(i));   % H2O release from weathering products, with reverser weathering 
    F_reLT(i) = r_LT(i) * F_LT(i);   % H2O release from low T altered OC at subduction zones
    F_reHT(i) = r_HT(i) * F_HT(i);   % H2O release from high T altered OC at subduction zones
    
    % formation and subduction of OC and sediments
    F_mw(i) = f_Krc(i) * F_mw_0;   % Formation rate of continental sediments
    F_sw(i) =  F_mw(i) + (1 - r_w(i)) * (F_w(i) + F_reverse(i));   % subduction rate of sediments
    F_mLT(i) = f_LT(i) * F_mLT_0;   % Formation rate of low temperature altered OC
    F_sLT(i) = F_mLT(i) + (1 - r_LT(i)) * F_LT(i);   % subduction rate of low temperature altered OC
    F_mHT(i) = f_HT(i) * F_mHT_0;   % Formation rate of high temperature altered OC
    F_sHT(i) = F_mHT(i) + (1 - r_HT(i)) * F_HT(i);   % subduction rate of high temperature altered OC
    
    % H218O uptake
    F_w18(i) = alpha_w * R_f(i) * F_w(i)/ (1 + alpha_w * R_f(i)) ;   % H218O uptake in continental weathering products
    F_LT18(i) = alpha_LT * R_f(i) * F_LT(i)/ (1 + alpha_LT * R_f(i)) ;   % H218O uptake in low T altered OC
    F_HT18(i) = alpha_HT * R_f(i) * F_HT(i)/ (1 + alpha_HT * R_f(i)) ;   % H218O uptake in high T altered OC
    F_reverse18(i) = alpha_rev * R_f(i) * F_reverse(i)/ (1 + alpha_rev * R_f(i)) ; % H218O uptake in reverse weathering
    
    % 18O loss during porewater formation and recycling
    F_pw18(i) = f_LT(i) * F_pw18_0;
    
    % 18O release
    phi_w(i) = M_w18(i-1) / M_w(i-1);
    phi_LT(i) = M_LT18(i-1) / M_LT(i-1);
    phi_HT(i) = M_HT18(i-1) / M_HT(i-1);
    F_rew18(i) = phi_w(i) * F_rew(i);  % H218O release from weathering products
    F_reLT18(i) = phi_LT(i) * F_reLT(i);   % H218O release from low T altered OC at subduction zones
    F_reHT18(i) = phi_HT(i) * F_reHT(i);   % H218O release from high T altered OC at subduction zones
    F_m18(i) = phi_m * F_m(i);   % H218O release through mantle degassing
    
    % formation of 18O
    F_mw18(i) = phi_si(i) * F_mw(i);   % formation of sedimentary rock 18O
    F_mLT18(i) = phi_oc * F_mLT(i);   % formation of low T altered OC 18O at spreading zones
    F_mHT18(i) = phi_oc * F_mHT(i);   % formation of high T altered OC 18O at spreading zones
    
    % 18O loss from rocks
    F_sw18(i) = phi_w(i) * F_sw(i);   % metamorphosis of sedimentary rock 18O
    F_sLT18(i) = phi_LT(i) * F_sLT(i);   % subduction of low T altered OC 18O
    F_sHT18(i) = phi_HT(i) * F_sHT(i);   % subduction of high T altered OC 18O
    
    % 18O exchange
    F_w_ex(i) = f_LT(i) * k_w * (omega_w(i) - 1);   % 18O exchange during weathering
    F_LT_ex(i) = f_LT(i) * k_LT * (omega_LT(i) - 1);   % 18O exchange between low T altered OC and seawater
    F_HT_ex(i) = f_HT(i) * k_HT * (omega_HT(i) - 1);   % 18O exchange between high T altered OC and seawater
    F_reverse_ex(i) = f_reverse(i) * f_LT(i) * k_rev * (omega_reverse(i) - 1);
    
    if (i <= 740)
        F_w(i) =  0;  % H2O fixation in weathering products
        F_rew(i) = 0;   % H2O release from weathering products, with reverser weathering
        F_mw(i) = 0;   % Formation rate of continental sediments
        F_sw(i) = 0;   % subduction rate of sediments
        F_w18(i) = 0;   % H218O uptake in continental weathering products
        F_rew18(i) = 0;  % H218O release from weathering products
        F_mw18(i) = 0;   % formation of sedimentary rock 18O
        F_sw18(i) = 0;   % metamorphosis of sedimentary rock 18O
        F_w_ex(i) = 0;   % 18O exchange during weathering
    end
    % change in free H2O mass, without reverse weathering
    dM_f = (F_rew(i) + F_reLT(i) + F_reHT(i) + F_m(i) - F_w(i) - F_reverse(i) - F_LT(i) - F_HT(i))*dt;
    % change in free H2O mass, without reverse weathering
    dM_w = (F_w(i) + F_mw(i) + F_reverse(i) - F_rew(i) - F_sw(i))*dt;
    % change in free H2O mass, without reverse weathering
    dM_LT = (F_LT(i) + F_mLT(i) - F_reLT(i) - F_sLT(i))*dt;
    % change in free H2O mass, without reverse weathering
    dM_HT = (F_HT(i) + F_mHT(i) - F_reHT(i) - F_sHT(i) )*dt;
    
    % change in free H218O mass
    dM_f18 = (F_rew18(i) + F_reLT18(i) + F_reHT18(i) + F_m18(i) + ...
        F_w_ex(i)  + F_LT_ex(i)  + F_reverse_ex(i) + F_HT_ex(i) - F_w18(i)  - F_LT18(i) - F_HT18(i) - F_reverse18(i) - F_pw18(i))*dt;
    % Change in 18O mass in sedimentary rocks
    dM_w18 = (F_mw18(i) + F_w18(i) + F_reverse18(i)+ F_pw18(i) - F_w_ex(i) - F_rew18(i) - F_sw18(i)- F_reverse_ex(i))*dt;
    % Change in 18O mass in low T altered OC 18O
    dM_LT18 = (F_mLT18(i) + F_LT18(i) - F_LT_ex(i) - F_reLT18(i) - F_sLT18(i))*dt;
    % Change in 18O mass in high T altered OC
    dM_HT18 = (F_mHT18(i) + F_HT18(i) - F_HT_ex(i) - F_reHT18(i) - F_sHT18(i))*dt;

    
    M_f(i) = M_f(i-1) + dM_f;
    M_w(i) = M_w(i-1) + dM_w;
    M_LT(i) = M_LT(i-1) + dM_LT;
    M_HT(i) = M_HT(i-1) + dM_HT;    
    M_f18(i) = M_f18(i-1) + dM_f18;
    M_w18(i) = M_w18(i-1) + dM_w18;
    M_LT18(i) = M_LT18(i-1) + dM_LT18;
    M_HT18(i) = M_HT18(i-1) + dM_HT18;

end