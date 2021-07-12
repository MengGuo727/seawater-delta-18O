function [Kmo,Kmc,Krc, Kmc_c] = ...
          Kmo_fun(t,Z,V,Kmo_tp,Mc,Mdd,Mud,Mcp,Mmp,...
          K40_BSE_tp,K40_CC_observe)
% This function simulates Ar degassing history of the Earth
% 
% Meng Guo, Yale University
% Summer 2019


% Input variables and their meanings:
% t: time, from 0 Gyr to 4.567 Gyr;
% V: plate velocity, in unit cm/yr, backward in time
% Z: the initial depth of mantle melting, in unit meter, backward in time
% Kmo_tp: presdent-day mass of the mantle processed 
%         by magmatism to produce oceanic crust, in unit kg/Gyr
% Mc: history of CC mass in unit kg
% Mdd: rate of continental crust recycling, in unit kg/Gyr
% Mud: rate of mass that is going from the mantle to create the crust, in unit kg/Gyr
% Mcp: present-day continental crust mass, in unit kg
% Mmp: present-day mantle mass, in unit kg

%% Time constant
nt = length(t);

%% use Z and V to calculate the changing of Kmo through time
Kmo = nan(size(t)); 
% Kmo(1) = Kmo_tp;
Z_tp = Z(nt);
V_tp = V(nt);
for i = 1:nt
    Kmo(i) = Kmo_tp * (Z(i)*V(i)/(Z_tp * V_tp));
end

%% set the other mass exchange rate
Kmc = Mud; % rate of mantle processed by magmatism to produce continental crust
Krc = Mdd; % rate of continental crust recycling

%% calculate resevoirs masses
M_BSE_tp = Mmp + Mcp ;% mass of BSE, in unit kg
M_BSE = nan(size(t)); 
for i = 1:nt
    M_BSE(i) = M_BSE_tp;
end
M_MM = nan(size(t));
M_CC = Mc;
for i = 1:nt
    M_MM(i) = M_BSE(i) - M_CC(i); 
end

%% decay coinstants
lamda_K = 5.305e-1;% Total decay constant of K, in unit Gyr-1

%% calcualte K40 in CC and MM 
K40_MM = nan(size(t));
K40_CC = nan(size(t));
K40_CC_t0 =  K40_CC_observe * exp(lamda_K*t(nt));
K40_BSE= nan(size(t));
K40_BSE_t0 =  K40_BSE_tp * exp(lamda_K*t(nt));
K40_BSE(1) = K40_BSE_t0;

for i = 2:nt
    K40_BSE(i) = K40_BSE(1)*exp(-lamda_K*t(i));
end

Mc_f = Mc/Mcp;

for i = 1:nt
    K40_CC(i) = K40_CC_t0 *exp(-lamda_K*t(i))* Mc_f(i);
    K40_MM(i) = K40_BSE(i)- K40_CC(i);
    if (K40_CC(i)<0)
        K40_CC(i)=0;
    end
end

%% set the crust growth and crust recycling rate, and the mantle processed rate
% calculate the mantle processed to generate continental crust
Kmc_c = nan(size(t)); % mantle processed rate to generate CC from mantle, in addition to secondary metling from OC

for i = 1:nt
    % mass balance between the mass transport to crust and the
    % mantle processed to create crust
    if M_CC(i) == 0
      Kmc_c(i) = 0;
    else
      Kmc_c(i) = (Kmc(i)/M_CC(i))*K40_CC(i)/K40_MM(i)*M_MM(i);
    end
end

for i = 1:nt
    Kmc(i) = max(Kmc_c(i)-Kmo(i),0);
end
