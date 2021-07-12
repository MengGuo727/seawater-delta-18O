%% 2019/08/05 Phase 1B of argon degassing -- Crustal growth model
% This file read "out_phase1a.dat" and select successful solutions
% according to the misfit values of your choice. The successful solutions
% are stored in file "out_phase1b.dat".
%
% Meng Guo, Yale University
% Summer, 2019

clear all;

% open output file
fileID = fopen('out_phase1b.dat','w');

%% read file "out_phase1a.dat"
data = load('out_phase1a.dat');
% the variables in 'out_phase1a.dat' are in following order:
% 1.Krw_factor, 2.kappa_r, 3.kappa_g, 4.Rs, 5.Rp, 6.ts, 7.RMSE_F, 8.RMSE_S

%% define the misfits for formation age and surface age distributions
RMSE_F = data(:,7);
RMSE_S = data(:,8);

%% select successful solutions for plotting results
[ndata,nparam] = size(data);% size of the data file
j = 0; % counter for the number of successful solutions

crit_misfit_FS = 0.2;% set the selection criterion for misfit
success_FS =[];% empty matrix for saving successful solutions

for i = 1:ndata
    if (RMSE_F(i)<crit_misfit_FS && RMSE_S(i)<crit_misfit_FS)
        j = j+1;
        for k = 1:nparam
            success_FS(j,k) = data(i,k);
        end % for k = 1:nparam
        
        % Save all the variables and RMSEs
        fprintf(fileID,['%6g %6g %6g ' ...
            '%6g %6g %6g ' ...
            '%6g  %6g\n'], ...
            success_FS(j,1), ...
            success_FS(j,2), success_FS(j,3), success_FS(j,4), ...
            success_FS(j,5), success_FS(j,6),...
            success_FS(j,7), success_FS(j,8) );
    end % if (RMSE_F(i)<0.2 && RMSE_S(i)<0.2)
    
end % for i = 1:ndata

fclose(fileID);

%% Set constants for further calculations
% Set the time period
tmax = 4.567;% the age of solar system, in unit Ga
dt = 0.001;% length of each timestep
t = 0:dt:tmax;
nt = length(t);% number of timesteps
t = t';

% constants for crustal growth pattern
Mcp = 2.09e22;% mass of continental crust at present-day, in unit kg

%% calculate the covariation and correlation matrixes, plot correlation
figure(1);
subplot(2,3,1);
histogram(success_FS(:,1)); xlabel('init Krw factor'); ylabel('count');
subplot(2,3,2);
histogram(success_FS(:,2)); xlabel('Kappa_r');
subplot(2,3,3);
histogram(success_FS(:,3)); xlabel('Kappa_g');
subplot(2,3,4);
histogram(success_FS(:,4)); xlabel('Rs'); ylabel('Count')
subplot(2,3,5);
histogram(success_FS(:,5)); xlabel('Rp');
subplot(2,3,6);
histogram(success_FS(:,6)); xlabel('Ts');

% calculate the covariation matrix and correlation matrix
correlation_FS = corrcoef(success_FS(:,1:6));

figure(2);
subplot(2,2,1);
plot(success_FS(:,1),success_FS(:,5),'b.');
xlabel('init Krw fac'); ylabel('Rp');
subplot(2,2,2);
plot(success_FS(:,2),success_FS(:,4),'b.');
xlabel('kappa_r'); ylabel('Rs');
subplot(2,2,3);
plot(success_FS(:,2),success_FS(:,5),'b.');
xlabel('kappa_r'); ylabel('Rp');


%% Calculate the data for the middle 50% and 90% of the selected solutions
% calculate formation age and surface age distributions
[size_FS,nparameter]= size(success_FS);
F = nan(nt,size_FS);
S = nan(nt,size_FS);
Mc = nan(nt,size_FS);
Mud = nan(nt,size_FS);
Mdd = nan(nt,size_FS);

for l=1:size_FS
    if mod(l,10)==0
        disp(['l=' num2str(l) ' of ' num2str(size_FS)]);% keep track of the calculation
    end
    % independent main variables
    Krw_factor_model = success_FS(l,1);
    kappa_r_model = success_FS(l,2);
    kappa_g_model = success_FS(l,3);
    Rs_model = success_FS(l,4);
    Rp_model = success_FS(l,5);
    ts_model = success_FS(l,6);
    
    % Calculate the dependent variable (crustal reworking rate)
    Krw_s_model = Rs_model * Krw_factor_model;% initial Krw_factor
    
    % Calculate the corresponding crustal growth pattern
    [Mc_model,Mdd_model,Mud_model,Krw_model_first] = CC_growth_fun1(t,ts_model,tmax,...
        Mcp,kappa_g_model,Rp_model,Rs_model,kappa_r_model,Krw_s_model);
    
    % Calculate formation age and surface age distributions and recalculate
    % crustal reworking rate using Formation_surface_age_fun function
    [F_model,S_model,m_tp,m,Krw_model] = ...
        Formation_surface_age_fun(t,Mud_model,Mdd_model,Mc_model,Krw_model_first);
    
    % save the results
    F(:,l) = F_model;
    S(:,l) = S_model;
    Mc(:,l) = Mc_model;
    Mud(:,l) = Mud_model;
    Mdd(:,l) = Mdd_model;
    
end % for l=1:itermax

%% Load in the observed formation age and surface age distributions
% Formation age distribution data from Korenaga (2018a)
data_formationage = load('korenaga18a_Tunmix_orig.dat');
% Surface age distribution data from Roberts & Spencer (2015)
data_zircon_surf = load('korenaga18a_T_U_Pb.dat');
% use load_FandS_fun function to set data in the same dimension as the time series
[F_Jun_same,S_Jun_same] = load_FandS_fun(t,nt,data_formationage,data_zircon_surf);

%% plot figures of successful runs from stage 1
%formation age distribution
[F_5,F_25,F_50,F_75,F_95] = calculate_percentile_fun(F,nt,t,size_FS);
figure(3);
subplot(2,1,1); hold off;
plot(t,F_25,'k--',t,F_75,'k--',t,F_5,'k:',t,F_95,'k:',t,F_50,'r-',...
    t, F_Jun_same,'c-');
xlabel('Time (Gyr)');
ylabel('Formation age distribution');
grid on;

% plot surface age distribution
[S_5,S_25,S_50,S_75,S_95] = calculate_percentile_fun(S,nt,t,size_FS);
subplot(2,1,2); hold off;
plot(t,S_25,'k--',t,S_75,'k--',t,S_5,'k:',t,S_95,'k:',t,S_Jun_same,'b-',...
    t,S_50,'r-');
xlabel('Time (Gyr)');
ylabel('Surface age distribution');
grid on;

%% plot net crustal growth(Mc), crustal generation rate(Mud), ...
%  and crustal recycling rate (Mdd)
[Mc_5,Mc_25,Mc_50,Mc_75,Mc_95] = calculate_percentile_fun(Mc,nt,t,size_FS);
[Mud_5,Mud_25,Mud_50,Mud_75,Mud_95] = calculate_percentile_fun(Mud,nt,t,size_FS);
[Mdd_5,Mdd_25,Mdd_50,Mdd_75,Mdd_95] = calculate_percentile_fun(Mdd,nt,t,size_FS);

figure(4);
subplot(2,2,1); hold off;
plot(t,Mc_25,'k--',t,Mc_75,'k--',t,Mc_5,'k:',t,Mc_95,'k:',t,Mc_50,'r-');
xlabel('Time (Gyr)'); ylabel('Net crustal growth'); grid on;

subplot(2,2,2); hold off;
semilogy(t,Mud_25,'k--',t,Mud_75,'k--',t,Mud_5,'k:',t,Mud_95,'k:',t,Mud_50,'r-'); hold on;
xlabel('Time (Gyr)'); ylabel('Crustal generation rate'); grid on;

subplot(2,2,3); hold off;
semilogy(t,Mdd_25,'k--',t,Mdd_75,'k--',t,Mdd_5,'k:',t,Mdd_95,'k:',t,Mdd_50,'b-');
xlabel('Time (Gyr)'); ylabel('Crustal recycling rate'); grid on;
