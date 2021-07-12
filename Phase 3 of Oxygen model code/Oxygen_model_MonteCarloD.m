%% 2020.03.17 Delta Oxygen Model D -- Plot the a posteriori distributions of parameters
% This file is used to plot histograms of parameters of successful solutions
%
% Meng Guo, Yale University
% Summer, 2019

clear all;
clc;

%% Read data into Matlab
% read file "out_phase1b.dat"
success_FS = load('out_phase1b.dat');
% the variables in 'out_phase1b.dat' are in following order:
% 1.Krw_factor, 2.kappa_r, 3.kappa_g, 4.Rs, 5.Rp, 6.ts,
% 7.RMSE_F, 8.RMSE_S,

% read file "out_phase2b.dat"
success_Tp = load('out_phase2b.dat');
% the variables in 'out_phase2b.dat' are in following order:
% 1.Krw_factor, 2.kappa_r, 3.kappa_g, 4.Rs, 5.Rp, ...
% 6.ts, 7.RMSE_F, 8.RMSE_S, 9.H_BSE_tp, 10.H_cc_tp,...
% 11.Q_total_tp, 12.Qc_tp, 13.d_Qc, 14.misfit_Tp

% read file "out_phase3c.dat"
success_Oxygen = load('out_oxygenmodel2_scenario2.dat');
% 1.Krw_factor, 2.kappa_r, 3.kappa_g, 4.Rs, 5.Rp, ...
% 6.ts, 7.H_BSE_tp, 8.H_cc_tp,...
% 9.Q_total_tp, 10.Qc_tp, 11.d_Qc, ...
% 12.f_reverse_factore, 13.r_factor, 14.alpha_w, 15.alpha_LT, 16.alpha_rev,...
% 17.F_HT_0, 18.M_f_ini, 19.delta_w18_ini, 20.delta_LT18_ini, 21.delta_HT18_ini,...
% 22.delta_f18_ini, 23,misfit_Mf, 24.misfit_deltaO_lateArchean,...
% 25.misfit_deltaO_phanerozoic1, 26.misfit_deltaO_phanerozoic2, 27.misfit_deltaO_phanerozoic3
% 28.n

%% Plot histograms of parameters from crustal growth model
figure(1);
subplot(2,3,1);
histogram(success_FS(:,4),'Normalization','probability','Displaystyle', 'stairs','LineWidth',2,'EdgeColor','g','BinWidth',1e22); ylabel('Probability');  hold on;
histogram(success_Tp(:,4),'Normalization','probability','Displaystyle', 'stairs','LineWidth',2,'EdgeColor','b','LineStyle',':','BinWidth',1e22);
histogram(success_Oxygen(:,4),'Normalization','probability','Displaystyle', 'stairs','LineWidth',2,'EdgeColor','r','BinWidth',1e22); xlabel('R_s (kg/Gyr)'); hold off;
set(gca,'FontSize',14);
subplot(2,3,2);
histogram(success_FS(:,5),'Normalization','probability','Displaystyle', 'stairs','LineWidth',2,'EdgeColor','g','BinWidth',2e21); hold on;
histogram(success_Tp(:,5),'Normalization','probability','Displaystyle', 'stairs','LineWidth',2,'EdgeColor','b','LineStyle',':','BinWidth',2e21);
histogram(success_Oxygen(:,5),'Normalization','probability','Displaystyle', 'stairs','LineWidth',2,'EdgeColor','r','BinWidth',2e21); xlabel('R_p (kg/Gyr)'); hold off;
set(gca,'FontSize',14);
subplot(2,3,3);
histogram(success_FS(:,2),'Normalization','probability','Displaystyle', 'stairs','LineWidth',2,'EdgeColor','g','BinWidth',0.5); hold on;
histogram(success_Tp(:,2),'Normalization','probability','Displaystyle', 'stairs','LineWidth',2,'EdgeColor','b','LineStyle',':','BinWidth',0.5);
histogram(success_Oxygen(:,2),'Normalization','probability','Displaystyle', 'stairs','LineWidth',2,'EdgeColor','r','BinWidth',0.5); xlabel('\kappa_r (Gyr^{-1})');hold off;
set(gca,'FontSize',14);
subplot(2,3,4);
histogram(success_FS(:,3),'Normalization','probability','Displaystyle', 'stairs','LineWidth',2,'EdgeColor','g','BinWidth',3); ylabel('Probability'); hold on;
histogram(success_Tp(:,3),'Normalization','probability','Displaystyle', 'stairs','LineWidth',2,'EdgeColor','b','LineStyle',':','BinWidth',3);
histogram(success_Oxygen(:,3),'Normalization','probability','Displaystyle', 'stairs','LineWidth',2,'EdgeColor','r','BinWidth',3); xlabel('\kappa_g (Gyr^{-1})'); hold off;
legend('Stage 1','Stage 2','Stage 3','Location','Northeast');
set(gca,'FontSize',14);
subplot(2,3,5);
histogram(success_FS(:,6),'Normalization','probability','Displaystyle', 'stairs','LineWidth',2,'EdgeColor','g','BinWidth',0.1); hold on;
histogram(success_Tp(:,6),'Normalization','probability','Displaystyle', 'stairs','LineWidth',2,'EdgeColor','b','LineStyle',':','BinWidth',0.1);
histogram(success_Oxygen(:,6),'Normalization','probability','Displaystyle', 'stairs','LineWidth',2,'EdgeColor','r','BinWidth',0.1); xlabel('t_s (Gyr)'); hold off;
set(gca,'FontSize',14);
subplot(2,3,6);
histogram(success_FS(:,1),'Normalization','probability','Displaystyle', 'stairs','LineWidth',2,'EdgeColor','g','BinWidth',0.1); hold on;
histogram(success_Tp(:,1),'Normalization','probability','Displaystyle', 'stairs','LineWidth',2,'EdgeColor','b','LineStyle',':','BinWidth',0.1);
histogram(success_Oxygen(:,1),'Normalization','probability','Displaystyle', 'stairs','LineWidth',2,'EdgeColor','r','BinWidth',0.1); xlabel('f_{rw}'); hold off;
set(gca,'FontSize',14);
hold off;

myfig = figure(1);
myfig.Renderer = 'Painters';

%% Plot histogram for n
figure(2);
histogram(success_Oxygen(:,28),'Normalization','probability','Displaystyle', 'stairs','LineWidth',2,'EdgeColor','g'); 
xlabel('n')
ylabel('Probability');
set(gca,'FontSize',14);

myfig = figure(2);
myfig.Renderer = 'Painters';
