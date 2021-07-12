%% 2020.03.17 Delta Oxygen Model B -- select successful runs
% This file reads the output from Oxygen_model_MonteCarloA and selects the
% successful runs according to observational constraints. The output file
% "out_oxygenmodel2.dat" contains the parameters for successful solutions
close all;
clear all;
clc;

clear all;

rng('shuffle');% to avoid have same random results everytime

% open output file
fileID = fopen('out_oxygenmodel21_noFandPHI.dat','w');

%% Read file "out_phase3a.dat"
data = load('out_oxygenmodel1_noFandPHI.dat');
% the variables in 'out_phase3a.dat' are in following order:
% 1.Krw_factor, 2.kappa_r, 3.kappa_g, 4.Rs, 5.Rp, ...
% 6.ts, 7.H_BSE_tp, 8.H_cc_tp,...
% 9.Q_total_tp, 10.Qc_tp, 11.d_Qc, ...
% 12.f_reverse_factore, 13.r_factor, 14.alpha_w, 15.alpha_LT, 16.alpha_rev,...
% 17.F_HT_0, 18.M_f_ini, 19.delta_w18_ini, 20.delta_LT18_ini, 21.delta_HT18_ini,...
% 22.delta_f18_ini, 23,misfit_Mf, 24.misfit_deltaO_lateArchean,...
% 25.misfit_deltaO_phanerozoic1, 26.misfit_deltaO_phanerozoic2, 27.misfit_deltaO_phanerozoic3
% 28.n

%% Define the misfits
misfit_Mf = data(:,23);
misfit_deltaO_lateArchean = data(:,24);
misfit_deltaO_phanerozoic1 = data(:,25);
misfit_deltaO_phanerozoic2 = data(:,26);
misfit_deltaO_phanerozoic3 = data(:,27);

%% Select successful solutions for plotting results
[ndata,nparam] = size(data);% the size of the data file
j = 0;% counter for the number of successful solutions

% set the selection criteria for misfits
crit_misfit1 = 1;
crit_misfit2 = 3;
crit_misfit3 = 2;

success_oxygen =[];% empty matrix for saving successful solutions


for i = 1:ndata
    if (misfit_Mf(i)<crit_misfit1  && misfit_deltaO_lateArchean(i)<crit_misfit3 && ...
            misfit_deltaO_phanerozoic1(i)<crit_misfit3 && ...
            misfit_deltaO_phanerozoic2(i)<crit_misfit3 && ...
            misfit_deltaO_phanerozoic3(i)<crit_misfit2)
 
        j = j+1;
        for k = 1:nparam
            success_oxygen(j,k) = data(i,k);
        end % for k = 1:nparam
        
        % save all the variables and misfits
        fprintf(fileID,['%6g %6g %6g %6g %6g ' ...
            '%6g %6g %6g %6g %6g '  ...
            '%6g %6g '  ...
            '%6g %6g %6g %6g %6g '  ...
            '%6g %6g %6g %6g '  ...
            '%6g %6g %6g %6g %6g %6g %6g\n'], ...
            success_oxygen(j,1), success_oxygen(j,2), success_oxygen(j,3), success_oxygen(j,4),success_oxygen(j,5),...
            success_oxygen(j,6), success_oxygen(j,7), success_oxygen(j,8),success_oxygen(j,9), success_oxygen(j,10),...
            success_oxygen(j,11), success_oxygen(j,12), success_oxygen(j,13), success_oxygen(j,14),...
            success_oxygen(j,15), success_oxygen(j,16),success_oxygen(j,17), success_oxygen(j,18),success_oxygen(j,19),...
            success_oxygen(j,20), success_oxygen(j,21), success_oxygen(j,22), success_oxygen(j,23),...
            success_oxygen(j,24),success_oxygen(j,25), success_oxygen(j,26), success_oxygen(j,27), success_oxygen(j,28));
    end % if (misfit_pujol(i)...
    
end % for i = 1:ndata

fclose(fileID);
