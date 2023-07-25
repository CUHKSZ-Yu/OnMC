clear all; clc; warning off;
addpath(genpath('Baseline'));
addpath(genpath('Our_Method'));

load('demo_data.mat');
dataset = 'Demo_data'; 

%% Parameter Setting
n = size(X, 2);           % number of all data points in the dataset
r_list = [0.2, 0.5, 0.8]; % values of missing ratio r
niter = 5;                % number of iterations / repeated experiments
noff = 1000;              % number of offline data points
non = 100;                % number of online data points
seed = 2023;              % random seed for reproducing results
metric = 'cosine';        % cosine similarity metric

fprintf('\nUAI-2023 "Online Estimation of Similarity Matrices with Incomplete Data"');
fprintf('\nDemo: online scenario with incomplete data in Section 4.2\n');

%% Online Scenario with Incomplete Data
for k = 1 : length(r_list)
    r = r_list(k);
    fprintf(['\n',dataset]); fprintf(': noff = %1.0f, non = %1.0f, r = %1.2f, niter = %1.0f', noff, non, r, niter);
    
    %% Online Scenario
    for i = 1 : niter
        fprintf('\nIter = %1.0f: ', i);
        rng(seed + i);
        
        %% Data Construction
        rp = randperm(n);
        Xoff = X(:, rp(1:noff));
        Xon = X(:, rp(noff+1 : noff+non));
        Xtrue = [Xoff, Xon];
        Strue = similarity(Xtrue, metric);

        Xon(rand(size(Xon)) < r) = NaN;  
        Xmiss = [Xoff, Xon];
        Smiss = similarity(Xmiss, metric, 'miss');
        Fnorm = norm(Smiss-Strue, 'fro')^2;
        
        %% Online Process
        % ====================== ZERO Imputation ==========================
        fprintf('ZERO, ');
        tic; Xzero = impute_zero(Xmiss); time{k}(i,1) = toc;
        Szero = similarity(Xzero, metric);
        rmse{k}(i,1) = norm(Szero-Strue, 'fro')^2 / Fnorm;
        
        % ====================== MEAN Imputation ==========================
        fprintf('MEAN, ');
        tic; Xmean = impute_mean(Xoff, Xon); time{k}(i,2) = toc;
        Smean = similarity(Xmean, metric);
        rmse{k}(i,2) = norm(Smean-Strue, 'fro')^2 / Fnorm;

        % ====================== kNN Imputation ===========================
        fprintf('kNN, '); 
        tic; Xknn = impute_knn(Xoff, Xon); time{k}(i,3) = toc;
        Sknn = similarity(Xknn, metric);
        rmse{k}(i,3) = norm(Sknn-Strue, 'fro')^2 / Fnorm;
        
        % ====================== OnMC-S Correction ========================
        fprintf('OnMC-S, '); 
        tic; Sonmcs = correct_onmc_s(Smiss, noff, non); time{k}(i,4) = toc;
        rmse{k}(i,4) = norm(Sonmcs-Strue, 'fro')^2 / Fnorm;

        % ====================== OnMC-B Correction ========================
        fprintf('OnMC-B, '); 
        tic; Sonmcb = correct_onmc_b(Smiss, noff, non); time{k}(i,5) = toc;
        rmse{k}(i,5) = norm(Sonmcb-Strue, 'fro')^2 / Fnorm;
        
        fprintf('Finish.');
    end

    %%
    fprintf(['\n',dataset]); fprintf(': noff=%1.0f, non=%1.0f, r=%1.2f, niter=%1.0f\n', noff, non, r, niter);

    stat{k} = [mean(rmse{k}); mean(time{k})];
    Stat{k} = roundn(stat{k}, -3);
    Table{k} = array2table(Stat{k}, 'VariableNames',{'ZERO','MEAN','kNN','OnMC-S','OnMC-B'}, 'RowNames',{'RMSE';'Time'});
    disp(Table{k});
end
