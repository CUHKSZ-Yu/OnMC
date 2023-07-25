clear all; clc;
addpath(genpath('Our_Method'));

load('demo_data.mat');
dataset = 'Demo_data'; 

%% Parameter Setting
n = size(X, 2);           % number of all data points in the dataset
r_list = [0.2, 0.5, 0.8]; % values of missing ratio p
niter = 5;                % number of iterations
noff = 2000;              % number of offline samples
non = 2000;               % number of online samples
seed = 2023;              % random seed for reproducing results
metric = 'cosine';        % cosine similarity metric

fprintf('\nUAI-2023 "Online Estimation of Similarity Matrices with Incomplete Data"');
fprintf('\nDemo: scalability analysis in Section 4.5\n');

%% Scalability Analysis
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

        %% Matrix Correction
        % ======================== OnMC-B Correction ======================
        fprintf('B, '); 
        tic; S_b1 = correct_onmc_b(Smiss, noff, non, 'svd', 'parfor'); 
        time{k}(i,1) = toc;
        rmse{k}(i,1) = norm(S_b1-Strue, 'fro')^2 / Fnorm;
        
        fprintf('B-RSVD, '); 
        tic; S_b2 = correct_onmc_b(Smiss, noff, non, 'rsvd', 'parfor'); 
        time{k}(i,2) = toc;
        rmse{k}(i,2) = norm(S_b2-Strue, 'fro')^2 / Fnorm;
        
        % ======================== OnMC-L Correction ======================
        fprintf('L, '); 
        tic; S_l1 = correct_onmc_l(Smiss, noff, non, 'svd', 'parfor'); 
        time{k}(i,3) = toc;
        rmse{k}(i,3) = norm(S_l1-Strue, 'fro')^2 / Fnorm;
        
        fprintf('L-RSVD, '); 
        tic; S_l2 = correct_onmc_l(Smiss, noff, non, 'rsvd', 'parfor'); 
        time{k}(i,4) = toc;
        rmse{k}(i,4) = norm(S_l2-Strue, 'fro')^2 / Fnorm;
        
        fprintf('Finish.');
    end
    
    %%
    fprintf(['\n',dataset]); fprintf(' Large-scale: noff=%1.0f, non=%1.0f, r=%1.2f, niter=%1.0f\n', noff, non, r, niter);

    stat{k} = [mean(rmse{k}); mean(time{k})];
    Stat{k} = roundn(stat{k}, -3);
    Table{k} = array2table(Stat{k}, 'VariableNames',{'B','B-RSVD','L','L-RSVD'}, 'RowNames',{'RMSE';'Time'});
    disp(Table{k});
end


