function [S] = correct_onmc_s(S0, noff, non, method)
% function [S] = correct_onmc_s(S0, noff, non, method)
%
% Algorithm - Online Similarity Matrix Correction for Sequential Data (OnMC-S)
% Correct a similarity matrix on sequential data by the OnMC-S method 
% (see Algorithm 2 in Section 3.2.1)
%
% @param  S0      Initial similarity matrix
% @param  noff    Number of offline samples
% @param  non     Number of online samples
% @param  method  Default SVD (RSVD is faster)
%
% @return S       Corrected similarity matrix
%
% <Reference>
% Fangchen Yu, Yicheng Zeng, Jianfeng Mao, and Wenye Li. "Online estimation 
% of similarity matrices with incomplete data." Uncertainty in Artificial 
% Intelligence. PMLR, 2023.

if (nargin < 4)
    method = 'svd';
end

Sobs = S0(1:noff, 1:noff); % observed offline similarity matrix
Simp = Sobs;

for i = 1 : non
    c = S0(noff+i, noff+i);                   % default c=1
    v0 = S0(1:noff+i-1, noff+i);              % online similarity vector
    vonl = correct_step(Simp, v0, c, method); % corrected similarity vector
    Simp = [Simp, vonl; vonl', c];            % corrected similarity matrix
end
S = Simp; % corrected similarity matrix
end
