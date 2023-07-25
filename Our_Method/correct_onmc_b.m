function [S] = correct_onmc_b(S0, noff, non, method, parallel, maxiter)
% function [S] = correct_onmc_b(S0, noff, non, method, parallel, maxiter)
%
% Algorithm - Online Similarity Matrix Correction for Batch Data (OnMC-B)
% Correct a similarity matrix on batch data by the OnMC-B method 
% (see Algorithm 3 in Section 3.2.2)
%
% @param  S0         Initial similarity matrix
% @param  noff       Number of offline samples
% @param  non        Number of online samples
% @param  method     Default SVD (RSVD is faster)
% @param  parallel   Default "parfor" for parallel computation
% @param  maxiter    Maximum number of iterations
%
% @return S          Corrected similarity matrix
%
% <Reference>
% Fangchen Yu, Yicheng Zeng, Jianfeng Mao, and Wenye Li. "Online estimation 
% of similarity matrices with incomplete data." Uncertainty in Artificial 
% Intelligence. PMLR, 2023.

if (nargin < 6)
    maxiter = 50;
end
if (nargin < 5)
    parallel = 'parfor';
end
if (nargin < 4)
    method = 'svd';
end

% use the OffMC algorithm to correct Son
Son_0 = S0(noff+1:noff+non, noff+1:noff+non);
Son = correct_offmc(Son_0, maxiter);

Soff_0 = S0(1:noff, 1:noff);
Soff = Soff_0; % if Soff_0 is already a PSD matrix

% if min(eig(Soff_0)) >= -tol
%     Soff = Soff_0;
% else
%     Soff = correct_offmc(Soff_0, maxiter);
% end

% to correct Spar_0 by parallel correction
Spar_0 = S0(1:noff, noff+1:noff+non);
Spar = correct_parallel(Soff, Spar_0, method, parallel);

S = [Soff, Spar; Spar', Son];

end


