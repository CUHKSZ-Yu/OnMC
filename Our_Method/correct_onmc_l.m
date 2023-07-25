function [S] = correct_onmc_l(S0, noff, non, method, parallel, maxiter, koff, kon)
% function [S] = correct_onmc_l(S0, noff, non, method, parallel, maxiter, koff, kon)
%
% Algorithm - Online Similarity Matrix Correction for Large-scale Data (OnMC-L)
% Correct a similarity matrix on large-scale data by the OnMC-L method 
% (see Algorithm 4 in Section 3.2.3)
%
% @param  S0         Initial similarity matrix
% @param  noff       Number of offline samples
% @param  non        Number of online samples
% @param  method     Default SVD (RSVD is faster)
% @param  parallel   Default "parfor" for parallel computation
% @param  maxiter    Maximum number of iterations
% @param  koff       Partition size of offline sub-matrix
% @param  kon        Partition size of online sub-matrix
%
% @return S          Corrected similarity matrix
%
% <Reference>
% Fangchen Yu, Yicheng Zeng, Jianfeng Mao, and Wenye Li. "Online estimation 
% of similarity matrices with incomplete data." Uncertainty in Artificial 
% Intelligence. PMLR, 2023.

if (nargin < 8)
    kon = 1000;
end
if (nargin < 7)
    koff = 1000;
end
if (nargin < 6)
    maxiter = 50;
end
if (nargin < 5)
    parallel = 'parfor';
end
if (nargin < 4)
    method = 'svd';
end

npar_off = noff / koff;
Soff_par = S0(1:noff, noff+1:noff+non);
% We can use "parfor" to execute parallel correction. 
if strcmp(parallel, 'parfor')
    parfor i = 1 : npar_off
        scale = [(i-1)*koff+1 : i*koff];
        Soff = S0(scale, scale);
        Spar_0 = S0(scale, noff+1:noff+non);
        Spar{i} = correct_parallel(Soff, Spar_0, method, 'parfor');
    end
    for i = 1 : npar_off
        scale = [(i-1)*koff+1 : i*koff];
        Soff_par(scale, :) = Spar{i};
    end
elseif strcmp(parallel, 'for')
    for i = 1 : npar_off
        scale = [(i-1)*koff+1 : i*koff];
        Soff = S0(scale, scale);
        Spar_0 = S0(scale, noff+1:noff+non);
        Spar = correct_parallel(Soff, Spar_0, method, 'for');
        Soff_par(scale, :) = Spar;
    end
end

npar_on = non / kon;
Son = zeros(non, non);
Son_0 = S0(noff+1:noff+non, noff+1:noff+non);
% We can use "parfor" to execute parallel correction. 
if strcmp(parallel, 'parfor')
    parfor j = 1 : npar_on
        scale = [(j-1)*kon+1 : j*kon];
        Son_ini{j} = Son_0(scale, scale);
        Son_new{j} = correct_offmc(Son_ini{j}, maxiter);
        Spar_ini = Son_0(scale, j*kon+1:non);
        Spar_new{j} = correct_parallel(Son_new{j}, Spar_ini, method, 'parfor');
    end
    for j = 1 : npar_on
        scale = [(j-1)*kon+1 : j*kon];
        Son(scale, (j-1)*kon+1:non) = [Son_new{j}, Spar_new{j}];
        Son(j*kon+1:non, scale) = Spar_new{j}';
    end
elseif strcmp(parallel, 'for')
    for j = 1 : npar_on
        scale = [(j-1)*kon+1 : j*kon];
        Son_ini = Son_0(scale, scale);
        Son_new = correct_offmc(Son_ini, maxiter);
        Spar_ini = Son_0(scale, j*kon+1:non);
        Spar_new = correct_parallel(Son_new, Spar_ini, method, 'for');
        Son(scale, (j-1)*kon+1:non) = [Son_new, Spar_new];
        Son(j*kon+1:non, scale) = Spar_new';
    end
end

S = [S0(1:noff, 1:noff), Soff_par; Soff_par', Son];

end


