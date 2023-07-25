function [Spar] = correct_parallel(S0, Spar_0, method, parallel)
% function [Spar] = correct_parallel(S0, Spar_0, method, parallel)
%
% Correct multiple similarity vectors in parallel (See Section 3.2.2)
%
% @param  S0         Pairwise similarity matrix
% @param  Spar_0     Intial similarity vectors
% @param  method     Default SVD (RSVD is faster)
% @param  parallel   Default "parfor" for parallel computation
%
% @return Spar       Corrected similarity vectors
%
% <Reference>
% Fangchen Yu, Yicheng Zeng, Jianfeng Mao, and Wenye Li. "Online estimation 
% of similarity matrices with incomplete data." Uncertainty in Artificial 
% Intelligence. PMLR, 2023.

if (nargin < 4)
    parallel = 'parfor';
end
if (nargin < 3)
    method = 'svd';
end

n = size(Spar_0, 2);
tol = 1e-4;

if isempty(Spar_0)
    Spar = Spar_0;
else
    if strcmp(method, 'svd')
        % Standard Singular Value Decomposition
        [V, S, ~] = svd(S0);
    elseif strcmp(method, 'rsvd')
        % Randomized Singular Value Decomposition
        [V, S, ~] = rsvd(S0, 100);
    end

    s = diag(S);
    C = V * diag(sqrt(s));
    Cinv = diag(1./s) * C';
    U0 = Cinv * Spar_0;
    U = U0;
    c = S0(end, end);

    % We can use "parfor" to execute parallel correction 
    if strcmp(parallel, 'parfor')
        parfor i = 1 : n 
            u0 = U0(:, i);
            U(:, i) = correct_step_reuse(s, u0, c, tol);
        end
    elseif strcmp(parallel, 'for')
        for i = 1 : n 
            u0 = U0(:, i);
            U(:, i) = correct_step_reuse(s, u0, c, tol);
        end
    end
    Spar = C * U;
end

end
