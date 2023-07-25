function [v] = correct_step(S0, v0, c, method)
% function [v] = correct_step(S0, v0, c, method)
%
% Correct one similarity vector by the One-Step Matrix Correction method (see Section 3.2.1)
%
% @param  S0      Intial similarity matrix
% @param  v0      Intial similarity vector
% @param  c       Similarity value of itself (default c=1)
% @param  method  Default SVD (RSVD is faster)
%
% @return v       Corrected similarity vector
%
% <Reference>
% Fangchen Yu, Yicheng Zeng, Jianfeng Mao, and Wenye Li. "Online estimation 
% of similarity matrices with incomplete data." Uncertainty in Artificial 
% Intelligence. PMLR, 2023.

if (nargin < 4)
    method = 'svd';
end
if (nargin < 3)
    c = 1;
end

n = size(S0, 1); % number of samples
tol = 1e-4;      % convergence tolerance

if strcmp(method, 'svd')
    % Standard Singular Value Decomposition
    [U, S, ~] = svd(S0);
elseif strcmp(method, 'rsvd')
    % Randomized Singular Value Decomposition
    [U, S, ~] = rsvd(S0, 100);
end

s = diag(S);
C = U * diag(sqrt(s));
Cinv = diag(1./s) * C';

y0 = Cinv * v0;
if norm(y0)^2 <= c
    % to test if y0 is a feasible solution
    ycal = y0;
else
    % to correct y0 by solving Eq.6 in Section 3.2.1
    lambda_min = 0;
    lambda_max = max(s) * norm(y0) / (2 * c^0.5);
    lambda = lambda_min;
    ycal = (s ./ (s+lambda)) .* y0;
    ylen = norm(ycal)^2;
    iter = 0;
    while (ylen > c) || (ylen < c-tol)
        iter = iter + 1;
        if iter > 50
            break
        end
        lambda = 0.5*(lambda_min + lambda_max);
        ycal = (s ./ (s+lambda)) .* y0;
        ylen = norm(ycal)^2;
        % use bisection search to find optimal lambda
        if ylen > c
            lambda_min = lambda;
        elseif ylen < c-tol
            lambda_max = lambda;
        end 
    end
end
v = C * ycal; % corrected similarity vector
end