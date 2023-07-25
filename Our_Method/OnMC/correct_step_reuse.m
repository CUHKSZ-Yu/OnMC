function [u] = correct_step_reuse(s, u0, c, tol)
% function [u] = correct_step_reuse(s, u0, c, tol)
% 
% Correct one similarity vector in the parallel correction method
%
% <Reference>
% Fangchen Yu, Yicheng Zeng, Jianfeng Mao, and Wenye Li. "Online estimation 
% of similarity matrices with incomplete data." Uncertainty in Artificial 
% Intelligence. PMLR, 2023.

if (nargin < 4)
    tol = 1e-4;
end
if (nargin < 3)
    c = 1;
end

if norm(u0)^2 <= c
    u = u0;
else
    lambda_min = 0;
    lambda_max = max(s) * norm(u0) / (2 * c^0.5);
    lambda = lambda_min;
    u = (s ./ (s+lambda)) .* u0;
    u_len = norm(u)^2;
    iter = 0;
    while (u_len > c) || (u_len < c-tol)
        iter = iter + 1;
        if iter > 50
            break
        end
        lambda = 0.5*(lambda_min + lambda_max);
        u = (s ./ (s+lambda)) .* u0;
        u_len = norm(u)^2;
        % use bisection search to find optimal lambda
        if u_len > c
            lambda_min = lambda;
        elseif u_len < c-tol
            lambda_max = lambda;
        end 
    end
end
end