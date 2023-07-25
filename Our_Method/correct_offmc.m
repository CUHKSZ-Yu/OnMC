function [S] = correct_offmc(S0, maxiter)
% function [S] = correct_offmc(S0, maxiter)
%
% Algorithm - Offline Similarity Matrix Correction (OffMC)
% Correct a similarity matrix by the OffMC method 
% (see Algorithm 1 in Section 3.1)
%
% @param  S0        Initial similarity matrix
% @param  maxiter   Maximum number of iterations
%
% @return S         Corrected similarity matrix
%
% <Reference>
% Wenye Li. "Estimating Jaccard index with missing observations: a matrix 
% calibration approach." In Advances in Neural Information Processing Systems,
% volume 28, pages 2620–2628, Canada, 2015.
% 
% Fangchen Yu, Yicheng Zeng, Jianfeng Mao, and Wenye Li. "Online estimation 
% of similarity matrices with incomplete data." Uncertainty in Artificial 
% Intelligence. PMLR, 2023.

if (nargin < 2)
    maxiter = 50;
end

S = nearpsd(S0, maxiter);

end

%%
function [Y] = nearpsd(A, maxiter, low, high)
% function [Y] = nearpsd(A, maxiter, low, high)
%
% Computes the nearest positive semi-definite matrix for a given square matrix.
%
% @param A        a square matrix to be corrected
% @param maxiter  max number of iterations, default 50
% @param low      default -1 for cosine similarity
% @param high     default +1 for cosine similarity
%
% @return Y       nearest psd matrix to A
% 
% <Reference>
% Wenye Li. "Estimating Jaccard index with missing observations: a matrix 
% calibration approach." In Advances in Neural Information Processing Systems,
% volume 28, pages 2620–2628, Canada, 2015.
% 
% Fangchen Yu, Yicheng Zeng, Jianfeng Mao, and Wenye Li. "Online estimation 
% of similarity matrices with incomplete data." Uncertainty in Artificial 
% Intelligence. PMLR, 2023.

if  ~isequal(A,A')
    A = (A + A') / 2;
end
if nargin < 4
    high = 1;
end
if nargin < 3
    low = -1; 
end
if nargin < 2
    maxiter = 50;
end

% threshold for convergence & eigs
tolconv = 1e-6;
toleigs = 1e-5;

n = size(A,1);

U = zeros(n);
Y = A;

[V, D] = eig(Y);
d = diag(D);

iter = 0;
while 1
    T = Y - U;

    % project onto psd matrices
    [Q, D] = eig(T);
    d = diag(D);
    p = d > toleigs*d(n);
    X = Q(:,p) * D(p,p) * Q(:,p)';

    % update correction
    U = X - T;

    % maximum iteration & convergence test
    iter = iter + 1;
    if iter == maxiter
        break; 
    end
    if norm(Y-X,'inf')/norm(Y,'inf') <= tolconv 
    	break;
    end
    
    % problem-dependent here
    Y = X;
    Y(Y<low) = low;
    Y(Y>high) = high;
end
Y(Y<low) = low;
Y(Y>high) = high;
end
