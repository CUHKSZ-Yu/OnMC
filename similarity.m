function [S] = similarity(X, metric, type)
% function [S] = similarity(X, metric, type)
%
% Approximate a similarity matrix for incomplete samples with NaN values.
% (see Eqs.1 & 2 in Section 2.1)
%
% @param  X        Data matrix of size d*n, each column is a sample
% @param  metric   Default "cosine" (cosine similarity)
% @param  type     Default "true"   (calculate the true similarity matrix)
%
% @return S        Similarity matrix of size n*n
%
% <Reference>
% Fangchen Yu, Yicheng Zeng, Jianfeng Mao, and Wenye Li. "Online estimation 
% of similarity matrices with incomplete data." Uncertainty in Artificial 
% Intelligence. PMLR, 2023.

if (nargin < 3)
    type = 'true';
end
if (nargin < 2)
    metric = 'cosine';
end

if strcmp(metric, 'cosine')
    %% Cosine Similarity
    if strcmp(type, 'true')
        % if X is complete
        Norm = sqrt(sum(X.^2, 1));
        S = (X' * X) ./ (Norm' * Norm);
    elseif strcmp(type, 'miss')
        % if X is incomplete
        [d, n] = size(X);
        Idx = isnan(X);
        Xzero = X; Xzero(Idx) = 0;
        XX = Xzero' * Xzero;
        Norm = zeros(n);
        Norm_ini = sqrt(sum(Xzero.^2, 1));
        for i = 1 : n
            idx = Idx(:, i);
            if sum(idx) == 0
                Norm(i, :) = Norm_ini;
            else
                idx = repmat(idx, 1, n);
                XI = Xzero; XI(idx) = 0;
                Norm(i, :) = sqrt(sum(XI.^2, 1));
            end
        end
        S = XX ./ (Norm .* Norm');
        S(isnan(S)) = 0;
    elseif strcmp(type, 'basic')
        [d, n] = size(X);
        O =~isnan(X);
        S = zeros(n);
        for i = 1 : n
            for j = i+1 : n
                k = O(:,i) & O(:,j);
                S(i,j) = X(k,i)'*X(k,j) / (norm(X(k,i))*norm(X(k,j)));
            end
        end
        S = S + S' + eye(n);
        S(isnan(S)) = 0;
    end
else
    fprintf('It can be applied to other similarity metrics.')
end

end