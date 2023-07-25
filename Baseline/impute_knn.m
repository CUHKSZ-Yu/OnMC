function [imputedX] = impute_knn(Xoff, Xon, model, k)
% function [imputedX] = impute_knn(Xoff, Xon, model, k)
%
% Impute a data matrix. Each column is a sample. Each NaN value is replaced
% by the mean of the sample's k-nearest neighbors with known values. If all 
% k-nearest samples' corresponding features are NaN, then replaced by zero.
%
% @param Xoff       Offline dataset, each column is a complete sample
% @param Xon        Online dataset, each column is an incomplete sample
% @param model      Default 'on' (online version)
% @param k          Default 10 (k-nearest neighbors)
% 
% @return imputedX  Imputed matrix with all data samples
%
% <Reference>
% Ki-Yeol Kim, Byoung-Jin Kim, and Gwan-Su Yi. "Reuse of imputed data in 
% microarray analysis increases imputation efficiency." BMC Bioinformatics, 
% 5(1):1–9, 2004.

if (nargin < 4)
    k = 10;
end
if (nargin < 3)
    model = 'on';
end
low = min(Xoff(:)); high = max(Xoff(:));

noff = size(Xoff, 2);
non = size(Xon, 2);
imputedX = [Xoff, Xon];

if strcmp(model, 'on')
    % online update
    for i = 1 : non
        X = [Xoff, Xon(:, i)];
        Ximp = knnimpute(X, k);
        imputedX(:, noff+i) = Ximp(:, end);
    end
elseif strcmp(model, 'seq')
    % sequential update
    for i = 1 : non
        X = [Xoff, Xon(:, i)];
        Ximp = knnimpute(X, k);
        imputedX(:, noff+i) = Ximp(:, end);
        Xoff = [Xoff, Ximp(:, end)];
    end    
elseif strcmp(model, 'off')
    % offline update
    X = [Xoff, Xon];
    imputedX = knnimpute(X, k);
end
imputedX(imputedX<low) = low;
imputedX(imputedX>high) = high;
end