function [imputedX] = impute_mean(Xoff, Xon, model)
% function [imputedX] = impute_mean(Xoff, Xon, model)
%
% Replace NaN values by row means.
%
% @param Xoff       Offline dataset, each column is a complete sample
% @param Xon        Online dataset, each column is an incomplete sample
% @param model      Default 'on' (online version)
% 
% @return imputedX  Imputed matrix with all data samples

if (nargin < 3)
    model = 'on';
end
low = min(Xoff(:)); high = max(Xoff(:));

noff = size(Xoff, 2);
non = size(Xon, 2);

if strcmp(model, 'on')
    % online update
    imputedX = [Xoff, Xon];
    rowmean = mean(Xoff, 2);
    for i = 1 : non
        x = Xon(:, i);
        idx = isnan(x);
        imputedX(idx, noff+i) = rowmean(idx);
    end
elseif strcmp(model, 'seq')
    % sequential update
    for i = 1 : non
        rowmean = mean(Xoff, 2);
        x = Xon(:, i);
        idx = isnan(x);
        x(idx) = rowmean(idx);
        Xoff = [Xoff, x];
    end
    imputedX = Xoff;
elseif strcmp(model, 'off')
    % offline update
    Ximp = [Xoff, Xon];
    idx = isnan(Ximp);
    imputedX = Ximp;
    imputedX(idx) = 0;
    imputedX = imputedX + nanmean(Ximp,2).*idx;
end
imputedX(imputedX<low) = low;
imputedX(imputedX>high) = high;
end