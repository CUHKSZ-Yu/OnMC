function [imputedX] = impute_zero(X)
% function [imputedX] = impute_zero(X)
%
% Impute a data matrix. Each NaN value in a sample is replaced by zero.
%
% @param X  Each column is a sample.
% @return imputedX  Imputed data matrix

imputedX = X;
imputedX(isnan(X)) = 0;

end