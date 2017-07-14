function [Y] = NotNaN(X)
% Keep only not NaN elements of a vector
Y = X(~isnan(X));
end
