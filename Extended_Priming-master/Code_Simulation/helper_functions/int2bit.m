function Y = int2bit(X,n)
% Summary:
%   Converts integers into n-digit binary. See also bit2int
% Inputs:
%   X: vector of non-negative integers, either a row or a column
%   n: positive integer
% Output:
%   Y: Matrix of size n-by-length of X. Each column corresponds to a
%   n-digit binary number that corresponds to each element of X.
if isrow(X)
    X = X';
end
Y = rem(floor(X./(2.^(n-1:-1:0))),2)';
end