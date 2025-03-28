function Y = bit2int(X,n)
%% Documentation
% Summary:
%   Converts n-digit binaries into integers. See also int2bit

% Inputs:
%   X: Matrix whose columns represent n-digit binary numbers
%   n: positive integer

% Output:
%   Y: Row vector equal to the number of columns of X, whose elements are
%   integers corresponding to the columns of X.
%%
[s1,s2] = size(X);
if n>s1
    error('n is greater than vector length')
end
Y = zeros(1,s2);
for i=1:s2
   Y(i) = sum(2.^(n-1:-1:0).*(X(1:n,i)')); 
end
end

