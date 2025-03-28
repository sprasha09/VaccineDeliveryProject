function arr = alignLeft(arr, trim)
%% Documentation
% Summary
% Helper function that takes in an array as input and then aligns non-zero 
% columns to the left.  

% Inputs
%  arr: 2-D array. Each column must be either all-zeros or non-zeros. 
%  trim: 0 or 1. If 0, the array is trimed at the last non-zero column

% Outputs
%  arr: Left-aligned and/or trimmed array

%%
temp = find(arr(1,:));
if trim
    arr = arr(:,temp);
else
    arr(:,1:length(temp)) = arr(:,temp);
    arr(:,length(temp)+1:end) = 0;
end
end